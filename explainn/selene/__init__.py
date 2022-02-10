"""
Adapted from Selene:
https://github.com/FunctionLab/selene/blob/master/selene_sdk/train_model.py
"""

import copy
from functools import partial
import logging
import numpy as np
import os
import sys
from time import time
import torch

class Trainer(object):
    """
    This class ties together the various objects and methods needed to
    train and validate a model.

    `Trainer` saves a checkpoint model (overwriting it after
    `save_checkpoint_every_n_steps`) as well as a best-performing model
    (overwriting it after `report_stats_every_n_steps` if the latest
    validation performance is better than the previous best-performing
    model) to `output_dir`.

    `Trainer` also outputs 2 files that can be used to monitor training
    as Selene runs: `selene_sdk.train_model.train.txt` (training loss) and
    `selene_sdk.train_model.validation.txt` (validation loss & average
    ROC AUC). The columns in these files can be used to quickly visualize
    training history (e.g. you can use `matplotlib`, `plt.plot(auc_list)`)
    and see, for example, whether the model is still improving, if there are
    signs of overfitting, etc.

    Parameters
    ----------
    model : torch.nn.Module
        The model to train.
    data_loaders : dict
        ...
    loss_criterion : torch.nn._Loss
        The loss function to optimize.
    metrics : dict
        ...
    optimizer : torch.optim.Optimizer
        The optimizer to minimize loss with.
    max_steps : int
        The maximum number of mini-batches to iterate over.
    report_stats_every_n_steps : int
        The frequency with which to report summary statistics. You can
        set this value to be equivalent to a training epoch
        (`n_steps * batch_size`) being the total number of samples
        seen by the model so far. Selene evaluates the model on the validation
        dataset every `report_stats_every_n_steps` and, if the model obtains
        the best performance so far (based on the user-specified loss function),
        Selene saves the model state to a file called `best_model.pth.tar` in
        `output_dir`.
    output_dir : str
        The output directory to save model checkpoints and logs in.
    save_checkpoint_every_n_steps : int or None, optional
        Default is 1000. If None, set to the same value as
        `report_stats_every_n_steps`
    save_new_checkpoints_after_n_steps : int or None, optional
        Default is None. The number of steps after which Selene will
        continually save new checkpoint model weights files
        (`checkpoint-<TIMESTAMP>.pth.tar`) every
        `save_checkpoint_every_n_steps`. Before this point,
        the file `checkpoint.pth.tar` is overwritten every
        `save_checkpoint_every_n_steps` to limit the memory requirements.
    n_validation_samples : int or None, optional
        Default is `None`. Specify the number of validation samples in the
        validation set. If `n_validation_samples` is `None` and the data sampler
        used is the `selene_sdk.samplers.IntervalsSampler` or
        `selene_sdk.samplers.RandomSampler`, we will retrieve 32000
        validation samples. If `None` and using
        `selene_sdk.samplers.MultiSampler`, we will use all
        available validation samples from the appropriate data file.
    n_test_samples : int or None, optional
        Default is `None`. Specify the number of test samples in the test set.
        If `n_test_samples` is `None` and
            - the sampler you specified has no test partition, you should not
              specify `evaluate` as one of the operations in the `ops` list.
              That is, Selene will not automatically evaluate your trained
              model on a test dataset, because the sampler you are using does
              not have any test data.
            - the sampler you use is of type `selene_sdk.samplers.OnlineSampler`
              (and the test partition exists), we will retrieve 640000 test
              samples.
            - the sampler you use is of type
              `selene_sdk.samplers.MultiSampler` (and the test partition
              exists), we will use all the test samples available in the
              appropriate data file.
    cpu_n_threads : int, optional
        Default is 1. Sets the number of OpenMP threads used for parallelizing
        CPU operations.
    use_cuda : bool, optional
        Default is `False`. Specify whether a CUDA-enabled GPU is available
        for torch to use during training.
    # data_parallel : bool, optional
    #     Default is `False`. Specify whether multiple GPUs are available
    #     for torch to use during training.
    logging_verbosity : {0, 1, 2}, optional
        Default is 2. Set the logging verbosity level.
            * 0 - Only warnings will be logged.
            * 1 - Information and warnings will be logged.
            * 2 - Debug messages, information, and warnings will all be
                  logged.
    checkpoint_resume : str or None, optional
        Default is `None`. If `checkpoint_resume` is not None, it should be the
        path to a model file generated by `torch.save` that can now be read
        using `torch.load`.
    # use_scheduler : bool, optional
    #     Default is `True`. If `True`, learning rate scheduler is used to
    #     reduce learning rate on plateau. PyTorch ReduceLROnPlateau scheduler 
    #     with patience=16 and factor=0.8 is used.

    Attributes
    ----------
    model : torch.nn.Module
        The model to train.
    sampler : selene_sdk.samplers.Sampler
        The example generator.
    criterion : torch.nn._Loss
        The loss function to optimize.
    optimizer : torch.optim.Optimizer
        The optimizer to minimize loss with.
    max_steps : int
        The maximum number of mini-batches to iterate over.
    nth_step_report_stats : int
        The frequency with which to report summary statistics.
    use_cuda : bool
        If `True`, use a CUDA-enabled GPU. If `False`, use the CPU.
    # data_parallel : bool
    #     Whether to use multiple GPUs or not.
    output_dir : str
        The directory to save model checkpoints and logs.
    """

    def __init__(self,
                 model,
                 data_loaders,
                 loss_criterion,
                 metrics,
                 optimizer,
                 max_steps=128000,
                 patience=32000,
                 report_stats_every_n_steps=1000,
                 output_dir="./",
                 cpu_n_threads=1,
                 use_cuda=False,
                 checkpoint_resume=None,
                 freeze_top_n_filters=0,
                 logging_verbosity=2,
                ):
        """
        Constructs a new `Trainer` object.
        """
        self.model = model
        self.data_loaders = data_loaders
        self.criterion = loss_criterion
        self.metrics = metrics
        self.optimizer = optimizer
        self.max_steps = max_steps
        self.patience = patience
        self.nth_step_report_stats = report_stats_every_n_steps

        torch.set_num_threads(cpu_n_threads)

        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir

        # Freeze filters for fine-tuning
        self.freeze_top_n_filters = freeze_top_n_filters

        self.logger = selene_logger(self.output_dir, logging_verbosity)

        self.use_cuda = use_cuda
        if self.use_cuda:
            self.model.cuda()
            self.criterion.cuda()
            self.logger.debug("Set modules to use CUDA")

        self._data_iterators = {}
        for k in self.data_loaders.keys():
            self._data_iterators.setdefault(k, [])

        self._init_train()
        self._init_validate()
        if checkpoint_resume is not None:
            self._load_checkpoint(checkpoint_resume)

    def _init_train(self):
        self._start_step = 1
        self._train_logger = metrics_logger("train", self.output_dir)
        self.logger.info("Training metrics: loss")
        self._train_logger.log(10, "loss")
        self._time_per_step = []
        self._train_loss = []

    def _init_validate(self):
        self._min_loss = float("inf")
        self._best_step = 1
        self._validation_logger = metrics_logger("validation", self.output_dir)
        metrics_str = ", ".join(["loss"] + [x for x in self.metrics.keys()])
        self.logger.info(f"Validation metrics: {metrics_str}")
        metrics_str = "\t".join(["loss"] + [x for x in self.metrics.keys()])
        self._validation_logger.log(10, metrics_str)

    def _load_checkpoint(self, checkpoint_resume):
        checkpoint = torch.load(checkpoint_resume)
        self.model.load_state_dict(checkpoint["state_dict"])
        self._start_step = checkpoint["step"]
        self._min_loss = checkpoint["min_loss"]
        self._best_step = checkpoint["step"]
        self.optimizer.load_state_dict(checkpoint["optimizer"])
        if self.use_cuda:
            for state in self.optimizer.state.values():
                for k, v in state.items():
                    if isinstance(v, torch.Tensor):
                        state[k] = v.cuda()
        self.logger.info(f"Resuming from checkpoint: step {self._start_step}" +
                         f", min loss {self._min_loss}")

    def _get_batch(self, which_data):
        """
        Fetches a mini-batch of examples

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            A tuple containing the examples and targets.
        """
        t_i_sampling = time()
        try:
            batch_sequences, batch_targets = \
                next(self._data_iterators[which_data])
        except:
            self._data_iterators[which_data] = \
                iter(self.data_loaders[which_data])
            batch_sequences, batch_targets = \
                next(self._data_iterators[which_data])
        t_f_sampling = time()

        t = t_f_sampling - t_i_sampling
        self.logger.debug(f"[BATCH] Time to sample batch: {t} s.")

        return(batch_sequences, batch_targets)

    def train_and_validate(self):
        """
        Trains the model and measures validation performance.
        """

        # Freeze filters for fine-tuning
        if self.freeze_top_n_filters > 0:
            p = partial(freeze_top_n_filters, use_cuda=self.use_cuda,
                freeze_top_n_filters=self.freeze_top_n_filters)
            self.model.linears[0].weight.register_hook(p)

        for step in range(self._start_step, self.max_steps + 1):
            self.step = step
            self.train()
            if self.step % self.nth_step_report_stats == 0:
                self.validate()
            if self.step >= self._best_step + self.patience:
                self.logger.info(f"Early stopping: stop training!")
                break

        # Reset loggers
        self.logger.handlers.clear()
        self._train_logger.handlers.clear()
        self._validation_logger.handlers.clear()

    def train(self):
        """
        Trains the model on a batch of data.

        Returns
        -------
        float
            The training loss.
        """
        t_i = time()
        self.model.train()
        inputs, targets = self._get_batch("train")
        if self.use_cuda:
            inputs = inputs.cuda()
            targets = targets.cuda()
        predictions = self.model(inputs)
        loss = self.criterion(predictions, targets)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        # if self.model._options["clamp_weights"]:
        #     self.model.final.weight.data.clamp_(0)
        self._train_loss.append(loss.item())
        t_f = time()

        self._time_per_step.append(t_f - t_i)
        if (self.step > 0 and \
            self.step % self.nth_step_report_stats == 0):
            self.logger.info(
                f"[STEP {self.step}] average number of steps per second" +
                f": {1. / np.average(self._time_per_step)}"
            )
            self.logger.info(f"Training loss: {np.average(self._train_loss)}")
            self._train_logger.log(10, np.average(self._train_loss))
            self._time_per_step = []
            self._train_loss = []

    def _evaluate_on_data(self, which_data):
        """
        Makes predictions for some labeled input data.

        Parameters
        ----------
        data_in_batches : list(tuple(numpy.ndarray, numpy.ndarray))
            A list of tuples of the data, where the first element is
            the example, and the second element is the label.

        Returns
        -------
        tuple(float, list(numpy.ndarray))
            Returns the average loss, and the list of all predictions.
        """
        self.model.eval()
        batch_losses = []
        all_predictions = []
        all_targets = []
        for inputs, targets in iter(self.data_loaders[which_data]):
            if self.use_cuda:
                inputs = inputs.cuda()
                targets = targets.cuda()
            with torch.no_grad():
                predictions = self.model(inputs)
                loss = self.criterion(predictions, targets)
                all_predictions.append(predictions.data.cpu().numpy())
                batch_losses.append(loss.item())
            all_targets.extend(targets.data.cpu().numpy())
        all_predictions = np.vstack(all_predictions)
        all_targets = np.vstack(all_targets)

        return(np.average(batch_losses), all_predictions, all_targets)

    def validate(self):
        """
        Measures model validation performance.

        Returns
        -------
        dict
            A dictionary, where keys are the names of the loss metrics,
            and the values are the average value for that metric over
            the validation set.
        """
        validation_loss, all_predictions, all_targets = \
            self._evaluate_on_data("validation")
        self.logger.info(f"Validation loss: {validation_loss}")
        valid_scores = {}
        for metric in self.metrics:
            y_true = all_targets.flatten()
            y_score = all_predictions.flatten()
            score = self.metrics[metric](y_true, y_score)
            if isinstance(score, float):
                valid_scores.setdefault(metric, score)
            else:
                valid_scores.setdefault(metric, score[0])
        for name, score in valid_scores.items():
            self.logger.info(f"Validation {name}: {score}")
        self._validation_metrics = valid_scores

        # Save best_model
        if validation_loss < self._min_loss:
            self._min_loss = validation_loss
            self._best_step = int(self.step)
            model = copy.deepcopy(self.model)
            checkpoint = {
                "step": self._best_step,
                "arch": model.__class__.__name__,
                "options": model._options,
                "state_dict": model.state_dict(),
                "min_loss": self._min_loss,
                "optimizer": copy.deepcopy(self.optimizer.state_dict()),
            }
            self._save_checkpoint(checkpoint)
            self.logger.info("Updating `best_model.pth.tar`")

        # Logging
        log_message = [validation_loss]
        for metric in valid_scores:
            log_message.append(valid_scores[metric])
        self._validation_logger.log(10, "\t".join(map(str, log_message)))

    def _save_checkpoint(self, state):
        """
        Saves snapshot of the model state to file. Will save a checkpoint
        with name `<filename>.pth.tar` and, if this is the model's best
        performance so far, will save the state to a `best_model.pth.tar`
        file as well.
        Models are saved in the state dictionary format. This is a more
        stable format compared to saving the whole model (which is another
        option supported by PyTorch). Note that we do save a number of
        additional, Selene-specific parameters in the dictionary
        and that the actual `model.state_dict()` is stored in the `state_dict`
        key of the dictionary loaded by `torch.load`.
        See: https://pytorch.org/docs/stable/notes/serialization.html for more
        information about how models are saved in PyTorch.

        Parameters
        ----------
        state : dict
            Information about the state of the model. Note that this is
            not `model.state_dict()`, but rather, a dictionary containing
            keys that can be used for continued training in Selene
            _in addition_ to a key `state_dict` that contains
            `model.state_dict()`.
        # filename : str, optional
        #     Default is "checkpoint". Specify the checkpoint filename. Will
        #     append a file extension to the end of the `filename`
        #     (e.g. `checkpoint.pth.tar`).

        Returns
        -------
        None
        """
        s = state["step"]
        self.logger.debug(f"[TRAIN] {s}: Saving model state to file.")
        f = os.path.join(self.output_dir, "best_model")
        torch.save(state, "{0}.pth.tar".format(f))

def selene_logger(logger_path="./", verbosity=2):
    """
    Initializes the logger for Selene. This function can only
    be called successfully once. If the logger has already been
    initialized with handlers, the function exits. Otherwise,
    it proceeds to set the logger configurations.

    Parameters
    ----------
    output_path : str
        The path to the output file where logs will be written.
    verbosity : int, {2, 1, 0}
        Default is 2. The level of logging verbosity to use.
            * 0 - Only log warnings
            * 1 - Log information and warnings
            * 2 - Log debug messages, information, and warnings
    """

    logger = logging.getLogger("selene")
    if len(logger.handlers):
        return

    if verbosity == 0:
        logger.setLevel(logging.WARNING)
    elif verbosity == 1:
        logger.setLevel(logging.INFO)
    elif verbosity == 2:
        logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handle = logging.FileHandler(os.path.join(logger_path, "selene.log"))
    handle.setFormatter(formatter)
    logger.addHandler(handle)

    formatter = logging.Formatter("%(asctime)s - %(message)s")
    stdout_handle = logging.StreamHandler(sys.stdout)
    stdout_handle.setFormatter(formatter)
    stdout_handle.setLevel(logging.INFO)
    logger.addHandler(stdout_handle)

    return(logger)

def metrics_logger(metric, logger_path="./", verbosity=2):
    """
    Parameters
    ----------
    output_path : str
        The path to the output file where logs will be written.
    verbosity : int, {2, 1, 0}
        Default is 2. The level of logging verbosity to use.
            * 0 - Only log warnings
            * 1 - Log information and warnings
            * 2 - Log debug messages, information, and warnings
    """

    logger = logging.getLogger(metric)
    if len(logger.handlers):
        return

    if verbosity == 0:
        logger.setLevel(logging.WARNING)
    elif verbosity == 1:
        logger.setLevel(logging.INFO)
    elif verbosity == 2:
        logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(message)s")
    handle = logging.FileHandler(os.path.join(logger_path, f"{metric}.txt"))
    handle.setFormatter(formatter)
    logger.addHandler(handle)

    formatter = logging.Formatter("%(asctime)s - %(message)s")
    stdout_handle = logging.StreamHandler(sys.stdout)
    stdout_handle.setFormatter(formatter)
    stdout_handle.setLevel(logging.INFO)
    logger.addHandler(stdout_handle)

    return(logger)

def freeze_top_n_filters(grad, freeze_top_n_filters=0, use_cuda=False):
    nulls = torch.zeros(grad.shape)
    if use_cuda:
        nulls = nulls.cuda()
    nulls[freeze_top_n_filters:, :,:] = 1

    return(grad*nulls)