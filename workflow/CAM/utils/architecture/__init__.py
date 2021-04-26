from collections import OrderedDict
import math
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings("ignore")

class ExpAct(nn.Module):
    """Exponential Activation."""

    def _init_(self):
        super(ExpAct, self)._init_()

    def forward(self, x):
        return(torch.exp(x))

class UnSqueeze(torch.nn.Module):
    def forward(self, x):
        return(x.unsqueeze(-1))

class __Model(nn.Module):

    def load_weights(self, weight_file):
        sd = torch.load(weight_file)
        ord_dict = OrderedDict()
        keys = list(self.state_dict().keys())
        values = list(sd.values())
        for i in range(len(values)):
            v = values[i]
            if v.dim() > 1:
                if v.shape[-1] == 1:
                    ord_dict[keys[i]] = v.squeeze(-1)
                    continue
            ord_dict[keys[i]] = v
        self.load_state_dict(ord_dict)

class CAM(__Model):
    """Convolutional Additive Model (CAM)."""

    def __init__(self, cnn_units, motif_length, seq_length, n_features=1,
        weights_file=None):
        """
        Parameters
        ----------
        cnn_units : int
            Total number of individual CNN units
        motif_length : int
            Max. 
        seq_length : int
            Input sequence length
        n_features : int
            Total number of individual CNN units
        weights_file : pass
            ...
        """
        super(CAM, self).__init__()

        self._cnn_units = cnn_units
        self._motif_length = motif_length
        self._padding = self._motif_length
        self._n_channels = math.floor(
            (seq_length + 2*self._padding - (self._motif_length-1)) / 7.
        )

        self.linears = nn.Sequential(
            nn.Conv1d(
                in_channels=4*self._cnn_units,
                out_channels=1*self._cnn_units,
                kernel_size=self._motif_length,
                padding=self._padding,
                groups=self._cnn_units,
            ),
            nn.BatchNorm1d(self._cnn_units),
            ExpAct(),
            nn.MaxPool1d(7, 7),
            nn.Flatten(), 
            UnSqueeze(),
            nn.Conv1d(
                in_channels=self._n_channels*self._cnn_units,
                out_channels=100*self._cnn_units,
                kernel_size=1,
                groups=self._cnn_units,
            ),
            nn.BatchNorm1d(100*self._cnn_units, 1e-05, 0.1, True),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Conv1d(
                in_channels=100*self._cnn_units,
                out_channels=1*self._cnn_units,
                kernel_size=1,
                groups=self._cnn_units),
            nn.BatchNorm1d(1*self._cnn_units, 1e-05, 0.1, True),
            nn.ReLU(),
            nn.Flatten(),
        )

        self.final = nn.Linear(self._cnn_units, n_features)

        if weights_file is not None:
            self.load_weights(weights_file)

    def forward(self, x):
        # ModuleList can act as an iterable, or be indexed using ints
        x = x.repeat(1, self._cnn_units, 1)
        outs = self.linears(x)
        out = self.final(outs)

        return(out)

def _flip(x, dim):
    """
    Adapted from Selene:
    https://github.com/FunctionLab/selene/blob/master/selene_sdk/utils/non_strand_specific_module.py

    Reverses the elements in a given dimension `dim` of the Tensor.
    source: https://github.com/pytorch/pytorch/issues/229
    """
    xsize = x.size()
    dim = x.dim() + dim if dim < 0 else dim
    x = x.contiguous()
    x = x.view(-1, *xsize[dim:])
    x = x.view(
        x.size(0), x.size(1), -1)[:, getattr(
            torch.arange(x.size(1)-1, -1, -1),
            ("cpu","cuda")[x.is_cuda])().long(), :]

    return x.view(xsize)

class NonStrandSpecific(nn.Module):
    """
    Adapted from Selene:
    https://github.com/FunctionLab/selene/blob/master/selene_sdk/utils/non_strand_specific_module.py

    A torch.nn.Module that wraps a user-specified model architecture if the
    architecture does not need to account for sequence strand-specificity.

    Parameters
    ----------
    model : torch.nn.Module
        The user-specified model architecture.
    mode : {'mean', 'max'}, optional
        Default is 'mean'. NonStrandSpecific will pass the input and the
        reverse-complement of the input into `model`. The mode specifies
        whether we should output the mean or max of the predictions as
        the non-strand specific prediction.

    Attributes
    ----------
    model : torch.nn.Module
        The user-specified model architecture.
    mode : {'mean', 'max'}
        How to handle outputting a non-strand specific prediction.
    """

    def __init__(self, model):
        super(NonStrandSpecific, self).__init__()

        self.model = model

    def forward(self, input):
        reverse_input = None
        reverse_input = _flip(_flip(input, 1), 2)

        output = self.model.forward(input)
        output_from_rev = self.model.forward(reverse_input)

        return((output + output_from_rev) / 2)

def get_loss_criterion(data="binary"):
    """
    Specify the appropriate loss function (criterion) for this model.

    Returns
    -------
    torch.nn._Loss
    """
    if data == "binary":
        return(nn.BCEWithLogitsLoss())
    return(nn.MSELoss())

def get_metrics(data="binary"):
    if data == "binary":
        return(dict(aucROC=roc_auc_score, aucPR=average_precision_score))
    return(dict(Pearson=pearsonr, Spearman=spearmanr))

def get_optimizer(params, lr=1e-03):
    return(torch.optim.Adam(params, lr=lr))