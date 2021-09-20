from collections import OrderedDict
import math
import numpy as np
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

class _Model(nn.Module):

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

class Chinook(_Model):
    """Chinook: a glass-box deep learning model for genomics."""

    def __init__(self, cnn_units, kernel_size, sequence_length, n_features=1,
        clamp_weights=False, no_padding=False, weights_file=None):
        """
        Parameters
        ----------
        cnn_units : int
            Total number of individual CNN units
        kernel_size : int
            Convolutional kernel size
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        weights_file : pass
            ...
        """
        super(Chinook, self).__init__()

        self._options = {
            "cnn_units": cnn_units,
            "kernel_size": kernel_size,
            "sequence_length": sequence_length,
            "n_features": n_features,
            "clamp_weights": clamp_weights,
            "no_padding": no_padding,
            "weights_file": weights_file,
        }

        if no_padding:
            self.__n_channels = math.floor((sequence_length-kernel_size+1)/7.)
            self.__padding = 0
        else:
            self.__n_channels = math.floor((sequence_length+kernel_size+1)/7.)
            self.__padding = kernel_size

        self.linears = nn.Sequential(
            nn.Conv1d(
                in_channels=4*cnn_units,
                out_channels=1*cnn_units,
                kernel_size=kernel_size,
                padding=self.__padding,
                groups=cnn_units,
            ),
            nn.BatchNorm1d(cnn_units),
            ExpAct(),
            nn.MaxPool1d(7, 7),
            nn.Flatten(), 
            UnSqueeze(),
            nn.Conv1d(
                in_channels=self.__n_channels*cnn_units,
                out_channels=100*cnn_units,
                kernel_size=1,
                groups=cnn_units,
            ),
            nn.BatchNorm1d(100*cnn_units, 1e-05, 0.1, True),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Conv1d(
                in_channels=100*cnn_units,
                out_channels=1*cnn_units,
                kernel_size=1,
                groups=cnn_units
            ),
            nn.BatchNorm1d(1*cnn_units, 1e-05, 0.1, True),
            nn.ReLU(),
            nn.Flatten(),
        )

        self.final = nn.Linear(cnn_units, n_features)

        if weights_file is not None:
            self.load_weights(weights_file)

    def forward(self, x):
        """Forward propagation of a batch."""
        x = x.repeat(1, self._options["cnn_units"], 1)
        outs = self.linears(x)
        outs = self.final(outs)

        return(outs)

# DanQ Pytorch implementation 
# From: https://github.com/PuYuQian/PyDanQ/blob/master/DanQ_train.py
class DanQ(_Model):
    """DanQ architecture (Quang & Xie, 2016)."""

    def __init__(self, sequence_length, n_features=1, weights_file=None):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        weights_file : pass
            ...
        """
        super(DanQ, self).__init__()

        self._options = {
            "sequence_length": sequence_length,
            "n_features": n_features,
            "weights_file": weights_file,
        }

        self.__n_channels = math.floor((sequence_length-25)/13.)
    
        self.Conv1 = nn.Conv1d(
            in_channels=4,
            out_channels=320,
            kernel_size=26,
            padding=0,
        )
        self.Maxpool = nn.MaxPool1d(kernel_size=13, stride=13)
        self.Drop1 = nn.Dropout(p=0.2)
        self.BiLSTM = nn.LSTM(
            input_size=320,
            hidden_size=320,
            num_layers=2,
            batch_first=True,
            dropout=0.5,
            bidirectional=True
        )
        self.Linear1 = nn.Linear(self.__n_channels*640, 925)
        self.Linear2 = nn.Linear(925, n_features)

        if weights_file is not None:
            self.load_weights(weights_file)

    def forward(self, x):
        """Forward propagation of a batch."""
        x = self.Conv1(x)
        x = nn.functional.relu(x)
        x = self.Maxpool(x)
        x = self.Drop1(x)
        x_x = torch.transpose(x, 1, 2)
        x, _ = self.BiLSTM(x_x)
        x = x.contiguous().view(-1, self.__n_channels*640)
        x = self.Linear1(x)
        x = nn.functional.relu(x)
        out = self.Linear2(x)

        return(out)

# def _flip(x, dim):
#     """
#     Adapted from Selene:
#     https://github.com/FunctionLab/selene/blob/master/selene_sdk/utils/non_strand_specific_module.py

#     Reverses the elements in a given dimension `dim` of the Tensor.
#     source: https://github.com/pytorch/pytorch/issues/229
#     """
#     xsize = x.size()
#     dim = x.dim() + dim if dim < 0 else dim
#     x = x.contiguous()
#     x = x.view(-1, *xsize[dim:])
#     x = x.view(
#         x.size(0), x.size(1), -1)[:, getattr(
#             torch.arange(x.size(1)-1, -1, -1),
#             ("cpu","cuda")[x.is_cuda])().long(), :]

#     return x.view(xsize)

# class NonStrandSpecific(nn.Module):
#     """
#     Adapted from Selene:
#     https://github.com/FunctionLab/selene/blob/master/selene_sdk/utils/non_strand_specific_module.py

#     A torch.nn.Module that wraps a user-specified model architecture if the
#     architecture does not need to account for sequence strand-specificity.

#     Parameters
#     ----------
#     model : torch.nn.Module
#         The user-specified model architecture.
#     mode : {'mean', 'max'}, optional
#         Default is 'mean'. NonStrandSpecific will pass the input and the
#         reverse-complement of the input into `model`. The mode specifies
#         whether we should output the mean or max of the predictions as
#         the non-strand specific prediction.

#     Attributes
#     ----------
#     model : torch.nn.Module
#         The user-specified model architecture.
#     mode : {'mean', 'max'}
#         How to handle outputting a non-strand specific prediction.
#     """

#     def __init__(self, model):
#         super(NonStrandSpecific, self).__init__()

#         self.model = model

#     def forward(self, input):
#         reverse_input = None
#         reverse_input = _flip(_flip(input, 1), 2)

#         output = self.model.forward(input)
#         output_from_rev = self.model.forward(reverse_input)

#         return((output + output_from_rev) / 2)

def get_loss(input_data="binary"):
    """
    Specify the appropriate loss function (criterion) for this model.

    Returns
    -------
    torch.nn._Loss
    """
    if input_data == "binary":
        return(nn.BCEWithLogitsLoss())
    return(nn.MSELoss())

def get_metrics(input_data="binary"):
    if input_data == "binary":
        return(dict(aucROC=roc_auc_score, aucPR=average_precision_score))
    return(dict(Pearson=pearsonr, Spearman=spearmanr))

def get_optimizer(params, lr=1e-03):
    return(torch.optim.Adam(params, lr=lr))