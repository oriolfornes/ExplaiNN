import math
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings("ignore")

# From: https://github.com/FunctionLab/selene/blob/master/models/danQ.py
class DanQ(nn.Module):
    """DanQ architecture (Quang & Xie, 2016)."""

    def __init__(self, sequence_length, n_features=1, output="binary"):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        """
        super(DanQ, self).__init__()

        self.nnet = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=26),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=13, stride=13),
            nn.Dropout(0.2)
        )

        self.bdlstm = nn.Sequential(
            nn.LSTM(320, 320, batch_first=True, bidirectional=True)
        )

        self._n_channels = math.floor((sequence_length - 25) / 13)
        
        self.classifier = nn.Sequential(
            nn.Dropout(0.5),
            nn.Linear(self._n_channels * 640, 925),
            nn.ReLU(inplace=True),
            nn.Linear(925, n_features)
        )
        # if output == "binary":
        #     self.classifier.add_module("Sigmoid", nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch."""
        o = self.nnet(x)
        reshaped = o.transpose(0, 1).transpose(0, 2)
        o, _ = self.bdlstm(reshaped)
        o = o.transpose(0, 1)
        reshaped = o.contiguous().view(o.size(0), 640 * self._n_channels)

        return(self.classifier(reshaped))

def get_criterion(output="binary"):
    """
    Specify the appropriate loss function (criterion) for this model.

    Returns
    -------
    torch.nn._Loss
    """
    if output == "binary":
        return(nn.BCELoss())
    return(nn.MSELoss())

def get_optimizer(params, lr):
    return(torch.optim.Adam(params, lr=lr))
    # return(torch.optim.SGD(params, lr=lr, momentum=0.9, weight_decay=1e-6))