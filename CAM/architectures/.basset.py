import math
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings("ignore")

class Basset(nn.Module):
    """Basset architecture (Kelley, Snoek & Rinn, 2016)."""

    def __init__(self, sequence_length, n_features=1, output="binary"):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        """
        super(Basset, self).__init__()

        padding = math.floor((200 - sequence_length) / 2.)

        self.blk1 = nn.Sequential(
            nn.Conv1d(4, 100, kernel_size=19, padding=padding),
            nn.BatchNorm1d(100),
            nn.ReLU(inplace=True)
        )
        self.max_pool = nn.MaxPool1d(kernel_size=3, stride=3)
    
        self.blk2 = nn.Sequential(
            nn.Conv1d(100, 200, kernel_size=7),
            nn.BatchNorm1d(200),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=3, stride=3)
        )

        self.blk3 = nn.Sequential(
            nn.Conv1d(200, 200, kernel_size=4),
            nn.BatchNorm1d(200),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=3, stride=3)
        )

        self.fc1 = nn.Sequential(
            nn.Linear(1000, 1000),
            nn.BatchNorm1d(1000, 1e-05, 0.1, True),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3)
        )

        self.fc2 = nn.Sequential(
            nn.Linear(1000, 1000),
            nn.BatchNorm1d(1000, 1e-05, 0.1, True),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3)
        )

        self.fc3 = nn.Sequential(
            nn.Linear(1000, n_features)
        )
        if output == "binary":
            self.fc3.add_module("Sigmoid", nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch."""
        o = self.blk1(x)
        # Save activations from 1st layer
        # (activations, act_index) = torch.max(o, dim=2)
        o = self.max_pool(o)
        o = self.blk2(o)
        o = self.blk3(o)
        o = torch.flatten(o, start_dim=1)
        o = self.fc1(o)
        o = self.fc2(o)

        return(self.fc3(o))

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
    return(torch.optim.SGD(params, lr=lr, momentum=0.9, weight_decay=1e-6))