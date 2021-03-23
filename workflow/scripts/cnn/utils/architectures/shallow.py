import math
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings("ignore")

class Model(nn.Module):
    """
    DeepSELEX architecture (Asif & Orenstein, 2020):
    We developed DeepSELEX to infer protein–DNA binding preferences from
    HT-SELEX data. DeepSELEX is based on a convolutional neural network
    comprised of a 1D-convolutional layer with 512 convolutional kernels
    (width 8, stride 1), a max pooling layer (pool-size 5, stride 5), a
    cascade of 3 fully connected layers with 64, 32, 32 nodes (ReLU
    activation), respectively, and an output layer with 5 nodes (sigmoid
    activation) (Fig. 2B). The input is a one-hot encoding of a DNA sequence
    in the length of a HT-SELEX oligonucleotide (14, 20, 30 or 40, depending
    on the experiment) plus 4 flanking nucleotides on each side. The labels
    are the corresponding enrichment cycles, i.e the network learns the
    enrichment cycle each DNA sequence belongs to. Model parameters were
    learned by Adam optimizer (learning rate α=10−3, β1=0.9, β2=0.999,
    decay=10−5⁠, batch size 64) using categorical cross-entropy loss.
    """

    def __init__(self, sequence_length, n_features=1, output="binary"):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        """
        super(Model, self).__init__()

        self.nnet = nn.Sequential(
            nn.Conv1d(4, 512, kernel_size=8),
            nn.BatchNorm1d(512),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=5, stride=5)
        )

        self._n_channels = math.floor((sequence_length - 7) / 5)

        self.classifier = nn.Sequential(
            nn.Linear(self._n_channels * 512, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 32),
            nn.ReLU(inplace=True),
            nn.Linear(32, 32),
            nn.ReLU(inplace=True),
            nn.Linear(32, n_features),
        )

        if output == "binary":
            self.classifier.add_module("Sigmoid", nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch."""
        out = self.nnet(x)
        out = torch.flatten(out, start_dim=1)
        out = self.classifier(out)

        return(out)

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