import math
import torch
import torch.nn as nn

# Adapted from:
# https://github.com/FunctionLab/selene/blob/master/models/danQ.py
class TBiNet(nn.Module):
    """TBiNet architecture (https://github.com/dmis-lab/tbinet)."""

    def __init__(self, sequence_length, n_features=1, output="binary"):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_features : int
            Total number of features to predict
        """
        super(TBiNet, self).__init__()

        # Convolutional Layer
        self.convolutional_layer = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=26),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=13, stride=13),
            nn.Dropout(0.2)
        )
        # # Convolutional Layer
        # output = Conv1D(320,kernel_size=26,padding="valid",activation="relu")(sequence_input)
        # output = MaxPooling1D(pool_size=13, strides=13)(output)
        # output = Dropout(0.2)(output)

        # Attention layer
        self.attention_layer = nn.Sequential(
            nn.Dense(1)
        )

        # Bidirectional LSTM layer
        self.bidirectional_lstm_layer = nn.Sequential(
            nn.LSTM(320, 320, batch_first=True, bidirectional=True),
            nn.Dropout(0.5)
        )
        # output = Bidirectional(LSTM(320,return_sequences=True))(output)
        # output = Dropout(0.5)(output)

        self._n_channels = max([math.floor((sequence_length - 25) / 13), 1])

        self.classifier = nn.Sequential(
            ,
            nn.Linear(self._n_channels * 640, 925),
            nn.ReLU(inplace=True),
            nn.Linear(925, n_features)
        )
        if output == "binary":
            self.classifier.add_module("Sigmoid", nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch."""
        out = self.nnet(x)
        print(out.shape)
        reshape_out = out.transpose(0, 1).transpose(0, 2)
        out, _ = self.bdlstm(reshape_out)
        out = out.transpose(0, 1)
        print(out.shape)
        reshape_out = out.contiguous().view(
            out.size(0), 640 * self._n_channels)
        print(reshape_out.shape)
        exit(0)
        predict = self.classifier(reshape_out)

        return(predict)

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

def get_optimizer(params, lr=0.0001):
    return(torch.optim.Adam(params, lr=lr))