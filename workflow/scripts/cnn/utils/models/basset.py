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

        # Block 1 :
        self.c1 = nn.Conv1d(4, 100, kernel_size=19, padding=padding)
        self.bn1 = nn.BatchNorm1d(100)
        self.rl1 = nn.modules.activation.ReLU()
        self.mp1 = nn.MaxPool1d(kernel_size=3, stride=3)

        # Block 2 :
        self.c2 = nn.Conv1d(100, 200, kernel_size=7)
        self.bn2 = nn.BatchNorm1d(200)
        self.rl2 = nn.modules.activation.ReLU()
        self.mp2 = nn.MaxPool1d(kernel_size=3, stride=3)

        # Block 3 :
        self.c3 = nn.Conv1d(200, 200, kernel_size=4)
        self.bn3 = nn.BatchNorm1d(200)
        self.rl3 = nn.modules.activation.ReLU()
        self.mp3 = nn.MaxPool1d(kernel_size=3, stride=3)

        # Block 4 : Fully Connected 1 :
        self.d4 = nn.Linear(1000, 1000) #1000 for 200 input size
        self.bn4 = nn.BatchNorm1d(1000, 1e-05, 0.1, True)
        self.rl4 = nn.modules.activation.ReLU() 
        self.dr4 =  nn.Dropout(0.3) 

        # Block 5 : Fully Connected 2 :
        self.d5 = nn.Linear(1000, 1000)
        self.bn5 = nn.BatchNorm1d(1000, 1e-05, 0.1, True)
        self.rl5 = nn.modules.activation.ReLU() 
        self.dr5 =  nn.Dropout(0.3) 

        # Block 6 : Fully connected 3
        self.d6 = nn.Linear(1000, n_features)

        # Output
        self.sigmoid = nn.modules.activation.Sigmoid()
        self._output = output        

    def forward(self, x):
        """Forward propagation of a batch."""
        # Block 1
        # x is of size - batch, 4, 200
        x = self.rl1(self.bn1(self.c1(x))) # output - batch, 100, 182
        # we save the activations of the first layer (interpretation)
        activations = x # batch, 100, 182
        x = self.mp1(x) # output - batch, 100, 60
        # Block 2
        # input is of size batch, 100, 60
        x = self.mp2(self.rl2(self.bn2(self.c2(x)))) #output - batch, 200, 18
        # Block 3
        # input is of size batch, 200, 18
        em = self.mp3(self.rl3(self.bn3(self.c3(x)))) #output - batch, 200, 5
        # Flatten
        o = torch.flatten(em, start_dim=1) #output - batch, 1000
        # FC1
        #input is of size - batch, 1000
        o = self.dr4(self.rl4(self.bn4(self.d4(o)))) #output - batch, 1000
        # FC2
        #input is of size - batch, 1000
        o = self.dr5(self.rl5(self.bn5(self.d5(o)))) #output - batch, 1000
        # FC3
        #input is of size - batch, 1000
        #o = self.sig(self.d6(o)) #output - batch, num_of_classes
        o = self.d6(o) #doing BCEWithLogits #output - batch, num_of_classes
        #maximum for every filter (and corresponding index)
        activations, act_index = torch.max(activations, dim=2)
        if self._output == "binary":
            return(self.sigmoid(o))
        return(o)

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