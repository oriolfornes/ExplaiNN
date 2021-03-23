from ignite.metrics import Metric
import torch

class PearsonR(Metric):
    """https://github.com/pytorch/pytorch/issues/1254"""

    def __init__(self, output_transform=lambda x: x):
        self._x = None
        self._y = None
        super(PearsonR, self).__init__(output_transform=output_transform)

    def reset(self):
        self._x = []
        self._y = []

    def update(self, output):
        y_pred, y = output[0].detach(), output[1].detach()
        y_pred = y_pred.clone().to(self._device)
        y = y.clone().to(self._device)
        self._x.append(y_pred)
        self._y.append(y)

    def compute(self):
        x = torch.flatten(torch.cat(self._x, dim=0))
        y = torch.flatten(torch.cat(self._y, dim=0))
        mean_x = torch.mean(x)
        mean_y = torch.mean(y)
        xm = x.sub(mean_x)
        ym = y.sub(mean_y)
        r_num = xm.dot(ym)
        r_den = torch.norm(xm, 2) * torch.norm(ym, 2)
        return(r_num / r_den)