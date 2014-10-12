import math

from Grid_pure import Grid_pure
from Params import Params


class Grid_uniform(Grid_pure):
    def __init__(self, data, param):
        Grid_pure.__init__(self, data, param)
        # compute the best grid size
        self.param.m = int(math.ceil((Params.NDATA * self.param.Eps / Params.c) ** (1.0 / 2)))
        self.param.maxHeightHTree = 2

    def getCountBudget(self):
        """ only leaf level get full budget """
        return [0, 0, self.param.Eps]