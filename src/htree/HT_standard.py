import sys
import logging

import numpy as np

from HT_true import HT_true
from Params import Params


class HT_standard(HT_true):
    """ standard private h-tree is a 2 level data dependent tree """

    def __init__(self, data, param):
        HT_true.__init__(self, data, param)
        self.param.maxHeightHTree = 2

    # overrided in HT_composite and HT_hybrid
    def getCountBudget(self):
        count_eps = self.param.Eps * (1 - self.param.PercentSplit)
        if Params.useLeafOnlyHTree:
            return [0, 0, count_eps]
        if self.param.geoBudget == 'optimal':
            # n1 = 2**self.maxSplit[0]
            #            n2 = n1*2**self.maxSplit[1]
            #            unit = count_eps / (n1**(1.0/3) + n2**(1.0/3))
            #            eps1 = unit*n1**(1.0/3)
            #            eps2 = unit*n2**(1.0/3)
            #            return [0, eps1, eps2]
            return [0, count_eps / 2, count_eps / 2]
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)

    # overrided in HT_composite and HT_hybrid
    def getSplitBudget(self, curr):
        dimP = curr.n_depth % Params.NDIM  # split dimension
        if curr.n_depth == 0:
            split_eps = self.param.Eps * self.param.PercentSplit * self.param.splitFrac
        else:
            split_eps = self.param.Eps * self.param.PercentSplit * (1 - self.param.splitFrac)
        split_no = self.maxSplit[dimP]
        ret = [split_eps / split_no for _ in range(split_no)]
        return ret

    def getNoisySlice(self, array, left, right, epsilon, partitions, slice):
        if self.param.splitScheme == 'expo':
            med = self.differ.getSplit_slice(array, left, right, epsilon, partitions, slice)
            pos = np.searchsorted(array, med)
            return pos, med
        else:
            logging.error('No such split scheme')
            sys.exit(1)

    def getNoisyMedian(self, array, left, right, epsilon):
        if self.param.splitScheme == 'expo':
            med = self.differ.getSplit_exp(array, left, right, epsilon)
            pos = np.searchsorted(array, med)
            return pos, med
        else:
            logging.error('No such split scheme')
            sys.exit(1)