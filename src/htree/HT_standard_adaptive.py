import math
from collections import deque

import numpy as np

from HT_standard import HT_standard
from Params import Params


class HT_standard_adaptive(HT_standard):
    def __init__(self, data, param):
        HT_standard.__init__(self, data, param)

        self.split_eps = self.param.Eps * self.param.PercentSplit / 2
        budget_x = self.split_eps / self.maxSplit[1]
        budget_y = self.split_eps / self.maxSplit[1]

        _data = data
        _idx = np.argsort(_data[0, :], kind='mergesort')
        _data[:, :] = _data[:, _idx]  # sorted by x
        split_no_x = self.getRecursiveMediansNo(_data[0, :], budget_x)

        _idx = np.argsort(_data[1, :], kind='mergesort')
        _data[:, :] = _data[:, _idx]  # sorted by x
        split_no_y = self.getRecursiveMediansNo(_data[1, :], budget_y)
        self.split_no = [split_no_x, split_no_y]
        print self.split_no

    def getSplitBudget(self, curr):
        dimP = curr.n_depth % Params.NDIM  # split dimension
        ret = [self.split_eps / self.split_no[dimP] for _ in range(self.split_no[dimP])]
        return ret

    def getCoordinates(self, curr):
        budget_s = self.getSplitBudget(curr)
        _data = curr.n_data
        _box = curr.n_box
        dimP = curr.n_depth % Params.NDIM  # split dimension

        _idx = np.argsort(_data[dimP, :], kind='mergesort')
        _data[:, :] = _data[:, _idx]  # sorted by dimP dimension

        split_arr = self.recursiveSlices(self.split_no[dimP] - 1, 2 ** self.split_no[dimP], _data[dimP, :], budget_s)

        split_arr.insert(0, _box[0, dimP])
        split_arr.append(_box[1, dimP])

        # get data points in these partitions
        n_data_arr = []
        for i in range(len(split_arr) - 1):
            posP1 = np.searchsorted(_data[dimP, :], split_arr[i])
            posP2 = np.searchsorted(_data[dimP, :], split_arr[i + 1])
            if i == 0:
                n_data = _data[:, :posP2]
            elif i == len(split_arr) - 2:
                n_data = _data[:, posP1:]
            else:
                n_data = _data[:, posP1:posP2]
            n_data_arr.append(n_data)

        return split_arr, n_data_arr

    def getRecursiveMediansNo(self, array, budget):
        queue = deque()
        queue.append(array)
        split = 0
        while len(queue) > 0:
            curr = queue.popleft()

            totalErr = 0
            testNo = 33
            for i in range(testNo):
                idP, splitP = self.getNoisyMedian(curr, curr[0], curr[-1], budget)
                idP_true = float(len(curr)) / 2
                relErr = abs(idP_true - idP) / idP_true
                totalErr += relErr
            avgErr = totalErr / testNo
            if avgErr > Params.maxRelativeError:
                return int(math.ceil(math.log(split, 2)))

            queue.append(curr[:idP + 1])
            queue.append(curr[idP:])
            split += 2