import math
from collections import deque

import numpy as np

from HT_hybrid import HT_hybrid
from Params import Params


class HT_hybrid_adaptive(HT_hybrid):
    def __init__(self, data, param):
        HT_hybrid.__init__(self, data, param)

    def getCoordinates(self, curr):
        budget_s = self.getSplitBudget(curr)
        _data = curr.n_data
        _box = curr.n_box
        curr_depth = curr.n_depth
        dimP = curr_depth % Params.NDIM  # split dimension

        isSorted = False
        if curr_depth <= 1:  # data dependent
            _idx = np.argsort(_data[dimP, :], kind='mergesort')
            _data[:, :] = _data[:, _idx]  # sorted by dimP dimension
            isSorted = True

            sw = self.switchSplit[dimP]  # switching level
            split_height = sw - 1
            swPartitions = self.param.switchPartitionsHTree[dimP]  # slicing

            split_arr = self.recursiveSlices(split_height, swPartitions, _data[dimP, :], budget_s)
            split_arr.insert(0, _box[0, dimP])
            split_arr.append(_box[1, dimP])
        elif curr_depth == 2:  # data independent (grid)
            # compute the best grid size at level 2
            N_prime = max(0, curr.n_count)
            m2 = int(math.floor((N_prime * self.gridBudget / Params.c2) ** (1.0 / 2)))
            curr.secondLevelPartitions = m2  # save for furture use
            split_arr = self.getEqualSplit(m2, _box[0, dimP], _box[1, dimP])
        else:  # get grid size stored in parent nodes
            split_arr = self.getEqualSplit(curr.secondLevelPartitions, _box[0, dimP], _box[1, dimP])

        # get data points in these partitions
        n_data_arr = [None for _ in range(len(split_arr) - 1)]
        if _data is not None and _data.shape[1] >= 1:
            if not isSorted:
                _idx = np.argsort(_data[dimP, :], kind='mergesort')
                _data[:, :] = _data[:, _idx]  # sorted by dimP dimension
            for i in range(len(split_arr) - 1):
                posP1 = np.searchsorted(_data[dimP, :], split_arr[i])
                posP2 = np.searchsorted(_data[dimP, :], split_arr[i + 1])
                if i == 0:
                    n_data = _data[:, :posP2]
                elif i == len(split_arr) - 2:
                    n_data = _data[:, posP1:]
                else:
                    n_data = _data[:, posP1:posP2]
                n_data_arr[i] = n_data

        return split_arr, n_data_arr


    def getRecursiveMediansNo(self, array, budget):
        queue = deque()
        queue.append(array)
        split = 0
        while len(queue) > 0:
            curr = queue.popleft()
            idP, splitP = self.getSplit(curr, curr[0], curr[-1], budget)
            r_k = np.searchsorted(curr, res[k])
            idP_true = float(len(curr)) / 2
            relErr = abs(idP_true - r_k) / idP_true
            cumulative_error += relErr
            if cumulative_error > Params.maxRelativeError:
                return split

            queue.append(curr[:idP + 1])
            queue.addpend(curr[idP:])
            split += 2