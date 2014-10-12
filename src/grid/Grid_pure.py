import logging

import numpy as np

from Grid import Grid
from Params import Params


class Grid_pure(Grid):
    """ non-private grid, use true median for split and exact node count """

    def __init__(self, data, param):
        Grid.__init__(self, data, param)

        self.param.m = self.param.partitionsHTree[1]
        self.param.maxHeightHTree = 2

        logging.debug("Grid_pure: size: %d" % self.param.m)

    def getCountBudget(self):
        """ zero count budget mean exact node count """
        return [0 for _ in range(self.param.maxHeightHTree + 1)]


    def getCoordinates(self, curr):
        """ 
        get corrdinates of the point which defines the subnodes: 
        return split_arr: split points
	n_data_arr: data in each partitions
        """
        _box = curr.n_box
        dimP = curr.n_depth % Params.NDIM  # split dimension

        split_arr = self.getEqualSplit(self.param.m, _box[0, dimP], _box[1, dimP])

        # get data points in these partitions
        n_data_arr = [None for _ in range(self.param.m)]
        _data = curr.n_data

        if _data is not None and _data.shape[1] >= 1:
            _idx = np.argsort(_data[dimP, :], kind='mergesort')
            _data[:, :] = _data[:, _idx]  # sorted by dimP dimension

            for i in range(self.param.m):
                posP1 = np.searchsorted(_data[dimP, :], split_arr[i])
                posP2 = np.searchsorted(_data[dimP, :], split_arr[i + 1])
                if i == 0:  # the first partition
                    n_data = _data[:, :posP2]
                elif i == len(split_arr) - 2:  # the last partition
                    n_data = _data[:, posP1:]
                else:
                    n_data = _data[:, posP1:posP2]
                n_data_arr[i] = n_data

        return split_arr, n_data_arr
