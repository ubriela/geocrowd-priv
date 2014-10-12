import math
import logging
from collections import deque

import numpy as np

from Grid_pure import Grid_pure
from Params import Params


class Grid_adaptive(Grid_pure):
    def __init__(self, data, param):
        """
        two levels grid
        """

        Grid_pure.__init__(self, data, param)
        self.param.maxHeightHTree = 4

        # compute the best grid size at level 1
        if Params.FIX_GRANULARITY:
            self.m = Params.PARTITION_AG[0]
        else:
            # print self.param.NDATA
            # print self.param.Eps
            # print Params.c
            self.m = int(max(10, int(0.25 * math.ceil((self.param.NDATA * self.param.Eps / param.c) ** (1.0 / 2)))))
        logging.debug("Grid_adaptive: Level 1 size: %d" % self.m)

    def testLeaf(self, curr):
        """test whether a node is a leaf node"""
        if curr.n_depth == self.param.maxHeightHTree:
            return True

        if (curr.n_depth == 2) and (curr.n_data is None or curr.n_data.shape[1] == 0):
            return True

        return False

    def getCountBudget(self):
        count_eps_1 = self.param.Eps * Params.PercentGrid
        count_eps_2 = self.param.Eps * (1 - Params.PercentGrid)
        return [0, 0, count_eps_1, 0, count_eps_2]

    def getCoordinates(self, curr):
        """
        get corrdinates of the point which defines the subnodes
        """
        _box = curr.n_box
        curr_depth = curr.n_depth
        dimP = curr.n_depth % Params.NDIM  # split dimension

        # find the number of partitions
        if curr_depth <= 1:
            _partitions = self.m
        elif curr_depth == 2:
            # compute the best grid size at level 2
            N_prime = max(0, curr.n_count)
            if Params.FIX_GRANULARITY:
                self.m2 = Params.PARTITION_AG[1]
            else:
                self.m2 = int(math.sqrt(N_prime * self.param.Eps * (1 - Params.PercentGrid) / self.param.c2) + 0.5)
                if Params.CUSTOMIZED_GRANULARITY:
                    self.m2 = int(math.sqrt(N_prime * self.param.Eps * (1 - Params.PercentGrid) / Params.c2_c) + 0.5)
            _partitions = curr.secondLevelPartitions = self.m2
            if _partitions <= 1:
                return None, None  # leaf node
        else:  # get grid size stored in parent nodes
            _partitions = curr.secondLevelPartitions

        split_arr = self.getEqualSplit(_partitions, _box[0, dimP], _box[1, dimP])

        # get data points in these partitions
        n_data_arr = [None for _ in range(len(split_arr) - 1)]
        _data = curr.n_data
        if _data is not None and _data.shape[1] >= 1:
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

    def adjustConsistency(self):  # used for adaptive grid only
        # similar to htree variants, adaptive grid do not apply constraint inference on root node
        for x in self.root.children:
            for curr1 in x.children:  # child1 is a first-level cell
                if curr1.n_isLeaf is False:
                    sum = 0
                    for y in curr1.children:
                        for curr2 in y.children:
                            sum += curr2.n_count
                    adjust = (curr1.n_count - sum) / len(curr1.children) ** 2

                    for y in curr1.children:
                        for curr2 in y.children:
                            curr2.n_count += adjust

    def adjustConsistency2(self):
        for x in self.root.children:
            for child1 in x.children:  # child1 is a first-level cell
                # find all child1's children
                m2 = len(child1.children)
                if m2 > 1:
                    alpha = Params.PercentGrid
                    child2s_sum = 0.0
                    for c_child2 in child1.children:
                        for child2 in c_child2.children:  # child2 is a second-level cell
                            child2s_sum += child2.n_count
                    coef = (alpha ** 2 * m2 ** 2) / ((1 - alpha) ** 2 + alpha ** 2 * m2 ** 2)
                    child1.n_count = coef * child1.n_count + (1 - coef) * child2s_sum

                    # update child2s
                    for c_child2 in child1.children:
                        for child2 in c_child2.children:  # child2 is a second-level cell
                            child2.n_count += (child1.n_count - child2s_sum)

    def getCoordinates2(self, curr):
        """ 
        get corrdinates of the point which defines the subnodes
        """
        _box = curr.n_box
        curr_depth = curr.n_depth
        dimP = curr.n_depth % Params.NDIM  # split dimension

        # find the number of partitions
        if curr_depth <= 1:
            _partitions = self.m
        elif curr_depth == 2:
            # compute the best grid size at level 2
            N_prime = max(0, curr.n_count)
            self.m2 = int(math.floor((N_prime * self.param.Eps * (1 - Params.PercentGrid) / Params.c2) ** (1.0 / 2)))
            _partitions = curr.secondLevelPartitions = self.m2
            if _partitions <= 1:
                return None, None  # leaf node
        else:  # get grid size stored in parent nodes
            _partitions = curr.secondLevelPartitions

        split_arr = self.getEqualSplit(_partitions, _box[0, dimP], _box[1, dimP])

        # get data points in these partitions
        _data = curr.n_data
        diff = _box[1, dimP] - _box[0, dimP]
        size = len(split_arr) - 1
        data_arr = [[] for _ in range(size)]
        if _data is not None and _data.shape[1] >= 1:
            for i in range(len(_data[dimP, :])):
                idx = min(size - 1, int((_data[dimP, :][i] - _box[0, dimP]) * size / diff))
                data_arr[idx].append(_data[:, i].tolist())

        n_data_arr = map(lambda data_item: np.array(data_item).T, data_arr)
        for i in range(len(n_data_arr)):
            if n_data_arr[i].size == 0:  # is empty?
                n_data_arr[i] = None
        return split_arr, n_data_arr


    # this function is not in use yet
    def getLeafGranularity(self, query):
        queue = deque()
        queue.append(self.root)
        granularity = 0
        while len(queue) > 0:
            curr = queue.popleft()
            _box = curr.n_box
            if curr.n_isLeaf is not True:  # if not leaf
                bool_matrix = np.zeros((2, query.shape[1]))
                bool_matrix[0, :] = _box[0, :] <= query[0, :]
                bool_matrix[1, :] = _box[1, :] >= query[1, :]

                if np.all(bool_matrix):  # if node range contains query range
                    if curr.n_depth == 2:
                        granularity = len(curr.children)

                    for node in curr.children:
                        queue.append(node)

        return granularity