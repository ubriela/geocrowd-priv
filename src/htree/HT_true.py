from collections import deque
import sys
import math
import logging

import numpy as np

from HT_pure import HT_pure
from Params import Params


class HT_true(HT_pure):
    """ non-private ht-tree, use true medians to split but noisy counts """

    def __init__(self, data, param):
        HT_pure.__init__(self, data, param)
        self.param.maxHeightHTree = 2

        m = int(math.ceil(math.log(self.param.partitionsHTree[0], 2)))
        n = int(math.ceil(math.log(self.param.partitionsHTree[1], 2)))
        self.maxSplit = [m, n]
        m2 = int(math.ceil(math.log(self.param.switchPartitionsHTree[0], 2)))
        n2 = int(math.ceil(math.log(self.param.switchPartitionsHTree[1], 2)))
        self.switchSplit = [m2, n2]

    def testLeaf(self, curr):
        """test whether a node is a leaf node"""
        if (curr.n_depth == self.param.maxHeightHTree) or \
                (curr.n_data is None or curr.n_data.shape[1] == 0) or \
                (curr.n_count < self.param.minPartSizeHTree):
            return True
        return False

    # overrided in HT_standard
    def getCountBudget(self):
        count_eps = self.param.Eps
        if Params.useLeafOnlyHTree:
            return [0, 0, count_eps]
        H = self.param.maxHeightHTree
        if self.param.geoBudget == 'none':
            ret = [count_eps / (H + 1) for _ in range(H + 1)]
            ret.insert(0, 0)
            return ret
        elif self.param.geoBudget == 'optimal':
            n1 = 2 ** self.maxSplit[0]
            n2 = n1 * 2 ** self.maxSplit[1]
            unit = count_eps / (n1 ** (1.0 / 3) + n2 ** (1.0 / 3))
            eps1 = unit * n1 ** (1.0 / 3)
            eps2 = unit * n2 ** (1.0 / 3)
            return [0, eps1, eps2]
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)

    def rangeCount(self, query):
        """
        Query answering function. Find the number of data points within a query rectangle.
        This function assume that the tree is contructed with noisy count for every node
        """

        queue = deque()
        queue.append(self.root)
        count = 0.0
        while len(queue) > 0:
            curr = queue.popleft()
            _box = curr.n_box
            if curr.n_isLeaf is True:
                frac = 1
                if self.intersect(_box, query):
                    for i in range(_box.shape[1]):
                        if _box[1, i] == _box[0, i] or Params.WorstCase == True:
                            frac *= 1
                        else:
                            frac *= (min(query[1, i], _box[1, i]) - max(query[0, i], _box[0, i])) / (
                                _box[1, i] - _box[0, i])
                    count += curr.n_count * frac
            else:  # if not leaf
                for node in curr.children:
                    bool_matrix = np.zeros((2, query.shape[1]))
                    bool_matrix[0, :] = query[0, :] <= _box[0, :]
                    bool_matrix[1, :] = query[1, :] >= _box[1, :]

                    if (not Params.useLeafOnlyHTree) and np.all(bool_matrix):  # if query range contains node range
                        count += node.n_count
                    if self.intersect(_box, query):
                        queue.append(node)
        return float(count)