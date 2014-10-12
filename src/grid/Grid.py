from collections import deque

import numpy as np

from Generic import Generic
from Params import Params


class Grid(Generic):
    """
    Generic grid
    """

    def __init__(self, data, param):
        Generic.__init__(self, data, param)


    def testLeaf(self, curr):
        """test whether a node is a leaf node"""
        if curr.n_depth == self.param.maxHeightHTree:
            return True
        return False


        # apply canonical range

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

                    if np.all(bool_matrix) and node.n_depth == 2:  # if query range contains node range
                        count += node.n_count
                    if self.intersect(_box, query):
                        queue.append(node)
        return float(count)