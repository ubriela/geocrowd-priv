import logging
from collections import deque
import math

import numpy as np

from Grid_adaptive import Grid_adaptive
from Params import Params
from Utils import distance


class Grid_adaptive_localness(Grid_adaptive):
    def __init__(self, data, param):
        Grid_adaptive.__init__(self, data, param)

        # compute the best grid size at level 1

        # grid size
        d1 = distance(Params.LOW[0], Params.LOW[1], Params.LOW[0], Params.HIGH[1])
        d2 = distance(Params.LOW[0], Params.LOW[1], Params.HIGH[0], Params.LOW[1])
        mind = min(d1, d2)
        self.m = int(math.floor(mind / (2 * Params.MTD)))
        logging.debug("Grid_adaptive_localness: Level 1 size: %d" % self.m)

    def getCountBudget(self):
        """ only leaf level get full budget """
        count_eps_1 = self.param.Eps * Params.PercentGridLocalness
        count_eps_2 = self.param.Eps * (1 - Params.PercentGridLocalness)
        return [0, 0, count_eps_1, 0, count_eps_2]


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
                if self.rect_intersect(_box, query):
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

                    if rect_intersectrsect(_box, query):
                        queue.append(node)
        return float(count)