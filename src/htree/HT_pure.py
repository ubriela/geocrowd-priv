import math

import numpy as np

from HTree import HTree
from Params import Params


class HT_pure(HTree):
    """ non-private h-tree, use true median for split and exact node count """

    def __init__(self, data, param):
        HTree.__init__(self, data, param)
        self.param.maxHeightHTree = 2

        if param.dynamicGranularity:
            count_eps = self.param.Eps * (1 - self.param.PercentSplit)
            self.param.partitionsHTree[0] = self.param.partitionsHTree[1] = int(
                math.ceil((Params.NDATA * count_eps / Params.c_htree) ** (1.0 / 2)))
            self.param.switchPartitionsHTree[0] = self.param.switchPartitionsHTree[1] = 2 ** (
                int(math.ceil(math.log(self.param.partitionsHTree[0], 2) / 2.0)))
            print self.param.partitionsHTree[0], self.param.switchPartitionsHTree[0]

        self.maxSplit = [8, 8]  # this is overrided in HT_true
        self.switchSplit = [4, 4]  # this is overrided in HT_true

    def testLeaf(self, curr):
        """test whether a node is a leaf node"""
        if (curr.n_depth == self.param.maxHeightHTree) or \
                (curr.n_data is None or curr.n_data.shape[1] == 0):
            return True
        return False

    def getSplitBudget(self, curr):
        """ zero split budget means true median """
        dimP = curr.n_depth % Params.NDIM  # split dimension
        return [0 for _ in range(self.maxSplit[dimP])]

    # overrided in HT_true
    def getCountBudget(self):
        """ zero count budget mean exact node count """
        return [0 for _ in range(self.param.maxHeightHTree + 1)]

    def getSplit(self, array, left, right, epsilon=0, partitions=0, slice=0):
        """ If epsilon is zero, return true split points """
        n = len(array)

        if epsilon < 10 ** (-6):
            if n % partitions == 1:
                return slice * n / partitions, array[slice * n / partitions]
            else:
                return slice * n / partitions, (array[slice * n / partitions] + array[slice * n / partitions - 1]) / 2.0
        else:  # use recursive slicing
            return self.getNoisySlice(array, left, right, epsilon, partitions, slice)

    def getCoordinates(self, curr):
        """ 
        get corrdinates of the point which defines the four subnodes: 
        split_arr, n_data_arr
        """
        budget_s = self.getSplitBudget(curr)
        _data = curr.n_data
        _box = curr.n_box
        dimP = curr.n_depth % Params.NDIM  # split dimension

        _idx = np.argsort(_data[dimP, :], kind='mergesort')
        _data[:, :] = _data[:, _idx]  # sorted by dimP dimension

        split_arr = self.recursiveSlices(len(budget_s) - 1, self.param.partitionsHTree[dimP], _data[dimP, :], budget_s)

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


    def recursiveSlices(self, height, partitions, array, budget_s, leftPriority=True):
        """recursively select slicing positions.
        partitions is the number of slices
        leftPriority is used to make the tree more balanced"""
        if len(array) < self.param.minPartSizeHTree or partitions == 1:
            return []

        if leftPriority:
            left = (partitions + 1) / 2
        else:
            left = partitions / 2
        right = partitions - left

        if height <= 0:
            idP, splitP = self.getSplit(array, array[0], array[-1], budget_s[height], partitions, left)
            if splitP <= array[0] or splitP >= array[-1]:
                splitP = float(self.differ.getUniform(array[0], array[-1]))

            return [splitP]

        else:
            idP, splitP = self.getSplit(array, array[0], array[-1], budget_s[height], partitions, left)

            # if splitP is not in range --> get a random value in the range
            if splitP <= array[0] or splitP >= array[-1]:
                splitP = float(self.differ.getUniform(array[0], array[-1]))
                idP = np.searchsorted(array, splitP)

            left_arrP = self.recursiveSlices(height - 1, left, array[:idP + 1], budget_s, not leftPriority)
            right_arrP = self.recursiveSlices(height - 1, right, array[idP:], budget_s, not leftPriority)
            return left_arrP + [splitP] + right_arrP


    def recursiveMedians(self, height, array, budget_s):
        """recursively select medians"""
        if len(array) < self.param.minPartSizeHTree:
            return []
        if height <= 0:
            idP, splitP = self.getSplit(array, array[0], array[-1], budget_s[height])
            if splitP <= array[0] or splitP >= array[-1]:
                splitP = float(self.differ.getUniform(array[0], array[-1]))

            return [splitP]

        else:  # height >=1
            idP, splitP = self.getSplit(array, array[0], array[-1], budget_s[height])

            # if splitP is not in range --> get a random value in the range
            if splitP <= array[0] or splitP >= array[-1]:
                splitP = float(self.differ.getUniform(array[0], array[-1]))
                idP = np.searchsorted(array, splitP)

            left_arrP = self.recursiveMedians(height - 1, array[:idP + 1], budget_s)
            right_arrP = self.recursiveMedians(height - 1, array[idP:], budget_s)
            return left_arrP + [splitP] + right_arrP

            # def __init__(self):
        # self.differ = Differential(1412)
        # return None
        #
        #    def getNoisyMedian(self, array, left, right, epsilon):
        #        med = self.differ.getSplit_exp(array, left, right, epsilon)
        #        pos = np.searchsorted(array, med)
        #        return pos, med
        #
        #    def getNoisySlice(self,array,left,right,epsilon,partitions,slice):
        #        med = self.differ.getSplit_slice(array, left, right, epsilon, partitions, slice)
        #        pos = np.searchsorted(array, med)
        #        return pos, med
        #
        #
        #    def getRecursiveMediansNo(self,array,budget):
        #	queue = deque()
        #        queue.append(array)
        #	split = 0;
        #	cumulative_error = 0
        #        while len(queue) > 0:
        #	    curr = queue.popleft()
        #            idP, splitP = self.getNoisyMedian(curr,curr[0],curr[-1],budget)
        #            idP_true = float(len(curr))/2
        #            relErr = abs(idP_true-idP)/idP_true
        #            cumulative_error = cumulative_error + relErr
        #            if cumulative_error > Params.maxRelativeError:
        #                return int(math.ceil(math.log(split, 2)));
        #
        #	    queue.append(curr[:idP+1])
        #	    queue.append(curr[idP:])
        #	    split = split + 2

        return 0

# debugging
if __name__ == '__main__':
    # ht_pure = HT_pure()

    # test recursiveMedians
    # height = 2
    # array = range(1000)
    # budget_s = [.1 for _ in range(3)]
    #    x = ht_pure.recursiveMedians(height,array,budget_s)
    #    print x
    #    print len(x)

    # test recursiveSlices
    #    height = 2
    #    array = range(1001)
    #    budget_s = [0 for _ in range(3)]
    #    x = ht_pure.recursiveSlices(height,8,array,budget_s)
    #    print x
    #    print len(x)

    # test getRecursiveMediansNo
    budget = .1
    array = range(1000000)
    x = ht_pure.getRecursiveMediansNo(array, budget)
    print x