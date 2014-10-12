from collections import deque

from Generic import Generic


class HTree(Generic):
    """
    Generic high fanout tree
    """

    def __init__(self, data, param):
        Generic.__init__(self, data, param)

    def getSplitBudget(self, curr):
        """return a list of budget values for split"""
        raise NotImplementedError

    def getNoisyMedian(self, array, left, right, epsilon):
        """return the split value of an array"""
        raise NotImplementedError

    def getSplit(self, array, left, right, epsilon=0, partitions=0, slice=0):
        """
        return the split point given an array, may be data-independent or
        true median or noisy median, depending on the type of the tree.
        partitions is the number of slices, used in slicing mechanism.
        slice is the splitting position      
        """
        raise NotImplementedError

    def adjustConsistency(self):  # used for htree variants
        queue = deque()
        for child in self.root.children:  # do not apply constraint inference on root node
            queue.append(child)
        while len(queue) > 0:
            curr = queue.popleft()
            # if curr is not None and curr.n_isLeaf is False:
            if curr.n_isLeaf is False:
                sum = 0
                for child in curr.children:
                    sum += child.n_count
                adjust = (curr.n_count - sum) / len(curr.children)

                for child in curr.children:
                    child.n_count += adjust
                    queue.append(child)