import numpy as np
from collections import deque
from HT_standard import HT_standard
from Params import Params

class HT_hybrid(HT_standard):
    """ Hybrid ht-tree, use private medians for the first h levels and (data-independent) 
    domain mid-points for the remaining levels as the split points. """
    
    def __init__(self, data, param):         
        HT_standard.__init__(self, data, param)
        self.param.maxHeightHTree = 4;
        self.gridBudget = 0;
        
        
    def testLeaf(self, curr):
        """test whether a node is a leaf node"""
        if curr.n_depth == self.param.maxHeightHTree:
            return True;
        
        if (curr.n_depth <= 2):
            if (curr.n_data is None or curr.n_data.shape[1] == 0) or \
                    (curr.n_count < self.param.minPartSizeHTree):
                    return True
        return False
    
    def getSplitBudget(self, curr):
        split_eps = self.param.Eps*self.param.PercentSplit/2
	curr_depth = curr.n_depth
        dimP = curr_depth % Params.NDIM # split dimension
	sw = self.switchSplit[dimP] # switching level
	split_height = self.maxSplit[dimP] - sw
        
	if (curr_depth <=1):
            return [split_eps/sw for _ in range(sw)]
	else:
            return [0 for _ in range(split_height)]
            
    def getCountBudget(self):
        count_eps = self.param.Eps*(1-self.param.PercentSplit)
        if Params.useLeafOnlyHTree:
            return [0, 0, 0, 0, count_eps]
	sw = [0 for _ in range(3)]
        sw[0] = self.switchSplit[0] # switching level
        sw[1] = self.switchSplit[1]
        sw[2] = (self.maxSplit[0]-sw[0])*(self.maxSplit[1]-sw[1])
            
        if self.param.geoBudget == 'optimal':
            eps = [0 for _ in range(3)]
            sum = 0
            for i in range(3):
                if i == 0:
                    eps[i] = 2**sw[i]
                else:
                    eps[i] = eps[i-1]*2**sw[i]
                sum += eps[i]**(1.0/3)
            unit = count_eps / sum
            ret = [unit*eps[j]**(1.0/3) for j in range(3)]
            self.gridBudget = ret[2]
            return [0, ret[0], ret[1], 0, ret[2]]
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)
	
    
    def getCoordinates(self, curr):
        """ 
        get corrdinates of the point which defines the four subnodes: 
        split_arr, n_data_arr. We defines the subnodes same as in HT_pure except 
	deal with the number of split in each level of the tree
        """
        budget_s = self.getSplitBudget(curr)
        _data = curr.n_data
        _box = curr.n_box
	curr_depth = curr.n_depth
        dimP = curr_depth % Params.NDIM # split dimension

        isSorted = False
	if curr_depth <= 1: # data dependent
	    _idx = np.argsort(_data[dimP,:], kind ='mergesort')
	    _data[:,:] = _data[:,_idx] # sorted by dimP dimension
            isSorted = True

            swPartitions = self.param.switchPartitionsHTree[dimP] # 1st htree layer
            split_arr = self.recursiveSlices(len(budget_s)-1, swPartitions, _data[dimP,:], budget_s)
	    split_arr.insert(0,_box[0,dimP])
	    split_arr.append(_box[1,dimP])
	else: # data independent (grid)
            swPartitions = self.param.partitionsHTree[dimP]/self.param.switchPartitionsHTree[dimP] # 2nd htree layer
            split_arr = self.getEqualSplit(len(budget_s)-1, _box[0,dimP], _box[1,dimP])

        # get data points in these partitions
	n_data_arr = [None for _ in range(len(split_arr)-1)]
	if _data is not None and _data.shape[1] >= 1:
            if not isSorted:
                _idx = np.argsort(_data[dimP,:], kind ='mergesort')
                _data[:,:] = _data[:,_idx] # sorted by dimP dimension
	    for i in range(len(split_arr)-1):
		posP1 = np.searchsorted(_data[dimP,:], split_arr[i])
		posP2 = np.searchsorted(_data[dimP,:], split_arr[i+1])
		if i == 0:
		    n_data = _data[:,:posP2]
		elif i == len(split_arr)-2:
		    n_data = _data[:,posP1:]
		else:
		    n_data = _data[:,posP1:posP2]
		n_data_arr[i] = n_data
                
        return (split_arr),(n_data_arr)
    
    
    
    # do not apply canonical range query on level 3
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
                if self.intersect(_box,query):
                    for i in range(_box.shape[1]):
                        if _box[1,i] == _box[0,i] or Params.WorstCase == True:
                            frac *= 1
                        else:
                            frac *= (min(query[1,i],_box[1,i])-max(query[0,i],_box[0,i])) / (_box[1,i]-_box[0,i])
                    count += curr.n_count*frac
            else: # if not leaf
                for node in curr.children:
                    bool_matrix = np.zeros((2,query.shape[1]))
                    bool_matrix[0,:] = query[0,:] <= _box[0,:]
                    bool_matrix[1,:] = query[1,:] >= _box[1,:]
                    
                    if not Params.useLeafOnlyHTree and np.all(bool_matrix) and node.n_depth != 3: # if query range contains node range
                        count += node.n_count
                    elif self.intersect(_box,query):
                        queue.append(node)
        return float(count)