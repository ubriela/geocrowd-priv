import numpy as np
from HT_standard import HT_standard
from Params import Params

class HT_composite(HT_standard):
    """
    4 layers htrees
    """
    
    def __init__(self, data, param):
        HT_standard.__init__(self, data, param)
        self.param.maxHeightHTree = 4
                    
    def getSplitBudget(self, curr):
	curr_depth = curr.n_depth
        if curr_depth == 0:
            split_eps = self.param.Eps*self.param.PercentSplit*self.param.splitFrac**2
        elif curr_depth == 1 or curr_depth == 2:
            split_eps = self.param.Eps*self.param.PercentSplit*(self.param.splitFrac)*(1-self.param.splitFrac)
        elif curr_depth == 3:
            split_eps = self.param.Eps*self.param.PercentSplit*(1-self.param.splitFrac)*(1-self.param.splitFrac)
        dimP = curr_depth % Params.NDIM # split dimension
	sw = self.switchSplit[dimP] # the number of recursive split
        split_2nd = self.maxSplit[dimP] - sw
        
	if (curr_depth <=1):
            return [split_eps/(2*sw) for _ in range(sw)]
	else:
            return [split_eps/(2*split_2nd) for _ in range(split_2nd)]

    def getCountBudget(self):
        count_eps = self.param.Eps*(1-self.param.PercentSplit)
        if Params.useLeafOnlyHTree:
            return [0, 0, 0, 0, count_eps]
	sw = [0 for _ in range(4)]
        sw[0] = self.switchSplit[0] # switching level
        sw[1] = self.switchSplit[1]
        sw[2] = self.maxSplit[0]-sw[0]
        sw[3] = self.maxSplit[1]-sw[1]
	
        if self.param.geoBudget == 'optimal':
            eps = [0 for _ in range(4)]
            sum = 0
            for i in range(4):
                if i == 0:
                    eps[i] = 2**sw[i]
                else:
                    eps[i] = eps[i-1]*2**sw[i]
                sum += eps[i]**(1.0/3)
            unit = count_eps / sum
	    ret = [unit*eps[i]**(1.0/3) for i in range(4)]
	    ret.insert(0,0)
            return ret
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
	
        _idx = np.argsort(_data[dimP,:], kind ='mergesort')
        _data[:,:] = _data[:,_idx] # sorted by dimP dimension

        if curr_depth <= 1:
            swPartitions = self.param.switchPartitionsHTree[dimP] # 1st htree layer
        else:
            swPartitions = self.param.partitionsHTree[dimP]/self.param.switchPartitionsHTree[dimP] # 2nd htree layer

        split_arr = self.recursiveSlices(len(budget_s)-1, swPartitions, _data[dimP,:], budget_s)
            
        split_arr.insert(0,_box[0,dimP])
        split_arr.append(_box[1,dimP])

        # get data points in these partitions
        n_data_arr = []
        for i in range(len(split_arr)-1):
	    posP1 = np.searchsorted(_data[dimP,:], split_arr[i])
	    posP2 = np.searchsorted(_data[dimP,:], split_arr[i+1])
            if i == 0:
                n_data = _data[:,:posP2]
            elif i == len(split_arr)-2:
                n_data = _data[:,posP1:]
            else:
                n_data = _data[:,posP1:posP2]
            n_data_arr.append(n_data)
                
        return (split_arr),(n_data_arr)