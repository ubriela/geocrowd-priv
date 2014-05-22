from HT_standard import HT_standard
from Params import Params

class HT_standard_skew(HT_standard):
    
    def __init__(self, data, param):
        HT_standard.__init__(self, data, param)
    
    def getSplitBudget(self, curr):
        split_eps = self.param.Eps*self.param.PercentSplit/2
	curr_depth = curr.n_depth
        dimP = curr_depth % Params.NDIM # split dimension
	split_height = self.maxSplit[dimP]
        
	budgets = self.getSkewMedianBudget(split_height,split_eps)
	return budgets
        
    def getSkewMedianBudget(self, medianNumber, totalBudget):
        budgets = []
        unitBudget = totalBudget/(2**medianNumber-1);
        for i in range(medianNumber):
            budgets.append(unitBudget*2**i)
        return budgets