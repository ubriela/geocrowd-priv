from HT_hybrid import HT_hybrid
from Params import Params

class HT_hybrid_skew(HT_hybrid):
    
    def __init__(self, data, param):
        HT_hybrid.__init__(self, data, param)
    
    def getSplitBudget(self, curr):
        split_eps = self.param.Eps*self.param.PercentSplit/2
	curr_depth = curr.n_depth
        dimP = curr_depth % Params.NDIM # split dimension
	sw = self.switchSplit[dimP] # switching level
	split_height = self.maxSplit[dimP] - sw
        
	if (curr_depth <=1):
            budgets = self.getSkewMedianBudget(sw,split_eps)
            return budgets
	else:
            return [0 for _ in range(split_height)]
        
    def getSkewMedianBudget(self, medianNumber, totalBudget):
        budgets = []
        unitBudget = totalBudget/(2**medianNumber-1);
        for i in range(medianNumber):
            budgets.append(unitBudget*2**i)
        return budgets