from HT_composite import HT_composite

class HT_composite_localness(HT_composite):
   
    def __init__(self, data, param):
        HT_composite.__init__(self, data, param)
        
        
    def getCountBudget(self):
	"""allocate budget to 2nd htree layer only"""
        count_eps = self.param.Eps*(1-self.param.PercentSplit)
        if self.param.geoBudget == 'optimal':
            n1 = 2**(self.maxSplit[0]-self.switchSplit[0])
            n2 = n1*2**(self.maxSplit[1]-self.switchSplit[1])
            unit = count_eps / (n1**(1.0/3) + n2**(1.0/3))
            eps1 = unit*n1**(1.0/3)
            eps2 = unit*n2**(1.0/3)
            return [0,0,0,eps1,eps2]
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)