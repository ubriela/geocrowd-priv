from HT_hybrid import HT_hybrid

class HT_hybrid_localness(HT_hybrid):
   
    def __init__(self, data, param):
        HT_hybrid.__init__(self, data, param)
            
    def getCountBudget(self):
        count_eps = self.param.Eps*(1-self.param.PercentSplit)
        if self.param.geoBudget == 'optimal':
            return [0, 0, 0, 0, count_eps]
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)