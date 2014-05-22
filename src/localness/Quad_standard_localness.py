import logging
from Quad_standard import Quad_standard
from Params import Params

class Quad_standard_localness(Quad_standard):
    
    def __init__(self, data, param):
        Quad_standard.__init__(self, data, param)


    def getCountBudget(self):
        count_eps = self.param.Eps
        H = Params.maxHeight
        Hx = Params.Hx
        
        if self.param.geoBudget == 'none':
            ret = [0 for i in range(Hx+1,H+1)]
            for i in range(0,Hx+1):
                ret.append(count_eps/(Hx+1))
            return ret
        elif self.param.geoBudget == 'optimal':
            ret = [0 for i in range(Hx+1,H+1)]
            unit = count_eps * ((2**(1.0/3))-1) / (2**((1.0/3)*(Hx+1))-1)
            for i in range(0,Hx+1):
                ret.append(unit*2**((1.0/3)*i))
            return ret
        else:
            logging.error('No such geoBudget scheme')
            sys.exit(1)
    