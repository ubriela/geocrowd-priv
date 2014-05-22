from Kd_hybrid import Kd_hybrid
from Params import Params

class Kd_hybrid_localness(Kd_hybrid):
    """ standard private kd-tree """
    def __init__(self, data, param):
        Kd_hybrid.__init__(self, data, param)
    
    
    def getCountBudget(self):
        count_eps = self.param.Eps
        H = Params.maxHeight
        Hx = Params.Hx
        
        ret = [0 for i in range(Hx+1,H+1)]
        unit = count_eps * ((2**(1.0/3))-1) / (2**((1.0/3)*(Hx+1))-1)
        for i in range(0,Hx+1):
            ret.append(unit*2**((1.0/3)*i))
        return ret