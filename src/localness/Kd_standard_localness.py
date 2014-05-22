from Kd_standard import Kd_standard
from Params import Params

class Kd_standard_localness(Kd_standard):
    """ similar to Kd-standard, except that noisy mean is used as the split point
    we only use leaves to answer queries in this case """
    def __init__(self, data, param):
        Kd_standard.__init__(self, data, param)
    
    def getCountBudget(self):
        count_eps = self.param.Eps*(1-self.param.Percent)
        H = Params.maxHeight
        Hx = Params.Hx
        
        ret = [0 for i in range(Hx+1,H+1)]
        unit = count_eps * ((2**(1.0/3))-1) / (2**((1.0/3)*(Hx+1))-1)
        for i in range(0,Hx+1):
            ret.append(unit*2**((1.0/3)*i))
        return ret
    