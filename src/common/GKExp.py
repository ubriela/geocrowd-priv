"""
Collection class to init and build all other actual indices
"""
import numpy as np
import time
import logging
from Params import Params
from Kd_true import Kd_true
from Kd_standard import Kd_standard
from Kd_hybrid import Kd_hybrid
from Quad_standard import Quad_standard
from Quad_standard_localness import Quad_standard_localness
from Kd_true_localness import Kd_true_localness
from Kd_standard_localness import Kd_standard_localness
from Kd_hybrid_localness import Kd_hybrid_localness

from HT_pure import HT_pure
from HT_true import HT_true
from HT_standard import HT_standard
from Grid_pure import Grid_pure
from Grid_uniform import Grid_uniform
from Grid_adaptive import Grid_adaptive
from Grid_adaptive_localness import Grid_adaptive_localness
from HT_hybrid import HT_hybrid
from HT_composite import HT_composite
from HT_hybrid_skew import HT_hybrid_skew
from HT_standard_skew import HT_standard_skew
from HT_standard_adaptive import HT_standard_adaptive
from HT_composite_localness import HT_composite_localness
from HT_hybrid_localness import HT_hybrid_localness

from Log import log

class GKExp(object):
    
    def __init__(self, data, query_list):
        self.data = data
        self.query_list = query_list
        logging.debug('Getting true query answers...')
        self.trueRes = np.array([self.getTrue(query) for query in query_list])
        self.selectivity()
        
        
    def getTrue(self, query):
        """Get true answer by linear search along each dimension"""
        _data = self.data.copy()
        _ndim = _data.shape[0]
        for dim in range(_ndim):
            if _data.shape[1] == 0:
                break
            idx = np.argsort(_data[dim,:],kind='mergesort')
            _data[:,:] = _data[:,idx]
            x = np.searchsorted(_data[dim,:], query[0,dim],side='left')
            y = np.searchsorted(_data[dim,:], query[1,dim],side='right')
            _data = _data[:,x:y+1]
#        print "true value" + str(_data.shape[1])
        return _data.shape[1]
    
    def selectivity(self):
        sel_array = np.sort(self.trueRes/float(Params.NDATA))
        logging.debug('selectivity min %.2f' % np.min(sel_array))
        logging.debug('selectivity max %.2f' % np.max(sel_array))
        logging.debug('selectivity avg %.2f' % np.average(sel_array))
        logging.debug('selectivity median %.2f' % sel_array[Params.nQuery/2])
        logging.debug('selectivity first quartile %.2f' % sel_array[Params.nQuery/4])

    def query(self, tree, method="None"):
        """ wrapper for query answering and computing query error """
        result = []
        for query in self.query_list:
            result.append(tree.rangeCount(query))
        Res = np.array(result)
        return self.computeError(Res, method)
        
    def computeError(self, Res, method = "None"):
        """ Compute median absolute and relative errors """
        absErr = np.abs(Res-self.trueRes)
        idx_nonzero = np.where(self.trueRes != 0)
        absErr_nonzero = absErr[idx_nonzero]
        true_nonzero = self.trueRes[idx_nonzero]
        relErr = absErr_nonzero/true_nonzero
	
#        log_str_rel = "\n".join(map(str, relErr))
#        log_str_abs = "\n".join(map(str, absErr))
            
        if Params.IS_LOGGING:
	    log_str = ""
            for i in range(len(self.query_list)):
#		area_str = 
                query_str = str(self.query_list[i][0][1]) + "\t" + str(self.query_list[i][0][1]) + "\t" + str(self.query_list[i][1][0]) + "\t" + str(self.query_list[i][1][1])
		err_str = str(Res[i]) + "\t" + str(self.trueRes[i]) + "\t" + str(relErr[i]) + "\t" +  str(absErr[i])
		log_str = log_str + query_str + err_str
            log(method , log_str_rel)
        
	absErr = np.sort(absErr)
        relErr = np.sort(relErr)
        n_abs = len(absErr)
        n_rel = len(relErr)
        return absErr[int(n_abs/2)], relErr[int(n_rel/2)]
    
    def run_HT_standard_adaptive(self, param):
        logging.debug('building HT_standard_adaptive...')
        tree = HT_standard_adaptive(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_standard_adaptive building time: %.2d ' %(end-start))
	return self.query(tree)
    
    def run_HT_standard_skew(self, param):
        logging.debug('building HT_standard_skew...')
        tree = HT_standard_skew(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_standard_skew building time: %.2d ' %(end-start))
	return self.query(tree)
    
    
    def run_HT_hybrid_skew(self, param):
        logging.debug('building HT_hybrid_skew...')
        tree = HT_hybrid_skew(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_hybrid_skew building time: %.2d ' %(end-start))
	return self.query(tree)
    
    def run_HT_hybrid_localness(self, param):
        logging.debug('building run_HT_hybrid_localness...')
        tree = HT_hybrid_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
        end = time.clock()
        logging.info('[T] run_HT_hybrid_localness building time: %.2d ' %(end-start))
	return self.query(tree)
    
    def run_HT_composite_localness(self, param):
        logging.debug('building HT_composite_localness...')
        tree = HT_composite_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_composite_localness building time: %.2d ' %(end-start))
	return self.query(tree)

    def run_HT_composite(self, param):
        logging.debug('building HT_composite...')
        tree = HT_composite(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_composite building time: %.2d ' %(end-start))
	return self.query(tree)
 
    def run_HT_hybrid(self, param):
        logging.debug('building HT_hybrid...')
        tree = HT_hybrid(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_hybrid building time: %.2d ' %(end-start))
	return self.query(tree)

    def run_Grid_adaptive(self, param):
        logging.debug('building Grid_adaptive...')
        tree = Grid_adaptive(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()        
        end = time.clock()
        logging.info('[T] Grid_adaptive building time: %.2d ' %(end-start))
        return self.query(tree, "Grid_adaptive")
    
    def run_Grid_adaptive_localness(self, param):
        logging.debug('building Grid_adaptive_localness...')
        tree = Grid_adaptive_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] Grid_adaptive_localness building time: %.2d ' %(end-start))
        return self.query(tree)
    
    def run_Grid_uniform(self, param):
        logging.debug('building Grid_uniform...')
        tree = Grid_uniform(self.data, param)
        start = time.clock()
        tree.buildIndex()
        end = time.clock()
        logging.info('[T] Grid_uniform building time: %.2d ' %(end-start))
        return self.query(tree, "Grid_uniform")
    
    def run_Grid_pure(self, param):
        logging.debug('building Grid_pure...')
        tree = Grid_pure(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] Grid_pure building time: %.2d ' %(end-start))
        return self.query(tree)
    
    def run_HT_standard(self, param):
        logging.debug('building HT_standard...')
        tree = HT_standard(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_standard building time: %.2d ' %(end-start))
        return self.query(tree, "HT_standard")
    
    def run_HT_true(self, param):
        logging.debug('building HT_true...')
        tree = HT_true(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_true building time: %.2d ' %(end-start))
        return self.query(tree)
    
    def run_HT_pure(self, param):
        logging.debug('building HT_pure...')
        tree = HT_pure(self.data, param)
        start = time.clock()
        tree.buildIndex()
        if Params.CONSTRAINT_INFERENCE:
                tree.adjustConsistency()         
        end = time.clock()
        logging.info('[T] HT_pure building time: %.2d ' %(end-start))
        return self.query(tree)
        
    def run_Kd_pure(self, param):
        logging.debug('building Kd_pure...')
        tree = Kd_pure(self.data, param)
        start = time.clock()
        tree.buildIndex()
        end = time.clock()
        logging.info('[T] Kd-pure building time: %.2f' % (end-start))
        return self.query(tree)
    
    def run_Kd_true(self, param):
        logging.debug('building Kd_true...')
        tree = Kd_true(self.data, param)
        start = time.clock()
        tree.buildIndex()
#        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd-true building time: %.2f' % (end-start))
        return self.query(tree)
    
    def run_Kd_standard(self, param):
        logging.debug('building Kd_standard...')
        tree = Kd_standard(self.data, param)
        start = time.clock()
        tree.buildIndex()
#        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd-standard building time: %.2f' % (end-start))
        return self.query(tree)
    
    def run_Kd_hybrid(self, param):
        logging.debug('building Kd_hybrid...')
        tree = Kd_hybrid(self.data, param)
        start = time.clock()
        tree.buildIndex()
        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd-hybrid building time: %.2f' % (end-start))
        return self.query(tree, "Kd_hybrid")

    def run_Kd_true_localness(self, param):
        logging.debug('building Kd_true_localness...')
        tree = Kd_true_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
#        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd_true_localness building time: %.2f' % (end-start))
        return self.query(tree)

    def run_Kd_standard_localness(self, param):
        logging.debug('building Kd_standard_localness...')
        tree = Kd_standard_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
#        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd_standard_localness building time: %.2f' % (end-start))
        return self.query(tree)

    def run_Kd_hybrid_localness(self, param):
        logging.debug('building Kd_hybrid_localness...')
        tree = Kd_hybrid_localness(self.data, param)
        start = time.clock()
        tree.buildIndex()
#        tree.adjustConsistency()
        end = time.clock()
        logging.info('[T] Kd_hybrid_localness building time: %.2f' % (end-start))
        return self.query(tree)
    
    def run_Quad_baseline(self, param):
        logging.debug('building Quad_baseline...')
        param.geoBudget = 'none'
        tree = Quad_standard(self.data, param)
        tree.buildIndex()
        return self.query(tree)
    
    def run_Quad_geo(self, param):
        logging.debug('building Quad_geo...')
        param.geoBudget = 'optimal'
        tree = Quad_standard(self.data, param)
        tree.buildIndex()
        return self.query(tree, "Quad_geo")
    
    def run_Quad_baseline_localness(self, param):
        logging.debug('building Quad_baseline_localness...')
        param.geoBudget = 'none'
        tree = Quad_standard_localness(self.data, param)
        tree.buildIndex()
        return self.query(tree)
    
    def run_Quad_geo_localness(self, param):
        logging.debug('building Quad_geo_localness...')
        param.geoBudget = 'optimal'
        tree = Quad_standard_localness(self.data, param)
        tree.buildIndex()
        return self.query(tree)