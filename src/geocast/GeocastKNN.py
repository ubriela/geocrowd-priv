"""
Geocast without privacy
"""

import numpy as np
from Params import Params
from Utils import performed_task
from Geocrowd import rect_query_points, hops_expansion

#from scipy import spatial
from Utils import distance, acc_rate, is_performed, performed_tasks

def geocast_knn(data, t):
    # find all workers in MTD
    
    # find all workers in the query
    MTD_RECT = np.array([[t[0]-Params.ONE_KM*Params.MTD,t[1]-Params.ONE_KM*Params.MTD],[t[0]+Params.ONE_KM*Params.MTD, t[1]+Params.ONE_KM*Params.MTD]])
    locs = rect_query_points(data, MTD_RECT).transpose()
    locs = sorted(locs, key = lambda loc: distance(loc[0], loc[1], t[0], t[1]))

    u, dist, found = 0, 0, False 
    workers = np.zeros(shape=(2,0))
    for loc in locs:
        workers = np.concatenate([workers, np.array([[loc[0]],[loc[1]]])], axis=1)
        _dist = distance(loc[0], loc[1], t[0], t[1])
        u_c = acc_rate(Params.MTD, _dist)
        u = 1 - (1-u)*(1-u_c)
        if is_performed(u_c):
            if not found:
                found = True
                dist = _dist
        if u >= Params.U:
            break
        
    # simulation
    isPerformed, worker, dist_fcfs  = performed_tasks(workers, Params.MTD, t, True)
    hops_count, coverage, hops_count2 = hops_expansion(t, workers.transpose(), Params.NETWORK_DIAMETER)
    
    if isPerformed:     # the task is performed 
        return workers.shape[1], True, dist, dist_fcfs, hops_count, coverage, hops_count2
    
    return workers.shape[1], False, None, None, hops_count, coverage, hops_count2