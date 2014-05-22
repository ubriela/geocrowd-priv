"""
Utility methods used in Geocast algorithm
"""
import numpy as np
import os
from Params import Params
from PSDExp import data_readin
import math
from sets import Set
#from scipy import spatial
from Utils import distance, distance_point

initQueryShape = [1,1]

def query_init(x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """
    Init a random query of some specific size within a rect [[x1,y1],[x2,y2]]
    """
    x_range = (x2-x1)*2**initQueryShape[0]/(2**Params.queryUnit[0])
    y_range = (y2-y1)*2**initQueryShape[1]/(2**Params.queryUnit[1])

    data = data_readin()
    ran_indices = np.random.randint(0,data.shape[1],1)
    ran_points = data[:,ran_indices]
    x_low = ran_points[0,:] - x_range/2 
    x_high = ran_points[0,:] + x_range/2
    y_low = ran_points[1,:] - y_range/2
    y_high = ran_points[1,:] + y_range/2

    query = [[max(x_low[0],x1),max(y_low[0],y1)],[min(x_high[0],x2),min(y_high[0],y2)]]
    return np.array(query)


def query_init_naive(t, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """
    init a query around location L
    """    
    x_range = (x2-x1)*2**initQueryShape[0]/(2**Params.queryUnit[0])
    y_range = (y2-y1)*2**initQueryShape[1]/(2**Params.queryUnit[1])

    x_low = t[0] - x_range/2
    x_high = t[0] + x_range/2
    y_low = t[1] - y_range/2
    y_high = t[1] + y_range/2

    query = [[max(x_low,x1),max(y_low,y1)],[min(x_high,x2),min(y_high,y2)]]
    return np.array(query)

def query_gen(queryShape, taskNo, seed, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """Generate query around a random data point"""
    
    np.random.seed(seed)
    querylist = []
    cell_size_x = (x2-x1)/(2**Params.queryUnit[0])
    cell_size_y = (y2-y1)/(2**Params.queryUnit[1])
    x_range, y_range = cell_size_x*2**queryShape[0], cell_size_y*2**queryShape[1]
    
    data = data_readin()
    ran_indices = np.random.randint(0,data.shape[1],taskNo)
    ran_points = data[:,ran_indices]
    x_low = ran_points[0,:] - x_range/2
    x_high = ran_points[0,:] + x_range/2
    y_low = ran_points[1,:] - y_range/2
    y_high = ran_points[1,:] + y_range/2
    for i in range(taskNo):
        query = [[max(x_low[i],x1),max(y_low[i],y1)],[min(x_high[i],x2),min(y_high[i],y2)]]
        querylist.append(np.array(query))

    return querylist

"""
Generate random points within data domain
"""
def task_locs_gen(taskNo, seed, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    np.random.seed(seed)
    x = np.random.random_sample(taskNo)*(x2-x1) + x1
    y = np.random.random_sample(taskNo)*(y2-y1) + y1
    return np.array([x,y]).transpose()


"""
Generate random points within dataset
"""
def worker_locs_gen(data, workerNo, seed, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    np.random.seed(seed)
    ran_indices = np.random.randint(0,data.shape[1],workerNo)
    ran_points = data[:,ran_indices]
    return ran_points.transpose()
    
def rect_query_count(data, query):
    """Get true answer by linear search along each dimension"""
    return rect_query_points(data, query).shape[1]


def rect_query_points(data, query):
    """Get true answer by linear search along each dimension"""
    _data = data.copy()
    _ndim = _data.shape[0] 
    for dim in range(_ndim):
        if _data.shape[1] == 0:
            break
        idx = np.argsort(_data[dim,:],kind='mergesort')
        _data[:,:] = _data[:,idx]
        x = np.searchsorted(_data[dim,:], query[0,dim],side='left')
        y = np.searchsorted(_data[dim,:], query[1,dim],side='right')
        _data = _data[:,x:y]

    return _data
    
def resize(query, scale, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """ 
    Resize the query by a scale factor
    
    @scale = 1 means the same query
    """
    center = [(query[0,0] + query[1,0])/2.0,(query[0,1] + query[1,1])/2.0]
    x_range = query[1,0] - query[0,0]
    y_range = query[1,1] - query[0,1]
    x_range = x_range * scale
    y_range = y_range * scale
    
    x_low = center[0] - x_range/2
    x_high = center[0] + x_range/2
    y_low = center[1] - y_range/2
    y_high = center[1] + y_range/2
    query = [[max(x_low,x1),max(y_low,y1)],[min(x_high,x2),min(y_high,y2)]]
    return np.array(query)

def query_tree(query_list, tree, data, aligned_queries = None, aligned_percent = None):
        """ wrapper for query answering and computing query error """
        result = []
        trueRes = []
        if aligned_queries == None:
            for query in query_list:
                result.append(tree.rangeCount(query))
            Res = np.array(result)
            trueRes = np.array([rect_query_count(data, query) for query in query_list])
            return compute_error(Res, trueRes) 
        else:
            for i in range(len(query_list)):
                result.append(tree.rangeCount(query_list[i]) * aligned_percent[i])
                trueRes.append(rect_query_count(data, aligned_queries[i]))

            return compute_error(np.array(result), np.array(trueRes)) 
        
def compute_error(Res, trueRes):
        """ Compute median absolute and relative errors """
        absErr = np.abs(Res-trueRes)
        idx_nonzero = np.where(trueRes != 0)
        absErr_nonzero = absErr[idx_nonzero]
        true_nonzero = trueRes[idx_nonzero]
        relErr = absErr_nonzero/true_nonzero    
        return np.average(relErr), np.average(absErr)


#x_min = sys.float_info.max
#y_min = sys.float_info.max
#x_max = -sys.float_info.max
#y_max = -sys.float_info.max

def generate_workers(seed, time_instance):
    """ Generate a set of workers per time instance"""
    global x_min, y_min, x_max, y_max
    data = data_readin()
    worker_locs = worker_locs_gen(data, Params.WorkerNo, seed,Params.x_min,Params.y_min,Params.x_max,Params.y_max)
    filename = "../dataset/taskworker/workers" + str(time_instance) + ".txt"
    if os.path.exists(filename):
        os.remove(filename)
    for loc in worker_locs:
        with open(filename, "a") as workers:
            workers.write(str(loc[0]) + ", " + str(loc[1])  + "\n")

def generate_tasks(seed, time_instance):
    """ Generate a set of tasks per time instance"""
    data = data_readin()
    task_locs = task_locs_gen(Params.TASK_NO, seed, Params.x_min,Params.y_min,Params.x_max,Params.y_max)
    filename = "../dataset/taskworker/tasks" + str(time_instance) + ".txt"
    if os.path.exists(filename):
        os.remove(filename)
    for loc in task_locs:
        with open(filename, "a") as tasks:
            tasks.write(str(loc[0]) + ", " + str(loc[1])  + "\n")


def hops_expansion2(task, workers, network_diameter):
    """
    Find the number of hops (hops count) in infra structureless communication
    @param task : task loction
    @param workers : worker locations
    @param diameter : network diameter (some constant)
    """
    hops_count2 = simple_hop_count(mbr(workers))
    
    queue = Set()  # queue
    traversed = Set([])
    hops_count = 0
    
    # create kdtree
    tree = spatial.KDTree(workers)
    queue = tree.query_ball_point(task, network_diameter)
    traversed = traversed | Set(queue)
    
    while len(queue) > 0:
        hops_count += 1
        candidates = Set([])
        for i in queue: 
            candidates = candidates | Set(tree.query_ball_point(workers[i], network_diameter))
        candidates = candidates - traversed
        queue = candidates
        traversed = traversed | candidates
    print hops_count, (len(traversed)+0.0)/len(workers), hops_count2
    return hops_count, (len(traversed)+0.0)/len(workers), hops_count2

def ball_query_points(seed, workers, network_diameter):
    points = Set()
    for point in workers:
        if distance(seed[0], seed[1], point[1], point[2]) <= network_diameter:
            points.add(point)
    return points
    
# don't use spatial.KDTree
def hops_expansion(task, workers, network_diameter):
    """
    Find the number of hops (hops count) in infra structureless communication
    @param task : task loction
    @param workers : worker locations
    @param diameter : network diameter (some constant)
    """
    hops_count2 = simple_hop_count(mbr(workers))
    
    queue = Set()  # queue
    traversed = Set([])
    remained = Set([])
    hops_count = 0
    
    # create kdtree
    remained = Set([])
    for i in range(len(workers)):
        remained.add((i, workers[i][0], workers[i][1]))
    queue.add((-1, task[0], task[1]))

    while len(queue) > 0:
        hops_count += 1
        candidates = Set([])
        for worker in queue:  
            candidates = candidates | ball_query_points((worker[1], worker[2]), remained, network_diameter)
        candidates = candidates - traversed
        queue = candidates
        traversed = traversed | candidates
        remained = remained - traversed

#    print hops_count, (len(traversed)+0.0)/len(workers), hops_count2
    return hops_count, math.ceil((len(traversed)+0.0)/len(workers)), hops_count2

def generate_worker_task(no_instance):
    """ Generate tasks and workers per time instance"""
    for i in range(no_instance):
        seed = random.randint(1,1000)
        generate_workers(seed, i);
        generate_tasks(seed, i);

def mbr(coords):
    min_x, min_y = np.min(coords[...,0]), np.min(coords[...,1])
    max_x, max_y = np.max(coords[...,0]), np.max(coords[...,1])
    return np.array([[min_x, min_y], [max_x, max_y]])

def simple_hop_count(mbr):
    return distance(mbr[0][0], mbr[0][1], mbr[1][0], mbr[1][1])/(2*Params.NETWORK_DIAMETER)

if __name__ == '__main__':
#    generate_worker_task()
#    print task_locs_gen(5,2,0,0,10,10)
#    query = [[0,0],[1,1]]
#    print resize(np.array(query), 2, Params.x_min,Params.y_min,Params.x_max,Params.y_max)

    workers = np.array(
    [[0, 0],
     [0, 1],
     [0, 2],
     [0, 3],
     [1, 0],
     [1, 1],
     [1, 2],
     [1, 3],
     [2, 0],
     [2, 1],
     [2, 2],
     [2, 3],
     [3, 0],
     [3, 1],
     [3, 2],
     [3, 3]])
    print hops_expansion([2,2], workers, 1)
    
    query = np.array([[1,1],[3,2.5]])
    print rect_query_points(workers.transpose(), query)
    print rect_query_count(workers.transpose(), query)
