"""
Contains all experimental methods
"""
import numpy as np
import time
import os
import multiprocessing as mult
import sys
import logging

sys.path.append('..')
sys.path.append('../common')
sys.path.append('../htree')
sys.path.append('../grid')
sys.path.append('../icde12')
sys.path.append('../localness')

from Params import Params
from GKExp import GKExp

seed_list = [2172]
shape_list = [(1,1)]
eps_list = [0.4]

#seed_list = [9924, 8253, 3375]
#shape_list = [(1,1),(2,2),(3,3)]
#eps_list = [0.1,0.4,0.7,1.0]

method_list = None
exp_name = None

def data_readin():
    """Read in spatial data and initialize global variables."""
    p = Params(0)
    p.select_dataset()
    data = np.genfromtxt(Params.dataset,unpack = True)
    Params.NDIM, Params.NDATA = data.shape[0], data.shape[1]
    Params.LOW, Params.HIGH = np.amin(data,axis=1), np.amax(data,axis=1)
    logging.debug(data.shape)
    logging.debug(Params.LOW)
    logging.debug(Params.HIGH)
    return data


    all_points = []
    if os.path.isfile(Params.TASKPATH):
        with open(Params.TASKPATH) as f:
            content = f.readlines()
        for i in range(len(seed_list)):
            ran_points = []
            for j in range(taskNo):
                ran_points.append(map(float, content[i*taskNo + j].split()))
            all_points.append(ran_points)
    else:
        tasks = ""
        logging.debug('tasks_gen: generating tasks...')
        
        boundary = np.array([[x1,y1],[x2,y2]])
        for seed in seed_list:
            ran_points = []
            np.random.seed(seed)
            count = 0
            while count < taskNo:
                idx = np.random.randint(0,data.shape[1])
                _ran_point = data[:,idx]
                if is_rect_cover(boundary, _ran_point):
                    ran_points.append(_ran_point)
                    count = count + 1
            all_points.append(ran_points)
            for item in ran_points:
                tasks = tasks + ("%s\n" % " " . join(map(str, item)))
        outfile = open(Params.TASKPATH, "w")
        outfile.write(tasks)
        outfile.close()
    return all_points

def query_gen(queryShape, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """Generate query around a random data point"""
    
    logging.debug('Generating queries...')
    all_queries = []
    for seed in seed_list:
	querylist = []
        np.random.seed(seed)
        cell_size_x = (x2-x1)/(2**Params.queryUnit[0])
        cell_size_y = (y2-y1)/(2**Params.queryUnit[1])
        x_range, y_range = cell_size_x*2**queryShape[0], cell_size_y*2**queryShape[1]

        data = data_readin()
        ran_indices = np.random.randint(0,data.shape[1],Params.nQuery)
        ran_points = data[:,ran_indices]
        x_low = ran_points[0,:] - x_range/2
        x_high = ran_points[0,:] + x_range/2
        y_low = ran_points[1,:] - y_range/2
        y_high = ran_points[1,:] + y_range/2
        for i in range(int(Params.nQuery)):
            #query = [[max(x_low[i],x1),max(y_low[i],y1)],[min(x_high[i],x2),min(y_high[i],y2)]]
            query = [[x_low[i],y_low[i]],[x_high[i],y_high[i]]]
            querylist.append(np.array(query))
	all_queries.append(querylist)

    return all_queries

def query_gen_mtd(queryShape, seed, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """Generate query around a random data point"""
    
    logging.debug('Generating queries...')
    np.random.seed(seed)
    querylist = []
    x_range, y_range = Params.ONE_KM*Params.MTD/2**queryShape[0], Params.ONE_KM*Params.MTD/2**queryShape[1]
    data = data_readin()
    ran_indices = np.random.randint(0,data.shape[1],Params.nQuery)
    ran_points = data[:,ran_indices]
    x_low = ran_points[0,:] - x_range/2
    x_high = ran_points[0,:] + x_range/2
    y_low = ran_points[1,:] - y_range/2
    y_high = ran_points[1,:] + y_range/2
    for i in range(int(Params.nQuery)):
        query = [[max(x_low[i],x1),max(y_low[i],y1)],[min(x_high[i],x2),min(y_high[i],y2)]]
        querylist.append(np.array(query))

    return querylist


#def query_gen(queryShape, seed, random=False, x1=-116.915680, y1=37.000293, x2=-109.050173, y2=45.543541):
#    """Query generation. Each of which has at least one corner point in data populated areas.
#    Due to the distribution of the spatial data set, we do not want to generate many queries in 
#    the blank area. Hence the query generation function here is highly related to the Washington-NewMexico
#    data we use here. x1,y1,x2,y2 are the inner boundaries of the two states respectively.
#    You are encouraged to write your own query generation function depending on the dataset you use."""
#    
#    logging.debug('Generating queries...')
#    np.random.seed(seed)
#    querylist = []
#
#    cell_size_x = (-103.0020+124.8193)/2**(Params.queryUnit[0])
#    cell_size_y = (49.0025-31.3322)/2**(Params.queryUnit[1])
#    x_range, y_range = cell_size_x*2**queryShape[0], cell_size_y*2**queryShape[1]
#    
#    point_x = np.random.uniform(Params.LOW[0],x1,Params.nQuery/2)
#    point_y = np.random.uniform(y2,Params.HIGH[1],Params.nQuery/2)
#
#    x_low = point_x 
#    x_high = point_x + x_range
#    y_low = point_y - y_range
#    y_high = point_y 
#    for i in range(int(Params.nQuery/2)):
#        querylist.append(np.array([[x_low[i],y_low[i]],[x_high[i],y_high[i]]]))
#
#    point_x = np.random.uniform(x2,Params.HIGH[0],Params.nQuery/2)
#    point_y = np.random.uniform(Params.LOW[1],y1,Params.nQuery/2)
#    x_low = point_x - x_range
#    x_high = point_x
#    y_low = point_y 
#    y_high = point_y + y_range
#    for i in range(int(Params.nQuery/2)):
#        querylist.append(np.array([[x_low[i],y_low[i]],[x_high[i],y_high[i]]]))
#
#    return querylist

#   gowalla_CA: x1=-124.3041, y1=32.1714, x2=-114.0043, y2=41.9984
#   mcdonald: x1=-159.5862, y1=32.1714, x2=-67.2804, y2=64.8490
#28.071099,-125.081591    46.618575,-70.14995
#def query_gen(queryShape, seed, x1=-125.081591 , y1=28.071099, x2=-70.14995, y2=46.618575):
#    """Query generation. Generate fix-size query."""
#    
#    logging.debug('Generating queries...')
#    np.random.seed(seed)
#    querylist = []
#
#    cell_size_x = (x2-x1)/(2**Params.queryUnit[0])
#    cell_size_y = (y2-y1)/(2**Params.queryUnit[1])
#    x_range, y_range = cell_size_x*2**queryShape[0], cell_size_y*2**queryShape[1]
#    
#    point_x = np.random.uniform(x1,x2-x_range,Params.nQuery)
#    point_y = np.random.uniform(y1,y2-y_range,Params.nQuery)
#    x_low = point_x
#    x_high = point_x + x_range
#    y_low = point_y
#    y_high = point_y + y_range
#    for i in range(int(Params.nQuery)):
#        querylist.append(np.array([[x_low[i],y_low[i]],[x_high[i],y_high[i]]]))
#
#    return querylist


def test_quadtreeOpt(data, queryShape, all_queries):
    global method_list, exp_name
    exp_name = 'quadtreeOpt'
    method_list = ['quad-geo']
    # method_list = ['quad-baseline', 'quad-geo', 'quad-baseline-localness', 'quad-geo-localness']
    res_cube_abs = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    res_cube_rel = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    
    for j in range(len(seed_list)):
        queryList = all_queries[j];
        kexp = GKExp(data, queryList)
        p = Params(seed_list[j])
	
        for i in range(len(eps_list)):
            p.Eps = eps_list[i]	    
            for k in range(len(method_list)):
		
                if method_list[k] == 'quad-baseline':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_baseline(p)
                elif method_list[k] == 'quad-baseline-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_baseline_localness(p)
                elif method_list[k] == 'quad-geo':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_geo(p)
                elif method_list[k] == 'quad-geo-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_geo_localness(p)
                else:
                    logging.error('No such index structure!')
                    sys.exit(1)
    
    res_abs_summary = np.average(res_cube_abs,axis=1)
    res_rel_summary = np.average(res_cube_rel,axis=1)
    #np.savetxt(Params.resdir+exp_name+'_abs_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_abs_summary, fmt='%.4f\t')
    np.savetxt(Params.resdir+exp_name+'_rel_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_rel_summary, fmt='%.4f\t')

def test_kdTrees(data, queryShape, all_queries):
    global method_list, exp_name
    exp_name = 'kdTrees'
    method_list = ['kd-hybrid']
    #'kd-true','kd-true-localness','kd-standard','kd-hybrid','kd-standard-localness','kd-hybrid-localness'
    res_cube_abs = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    res_cube_rel = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    
    for j in range(len(seed_list)):
        queryList = all_queries[j]
        kexp = GKExp(data, queryList)
	p = Params(seed_list[j])
	
        for i in range(len(eps_list)):
	    p.Eps = eps_list[i]	    
            for k in range(len(method_list)):
		
                if method_list[k] == 'kd-pure':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_pure(p)
                elif method_list[k] == 'kd-true':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_true(p)
		elif method_list[k] == 'kd-true-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_true_localness(p)
                elif method_list[k] == 'kd-standard':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_standard(p)
		elif method_list[k] == 'kd-standard-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_standard_localness(p)
		elif method_list[k] == 'kd-hybrid':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_hybrid(p)    
                elif method_list[k] == 'kd-hybrid-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_hybrid_localness(p)
                else:
                    logging.error('No such index structure!')
                    sys.exit(1)
    
    res_abs_summary = np.average(res_cube_abs,axis=1)
    res_rel_summary = np.average(res_cube_rel,axis=1)
    #np.savetxt(Params.resdir+exp_name+'_abs_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_abs_summary, fmt='%.4f\t')
    np.savetxt(Params.resdir+exp_name+'_rel_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_rel_summary, fmt='%.4f\t')

def test_grids(data, queryShape, all_queries):
    global method_list, exp_name
    exp_name = 'grids'
    method_list = ['grid-uniform','grid-adaptive']
    #'grid-pure','grid-uniform','grid-adaptive','grid-adaptive-localness'
    res_cube_abs = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    res_cube_rel = np.zeros((len(eps_list),len(seed_list),len(method_list)))
	    
    for j in range(len(seed_list)):
	queryList = all_queries[j]
	kexp = GKExp(data, queryList)
	p = Params(seed_list[j])	
        
        for i in range(len(eps_list)):
            p.Eps = eps_list[i]
            for k in range(len(method_list)):		
                if method_list[k] == 'grid-pure':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Grid_pure(p)
		elif method_list[k] == 'grid-uniform':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Grid_uniform(p)
                elif method_list[k] == 'grid-adaptive':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Grid_adaptive(p)
                elif method_list[k] == 'grid-adaptive-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Grid_adaptive_localness(p)
                else:
                    logging.error('No such index structure!')
                    sys.exit(1)
    
    res_abs_summary = np.average(res_cube_abs,axis=1)
    res_rel_summary = np.average(res_cube_rel,axis=1)
    #np.savetxt(Params.resdir+exp_name+'_abs_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_abs_summary, fmt='%.4f\t')
    np.savetxt(Params.resdir+exp_name+'_rel_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_rel_summary, fmt='%.4f\t')
    
def test_htrees(data, queryShape, all_queries):
    global method_list, exp_name
    exp_name = 'htrees'
    method_list = ['ht-standard']
#    method_list = ['ht-standard','ht-composite']
    #'ht-pure','ht-true','ht-standard','ht-composite','ht-hybrid','ht-hybrid-skew','ht-composite-localness','ht-hybrid-localness'
    res_cube_abs = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    res_cube_rel = np.zeros((len(eps_list),len(seed_list),len(method_list)))
    
    for j in range(len(seed_list)):
        queryList = all_queries[j]
        kexp = GKExp(data, queryList)
	p = Params(seed_list[j])
        
	for i in range(len(eps_list)):
	    p.Eps = eps_list[i]
            for k in range(len(method_list)):

                if method_list[k] == 'ht-pure':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_pure(p)
                elif method_list[k] == 'ht-true':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_true(p)
                elif method_list[k] == 'ht-standard':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_standard(p)
                elif method_list[k] == 'ht-composite':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_composite(p)
                elif method_list[k] == 'ht-composite-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_composite_localness(p)
                elif method_list[k] == 'ht-hybrid':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid(p)
                elif method_list[k] == 'ht-standard-skew':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_standard_skew(p)    
		elif method_list[k] == 'ht-hybrid-skew': 
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid_skew(p)     
		elif method_list[k] == 'ht-standard-adaptive':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_standard_adaptive(p)     
                elif method_list[k] == 'ht-hybrid-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid_localness(p)
                else:
                    logging.error('No such index structure!')
                    sys.exit(1)
    res_abs_summary = np.average(res_cube_abs,axis=1)
    res_rel_summary = np.average(res_cube_rel,axis=1)
    #np.savetxt(Params.resdir+exp_name+'_abs_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_abs_summary, fmt='%.4f\t')
    np.savetxt(Params.resdir+exp_name+'_rel_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`, res_rel_summary, fmt='%.4f\t')

def test_height(data, queryShape, all_queries):
    heightList = [5,6,7,8,9]
    method_list = ['grid-uniform','kd-hybrid','ht-hybrid','kd-hybrid-localness'];
    #['grid-uniform','quad-geo','kd-hybrid','ht-composite','ht-hybrid','ht-hybrid-skew','quad-geo-localness','kd-hybrid-localness','ht-composite-localness','ht-hybrid-localness'];
    res_cube_abs = np.zeros((len(heightList),len(seed_list),len(method_list)))
    res_cube_rel = np.zeros((len(heightList),len(seed_list),len(method_list)))

    for j in range(len(seed_list)):
        queryList = all_queries[j]
        kexp = GKExp(data, queryList)
        p = Params(seed_list[j])
        for i in range(len(heightList)):
            Params.maxHeight = heightList[i]
            Params.switchLevel = heightList[i]/2
            Params.maxSplit = [heightList[i],heightList[i]]
            Params.switchLevelHTree = [heightList[i]/2,heightList[i]/2]
	    Params.partitions = [2**heightList[i],2**heightList[i]]
            Params.switchPartitions = [2**heightList[i]/2,2**heightList[i]/2]
            for k in range(len(method_list)):
                if method_list[k] == 'grid-uniform':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Grid_uniform(p)
		elif method_list[k] == 'quad-geo':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_geo(p)
                elif method_list[k] == 'ht-standard':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_standard(p)
                elif method_list[k] == 'kd-standard':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_standard(p)
		elif method_list[k] == 'kd-hybrid':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_hybrid(p)
                elif method_list[k] == 'ht-composite':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_composite(p)                
                elif method_list[k] == 'ht-hybrid':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid(p)
		elif method_list[k] == 'ht-standard-skew':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_standard_skew(p)   
                elif method_list[k] == 'ht-hybrid-skew':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid_skew(p)
		elif method_list[k] == 'quad-geo-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Quad_geo_localness(p)
		elif method_list[k] == 'kd-hybrid-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_Kd_hybrid_localness(p)               
                elif method_list[k] == 'ht-composite-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_composite_localness(p)
		elif method_list[k] == 'ht-hybrid-localness':
                    res_cube_abs[i,j,k], res_cube_rel[i,j,k] = kexp.run_HT_hybrid_localness(p)  
                else:
                    logging.error('No such index structure!')
                    sys.exit(1)
    res_abs_summary = np.average(res_cube_abs,axis=1)
    res_rel_summary = np.average(res_cube_rel,axis=1)
    
    for str in ['abs','rel']:
        summary = eval('res_'+str+'_summary')
        outName = Params.resdir+'height_'+`int(queryShape[0]*10)`+'_'+`int(queryShape[1]*10)`+str
        outFile = open(outName,'w')
        outFile.write('#; ' + '\t'.join(method_list) + '\n')
        for i in range(len(heightList)):
            outFile.write(`heightList[i]`)
            for j in range(len(method_list)):
                outFile.write('\t'+`summary[i,j]`)
            outFile.write('\n')
        outFile.close()

    
def createGnuData():
    """
    Post-processing output files to generate Gnuplot-friendly outcomes
    """
    line = 0
    for eps in eps_list:
        for type in ['abs','rel']:
            out = open(Params.resdir + exp_name + '_eps' + str(int(eps*10)) + '_' + type, 'w')
            out.write('#; ' + '\t'.join(method_list) + '\n')
            q_num = 1
            for queryShape in shape_list:
                fileName = Params.resdir + exp_name + '_' + type + '_' + `int(queryShape[0]*10)` + '_' + `int(queryShape[1]*10)`
                try:
                    thisfile = open(fileName, 'r')
                except:
                    sys.exit('no input result file!')
                out.write(str(q_num) + '\t' + thisfile.readlines()[line])
                thisfile.close()
                q_num += 1
            out.close()
        line += 1
    
    for type in ['abs','rel']:
        for queryShape in shape_list:
            fileName = Params.resdir + exp_name + '_' + type + '_' + `int(queryShape[0]*10)` + '_' + `int(queryShape[1]*10)`
            os.remove(fileName)
            
            
            
if __name__ == '__main__':
#    sys.path.append('Spatial_DP_v1')
#    sys.path.append('geocrowd')
    
    logging.basicConfig(level=logging.DEBUG, filename='debug.log')
    logging.info(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "  START") 
    
    data = data_readin()
    
    for q_shape in shape_list:
	all_queries = query_gen(q_shape,Params.x_min,Params.y_min,Params.x_max,Params.y_max)
        
        test_htrees(data, q_shape, all_queries)
        test_grids(data, q_shape, all_queries)
        test_kdTrees(data, q_shape, all_queries)  
        test_quadtreeOpt(data, q_shape, all_queries)
        
    
    filenames = ["HT_standard", "Grid_uniform", "Grid_adaptive", "Kd_hybrid", "Quad_geo"]
    file_rel = open("../log/dump_rel.log", "a")
    file_abs = open("../log/dump_abs.log", "a")
    files_contents_rel = []
    files_contents_abs = []
    for filename in filenames:
        content_rel = open("../log/" + filename + "_rel.log", "r").readlines()
        files_contents_rel.append(content_rel)
        content_abs = open("../log/" + filename + "_abs.log", "r").readlines()
        files_contents_abs.append(content_abs)
    
    for i in range(len(files_contents_rel[0])):
        line_rel = ""
        line_abs = ""
        for j in range(len(files_contents_rel)):
            line_rel = line_rel + str(float("%.3f" % float(files_contents_rel[j][i].strip()))) + "\t"
            line_abs = line_abs + str(float("%.3f" % float(files_contents_abs[j][i].strip()))) + "\t"
        file_rel.write(line_rel + "\n")
        file_abs.write(line_abs + "\n")
    file_rel.close()
    file_abs.close()
    
    # Experiment 5: test the impact of max tree height ###
#    for q_shape in shape_list:
#        test_height(data, q_shape)
         
        
    ######### run with multi-core support, number of cores should be at least len(shape_list) #########    
#    pool = mult.Pool(processes=len(shape_list))
#    pool.map(test_htrees,shape_list)
#    pool.map(test_grids,shape_list)
#    pool.map(test_kdTrees,shape_list)
    #    pool.map(test_quadtreeOpt,shape_list)
    
    #    pool.map(test_height,shape_list)
#    pool.close()
#    pool.join()

    ######### 
    
    
    logging.info(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "  END")