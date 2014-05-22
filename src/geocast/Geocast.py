"""
Contains all methods used in geocast selection algorithm
"""
import math 

import logging
import numpy as np
import heapq
from Params import Params
import copy
from Utils import utility, rect_intersect, is_range_overlap, rect_area, is_intersect, performed_tasks, acc_rate, is_intersect_segment, distance, distance_to_rect, rect_center, rect_vertex_set, is_rect_cover, is_rect_cover_rect, distance_point
from sets import Set
from minball.smallestenclosingcircle import make_circle, make_circle_one_point
from Geocrowd import rect_query_count, rect_query_points, hops_expansion

def neighbors_dir(node, dir, neighbor):
    """
    Find all neighbors (leaf nodes) of a node based on a direction
    """
    neighbors = []
    if neighbor.n_isLeaf == True:
        neighbors.append(neighbor)
    else:
        if dir == "N":
            range = [node.n_box[0,0], node.n_box[1,0]]
            for n in neighbor.children:
                m = n.children[len(n.children)-1]   # a north neighbor
                range_neighbor = [m.n_box[0,0], m.n_box[1,0]]
                if is_range_overlap(range, range_neighbor):
                    neighbors.append(m)
        elif dir == "S":
            range = [node.n_box[0,0], node.n_box[1,0]]
            for n in neighbor.children:
                m = n.children[0]   # a south neighbor
                range_neighbor = [m.n_box[0,0], m.n_box[1,0]]
                if is_range_overlap(range, range_neighbor):
                    neighbors.append(m)
        elif dir == "W":
            range = [node.n_box[0,1], node.n_box[1,1]]
            for m in neighbor.children[len(neighbor.children)-1].children:
                range_neighbor = [m.n_box[0,1], m.n_box[1,1]]
                if is_range_overlap(range, range_neighbor):
                    neighbors.append(m)
        elif dir == "E":
            range = [node.n_box[0,1], node.n_box[1,1]]
            for m in neighbor.children[0].children:
                range_neighbor = [m.n_box[0,1], m.n_box[1,1]]
                if is_range_overlap(range, range_neighbor):
                    neighbors.append(m)
    return neighbors

def log_cell_info(node):
    """
    Find the node's info and its statistic
    """
    parent = node.parent
    i, j = node.index, parent.index
    workers = 0
    if node.n_data is not None:
        workers = node.n_data.shape[1]
    if node.n_depth == 4:  # 2nd level
        k, l = parent.parent.index, parent.parent.parent.index
        return k, l, i, j, workers, float("%.1f" % node.n_count)
    
    return i, j, -1, -1, workers, float("%.1f" % node.n_count)
        
        
def neighbors(tree, node):
    """
    Find all neighbors (leaf nodes) of a node
    """
    q_neighbors = []
    gran_1st = len(tree.root.children)
    parent = node.parent
    i, j = parent.index, node.index
    grand_parent = parent.parent
    if node.n_depth == 4:  # 2nd level
        grand_parent_parent = grand_parent.parent
        k, l = grand_parent_parent.index, grand_parent.index
        gran_2nd = len(parent.children)
        if j==0:
            q_neighbors.append(parent.children[j+1])
            if i==0: # [0,0]
                q_neighbors.append(grand_parent.children[i+1].children[j])
                if l >= 1:
                    q_neighbors.extend(neighbors_dir(node, "N", grand_parent_parent.children[l-1]))
                if k >= 1:
                    q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[k-1].children[l]))
            elif i==gran_2nd-1: # [2,0]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                if l >= 1:
                    q_neighbors.extend(neighbors_dir(node, "N", grand_parent_parent.children[l-1]))
                if k <= gran_1st - 2:
                    q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[k+1].children[l]))
            else: # [1,0]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                q_neighbors.append(grand_parent.children[i+1].children[j])
                if l >= 1:
                    q_neighbors.extend(neighbors_dir(node, "N", grand_parent_parent.children[l-1]))
                    
        elif j==gran_2nd-1:
            q_neighbors.append(parent.children[j-1])
            if i==0: # [0,2]'
                q_neighbors.append(grand_parent.children[i+1].children[j])
                if l <= gran_1st-2:
                    q_neighbors.extend(neighbors_dir(node, "S", grand_parent_parent.children[l+1]))
                if k >= 1:
                    q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[k-1].children[l]))
            elif i==gran_2nd-1: # [2,2]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                if l <= gran_1st-2:
                    q_neighbors.extend(neighbors_dir(node, "S", grand_parent_parent.children[l+1]))
                if k <= gran_1st-2:
                    q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[k+1].children[l]))
            else: # [1,2]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                q_neighbors.append(grand_parent.children[i+1].children[j])
                if l <= gran_1st-2:
                    q_neighbors.extend(neighbors_dir(node, "S", grand_parent_parent.children[l+1]))
        else: 
            q_neighbors.append(parent.children[j-1])
            q_neighbors.append(parent.children[j+1])
            
            if i==0: # [0,1]
                q_neighbors.append(grand_parent.children[i+1].children[j])
                if k >= 1:
                    q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[k-1].children[l]))
            elif i==gran_2nd-1: # [2,1]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                if k <= gran_1st-2:
                    q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[k+1].children[l]))
            else: # [1,1]
                q_neighbors.append(grand_parent.children[i-1].children[j])
                q_neighbors.append(grand_parent.children[i+1].children[j])
    elif node.n_depth == 2: # 1st level    
        if j==0:
            q_neighbors.extend(neighbors_dir(node, "S", parent.children[j+1])) 
            if i==0: # [0,0]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
            elif i==gran_1st-1: # [2,0]
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
            else: # [1,0]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
        elif j==gran_1st-1:
            q_neighbors.extend(neighbors_dir(node, "N", parent.children[j-1]))
            if i==0: # [0,2]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
            elif i==gran_1st-1: # [2,2]
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
            else: # [1,2]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
        else:
            q_neighbors.extend(neighbors_dir(node, "N", parent.children[j-1]))
            q_neighbors.extend(neighbors_dir(node, "S", parent.children[j+1]))
            if i==0: # [0,1]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
            elif i==gran_1st-1: # [2,1]
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
            else: # [1,1]
                q_neighbors.extend(neighbors_dir(node, "E", tree.root.children[i+1].children[j]))
                q_neighbors.extend(neighbors_dir(node, "W", tree.root.children[i-1].children[j]))
                        
    return q_neighbors

def is_level_expansion(parent, childrens):
    granularity = len(parent.children)
    
    percent = len(childrens)/granularity**2
    if len(childrens) >= 2 and percent > Params.GAMMA:    # apply level expansion
        return percent
    return -1
#
## check if the cell is covered by a node
#    id_xmin_1st = int(math.floor((cell[0][0] - Params.LOW[0])/x_width_1st))
#    id_ymin_1st = int(math.floor((cell[0][1] - Params.LOW[1])/y_width_1st))
#    id_xmax_1st = int(math.ceil((cell[1][0] - Params.LOW[0])/x_width_1st))
#    id_ymax_1st = int(math.ceil((cell[1][1] - Params.LOW[1])/y_width_1st))
#    if id_xmax_1st - id_xmin_1st == 1 and id_ymax_1st - id_ymin_1st == 1:   # contain one 2nd grid only
#        # one grid
#        # compute 2nd granurality
#        granularity_2nd = tree.getLeafGranularity(cell)
#        if granularity_2nd >= 2:
#            x_width_2nd = x_width_1st/granularity_2nd
#            y_width_2nd = y_width_1st/granularity_2nd
#            x_min = id_xmin_1st*x_width_1st + Params.LOW[0]
#            y_min = id_ymin_1st*y_width_1st + Params.LOW[1]
#            x_max = id_xmax_1st*x_width_1st + Params.LOW[0]
#            y_max = id_ymax_1st*y_width_1st + Params.LOW[1]
#            id_xmin_2nd = int(math.floor((cell[0][0] - x_min)/x_width_2nd))
#            id_ymin_2nd = int(math.floor((cell[0][1] - y_min)/y_width_2nd))
#            id_xmax_2nd = int(math.ceil((cell[1][0] - x_min)/x_width_2nd))
#            id_ymax_2nd = int(math.ceil((cell[1][1] - y_min)/y_width_2nd))
#
#            cell_aligned = [[id_xmin_2nd*x_width_2nd + x_min, id_ymin_2nd*y_width_2nd + y_min],[id_xmax_2nd*x_width_2nd + x_min, id_ymax_2nd*y_width_2nd + y_min]]
#            cell_list_aligned.append(np.array(query_aligned))
#
#    #                             find expanded cell
#            N = granularity_2nd
#            n_x = id_xmax_2nd - id_xmin_2nd
#            n_y = id_ymax_2nd - id_ymin_2nd
#            percent = (n_x*n_y + 0.0)/N**2
#            if n_x*n_y >= 2 and percent > gamma:    # apply level expansion
#                cell_expanded = [[x_min, y_min],[x_max, y_max]]
#                cell_list_expanded_aligned.append(np.array(query_aligned))
#                cell_list_expanded_percent.append(percent)
#                cell_list_expanded.append(np.array(query_expanded))

def parent_index(tree, leaf_node):
    """
    get parent index (i,j) of a leaf node
    """
    granularity_1st = len(tree.root.children)
    x_width_1st, y_width_1st = (Params.HIGH[0] - Params.LOW[0])/granularity_1st, (Params.HIGH[1] - Params.LOW[1])/granularity_1st
   
    i, j = int(((leaf_node.n_box[0][0] + leaf_node.n_box[1][0])/2 - leaf_node.n_box[0][0])/x_width_1st), int(((leaf_node.n_box[0][1] + leaf_node.n_box[1][1])/2 - leaf_node.n_box[0][1])/y_width_1st)
    return (i,j)

def leaf_granularity(tree, index):
    return len(tree.root.children[index[0]].children[index[1]].children)
        
#def le_heuristic():
#	d = dict()
#	d_anchor = dict()
#        # LE heuristic
#        index = parent_index(tree, cell)
#        leafGran = leaf_granularity(tree, index)
#        if Params.LEVEL_EXPANSION and leafGran >= 2:
#            parent = tree.root.children[index[0]].children[index[1]]
#            if d.has_key(index):
#                d[index].add(cell)
##                print index, len(d[index])
#            else:
#                u_parent, dist = utility(parent, Params.MTD, L, UNDER_EST)
#                d[index] = set([cell])
#            
#            # check the parent node
#            percent = is_level_expansion(parent, d[index])
#            if percent > 0:     # LE is applied
#                # update utility
#		if not d_anchor.has_key(index): # first time
#		    d_anchor[index] = percent
#                    for child in d[index]:
#                        u_child, dist = utility(child, Params.MTD, L, UNDER_EST)
#                        u = max(0, (u - u_child)/(1-u_child))
##			print u_child, u
##                        print u_child, u
#                else:	# update
#                    last_percent = d_anchor[index]
#		    last_u_parent, dist = utility(parent, Params.MTD, L, UNDER_EST, last_percent)
#                    u = max(0, (u - last_u_parent)/(1-last_u_parent))
##		    print last_percent, u
#                
#                u_parent, dist = utility(parent, Params.MTD, L, UNDER_EST, percent)
#                u = 1 - (1-u)*(1-u_parent)
#	    else:
#		# update utility
#		u_c = c[1][0]
#		u = 1 - (1-u)*(1-u_c)
##		print "U: " + str(u)
#		if len(q) > 200:
#		    print "Cannot resolve this task"
#		    break
#		if u >= Params.U:
#		    break
    
def select_partial_cell(t, cell, u_prev):
    """
    Partial selectio on the last grid cell
    
    @param L: task location
    @param cell
    @u_prev: previous utility
    """
    sub_cell = copy.copy(cell)
    [[min_x, min_y],[max_x, max_y]] = cell.n_box
    dist = distance_to_rect(t[0],t[1],cell.n_box)

    # find a sub-region (i.e., sub-cell) whose utility is just enough
    ar = acc_rate(Params.MTD, dist)
    u_need = (Params.U-u_prev)/(1-u_prev)
    sub_cell.n_count = math.floor(math.log(1-u_need, 1-ar))
    percentile = sub_cell.n_count/cell.n_count    # area fraction
        
    # compute the partial region
    if is_rect_cover(cell.n_box, t):    # there is only one cell in geocast query
        # find a sub rectangle that is similar to a larger one, but smaller by a percentile
        x_len = math.sqrt(percentile)*(max_x - min_x)
        y_len = math.sqrt(percentile)*(max_y - min_y)
        small_rect = np.array([[min_x+x_len/2, min_y+y_len/2], [max_x-x_len/2, max_y-y_len/2]])
    
        if (t[0] < small_rect[0][0]):
            if (t[1] < small_rect[0][1]):   # [0,0]
                c_x, c_y = small_rect[0][0], small_rect[0][1]
            elif (small_rect[0][1] <= t[1] <= small_rect[1][1]):    # [0,1]
                c_x, c_y = small_rect[0][0], t[1]
            elif (small_rect[1][1] < t[1]): # [0,2]
                c_x, c_y = small_rect[0][0], small_rect[1][1]
        elif (small_rect[0][0] <= t[0] <= small_rect[1][0]):  
            if (t[1] < small_rect[0][1]):   # [1,0]  
                c_x, c_y = t[0], small_rect[0][1]
            elif (small_rect[0][1] <= t[1] <= small_rect[1][1]):    # [1,1]
                c_x, c_y = t[0], t[1]
            elif (small_rect[1][1] < t[1]): # [1,2]
                c_x, c_y = t[0], small_rect[1][1]
        elif (small_rect[1][0] < t[0]):
            if (t[1] < small_rect[0][1]):   # [2,0]  
                c_x, c_y = small_rect[1][0], small_rect[0][1]
            elif (small_rect[0][1] <= t[1] <= small_rect[1][1]):    # [2,1]
                c_x, c_y = small_rect[1][0], t[1]
            elif (small_rect[1][1] < t[1]): # [2,2]
                c_x, c_y = small_rect[1][0], small_rect[1][1]
        sub_cell.n_box = np.array([[c_x-x_len/2, c_y-y_len/2], [c_x+x_len/2, c_y+y_len/2]])
    else:   # the task is outside the cell (i.e., the geocast query includes multiple cells)
        [mid_x, mid_y] = rect_center(cell.n_box)
        mid_neighbor = rect_center(cell.neighbor.n_box)
        if is_intersect_segment(mid_neighbor, [mid_x, mid_y], [min_x, min_y], [max_x, min_y]):  # bottom
            sub_cell.n_box = np.array([[min_x, min_y],[max_x, min_y+(max_y - min_y)*percentile]])
        elif is_intersect_segment(mid_neighbor, [mid_x, mid_y], [max_x, min_y], [max_x, max_y]):    # right
            sub_cell.n_box = np.array([[max_x-(max_x - min_x)*percentile, min_y], [max_x, max_y]])
        elif is_intersect_segment(mid_neighbor, [mid_x, mid_y], [min_x, max_y], [max_x, max_y]):    # ceil
            sub_cell.n_box = np.array([[min_x, max_y-(max_y - min_y)*percentile], [max_x, max_y]])
        elif is_intersect_segment(mid_neighbor, [mid_x, mid_y], [min_x, min_y], [min_x, max_y]):  # left
            sub_cell.n_box = np.array([[min_x, min_y], [min_x+(max_x - min_x)*percentile, max_y]])
    return sub_cell

def new_neighbors(tree, centers, cell, MTD_RECT):
    """
    Find neighbors of a cell that are not added into geocast query yet
    
    @param tree : WorkerPSD
    @param q : geocast query
    @param cell : a cell
    @param MTD_RECT : MTD
    """
    q_valid_neighbors = []
    for n in neighbors(tree, cell):
        flag = True
        c = rect_center(n.n_box)
        if (c[0], c[1]) in centers:
            flag = False
            break
        if flag:
            _q = rect_intersect(MTD_RECT, np.array(n.n_box))    # MTD filter
            if _q is not None:
                # update n_box and n_count of neighbor node n
                n_copy = copy.copy(n)   # this function may cause error
                n_copy.n_count = math.floor(n_copy.n_count * rect_area(_q)/rect_area(n.n_box))
                n_copy.n_box = _q
                n_copy.neighbor = cell
                q_valid_neighbors.append(n_copy)
    return q_valid_neighbors

                    
def cost_function(utility, compactness, eps, alpha = Params.ALPHA):
    """
    Hybrid cost function
    """
    return (1-eps)*utility*(1-alpha) + eps*compactness*alpha

def geocast(tree, t, eps):
    """
    Given WorkerPSD and a task, find the region around the task that covers enough workers so that
    if we geocast the task to the region (i.e., geocast region/query), the task would be performed 
    with high probability
    
    @param tree : WorkerPSD
    @param L : task location
    """
    centers = Set([])
    MTD_RECT = np.array([[t[0]-Params.ONE_KM*Params.MTD,t[1]-Params.ONE_KM*Params.MTD],[t[0]+Params.ONE_KM*Params.MTD, t[1]+Params.ONE_KM*Params.MTD]])
    
    # initialize q as the cell that covers L
    q_init = tree.leafCover(t)  # get the leaf node that covers L
    if q_init == None:
        print t
        return None, None
    q_init.neighbor = q_init
    q_init.n_count = math.floor(q_init.n_count)
    u_q, dist = utility(q_init, Params.MTD, t)
    f_q = u_q	# cost function
    Q = [(f_q, [u_q, q_init, 2/math.pi, -dist])]
    q = []  # final cell {a set of cells}
    q_log = []	# [(lambda, cell[i]'s utility, utility, compactness, distance, cx, cy, r),..] of the cells in q
    u = 0   # current utility
    A = 0   # current geocast area
    corner_points = Set([])    # all corner points of geocast cell, used as an input to compute minimum bounding circle
    
    while (True):
        if len(Q) == 0: # Empty queue, no more expansion
            break
        
        subQ = list(filter(lambda p: p[1][0] > 0, Q))
        if len(subQ) > 0:
            c = max(subQ, key = lambda	p : p[0])    # (lambda, [utility, cell, compactness, distance, cx, cy, r])
            print c[0]
        else:
            c = max(Q, key = lambda	p : p[0])    # (lambda, [utility, cell, compactness, distance, cx, cy, r])
        Q.remove(c)
        
        # update cell
        cell = c[1][1]
        cen = rect_center(cell.n_box)
        if (cen[0], cen[1]) in centers:
            continue
            
	# update utility
	u_c = c[1][0]
        u_curr = u
	if u_c >= 0:
	    u = 1 - (1-u)*(1-u_c)
	elif Params.NEGATIVE_CELL:
	    u = (u - u_c)/(1-u_c)

	q_log.append((float("%.1f" % c[0]), float("%.3f" % c[1][0]), float("%.3f" % u), float("%.3f" % c[1][2]), float("%.3f" % c[1][3])))
        
        # if the new cell is too much (over-provision), then we may want to take part of the cell
	if u >= Params.U:
            if Params.PARTIAL_CELL_SELECTION and u >= Params.U + 0.05:
                sub_cell = select_partial_cell(t, cell, u_curr)
                q.append(sub_cell)
            else:
                q.append(cell)
            break
            
        q.append(cell)   
        centers.add((cen[0], cen[1]))
            
        corner_points = corner_points | rect_vertex_set(cell.n_box)
        A += rect_area(cell.n_box)

        # filter out the neighbors that already in q, then cut off region that is outside MTD
        q_valid_neighbors = new_neighbors(tree, centers, cell, MTD_RECT)

	# update Q
        for nb in q_valid_neighbors:
            u_q, dist = utility(nb, Params.MTD, t)
            Q.append((u_q, [u_q, nb, 0, -dist, 0.0, 0.0, 0.0]))
            
	for i in range(len(Q)):
	    # update key
            points = list(corner_points | rect_vertex_set(Q[i][1][1].n_box))
            x = make_circle(points)
            if x == None:
                print x
                return None, None
            cx,cy,r = x
            compactness = min(1.0, (A + rect_area(Q[i][1][1].n_box))/(math.pi*(r/Params.ONE_KM)**2))  # the area of the shape divided by the area of the minimum bounding circle                
            if Params.COST_FUNCTION == "hybrid":
                f_q = cost_function(Q[i][1][0], compactness, eps)
            elif Params.COST_FUNCTION == "compactness":
                f_q = compactness
            elif Params.COST_FUNCTION == "distance":
                f_q = Q[i][1][3]
            elif Params.COST_FUNCTION == "utility":
                f_q = Q[i][1][0]
                
            Q[i] = (f_q, [Q[i][1][0], Q[i][1][1], compactness, Q[i][1][3]])
        
    return q, q_log


def post_geocast(t, q, q_log):
    """
    Compute actual utility & average travel cost in simulation
    """
    if q == None:
        return [None for k in range(6)]
    cells = []
    no_workers = 0      
    workers = np.zeros(shape=(2,0)) # worker locations
    for i in range(len(q)):
        cell = q[i]
        cells.append([cell, log_cell_info(cell), q_log[i]])
        if cell.n_data != None:
            if Params.PARTIAL_CELL_SELECTION and i == len(q) - 1:
                _workers = rect_query_points(cell.n_data, cell.n_box)
            else:
                _workers = cell.n_data
            no_workers += _workers.shape[1] 
            workers = np.concatenate([workers, _workers], axis=1)

    hops_count, coverage, hops_count2 = 0, 0, 0
    if workers.shape[1] > 0:
        hops_count, coverage, hops_count2 = hops_expansion(t, workers.transpose(), Params.NETWORK_DIAMETER)
    
    return no_workers, workers, cells, hops_count, coverage, hops_count2

def simulation_geocast(t, q, FCFS):    
    for i in range(len(q)):
	cell = q[i]
        if Params.PARTIAL_CELL_SELECTION and i == len(q) - 1 and cell.n_data != None:
            _workers = rect_query_points(cell.n_data, cell.n_box)
        else:
            _workers = cell.n_data
        performed, worker, dist  = performed_tasks(_workers, Params.MTD, t, FCFS)
        if performed: # the task is performed
            return True, worker, dist
    return False, None, dist