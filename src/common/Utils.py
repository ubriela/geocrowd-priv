import math
import random
import bisect
import itertools as it
# noinspection PyDeprecation
from sets import Set

import numpy as np

from Params import Params

# from collections import Counter
# from operator import itemgetter


def distance(lat1, lon1, lat2, lon2):
    """
    Distance between two geographical location
    """
    R = 6371  # km
    dLat = math.radians(abs(lat2 - lat1))
    dLon = math.radians(abs(lon2 - lon1))
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.sin(dLon / 2) * math.sin(dLon / 2) * math.cos(lat1) * math.cos(
        lat2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c
    return d


def distance_point(lat1, lon1, lat2, lon2):
    """
    distance between two point
    """
    return math.sqrt((lat2 - lat1) ** 2 + (lon2 - lon1) ** 2)


_distance = distance


def distance_to_rect(lat, lon, rect):
    d1 = _distance(lat, lon, rect[0][0], rect[0][1])
    d2 = _distance(lat, lon, rect[1][0], rect[0][1])
    d3 = _distance(lat, lon, rect[1][0], rect[1][1])
    d4 = _distance(lat, lon, rect[0][0], rect[1][1])
    return sum([d1, d2, d3, d4]) / 4


"""
Step function/starcase function
"""


def _step_function(max_distance):
    steps = Params.STEPS
    max_y = Params.MAR

    step_x = (max_distance + 0.0) / steps
    x = [step_x]
    for _i in range(int(max_distance / step_x) - 1):
        x.append(x[_i] + step_x)

    step_y = max_y / len(x)
    y = [max_y]

    for _i in range(len(x)):
        y.append(abs(y[_i] - step_y))
    return x, y


def acc_rate(max_distance, dist):
    if Params.AR_FUNCTION == "zipf":
        k = max(1, int(dist * Params.ZIPF_STEPS / max_distance))  # rank
        return _zipf_pmf(k, Params.s, Params.ZIPF_STEPS) * Params.MAR
    elif Params.AR_FUNCTION == "step":
        x, y = _step_function(max_distance)
        pos = np.searchsorted(x, dist)
        return max(0, y[pos])
    elif Params.AR_FUNCTION == "linear":
        return max(0, (1 - dist / max_distance) * Params.MAR)
    elif Params.AR_FUNCTION == "const":
        return Params.MAR


_acc_rate = acc_rate


"""
acceptance rate
"""
def is_performed(ar):
    """
    Simulate whether a task is performed or not 0<=acc_rate<=1
    """
    ran = random.random()
    if ran <= ar:
        return True
    return False


# This function may be slow
def performed_tasks(workers, max_dist, t, FCFS):
    """
    find the performed task, given the workers being geocasted and their acceptance rates
    
    @param locs : a list of worker locations
    @param max_dist : MTD, accepatance rate is zero at MTD
    @param t : task location    
    @param FCFS : first-come-first-serve mode
    """
    if workers is None:  # double check
        return False, None, None
    workers_copy = workers.transpose()
    if FCFS:
        ar_weights = [_acc_rate(max_dist, distance(t[0], t[1], w[0], w[1])) for w in workers_copy]
        while len(workers_copy) > 0:
            idx = min(list(it.islice(_wrg(ar_weights), 1))[0], len(ar_weights) - 1)
            ar = ar_weights[idx]
            if is_performed(ar):
                worker = workers_copy[idx]
                return True, worker, distance(t[0], t[1], worker[0], worker[1])
            del ar_weights[idx]
            workers_copy = np.delete(workers_copy, idx, 0)
    else:
        workers_copy = sorted(workers_copy, key=lambda loc: distance(loc[0], loc[1], t[0], t[1]))
        for worker in workers_copy:
            dist = distance(t[0], t[1], worker[0], worker[1])
            ar = _acc_rate(max_dist, dist)
            if is_performed(ar):
                return True, worker, dist
    return False, None, None


def performed_tasks_naive(locs, max_dist, t, FCFS, seed):
    """
    compute performed task, given the number of workers being geocasted and their acceptance rate
    
    @param locs : a list of worker locations
    @param max_dist : MTD, acceptance rate is zero at MTD
    @param t : task location
    """
    performed, dist = False, 0
    if FCFS:
        purmuted_locs = range(len(locs))
        random.shuffle(purmuted_locs)
        for i in purmuted_locs:
            performed, dist = performed_task(locs[i], max_dist, t)
            if performed:
                return True, dist
    else:
        locs = sorted(locs, key=lambda loc: distance(loc[0], loc[1], t[0], t[1]))
        for loc in locs:
            performed, dist = performed_task(loc, max_dist, t)
            if performed:
                return True, dist

    return False, None


def performed_task(loc, max_dist, t):
    """
    Simulate whether a task is perform given location of the worker
    
    @param loc : worker location
    @param max_dist : MTD, accepatance rate is zero at MTD
    @param t : task location
    """
    dist = distance(t[0], t[1], loc[0], loc[1])
    ar = _acc_rate(max_dist, dist)
    if is_performed(ar):
        return True, dist
    return False, None


def utility(node, max_dist, t):
    """
    Compute utility of a cell with respect to location of a task
    
    @param node : node
    @param max_dist : MTD, utility = 0 at  MTD
    @param t : location of the task
    """
    dist = distance_to_rect(t[0], t[1], node.n_box)
    ar = _acc_rate(max_dist, dist)
    # print ar, node.n_count
    return np.sign(node.n_count) * (1 - (1 - ar) ** abs(node.n_count)), dist


def utility_naive(query, w, max_dist):
    dist = math.sqrt((query[1][0] - query[0][0]) ** 2 + (query[1][1] - query[0][1]) ** 2) / 2
    ar = Params.MAR
    return np.sign(w) * (1 - (1 - ar) ** abs(w)), dist


def is_intersect(rec, query):
    bool_m1 = query[0, :] >= rec[1, :]
    bool_m2 = query[1, :] <= rec[0, :]
    bool_m = np.logical_or(bool_m1, bool_m2)
    if np.any(bool_m):
        return False
    else:
        return True


__is_intersect = is_intersect


def rect_intersect(rec, query):
    if __is_intersect(rec, query):
        min_x = max(rec[0][0], query[0, 0])
        min_y = max(rec[0][1], query[0, 1])

        max_x = min(rec[1][0], query[1, 0])
        max_y = min(rec[1][1], query[1, 1])
        return np.array([[min_x, min_y], [max_x, max_y]])
    else:
        return None


def rect_area(rect):
    """
    Geographical coordinates
    """
    return distance(rect[0][0], rect[0][1], rect[0][0], rect[1][1]) * distance(rect[0][0], rect[0][1], rect[1][0],
                                                                               rect[0][1])


def rect_center(rect):
    return [(rect[0][0] + rect[1][0]) / 2, (rect[0][1] + rect[1][1]) / 2]


def rect_vertex_set(rect):
    return Set([(rect[0][0], rect[0][1]), (rect[0][0], rect[1][1]), (rect[1][0], rect[0][1]), (rect[1][0], rect[1][1])])


def is_rect_cover_rect(rect, query):
    bool_matrix = np.zeros((2, rect.shape[1]))
    bool_matrix[0, :] = rect[0, :] <= query[0, :]
    bool_matrix[1, :] = rect[1, :] >= query[1, :]

    if np.all(bool_matrix):  # if query range contains node range
        return True
    return False


def is_rect_cover(rect, loc):
    """
    checks if the rectangle covers a point
    """
    bool_m1 = rect[0, 0] <= loc[0] <= rect[1, 0]
    bool_m2 = rect[0, 1] <= loc[1] <= rect[1, 1]
    bool_m = np.logical_and(bool_m1, bool_m2)
    if bool_m:
        return True
    else:
        return False


def is_cover_or_intersect(rect, query):
    bool_matrix = np.zeros((2, rect.shape[1]))
    bool_matrix[0, :] = rect[0, :] <= query[0, :]
    bool_matrix[1, :] = rect[1, :] >= query[1, :]

    if np.all(bool_matrix):  # if query range contains node range
        return True
    elif is_intersect(rect, query):
        return True
    return False


__is_cover_or_intersect = is_cover_or_intersect


def is_range_overlap(range1, range2):
    """
    check if two ranges overlap each others
    """
    if range2[0] <= range1[0] <= range2[1] or range2[0] <= range1[1] <= range2[1] or range1[0] <= range2[0] <= range1[
        1] or range1[0] <= range2[1] <= range1[1]:
        return True

    return False


# http://en.wikipedia.org/wiki/Zipf's_law
def _zipf_pmf(k, s, N):
    return (1.0 / k ** s) / np.sum([float(n) ** -s for n in range(1, N + 1)])


def _zipf_cdf(k, s, N):
    return np.sum([float(n) ** -s for n in range(1, k + 1)]) / np.sum([float(n) ** -s for n in range(1, N + 1)])


# check if three points are counterclock wise order
def _ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) >= (B[1] - A[1]) * (C[0] - A[0])


def is_intersect_segment(A, B, C, D):
    """
    check if two segments intersect each others
    http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return _ccw(A, C, D) != _ccw(B, C, D) and _ccw(A, B, C) != _ccw(A, B, D)


# http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
def _wrg(wgts):
    totals = np.cumsum(wgts)
    wgtSum = totals[-1]

    # speed up namespace lookups
    ru01 = random.random
    bi_r = bisect.bisect_right

    while True:
        yield bi_r(totals, ru01() * wgtSum)

        # def displayCounts(container):

# cts = Counter(container)
# for c in cts:
# print c, cts[c]

if __name__ == "__main__":
    # using numpy indexing on numpy array
    colors = np.array(["yellow", "red", "green"])
    weights = [.19, .01, .8]
    print list(it.islice(_wrg(weights), 1))

#print is_rect_intersect_segment([1,2],[3,1],[1,0],[2,2])
#rec = np.array([[2,2],[5,5]])
#query = np.array([[0,0],[3,10]])
#print isIntersect(rec, query)


#lat1, lon1, lat2, lon2 = 39.436140 - Params.ONE_KM, -77.094491 - Params.ONE_KM, 39.436140 + Params.ONE_KM, -77.094491 + Params.ONE_KM
#print distance(lat2, lon2, lat1, lon1)
