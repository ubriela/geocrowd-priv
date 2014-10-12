"""
Geocast with scaling algorithm (scale the size of the geocast query by a factor)
"""

import sys
import math
import random

import numpy as np

from Params import Params
from Geocrowd import rect_query_points, query_init_naive
from Utils import performed_tasks_naive, utility_naive


def resize(query, scale, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    center = [(query[0, 0] + query[1, 0]) / 2.0, (query[0, 1] + query[1, 1]) / 2.0]
    x_range, y_range = query[1, 0] - query[0, 0], query[1, 1] - query[0, 1]
    x_range, y_range = x_range * scale, y_range * scale

    x_low, x_high = center[0] - x_range / 2, center[0] + x_range / 2
    y_low, y_high = center[1] - y_range / 2, center[1] + y_range / 2
    query = [[max(x_low, x1), max(y_low, y1)], [min(x_high, x2), min(y_high, y2)]]
    return np.array(query)


# simple scaling algorithm
def scaling_query2(tree, w, L):
    eps = 10
    c = 5
    iter = 0
    query = query_init_naive(L, Params.x_min, Params.y_min, Params.x_max, Params.y_max)

    last_est, last_last_est = sys.maxint, sys.maxint
    while True:
        est = tree.rangeCount(query)
        if math.fabs(est - w) <= c:
            break
        else:
            if est == 0:
                scale = random.random() + 0.5
            else:
                scale = math.sqrt(math.fabs(( w + 0.0) / est))
            query = resize(query, scale, Params.x_min, Params.y_min, Params.x_max, Params.y_max)
            iter += 1
        if iter >= 10 or math.fabs(est - last_est) < eps or math.fabs(est - last_last_est) < eps:
            break

        last_last_est = last_est
        last_est = est
    return query


def scaling_query(tree, L, U):
    iter = 0
    query = query_init_naive(L, Params.x_min, Params.y_min, Params.x_max, Params.y_max)

    while True:
        est = tree.rangeCount(query)
        u, dist = utility_naive(query, est, Params.MTD)
        if u >= U:
            break
        else:
            scale = 1.1
            query = resize(query, scale, Params.x_min, Params.y_min, Params.x_max, Params.y_max)
            iter += 1
            if iter >= 50:
                break
    return query


def geocast_naive(tree, data, L, FCFS, U, seed):
    query = scaling_query(tree, L, U)

    # find all workers in the query
    locs = rect_query_points(data, query).transpose()
    # actual utility & average travel cost
    isPerformed, dist = performed_tasks_naive(locs, Params.MTD, L, FCFS, seed)

    if isPerformed:  # the task is performed
        return len(locs), True, dist
    return len(locs), False, dist