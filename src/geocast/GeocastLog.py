"""
This file dump the geocast queries into a log file.

Each geocast query includes a list of cells [Cell].
Each cell has a unique index in the the tree structure <level, workers,(k,l), [(i,j)]>
    + level:  the level, either 1 or 2 for adaptive grid
    + workers: the number of workers in the node
    + noisy count
    + (k, l): index of the first level grid
    + (i,j): index of the second level grid (in case level = 2)
    
Output file: geocast.log
Each line contains <lat lng, level1 k1 l1 i1 j1, level2 k2 l2....>
"""


def geocast_log(prefix, line, eps):
    f = open("../" + prefix + "_" + str(eps) + ".log", "a")
    f.write(line)
    f.close()
