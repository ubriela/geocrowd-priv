*** This repository is the implementation of the following papers: ***

1) Hien To, Gabriel Ghinita, Cyrus Shahabi. A Framework for Protecting Worker Location Privacy in Spatial Crowdsourcing. In Proceedings of the 40th International Conference on Very Large Data Bases (VLDB 2014)

2) Hien To, Gabriel Ghinita, Cyrus Shahabi. PrivGeoCrowd: A Toolbox for Studying Private Spatial Crowdsourcing (demo). 40th International Conference on Very Large Data Bases (ICDE 2015)

Related studies:

https://bitbucket.org/hto/geocast/

https://bitbucket.org/hto/geocrowd

https://bitbucket.org/hto/maximumcomplextask/

----------- VERSION ---------------------------------

1.1

----------- REQUIRED LIBRARIES ----------------------

numpy

----------- PACKAGES --------------------------------

    exp
        : experimentation
    common
        common classes
    geocast
        : geocast classes
    grid
        : grid-based PSDs
    htree
        : htree-base PSDs
    idce12
        : quadtree-based and kdtree-based PSDs
    localness
        : localness approaches
    tornado
        : a package to provide services for VLDB14's demo
    minball
        : library to find smallest enclosing circle
    log
        : PSD logs and geocast logs

----------- USAGES -----------------------------------

A) PSD Experiments

Parameters in htree/Params.py
    
    DATASET: name of the dataset, e.g., "yelp"
    + PercentSplit: budget allocated for split
    + splitFrac: splitting budget percentage for 1st level, (1-splitFrac) for 2nd level
    + minPartSizeHTree: maximum number of data points in a leaf node
    + dynamicGranularity: compute granularity of the htree automatically
    + c_htree: a constant for computing granularity automatically
    + CONSTRAINT_INFERENCE: applying constraint inference (for both htrees and grids)
    + partitionsHTree: maximum size of the htree
    + switchPartitionsHTree: determine the size of first level tree (used for both h-trees and grids)
    + nQuery: the number of queries
    + queryUnit: the maximum query size
    + IS_LOGGING: enable logging individual query result

Parameters in PSDExp.py (starting point is at the end)
    seed_list: an array of seeds
    shape_list: an array of query size (the maximum query size is queryUnit in htree/Params.py)
    eps_list: an array of budgets

Run experiments with the following command
>> python PSDExp.py

The aggregated estimation error is in the corresponding folder in ../output/ (see the PSD RESULT)

The core dump estimation error is in folder ../log/. Each file is associated with a method. (see the CORE DUMP)

The debug file is in exp/debug.log

B) Geocast Experiments
>> python GeocastExp.py

----------- PSD RESULT --------------------------------------

Each file contains aggregated results of nQuery*size(seed_list) queries.

Data in each file is in tabular format. Rows are corresponding to budgets (i.e., eps_list) 
while columns are corresponding to PSD variants (method_list in PSDExp.py)

----------- PSD CORE DUMP -----------------------------------

Each file in ../log/ dumps the following information: query coordinates, query area (square km), values,  errors

    Format

    minLat    minLng    maxLat    maxLng    area    real_val    est_val    abs_err    rel_err

----------- GEOCAST CORE DUMP -------------------------------

    File name: geocast_{eps}.log
    Each line format: <is_assigned, #cells, lat lon, cell_1 info, cell_2 info....>
    1, 1, -120.081499 47.863793, 1 4 4 33 678 674.4 3.2 1.0 1.0
    1, 2, -108.865274 32.231816, 3 0 23 9 75 47.1 2.1 0.906 0.906, 3 0 23 8 154 190.7 2.9 1.0 1.0

    Info of each cell includes: <k l i j  #workers noisy_count cost u_c u  >
        >> (k, l): index of the first level grid
        >> (i, j): index of the second level grid (in case level = 2), i,j=1 means this is the first level cell
        >> cost: distance between the task and the cell's center
        >> u_c: utility of the cell
        >> u: updated utility