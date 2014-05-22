Dataset         size 
-----------------------
gowalla_CA      736,724
landmark        870,051
mcdonald        13,913
parkrec         48,588
restrnts        449,853
shopping        15,046
tiger_NMWA      1,634,167
zipcode         33,178
gowalla_SA	22,373


----------- GEOCAST LOG -----------
File name: geocast_{eps}.log
Each line format: <is_assigned, #cells, lat lon, cell_1 info, cell_2 info....>
1, 1, -120.081499 47.863793, 1 4 4 33 678 674.4 3.2 1.0 1.0
1, 2, -108.865274 32.231816, 3 0 23 9 75 47.1 2.1 0.906 0.906, 3 0 23 8 154 190.7 2.9 1.0 1.0

Info of each cell includes: <k l i j  #workers noisy_count cost u_c u  >
    + (k, l): index of the first level grid
    + (i, j): index of the second level grid (in case level = 2), i,j=1 means this is the first level cell
    + cost: distance between the task and the cell's center
    + u_c: utility of the cell
    + u: updated utility