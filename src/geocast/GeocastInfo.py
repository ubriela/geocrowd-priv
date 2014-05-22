"""
Stored information of a geocast query
"""
class GeocastInfo:
    loc = None # task location
    cells = None   # a list of geocast cells
    
    ANW = None
    ATD = None
    DPPT = None
    
    def __init__(self, assigned, loc, cells):
	self.assigned = assigned
        self.loc = loc
        self.cells = cells
        
    def logging(self):
        if self.cells == None:
            return None
        str_loc = " ".join(map(str, [self.loc[0], self.loc[1]]))
        str_cells = []
        for i in range(len(self.cells)):
            str_cells.append(" ".join(map(str, self.cells[i][1])))
        return ", ".join([str(self.assigned), str(len(self.cells)), str_loc, ", ".join(str_cells)])