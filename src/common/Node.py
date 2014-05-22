class Node(object):
    """ generic node class """
    
    def __init__(self):
        self.n_data = None # list of data points
        self.n_box = None # 2x2 matrix [[x_min,y_min],[x_max,y_max]]
        self.children = [] # list of its children
        self.n_count = 0 # noisy count of this node
        self.n_depth = 0
        self.n_budget = 0 # represented by height of the tree
        self.n_isLeaf = False
	self.secondLevelPartitions = None # for adaptive grid only
        self.parent = None; # parent node
        self.index = None; # the order of the node in the parent's children
        self.eps = None
        self.neighbor = None
        
    def debug(self):
        print "Actual count " + str(len(self.n_data))
        print "Noisy count " + str(self.n_count)
        print "Depth " + str(self.n_depth)
        print "Leaf? " + str(self.n_isLeaf)
        print "Second level partition: " + str(self.secondLevelPartitions)
        print "Boundary " + str(self.n_box)
        print "Index " + str(self.index)
        
        