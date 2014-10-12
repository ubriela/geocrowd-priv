# # Basic parameters
class Params(object):
    #
    # dataset = '../../dataset/gowalla_CA.dat'; resdir = '../../output/gowalla_ca/'
    # x_min=-124.3041; y_min=32.1714; x_max=-114.0043; y_max=41.9984 #gowalla_CA

    dataset = '../../dataset/mcdonald.dat'
    resdir = '../../output/mcdonald/'
    x_min = -159.5862
    y_min = 20.7507
    x_max = -67.2804
    y_max = 64.8490  # mcdonald

    # dataset = '../../dataset/tiger_NMWA.dat'; resdir = '../../output/tiger/'
    # x_min=-124.8193; y_min=31.3322; x_max=-103.0020; y_max=49.0025 # tiger

    # dataset = '../../dataset/landmark.dat'; resdir = '../../output/landmark/'
    #    x_min=-124.4384; y_min=24.5526; x_max=-67.0255; y_max=49.0016 # landmark

    #    dataset = '../../dataset/restrnts.dat'; resdir = '../../output/restrnts/'
    #    x_min=-124.4972; y_min=24.5473; x_max=-66.9844; y_max=48.9999 # restrnts

    #    dataset = '../../dataset/shopping.dat'; resdir = '../../output/shopping/'
    #    x_min=-124.2640; y_min=24.5515; x_max=-68.2106; y_max=48.9939 # shopping

    #    dataset = '../../dataset/parkrec.dat'; resdir = '../../output/parkrec/'
    #    x_min=-124.5249; y_min=24.5510; x_max=-66.9687; y_max=49.0010 # parkrec

    #    dataset = '../../dataset/zipcode.dat'; resdir = '../../output/zipcode/'
    #    x_min=-176.6368; y_min=17.9622; x_max=-65.2926; y_max=71.2995 # zipcode

    #    dataset = '../../dataset/truck.dat'; resdir = '../../output/truck/'
    #    x_min=23.5100; y_min=37.8103; x_max=24.0178; y_max=38.2966 # truck

    #    dataset = '../../dataset/ne.dat'; resdir = '../../output/ne/'
    #    x_min=0.0470; y_min=0.0543; x_max=0.9530; y_max=0.9457 # ne

    #    dataset = '../../dataset/na.dat'; resdir = '../../output/na/'
    #    x_min=-174.2057; y_min=14.6805; x_max=-52.7316; y_max=79.9842 # na

    #    dataset = '../../dataset/buses.dat'; resdir = '../../output/buses/'
    #    x_min=22.3331; y_min=37.8329; x_max=24.0203; y_max=38.7417 # buses

    #    dataset = '../../dataset/gowalla_SA.dat'; resdir = '../../output/gowalla_sa/'
    #    x_min=-63.3209; y_min=-176.3086; x_max=13.5330; y_max=-32.4150 # gowalla_SA

    #    dataset = '../../dataset/brightkite.dat'; resdir = '../../output/brightkite/'
    #    x_min=-94.5786; y_min=-179.1541; x_max=90; y_max=178 # brightkite

    NDATA = None
    NDIM = None
    LOW = None
    HIGH = None
    SRT = 0.01  # sampling rate
    nQuery = 50  # number of queries
    maxHeight = 8  # maximum tree height, for kdtrees and quadtrees
    unitGrid = 0.01  # cell unit in kd-cell


    # Common parameters
    WorstCase = False  # worst case in selectivity estimation
    Hx = 4  # max query height
    queryUnit = [8, 8]
    GridAlignment = True
    LevelExpansion = True

    # HTree
    maxHeightHTree = 2  # = 4 for ht-composite and ht-hybrid
    partitionsHTree = [64, 64]  # m,m
    switchPartitionsHTree = [8, 8]  # m2,m2
    PercentSplit = 0.4  # budget allocated for split
    minPartSizeHTree = 2 ** 5  # maximum number of data points in a leaf node
    maxRelativeError = .95  # maximum error allowed in median selection
    dynamicGranularity = False  # automatically compute granularity for htree
    c_htree = 1
    useLeafOnlyHTree = True

    # Grid (the parameter setting is recommended in ICDE'13)
    PercentGrid = 0.5  # budget allocated for first level, used for adaptive grid
    m = partitionsHTree[1]  # grid size, computed in Grid_uniform
    c = 10
    c2 = 5

    ONE_KM = 0.0089982311916  # convert km to degree

    ### Geocast algorithm
    # Variables
    MAR = 0.05
    MTD = 100  # kms
    U = 0.90  # Utility

    TaskNo = 200  # tasks per instance
    NEGATIVE_CELL = False  # couting cells with negative count into geocast query
    AR_FUNCTION = 1  # 1 means step function, 2 means linear function
    STEPS = 100
    GAMMA = 0.9
    LEVEL_EXPANSION = False
    UNDER_EST = .5
    GEOCAST_LOG = True

    # This parameter setting is of the same as ICDE'12
    def __init__(self, seed):
        self.Eps = 0.5  # epsilon
        self.Seed = seed
        self.minPartSize = 2 ** 5  # maximum number of data points in a leaf node
        self.Percent = 0.3  # budget allocated for split
        self.switchLevel = 3  # switch level for hybrid tree
        self.Res = 18  # Hilbert-R tree resolution
        self.useLeafOnly = False  # only True for kd-cell and kd-noisymean
        self.cellDistance = 40  # threshold to test uniformity in kd-cell

        #        self.geoBudget = 'none' # uniform budgeting
        #        self.geoBudget = 'aggressive' # geometric exponent i
        #        self.geoBudget = 'quadratic' # geometric exponent i/2
        self.geoBudget = 'optimal'  # geometric exponent i/3
        #        self.geoBudget = 'quartic' # geometric exponent i/4

        self.splitScheme = 'expo'  # exponential mechanism
        #        self.splitScheme = 'noisyMean' # noisy mean approximation