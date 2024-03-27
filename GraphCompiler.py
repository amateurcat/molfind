import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import ase.io, sys

#if not insert this path to the front, jupyter will use scipy in ~/.local/python3.8/site-packages instead
#sys.path.insert(0, "/home/shuhao/softwares/miniconda3/envs/autoani/lib/python3.8/site-packages")
from scipy.spatial.distance import pdist   #scipy v1.8.0

from primePy import primes

PRIME_LIST = primes.first(119)[1:]

#all of these bond length data are from google search, need a better data source
BOND_LENGTH_LIB = {(1,1):(0.737,),
                   (1,6):(1.094,), #C-H bondlength in ethane
                   (1,7):(1.017,),
                   (1,8):(0.969,), #O-H bond in WATER, this might be different in RO-H, need further check!!!
                   (6,6):(1.535,1.339,1.203,), #C-C, C=C, C#C
                   (6,7):(1.469,1.279,1.154,), #R-NH2, RC=N-R, R-N#C
                   (6,8):(1.430,1.230,1.128,), #C-O, C=O, C#O+
                   (7,8):(1.360,1.150,1.060,), #N-O, N=O, N#O+
                   (8,8):(1.208,), #O=O
                  }

BOND_LENGTH_LIB = {PRIME_LIST[k[0]-1]*PRIME_LIST[k[1]-1] : v for k,v in BOND_LENGTH_LIB.items()}
BUFFER = 1.05

def numbers2prime(numbers):
    return np.array([PRIME_LIST[n-1] for n in numbers])


###TODO: This procedure may be accelerated by using Pytorch 
# if we can find a reliable pytorch implementation of pairwise distance of point cloud
# with periodic boundary condition with arbitrary cell shape
# Simply calculating the adjacency matrix by batch with tensor operation


def get_CompositeNumberMatrix(atoms, buffer=BUFFER, bond_length_lib=BOND_LENGTH_LIB):
    #assuming pbc and cell info is stored into the input Atoms instance
    #considering pbc by defalut, turn off pbc by switching Atoms.pbc to False
    #bond length cutoff will be multiple by bufffer
    #may be better to use different buffer for different type of bond
    
    numbers = atoms.get_atomic_numbers()
    N = len(numbers)
    
    pl = numbers2prime(numbers).reshape(1,-1)
    temp_cnm = np.dot(pl.T, pl)   #O(n^2), but very easy to accelerate
    
    dm = atoms.get_all_distances(mic=True)   #O(n^2), but very easy to accelerate   ###TODO: check if this consider pbc
    
    cnm = np.ones(temp_cnm.shape)   #need to find a better abbreviation
    
    #O(n^2), can accelerate by using the fact that temp_cnm and dm are symmetric matrices 
    for x in range(N):
        for y in range(N):
            c = 0
            if x!=y:
                t = temp_cnm[x][y]
                if t in bond_length_lib.keys():
                    for o,bl in enumerate(bond_length_lib[t]):
                        if dm[x][y] <= bl * buffer:
                            c = t * 2**o
            
            cnm[x][y] = c
    
    cnm = cnm.astype(int)        
    return numbers, cnm

def cnm2labeledgraph(cnm):
    N = len(cnm)
    G = nx.Graph()
    for x in range(N):
        for y in range(x,N):
            if cnm[x][y] > 0:
                G.add_edge(x, y, type=cnm[x][y])
    
    return G

def print_graph(G):
    #for visulization
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True)
    edge_labels = dict([((n1, n2), d['type']) for n1, n2, d in G.edges(data=True)])
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, label_pos=0.5,)
    plt.show()
    
def atoms2graph(atoms):
    # This function assume the input Atoms instance contains multiple disjoint
    # fragments, and return a list of info of fragments
    # By doing so I can reuse this function when building the lib 
    # where I know input is an intact molecule and therefore the output 
    # is a length 1 list 
    numbers, cnm = get_CompositeNumberMatrix(atoms)
    G = cnm2labeledgraph(cnm)
    
    #I planned to use product of all prime indexes instead
    #but large fragments would have huge product
    #need to check if it's still faster than comparing array
    ret = []
    for subgraph_indices in nx.connected_components(G):
        subgraph_numbers = numbers[list(subgraph_indices)]
        subgraph_N = len(subgraph_numbers)
        subgraph_tfp = tuple(np.sort(subgraph_numbers.copy()))
        subgraph = G.subgraph(subgraph_indices).copy()
        ret.append( (subgraph_N, subgraph_tfp, subgraph) )
    
    return ret

class GraphFind():
    def __init__(self, fragment_lib="./FragmentLib.pkl", buffer=1.05):
        with open(fragment_lib, 'rb') as fr:
            self.fragment_lib = pickle.load(fr)
        self.buffer = buffer

    def __call__(self, atoms):
        nl = []
        Gs = atoms2graph(atoms)
        for info in Gs:
            n = FRAGMENT_LIB.search(*info)
            nl.append(n)
        ret = Counter(nl)
        
        return (ret, Gs)