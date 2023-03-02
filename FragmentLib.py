from GraphCompiler import atoms2graph
from pathlib import Path
from networkx.algorithms import isomorphism
import pickle, time, ase.io

class HashTable():
    
    def __init__(self):
        self.table = {}
        
    def append(self, x, add_to=None):
        assert isinstance(x,tuple) and (len(x)>1), "You must have at least one key and one value"
        if add_to==None:
            add_to = self.table
        if len(x) > 2:
            try:
                add_to = add_to[x[0]]
            except KeyError:
                add_to[x[0]] = {}
                add_to = add_to[x[0]]
                
            return self.append(x[1:],add_to)
        else:
            add_to[x[0]] = x[1]
            
    def get(self,x,find_in=None):
        if find_in==None:
            find_in = self.table
            
        if len(x) == 1:
            return find_in[x[0]]   # Must be tail recursion here
        else:
            find_in = find_in[x[0]]
            return self.get(x[1:], find_in)   # Must be tail recursion here
                

class FragmentLib():
    # molecule structure lib that save a molecule's(fragment) name and it's topological
    # info in nx.Graph, save total number of atoms and sorted atom type array in addition
    # and make them a hashtable to accelerate the query process
    # may need to add delete/edit function, but since we don't have too many candidate fragments
    # simply remake the FragmentLib isntance when you want to make changes
    
    def __init__(self):
        self.lib = HashTable()
        
    def append(self, name, atoms):
        gs = atoms2graph(atoms)
        assert len(gs) == 1, "The input structure contain disjoint fragments"
        self.lib.append( gs[0] + (name,) )
    
    def search(self, n, tfp, g):
        
        ret = "Not Found"
        try:
            possibleG = self.lib.get( (n,tfp) )
            for G, name in possibleG.items():
                # need to check if this actually used edge info or not
                GM = isomorphism.GraphMatcher(G, g)   
                if GM.is_isomorphic():
                    ret = name
                    
        except KeyError:
            # May need to print which keyword made the query failed here
            pass
                
        return ret
    

if __name__ == "__main__":
    
    # I stored some xyz files of small molecules in ./StructureFiles/
    # named by their chemical name, so just traverse that folder to make a test FragmentLib
    # may need to add more or import from PuBChem database
    
    start = time.time()
    
    p = Path('StructureFiles/')
    FRAGMENT_LIB = FragmentLib()
    for f in p.glob('*.xyz'):
        name = f.stem
        atoms = ase.io.read(f)
        FRAGMENT_LIB.append(name,atoms)
        
    with open("FragmentLib.pkl", 'wb') as fw:
        pickle.dump(FRAGMENT_LIB, fw)
    
    print("Finish building the lib, time cost: " +  str(time.time()-start))