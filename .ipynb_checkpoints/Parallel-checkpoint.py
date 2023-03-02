import numpy as np
import ase.io, pickle, io, sys
from FragmentLib import *
from GraphCompiler import atoms2graph
from multiprocessing import JoinableQueue, Process, Manager
from collections import Counter
import numpy as np

sys.path.append('/home/shuhao/softwares/miniconda3/envs/ani_3.6/lib/python3.6/site-packages/openbabel')
import openbabel

with open("FragmentLib.pkl", 'rb') as fr:
    FRAGMENT_LIB = pickle.load(fr)
    
def OBfind(atoms):
    # Based on openbabel python interface v3.1.1
    # load xyz file and use ConnectTheDots() to build bond connections
    # then output SMILES and count how many fragments in SMILES
    # May contain different SMILES that actually represent a same species
    cmatrix = atoms.cell.array
    #print(cmatrix)
    ucell = openbabel.OBUnitCell()
    ucell.SetData(openbabel.vector3(*cmatrix[0]), openbabel.vector3(*cmatrix[1]), openbabel.vector3(*cmatrix[2]))
    ucell.SetSpaceGroup('P 1')
    
    s = io.StringIO()
    ase.io.write(s, atoms, format='xyz')
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, s.getvalue())
    mol.CloneData(ucell)
    mol.SetPeriodicMol(True)

    openbabel.OBUnitCell.FillUnitCell(ucell,mol)
    mol.ConnectTheDots()
    obConversion.WriteString(mol)

    smile = obConversion.WriteString(mol)
    #print(smile)
    ret = Counter(smile.split("\t")[0].split('.'))
    #print(ret)
    
    return ret

def Graphfind(atoms):
    nl = []
    for info in atoms2graph(atoms):
        n = FRAGMENT_LIB.search(*info)
        nl.append(n)
    ret = Counter(nl)
    
    return ret
    

class MolFinder():
    def __init__(self, queue, add_to, finder=OBfind):
        self.queue = queue
        self.p = Process(target=self.wrapper)

        #define how to treat returns from single process
        #and how it interact to the collector
        self.add_to = add_to
        self.finder = finder

    def wrapper(self):
        #print('calling wrapper\n')
        while True:
            args = self.queue.get()   ###!!! This will NOT raise error and will NOT stop when the queue is empty
            if "STOP CONSUMER" in args:   ###!!! args is a tuple! Do NOT use "==" here!
                print("STOP CONSUMER")
                self.queue.task_done()
                break
            else:
                index, atoms = args
                ret = self.finder(atoms)
                self.add_to.append((index,ret))
                self.queue.task_done()

    def start(self):
        #print('consumer start!\n')
        self.p.start()

    def join(self):
        self.p.join()

class AtomsQueueGenerator():
    def __init__(self, xyzfile, Nconsumer, queue=None, cell=np.array([0,0,0])):
        self.f = xyzfile
        self.Nconsumer = Nconsumer
        self.cell = cell
        if queue:
            self.queue = queue
        else:
            print('create new queue')
            self.queue = JoinableQueue()
            
        self.p = Process(target=self.wrapper)
        
    def wrapper(self):
        atoms = ase.io.read(self.f,":")
        pbc = self.cell.astype(bool)
        for i,a in enumerate(atoms):
            # it's strange that ASE cannot read cell info from a xyz file 
            # that contain multiple conformers, so just add it each time
            a.pbc = pbc 
            a.cell = self.cell
            self.queue.put((i,a))
            
        for _ in range(self.Nconsumer):
            self.queue.put(("STOP CONSUMER",))
                
    def start(self):
        self.p.start()

    def join(self):
        self.p.join()