import numpy as np
import ase.io, pickle, io, sys
from FragmentLib import *
from GraphCompiler import atoms2graph
from multiprocessing import JoinableQueue, Process, Manager
from collections import Counter
import numpy as np

with open("FragmentLib.pkl", 'rb') as fr:
    FRAGMENT_LIB = pickle.load(fr)

def Graphfind(atoms):
    nl = []
    Gs = atoms2graph(atoms)
    for info in Gs:
        n = FRAGMENT_LIB.search(*info)
        nl.append(n)
    ret = Counter(nl)
    
    return (ret, Gs)
    

class MolFinder():
    def __init__(self, queue, collector, finder=Graphfind):
        self.queue = queue
        self.p = Process(target=self.wrapper)

        #define how to treat returns from single process
        #and how it interact to the collector
        self.collector = collector
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
                self.collector((index,ret))
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