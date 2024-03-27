import numpy as np
from multiprocessing import JoinableQueue, Process
import ase.io

class MolFinder():
    def __init__(self, queue, collector, finder, postprocess=None):
        self.queue = queue
        self.p = Process(target=self.wrapper)

        #define how to treat returns from single process
        #and how it interact to the collector
        self.collector = collector
        self.finder = finder
        self.postprocess = postprocess

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

                if self.postprocess is not None:
                    ret = self.postprocess(ret)

                self.collector.append((index,ret))
                self.queue.task_done()

    def start(self):
        #print('consumer start!\n')
        self.p.start()

    def join(self):
        self.p.join()

class AtomsQueueGenerator():
    def __init__(self, xyzfile, Nconsumer, queue=None, cell=np.array([0,0,0]), index_range=None, add_stop_signal=True):
        self.f = xyzfile
        self.Nconsumer = Nconsumer
        self.cell = cell
        self.index_range = index_range
        if queue:
            self.queue = queue
        else:
            print('create new queue')
            self.queue = JoinableQueue()
            
        self.p = Process(target=self.wrapper)
        self.add_stop_signal = add_stop_signal
        
    def wrapper(self):
        atoms = ase.io.read(self.f,":" if self.index_range is None else "%d:%d"%(self.index_range[0],self.index_range[1]))
        pbc = self.cell.astype(bool)
        for i,a in enumerate(atoms):
            # it's strange that ASE cannot read cell info from a xyz file 
            # that contain multiple conformers, so just add it each time
            a.pbc = pbc 
            a.cell = self.cell
            self.queue.put((i,a))
        
        if self.add_stop_signal:
            for _ in range(self.Nconsumer):
                self.queue.put(("STOP CONSUMER",))
                
    def start(self):
        self.p.start()

    def join(self):
        self.p.join()