from FragmentLib import *
from OBFinder import *
from Collector import RedisWrapper
import pickle, json
import ase.io

sys.path.append('/home/shuhao/')
from ParallelTools import *

N = 4
cell = np.array([22.0,36.0,41.52])
'''
with Manager() as cl:
    collector = cl.list()        
    producer = AtomsQueueGenerator('../cellulose_trajectory.xyz',N,None,cell)
    tasks = producer.queue
    consumer_list = [MolFinder(tasks, collector) for i in range(N)]
    
    producer.start()
    for c in consumer_list:
        c.start()
        
    producer.join()

    for c in consumer_list:
        c.join()
        
    results = list(collector)
    
results = sorted(results, key=lambda x:x[0])

with open("OBsearch_sorted_molfind_result_1000Ksecondhalf.pkl", 'wb') as fw:
    pickle.dump(results, fw)

'''

with open('AUTH.json') as f:
    AUTH = json.load(f)
    
RedisCollector = RedisWrapper(prefix='test:', host=AUTH['host'], port=AUTH['port'], db=AUTH['db'], password=AUTH['password'])

def read_conformers_xyz(filename):
    l = ase.io.read(filename, index=":")
    for atoms in l:
        atoms.set_cell(cell)
        yield atoms


AtomProducer = Producer(read_conformers_xyz('test/old_cellulose_0-7.xyz'), queue=None, add_stop_signal=N)
tasks = AtomProducer.queue

consumer_list = [Consumer(tasks, worker=OBfind, collector=RedisCollector) for i in range(N)] 

AtomProducer.start()
for c in consumer_list:
    c.start()

AtomProducer.join()

for c in consumer_list:
    c.join()