from FragmentLib import *
from Parallel import *
from multiprocessing import Manager
import pickle

N = 25
cell = np.array([22.0,36.0,41.52])
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