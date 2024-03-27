from Parallel import AtomsQueueGenerator, MolFinder
from GraphCompiler import GraphFind
from OBFinder import OBfind
from Collector import RedisWrapper

import pickle
from multiprocessing import Manager
import numpy as np


def molfind(traj_file, n_cpu=8, index_range=None, cell=None, engine=OBfind, collector=None, output=None):
    #separate the main funtion out so that you can import and call it elsewhere for testing
    #also separate collector so that you can use existing collector in an interactive session
    ###TODO: what if the coordinates contain negative values? how to redefine the cell?

    if (output is None) and (collector is None):
        print("No output method specified, will still run but the result will not be saved")        

    with Manager() as cl:   ###TODO: seems this is not necessary if not using cl.list()
        if collector is None:
            collector = cl.list()

        producer = AtomsQueueGenerator(traj_file, n_cpu, None, cell, index_range)
        tasks = producer.queue

        consumer_list = [MolFinder(tasks, collector, engine) for i in range(n_cpu)]
        
        producer.start()
        for c in consumer_list:
            c.start()

        for c in consumer_list:
            c.join()
            
        producer.join()   ###TODO: Check if this should be put after consumers joined

        if output is not None:
            results = list(collector)   ###TODO: think how to get results from other collectors
            results = sorted(results, key=lambda x:x[0])
            with open(output, 'wb') as fw:
                pickle.dump(results, fw)


if __name__ == "__main__":
    import argparse, json
    from pathlib import Path
    import numpy as np

    parser = argparse.ArgumentParser()
    parser.add_argument('--traj_file', type=str, help='path to the trajectory file')
    parser.add_argument('--index_range', type=str, default=None, help='index range in the form of "(start_index,end_index)"')
    parser.add_argument('--n_cpu', type=int, default=8)
    parser.add_argument('--cell', type=str, default=None, help='cell parameters in the form of "x,y,z"')
    parser.add_argument('--engine', type=str, default='OBFinder', help='engine to search for fragments, either GraphFind or OBFind')
    parser.add_argument('--to_redis', type=bool, default=False, help='if True, upload the results to a Redis server')
    parser.add_argument('--redis_login', type=str, default='AUTH.json', help='json file that store the login information for the Redis server')
    parser.add_argument('--redis_prefix', type=str, default=None, help='prefix for the keys in the Redis server')
    parser.add_argument('--output', type=str, default=None, help='output file name')

    
    args = parser.parse_args()
    traj_file = Path(args.traj_file)
    index_range = eval(args.index_range) if args.index_range is not None else None
    n_cpu = args.n_cpu

    cell = np.array([float(i) for i in args.cell.split(',')]) if args.cell is not None else None
    engine = OBfind if args.engine == 'OBFinder' else GraphFind
    to_redis = args.to_redis
    redis_login = args.redis_login
    redis_prefix = args.redis_prefix if args.redis_prefix is not None else traj_file.stem + ':'
    output = args.output

    collector = None
    if to_redis:
        with open(redis_login,'r') as f:
            AUTH = json.load(f)
        collector = RedisWrapper(prefix=redis_prefix, host=AUTH['host'], port=AUTH['port'], db=AUTH['db'], password=AUTH['password'])


    molfind(traj_file, n_cpu, index_range, cell, engine, collector, output)



'''
N = 8
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

'''