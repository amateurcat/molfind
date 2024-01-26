from Parallel import AtomsQueueGenerator, MolFinder, GraphFind
from OBFinder import OBfind
from Collector import RedisWrapper



if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from multiprocessing import Manager
    import pickle, json
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument('--traj_file', type=str, help='path to the trajectory file')
    parser.add_argument('--index_range', type=str, default=None, help='index range in the form of "(start_index,end_index)"')
    parser.add_argument('--n_cpu', type=int, default=8)
    parser.add_argument('--cell', type=str, default=None, help='cell parameters in the form of "x,y,z"')
    parser.add_argument('--engine', type=str, default='OBFinder', help='engine to search for fragments, either GraphFind or OBFind')
    parser.add_argument('--to_redis', type=bool, default=False, help='if True, upload the results to a Redis server')
    parser.add_argument('--redis_login', type=str, default='AUTH.json', help='json file that store the login information for the Redis server')
    parser.add_argument('--output', type=str, default=None, help='output file name')

    ###TODO: what if the coordinates contain negative values? how to redefine the cell?
    args = parser.parse_args()
    traj_file = Path(args.traj_file)
    index_range = eval(args.index_range) if args.index_range is not None else None
    n_cpu = args.n_cpu

    cell = np.array([float(i) for i in args.cell.split(',')]) if args.cell is not None else None
    engine = OBfind if args.engine == 'OBFinder' else GraphFind
    to_redis = args.to_redis
    redis_login = args.redis_login
    output = args.output

    if to_redis:
        with open(redis_login,'r') as f:
            AUTH = json.load(f)
        collector = RedisWrapper(prefix='%s:'%(traj_file.stem), host=AUTH['host'], port=AUTH['port'], db=AUTH['db'], password=AUTH['password'])
    else:
        if output is None:
            print("No output method specified, will still run but the result will not be saved")
            

    with Manager() as cl:
        if not to_redis:
            collector = cl.list()

        producer = AtomsQueueGenerator(traj_file, n_cpu, None, cell, index_range)
        tasks = producer.queue

        consumer_list = [MolFinder(tasks, collector, engine) for i in range(n_cpu)]
        
        producer.start()
        for c in consumer_list:
            c.start()
            
        producer.join()

        for c in consumer_list:
            c.join()
        
        if not to_redis:
            results = list(collector)
        
    results = sorted(results, key=lambda x:x[0])

    if output is not None:
        with open(output, 'wb') as fw:
            pickle.dump(results, fw)



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