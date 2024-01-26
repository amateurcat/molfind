import redis, pickle

class RedisWrapper():
    def __init__(self, prefix, host, port, db, password):
        self.prefix=prefix
        self.r = redis.Redis(host=host, port=port, db=db, password=password)
        self.r.auth(password)
    def append(self, ret):
        # named it append instead of __call__ in accordance of the form when using multiprocessing.Manager().list() as collector
        idx,v = ret
        v = pickle.dumps(v)
        self.r.set(self.prefix+str(idx), v)