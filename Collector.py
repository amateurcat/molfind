import redis, pickle

class RedisWrapper():
    def __init__(self, prefix, host, port, db, password):
        self.prefix=prefix
        self.r = redis.Redis(host=host, port=port, db=db, password=password)
        self.r.auth(password)
    def __call__(self, ret):
        idx,v = ret
        v = pickle.dumps(v)
        self.r.set(self.prefix+str(idx), v)