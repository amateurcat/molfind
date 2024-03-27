import redis, pickle

class RedisWrapper():
    def __init__(self, prefix, host, port, db, password):
        self.prefix=prefix
        self.r = redis.Redis(host=host, port=port, db=db, password=password)
        self.r.auth(password)
    def append(self, ret):
        # named it append instead of __call__ to fit the form when using multiprocessing.Manager().list() as collector
        idx,v = ret
        v = pickle.dumps(v)
        self.r.set(self.prefix+str(idx), v)
    def __iter__(self, prefix=None, sort=False):
        # assuming the keys are in the form of ${prefix}:index
        # this will return all keys with the given prefix in the current db, and their values
        # if you are using an existing Wrapper in the interactive model, 
        # make sure it only contains the keys you want to retrieve
        if prefix is None:
            prefix = self.prefix
        keys = self.r.keys(prefix+'*')
        if sort:
            keys = sorted(keys, key=lambda x:int(x.decode("utf-8").split(":")[-1]))

        for k in keys:
            yield (int(k.decode("utf-8").split(":")[-1]), pickle.loads(self.r.get(k)))