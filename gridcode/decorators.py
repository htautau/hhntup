import functools

class memoize(object):
    """
    Memoization for functions and instance methods
    """
    
    def __init__(self, func):
        
        self.func = func
        self.cache = {}
    
    def __call__(self, *args):
        
        return self.cache_get(args, lambda: self.func(*args))
    
    def __get__(self, obj, objtype):
        
        return self.cache_get(obj, lambda: self.__class__(functools.partial(self.func, obj)))
    
    def cache_get(self, key, func):
        
        try:
            return self.cache[key]
        except KeyError:
            
            self.cache[key] = func()
            return self.cache[key]
