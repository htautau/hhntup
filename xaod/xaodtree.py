import time
from rootpy.tree.filtering import EventFilterList
from rootpy.tree import Tree
from rootpy import log; log = log[__name__]
import ROOT

from .xaodobject import xAODTreeCollection

class xAODTree(object):
    
    def __init__(self, chain, filters=None, events=-1):
        self._chain = chain
        self._tree = ROOT.xAOD.MakeTransientTree(self._chain)
        # Create the TStore that hold the shallow copies
        self._store = ROOT.xAOD.TStore()
        log.info(self._tree)
        self._collections = {}
        self._events = events
        self._init_time = time.time()
        if filters is None:
            self._filters = EventFilterList([])
        else:
            self._filters = filters

    def define_collection(self, name, collection_name, mix=None, decorate_func=None):
        self._collections[name] =  (collection_name, mix, decorate_func)
        # return coll

    def reset_collections(self):
        for coll in self._collections.iterkeys():
            coll.reset()
    
    def __iter__(self):
        passed_events = 0
        entries = 0
        total_entries = float(self._tree.GetEntries())
        t2 = self._init_time
        for i in xrange(self._tree.GetEntries()):
            entries += 1
            self._tree.GetEntry(i)
            for name, (coll_name, mix, decorate_func) in self._collections.items():
                coll = xAODTreeCollection(
                    self._tree, name, coll_name, 
                    mix=mix, decorate_func=decorate_func)
                object.__setattr__(self._tree, name, coll)
            if self._filters(self._tree):
                yield self._tree
                passed_events +=1
                if self._events == passed_events:
                    break
            if time.time() - t2 > 2:
                entry_rate = int(entries / (time.time() - self._init_time))
                log.info(
                    "{0:d} entries per second. "
                    "{1:.0f}% done current tree".format(
                        entry_rate, 
                        100 * entries / total_entries))
                t2 = time.time()
            self._filters.finalize()
            self._store.clear()

    def __len__(self):
        return self._tree.GetEntries()

