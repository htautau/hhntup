# Copyright 2012 the rootpy developers
# distributed under the terms of the GNU General Public License
from __future__ import absolute_import
from rootpy import log; log = log[__name__]

__all__ = [
    'xAODTreeCollection',
]

__MIXINS__ = {}



class xAODTreeCollection(object):

    def __init__(self, tree, name, collection_name, cache=True, mix=None):

        self.tree = tree
        self.name = name
        self.collection = getattr(self.tree, collection_name)
        self.selection = None
        self.mix = mix

        self.__cache_objects = cache
        self.__cache = {}


    def __nonzero__(self):

        return len(self) > 0

    def reset(self):

        self.reset_selection()
        self.reset_cache()

    def reset_selection(self):

        self.selection = None

    def reset_cache(self):

        self.__cache = {}

    def remove(self, thing):

        if self.selection is None:
            self.selection = range(len(self))
        for i, other in enumerate(self):
            if thing == other:
                self.selection.pop(i)
                break

    def pop(self, index):

        if self.selection is None:
            self.selection = range(len(self))
        thing = self[index]
        self.selection.pop(index)
        return thing

    def select(self, func):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection = [
            i for i, thing in zip(self.selection, self)
            if func(thing)]

    def select_indices(self, indices):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection = [self.selection[i] for i in indices]

    def mask(self, func):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection = [
            i for i, thing in zip(self.selection, self)
            if not func(thing)]

    def mask_indices(self, indices):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection = [
            j for i, j in enumerate(self.selection)
            if i not in indices]

    def _wrap_sort_key(self, key):

        def wrapped_key(index):
            return key(self.getitem(index))
        return wrapped_key

    def sort(self, key, **kwargs):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection.sort(key=self._wrap_sort_key(key), **kwargs)

    def slice(self, start=0, stop=None, step=1):

        if self.selection is None:
            self.selection = range(len(self))
        self.selection = self.selection[slice(start, stop, step)]

    def make_persistent(self):
        """
        Perform actual selection and sorting on underlying
        attribute vectors
        """
        pass

    def getitem(self, index):
        """
        direct access without going through self.selection
        """
        if index >= len(self.collection):
            raise IndexError(index)
        if self.__cache_objects and index in self.__cache:
            return self.__cache[index]
        if self.mix is not None:
            obj = self.mix(self.collection[index])
        else:
            obj = self.collection[index]
        if self.__cache_objects:
            self.__cache[index] = obj
        return obj

    def __getitem__(self, index):

        if isinstance(index, slice):
            return [self.collection[i] for i in xrange(*index.indices(len(self)))]
        if index >= len(self):
            raise IndexError(index)
        if self.selection is not None:
            index = self.selection[index]
        if self.__cache_objects and index in self.__cache:
            return self.__cache[index]
        if self.mix is not None:
            obj = self.mix(self.collection[index])
        else:
            obj = self.collection[index]
        if self.__cache_objects:
            self.__cache[index] = obj
        return obj

    def len(self):
        """
        length of original collection
        """
        return len(self.collection)

    def __len__(self):

        if self.selection is not None:
            return len(self.selection)
        return len(self.collection)

    def __iter__(self):

        for index in xrange(len(self)):
            yield self.__getitem__(index)


