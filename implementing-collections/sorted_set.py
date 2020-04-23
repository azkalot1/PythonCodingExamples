from collections.abc import Sequence
from bisect import bisect_left
from itertools import chain


class SortedSet(Sequence):

    def __init__(self, items=None):
        self._items = sorted(set(items)) if items is not None else []

    def __contains__(self, item):  # allows to use %in%
        index = bisect_left(self._items, item)
        found = (index != len(self._items)) and (self._items[index] == item)
        return found

    def __len__(self):  # allows to use len()
        return len(self._items)

    def __iter__(self):  # allows to use iter(iterable) function
        return (iter(self._items))

    def __getitem__(self, index):
        result = self._items[index]
        return SortedSet(result) if isinstance(index, slice) else result
        # We want to return SortedSet if we use slice
        # If we use index we want to return just value

    def __repr__(self):  # allows to use print(obj) and have representation
        return "SortedSet({})".format(
            repr(self._items) if self._items else ''
        )

    def __eq__(self, rhs):  # How to compare? allows to use ==
        if not isinstance(rhs, SortedSet):
            return NotImplemented
        return self._items == rhs._items

    def __ne__(self, rhs):  # How to compare? allows to use ==
        if not isinstance(rhs, SortedSet):
            return NotImplemented
        return self._items != rhs._items

    # Now implement sequence protocol
    # We don't need __reversed__, since we have __getitem__ and __len__
    # Easy way - just inherit from Sequence
    # We write custom index, since it will be faster

    def index(self, item):
        index = bisect_left(self._items, item)
        if (index != len(self._items)) and (self._items[index] == item):
            return index
        raise ValueError('{} not found'.format(repr(item)))

    # Vanilla count works O(n), however, since
    # we alsways return 0 or 1 we can just run
    # binary search which runs O(logn)
    def count(self, item):
        return int(item in self)

    # To concat sets we impleement __add__
    def __add__(self, rhs):
        return SortedSet(chain(self._items, rhs._items)) 