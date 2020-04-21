class SimpleList:
    def __init__(self, items):
        self._items = list(items)

    def add(self, item):
        self._items.append(item)

    def __getitem__(self, index):
        return self._items[index]

    def sort(self):
        self._items.sort()

    def __len__(self):
        return len(self._items)

    def __repr__(self):
        return "SimpleList({!r})".format(self._items)


class SortedList(SimpleList):
    def __init__(self, items=()):
        super().__init__(items)
        self.sort()

    def add(self, item):
        super().add(item)
        self.sort()

    def __repr__(self):
        return "SortedList({!r})".format(self._items)


class IntList(SimpleList):
    def __init__(self, items=()):
        for x in items:
            self._validate(x)
        super().__init__(items)

    @staticmethod
    def _validate(x):
        if not isinstance(x, int):
            raise TypeError('Only ints')

    def add(self, item):
        self._validate(item)
        super().add(item)

    def __repr__(self):
        return "IntList({!r})".format(self._items)

# If a class
# has multiple base classes + defines no initializer
# the initializer of the first base class is called

# To access the base classes use cls.__bases__
# method resolution order : MRO : cls.__mro__ \ cls.mro()
# ordering inheritance graph


class SortedIntList(IntList, SortedList):
    def __repr__(self):
        return "SortedIntList({!r})".format(self._items)
    # Will sort and reject non int values
