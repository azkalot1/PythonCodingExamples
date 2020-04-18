class Trace:
    def __init__(self):
        self.enabled = True

    def __call__(self, f):
        def wrap(*args, **kwargs):
            if self.enabled:
                print('Calling {}'.format(f))
            return f(*args, **kwargs)
        return wrap


def escape_unicode(f):
    # This is a decorator
    # It takes a function, returns a function
    # Just modifies it
    def wrap(*args, **kwargs):
        x = f(*args, **kwargs)
        return ascii(x)
        # How function should be decorates is defined here
    return wrap


tracer = Trace()


@tracer
@escape_unicode
def norwegian_island_maker(name):
    return name + 'kjøpstad'


class IslandMaker:
    # This is a class version of example
    # where we decorate class methods
    def __init__(self, suffix):
        self.suffix = suffix

    @tracer
    @escape_unicode
    def make_island(self, name):
        return name + self.suffix
