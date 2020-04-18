import socket


class Resolver:
    """
    This is an example of callabale instance
    """
    def __init__(self):
        self._cache = {}

    def __call__(self, host):
        if host not in self._cache:
            self._cache[host] = socket.gethostbyname(host)
        return self._cache[host]

    def clear(self):
        self._cache.clear()

    def has_host(self, host):
        return host in self._cache


def sequence_class(immutable):
    """
    Return not an instance of tulpe/list
    for class itself,
    which can be called
    Example:
    seq = sequence_class(immutable=True)
    t = seq('Azaza')
    Note that seq can be called
     ecause it is an object of a class
    Calling a class invokes a sontructor
    """
    if immutable:
        cls = tuple
    else:
        cls = list
    return(cls)


def sequence_class_alterntaiv(immutable):
    """
    Return not an instance of tulpe/list
    for class itself,
    which can be called
    Example:
    seq = sequence_class(immutable=True)
    t = seq('Azaza')
    Note that seq can be called
     ecause it is an object of a class
    """
    return tuple if immutable else list


# Extended argument syntax

def hypervolume(length, *lengths):
    v = length
    for item in lengths:
        v *= item
    return v


def tag(name, **attributes):
    result = '<' + name
    for key, value in attributes.itms():
        result += ' {k}="{v}"'.format(k=key, v=str(value))
    result += '>'
    return result


def color(red, green, blue, **kwargs):
    print(red)
    print(green)
    print(blue)
    print(kwargs)

# we can use this construct : k = {'red': 25, .... 'alpha': 255}
# color(**k)
# We can use it to forward arguments to the underlying functions


def trace(f, *args, **kwargs):
    print(args)
    print(kwargs)
    result = f(*args, **kwargs)
    return result


def outer():
    x = 3
    # How does the a local uses binding to x?
    # Enclosing scope is gone after calling outer
    # Local functions forme closures
    # i = outer()
    # i(1) <- 4
    # i(10) <- 13

    def inner(y):
        return x + y
    return inner

# Usefull example of that: function factories
# This function will create functon to raise x to the
# power of exp


def raise_to(exp):
    def raise_to_exp(x):
        return pow(x, exp)
    return raise_to_exp


