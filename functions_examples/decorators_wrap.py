import functools
# examples of preserving metadata of functions
# after using decorators


def noop(f):
    @functools.wraps(f)
    def noop_wrapper(f):
        return f()
    return noop_wrapper


@noop
def hello():
    "Print a hello message"
    print('Hello!')
