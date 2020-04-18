def escape_unicode(f):
    # This is a decorator
    # It takes a function, returns a function
    # Just modifies it
    def wrap(*args, **kwargs):
        x = f(*args, **kwargs)
        return ascii(x)
        # How function should be decorates is defined here

    return wrap


@escape_unicode
def northern_city():
    return 'Nord-Jæren'
