# Example of enclosure
message = 'global'


def enclosing():
    message = 'enclosing'

    def local():
        nonlocal message
        # nonlocal no_such_name
        # will fail, there is no such variable
        message = 'local'
    print('enclosing message', message)
    local()
    print('enclosing message', message)


print('enclosing message', message)
enclosing()
print('enclosing message', message)
