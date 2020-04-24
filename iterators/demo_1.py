import timeit
import numpy as np

"""
Contains examples of iterators
and generators and their comparison
"""


def fibonacci(n):
    """
    Classic recursive fibonacci
    """
    assert n > 0, 'Negative n is not supported'
    if n == 1:
        return 1
    elif n == 2:
        return 2
    else:
        return fibonacci(n-1) + fibonacci(n-2)


def fibonacci_optimized(n):
    a = 1
    b = 2
    assert n > 0, 'Negative n is not supported'
    if n == 1:
        return a
    elif n == 2:
        return b
    else:
        for i in range(2, n):
            c = a + b
            a = b
            b = c
        return b


class FibonacciIterable:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        return(
            'FibonacciIterable({},{})'.format(
                self.start,
                self.end))

    def __iter__(self):
        return FibonacciIterator(self.start, self.end)


class FibonacciIterator:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.increment = 0

    def __next__(self):
        if self.start + self.increment > self.end:
            raise StopIteration()
        fibonacci_number = fibonacci(self.start + self.increment)
        self.increment += 1
        return fibonacci_number

    def __iter__(self):
        return self


class FibonacciGenerator:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        return(
            'FibonacciGenerator({},{})'.format(
                self.start,
                self.end))

    def __iter__(self):
        for i in range(self.start, self.end + 1):
            yield fibonacci(i)


class FibonacciGeneratorOptimized:
    """
    Here we use generator and optimized fidonacci
    yeild nmber by computing them
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        return(
            'FibonacciGeneratorOptimized({},{})'.format(
                                                self.start,
                                                self.end))

    def __iter__(self):
        for i in range(self.start, self.end + 1):
            yield fibonacci_optimized(i)
            # We use yeld so no need to worry about state


class FibonacciGeneratorLazyOptimized:
    """
    Here we precompute all fibonacci numbers
    and just yield them, but looking from dict
    Precomputing may take more time
    So it is unclear if fast lookup from dict
    worth it
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.fibonnaci_dict = dict({i: fibonacci_optimized(i)
                                    for i in range(self.start, self.end+1)})

    def __repr__(self):
        return(
            'FibonacciGeneratorLazyOptimized({},{})'.format(
                                                self.start,
                                                self.end))

    def __iter__(self):
        for i in range(self.start, self.end + 1):
            yield self.fibonnaci_dict[i]


class FibonacciIterableLazyOptimized:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.increment = 0

    def __repr__(self):
        return(
            'FibonacciIterableLazyOptimized({},{})'.format(
                                                self.start,
                                                self.end))

    def __iter__(self):
        return FibonacciIteratorLazyOptimized(self.start, self.end)


class FibonacciIteratorLazyOptimized:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.increment = 0
        self.fibonnaci_dict = dict({i: fibonacci_optimized(i)
                                    for i in range(self.start, self.end+1)})

    def __next__(self):
        if self.start + self.increment > self.end:
            raise StopIteration()
        fibonacci_number = self.fibonnaci_dict[self.start + self.increment]
        self.increment += 1
        # Here we use increment to check which value we need to return
        return fibonacci_number

    def __iter__(self):
        return self


class Call:
    def __init__(self, it, start, end):
        self.it = it
        self.start = start
        self.end = end

    def __call__(self, *args, **kwargs):
        for _ in self.it(self.start, end=self.end):
            pass


if __name__ == "__main__":
    start_fib = 2
    end_fib = 30
    fiter = FibonacciIterable
    fgen = FibonacciGenerator
    fgen_optim = FibonacciGeneratorOptimized
    fgen_optim_lazy = FibonacciGeneratorLazyOptimized
    fiter_optim_lazy = FibonacciIterableLazyOptimized
    iterable_list = [
        fiter,
        fgen,
        fgen_optim,
        fgen_optim_lazy,
        fiter_optim_lazy
        ]
    for iterable in iterable_list:
        print('Iterable {}:\n'.format(iterable(start_fib, end_fib)))
        times = timeit.repeat(
            Call(iterable, start_fib, end_fib),
            number=5,
            repeat=5)
        mean_time = np.mean(times)
        print(f'{mean_time}\n')

# FibonacciIterable(2,30): 2.28
# FibonacciGenerator(2,30): 2.291077638696879
# FibonacciGeneratorOptimized(2,30) 0.000121
# FibonacciGeneratorLazyOptimized(2,30) 0.000134
# FibonacciIterableLazyOptimized(2,30) 0.000179
