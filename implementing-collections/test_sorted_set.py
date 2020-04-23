import unittest
from sorted_set import SortedSet
from collections.abc import (Container, Sized, Iterable, Sequence)
# Tests firts!


class TestConstruction(unittest.TestCase):
    def test_empty(self):
        s = SortedSet([])

    def test_from_sequence(self):
        s = SortedSet([7, 8, 3, 1])

    def test_with_duplicated(self):
        s = SortedSet([7, 8, 3, 1, 8, 7])

    def test_from_iterable(self):
        def gen_small_seq():
            yield 6
            yield 8
            yield 4
            yield 2
            yield 3

        g = gen_small_seq()
        s = SortedSet(g)

    def test_default_empty(self):
        s = SortedSet()


class TestContainerProtocol(unittest.TestCase):
    def setUp(self):
        self.s = SortedSet([6, 7, 3, 9])

    def test_positive_contained(self):
        self.assertTrue(6 in self.s)

    def test_negative_contained(self):
        self.assertFalse(2 in self.s)

    def test_positive_not_contained(self):
        self.assertTrue(5 not in self.s)

    def test_negative_not_contained(self):
        self.assertFalse(9 not in self.s)

    def test_protocol(self):
        self.assertTrue(issubclass(SortedSet, Container))


class TestSizedProtocol(unittest.TestCase):
    def test_empty(self):
        s = SortedSet()
        self.assertEqual(len(s), 0)

    def test_one(self):
        s = SortedSet([42])
        self.assertEqual(len(s), 1)

    def test_more_one(self):
        s = SortedSet([42, 13])
        self.assertEqual(len(s), 2)

    def test_with_dup(self):
        s = SortedSet([42, 42, 42])
        self.assertEqual(len(s), 1)

    def test_protocol(self):
        self.assertTrue(issubclass(SortedSet, Sized))


class TestIterableProtocol(unittest.TestCase):
    def setUp(self):
        self.s = SortedSet([7, 2, 1, 1, 9])

    def test_iter(self):
        i = iter(self.s)
        self.assertEqual(next(i), 1)
        self.assertEqual(next(i), 2)
        self.assertEqual(next(i), 7)
        self.assertEqual(next(i), 9)  # Because sorted
        self.assertRaises(StopIteration, lambda: next(i))

    def test_for_loop(self):
        index = 0
        expected = [1, 2, 7, 9]
        for item in self.s:
            self.assertEqual(item, expected[index])
            index += 1

    def test_for_loop_zip(self):
        expected = [1, 2, 7, 9]
        [self.assertEqual(exp, val) for exp, val in zip(expected, iter(self.s))]

    def test_protocol(self):
        self.assertTrue(issubclass(SortedSet, Iterable))


class TestSequenceProtocol(unittest.TestCase):
    def setUp(self):
        self.s = SortedSet([1, 4, 9, 13, 15])

    def test_index_zero(self):
        self.assertEqual(self.s[0], 1)

    def test_index_four(self):
        self.assertEqual(self.s[4], 15)

    def test_index_beyond(self):
        with self.assertRaises(IndexError):
            self.s[5]

    def test_index_minus_one(self):
        self.assertEqual(self.s[-1], 15)

    def test_index_minus_five(self):
        self.assertEqual(self.s[-5], 1)

    def test_index_beyond_beggining(self):
        with self.assertRaises(IndexError):
            self.s[-6]

    # The following set will fail in vanilla representation,
    # because base equality is for reference equality,
    # not for values. We need to overwrite the default
    # equality
    def test_slice_from_start(self):
        self.assertEqual(self.s[:3], SortedSet([1, 4, 9]))

    def test_slice_to_end(self):
        self.assertEqual(self.s[3:], SortedSet([13, 15]))

    def test_slice_empty(self):
        self.assertEqual(self.s[10:], SortedSet())

    def test_slice_middle(self):
        self.assertEqual(self.s[2:4], SortedSet([9, 13]))

    def test_slice_full(self):
        self.assertEqual(self.s[:], self.s)

    def test_reversed(self):
        s = SortedSet([1, 3, 5, 7])
        r = reversed(s)
        self.assertEqual(next(r), 7)
        self.assertEqual(next(r), 5)
        self.assertEqual(next(r), 3)
        self.assertEqual(next(r), 1)
        with self.assertRaises(StopIteration):
            next(r)

    def test_index_positive(self):
        s = SortedSet([1, 2, 3])
        self.assertEqual(s.index(1), 0)

    def test_index_negative(self):
        s = SortedSet([1, 2, 3])
        with self.assertRaises(ValueError):
            s.index(4)

    def test_count_one(self):
        s = SortedSet([1, 2, 3])
        self.assertEqual(s.count(1), 1)

    def test_count_zero(self):
        s = SortedSet([1, 2, 3])
        self.assertEqual(s.count(5), 0)

    def test_protocol(self):
        self.assertTrue(issubclass(SortedSet, Sequence))

    # Additional test: concat disjoint, interescting and equal sets
    # Implemented only one
    def test_concat_intersecting(self):
        s = SortedSet([1, 2, 3])
        t = SortedSet([3, 4, 5])
        self.assertEqual(s + t, SortedSet([1, 2, 3, 4, 5]))


class TestReplProtocol(unittest.TestCase):
    def test_repr_empty(self):
        s = SortedSet()
        self.assertEqual(repr(s), "SortedSet()")

    def test_repr_some(self):
        s = SortedSet([42, 40, 19])
        self.assertEqual(repr(s), "SortedSet([19, 40, 42])")


class TestEqualityProtocol(unittest.TestCase):
    def test_positive_equal(self):
        self.assertTrue(SortedSet([4, 5, 6]) == SortedSet([4, 5, 6]))

    def test_negative_equal(self):
        self.assertFalse(SortedSet([4, 5, 6]) == SortedSet([8, 5, 6]))

    def test_type_mismatch(self):
        self.assertFalse(SortedSet([4, 5, 6]) == [8, 5, 6])

    def test_identical(self):
        s = SortedSet([10, 11, 12])
        self.assertTrue(s == s)


class TestUnequalityProtocol(unittest.TestCase):
    def test_positive_unequal(self):
        self.assertTrue(SortedSet([4, 5, 6]) != SortedSet([1, 2, 3]))

    def test_negative_unequal(self):
        self.assertFalse(SortedSet([4, 5, 6]) != SortedSet([4, 5, 6]))

    def test_type_mismatch(self):
        self.assertTrue(SortedSet([4, 5, 6]) != [8, 5, 6])

    def test_identical(self):
        s = SortedSet([10, 11, 12])
        self.assertFalse(s != s)


if __name__ == '__main__':
    unittest.main()
