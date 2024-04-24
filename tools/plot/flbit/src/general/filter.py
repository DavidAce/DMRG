from collections import Iterable
from itertools import groupby
from numpy import ndarray


def all_equal(iterable):
    if any(isinstance(el, list) or isinstance(el, ndarray) for el in iterable):
        iterator = iter(iterable)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all((first == x).all() for x in iterator)
    else:
        g = groupby(iterable)
        return next(g, True) and not next(g, False)


def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def contains(key, filter, mode='any'):
    if not filter:
        return True
    if isinstance(filter, list) or isinstance(filter, dict):
        matchlist = []
        for fitem in filter:
            if isinstance(fitem, list) or isinstance(filter, dict):
                matchlist.append(any(f in key for f in flatten(fitem)))
            else:
                matchlist.append(fitem in key)
        if mode == 'all':
            return all(matchlist)
        elif mode == 'any':
            return any(matchlist)
        else:
            raise AssertionError("Unexpected argument [" + str(mode) + "]. Must be [all] or [any]")
    else:
        return filter in key
