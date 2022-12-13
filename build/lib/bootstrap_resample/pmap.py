"""Parallel (multi-threaded) map function for python. 

Uses multiprocessing.Pool with error-resistant importing. There are two map
functions:

1) pmap(function, iterable) -> rapid fork-based multi-threaded map function.

2) low_memory_pmap(function, iterable) -> a more memory-efficient version
    intended for function calls that are individually long & memory-intensive.

"""
import os
import progressbar as pb
from time import sleep
import numpy as np
import pandas as pd
import traceback
import multiprocessing
from multiprocessing.pool import Pool

CPUs = multiprocessing.cpu_count()


def simple_pmap(func, Iter, processes=CPUs - 1):
    with Pool(processes=processes) as P:
        return P.map(func, Iter)


def low_memory_pmap(func, Iter, processes=int(round(CPUs / 2)), chunksize=1):
    with Pool(processes=processes) as P:
        return [result for result in P.imap(func, Iter)]


def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)
        except Exception:
            error(traceback.format_exc())
            raise
        return result


class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)


def pmap(func, Iter, processes=CPUs, nice=10, *args, **kwargs):
    os.nice(nice)
    with LoggingPool(processes=processes) as P:
        results = [P.apply_async(func, (it,) + args, kwargs) for it in Iter]
        maxval = len(results)
        bar = pb.ProgressBar(
            maxval=maxval,
            widgets=[pb.Percentage(), " ", pb.Bar("=", "[", "]"), " ", pb.ETA()],
        )
        while P._cache:
            bar.update(maxval - len(P._cache))
            sleep(2)
        bar.finish()
        P.close()
    return [r.get() for r in results]


def groupby_apply(gb, func):
    keys, values = zip(*gb)
    out = pmap(func, values)
    index = gb.size().index
    if not type(out[0]) in [pd.Series, pd.DataFrame, pd.Panel]:
        return pd.DataFrame(out, index=index)
    elif len(set(index.names) - set(out[0].index.names)) == 0:
        return pd.concat(out)
    return pd.concat(out, keys=keys, names=index.names)

