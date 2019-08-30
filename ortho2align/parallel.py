from concurrent.futures import ProcessPoolExecutor
from collections import deque
from tqdm import tqdm


def apply(func, arg):
    return func(arg)


def starapply(func, args):
    return func(*args)


class NonExceptionalProcessPool:

    def __init__(self, max_workers, verbose=True, suppress_exceptions=False):
        self.executor = ProcessPoolExecutor(max_workers=max_workers)
        self.max_workers = max_workers
        self.verbose = verbose
        self.suppress_exceptions = suppress_exceptions

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.executor.shutdown(wait=True)
        return False

    def shutdown(self, wait=True):
        self.executor.shutdown(wait=wait)

    def _map_async(self, func, iterable, asynchronous=True, applier=apply):
        if self.verbose:
            pbar = tqdm()
        futures = deque()
        try:
            collection = (item for item in iterable)
        except TypeError:
            raise ValueError("Provided iterable is not iterable.")
        for _ in range(self.max_workers):
            try:
                futures.append(self.executor.submit(applier, func, next(collection)))
            except StopIteration:
                break
        results = list()
        exceptions = list()
        while futures:
            future = futures.popleft()
            if future.done():
                pbar.update()
                try:
                    futures.append(self.executor.submit(applier, func, next(collection)))
                except StopIteration:
                    pass
                if future.exception():
                    if not self.suppress_exceptions:
                        exceptions.append(future.exception())
                    continue
                results.append(future.result())
            else:
                if asynchronous:
                    futures.append(future)
                else:
                    futures.appendleft(future)
        return results, exceptions

    def map(self, func, iterable):
        return self._map_async(func, iterable,
                               asynchronous=False, applier=apply)

    def starmap(self, func, iterable):
        return self._map_async(func, iterable,
                               asynchronous=False, applier=starapply)

    def map_async(self, func, iterable):
        return self._map_async(func, iterable,
                               asynchronous=True, applier=apply)

    def starmap_async(self, func, iterable):
        return self._map_async(func, iterable,
                               asynchronous=True, applier=starapply)
