from concurrent.futures import ProcessPoolExecutor
from collections import deque
from tqdm import tqdm


def apply(func, arg):
    """Applies one arg to the func."""
    return func(arg)


def starapply(func, args):
    """Applies starred args to the func."""
    return func(*args)


class NonExceptionalProcessPool:
    """Provides multiprocessing with no exceptions invoked.

    Based on concurrent.futures.ProcessPoolExecutor
    mimics multiprocessing.Pool API for mapping
    function to an iterable containing arguments.

    Attributes:
        executor (ProcessPoolExecutor): multiprocessing instance.
        max_workers (int): number of max_workers for executor
            (default: None, i.e. os.cpu_count()).
        verbose (bool): if True will show progress via tqdm
            (default: True).
        suppress_exceptions: if True, exceptions won't be
            recordered (default: False).

    Notes:
        There are 4 methods for multiprocessing: `map`, `starmap`,
        `map_async`, `starmap_async`. `map` and `map_async` allows
        multiprocessing of one-argument function execution.
        `starmap` and `starmap_async` allows multiprocessing of
        many-arguments function execution. Each method takes
        function and iterable containing arguments one argument
        (or group of arguments) per item. Iterable can be a collection,
        generator and coroutine as well. `map` and `starmap`
        return results and exceptions in the same order as arguments
        were passed by an iterable. Async versions of these functions
        do not guarantee that.

    Usage:
        def foo(x):
            for _ in range(2000000):
                continue
            if x % 10 == 0:
                raise ValueError(f"{x} is divided by ten.")
            return x

        with NonExceptionalProcessPool(max_workers=5) as p:
            results, exceptions = p.map(foo, range(100))

    Explanation of usage:
        All results returned by the functions are stored in
        the results list. All the exceptions are stored in
        the exceptions list unless `suppress_exceptions` is
        True.
    """

    def __init__(self, max_workers=None, verbose=True,
                 suppress_exceptions=False):
        """Initializes NonExceptionalProcessPool instance.

        Args:
            max_workers (int): number of max_workers for executor
                (default: None, i.e. os.cpu_count()).
            verbose (bool): if True will show progress via tqdm
                (default: True).
            suppress_exceptions: if True, exceptions won't be
                recordered (default: False).
        """
        self.executor = ProcessPoolExecutor(max_workers=max_workers)
        self.max_workers = max_workers
        self.verbose = verbose
        self.suppress_exceptions = suppress_exceptions

    def __enter__(self):
        """Context manager enter."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.shutdown(wait=True)
        return False

    def shutdown(self, wait=True):
        """Explicitly shuts down self.executor."""
        self.executor.shutdown(wait=wait)

    def _map_async(self, func, iterable, asynchronous=True, applier=apply):
        """Helper function to sustain API.

        If `self.verbose`is True, than the progress will be
        reported with tqdm.tqdm progress bar.

        Args:
            func (function): function object to multiprocess.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.
            asynchronous (bool): whether to process calls
                asynchronously or not (default: True).
            applier (function): how to apply arguments in
                iterable to func. One of `apply`, `starapply`
                (default: `apply`).

        Returns:
            (list, list) two lists. The first contains results
            of successful call, the second contains exceptions
            if `self.suppress_exceptions` is False.

        Raises:
            ValueError in case provided iterable is not iterable.
        """
        if self.verbose:
            pbar = tqdm()
        futures = deque()
        try:
            collection = (item for item in iterable)
        except TypeError:
            raise ValueError("Provided iterable is not iterable.")
        for _ in range(self.max_workers):
            try:
                futures.append(self.executor.submit(applier,
                                                    func,
                                                    next(collection)))
            except StopIteration:
                break
        results = list()
        exceptions = list()
        while futures:
            future = futures.popleft()
            if future.done():
                pbar.update()
                try:
                    futures.append(self.executor.submit(applier,
                                                        func,
                                                        next(collection)))
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
        """Maps one-argument function to an iterable.

        Args:
            func (function): function object to multiprocess.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            (list, list) two lists. The first contains results
            of successful call, the second contains exceptions
            if `self.suppress_exceptions` is False. Results
            and exceptions are returned in the same order as
            corresponding arguments were passed via iterable.

        Raises:
            ValueError in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=False, applier=apply)

    def starmap(self, func, iterable):
        """Maps many-argument function to an iterable.

        Args:
            func (function): function object to multiprocess.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            (list, list) two lists. The first contains results
            of successful call, the second contains exceptions
            if `self.suppress_exceptions` is False. Results
            and exceptions are returned in the same order as
            corresponding arguments were passed via iterable.

        Raises:
            ValueError in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=False, applier=starapply)

    def map_async(self, func, iterable):
        """Maps one-argument function to an iterable in asynchronous way.

        Args:
            func (function): function object to multiprocess.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            (list, list) two lists. The first contains results
            of successful call, the second contains exceptions
            if `self.suppress_exceptions` is False. Results
            and exceptions are not guaranteed to be returned
            in the same order as corresponding arguments were
            passed via iterable.

        Raises:
            ValueError in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=True, applier=apply)

    def starmap_async(self, func, iterable):
        """Maps multi-argument function to an iterable in asynchronous way.

        Args:
            func (function): function object to multiprocess.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            (list, list) two lists. The first contains results
            of successful call, the second contains exceptions
            if `self.suppress_exceptions` is False. Results
            and exceptions are not guaranteed to be returned
            in the same order as corresponding arguments were
            passed via iterable.

        Raises:
            ValueError in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=True, applier=starapply)
