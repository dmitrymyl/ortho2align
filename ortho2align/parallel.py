from concurrent.futures import ProcessPoolExecutor
from collections import deque
from functools import partial
from tqdm import tqdm
from pebble import ProcessPool


def apply(func, arg):
    """Applies one arg to the func.

    Args:
        func (function): a function to apply an arg to.
        arg (any): the argument to apply to the function.

    Returns:
        Result of the function with the applied arg.
    """
    return func(arg)


def starapply(func, args):
    """Applies starred args to the func.

    Args:
        func (function): a function to apply args to.
        args (collection): a collection of args.

    Returns:
        Result of the function with args applied.
    """
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
            suppress_exceptions (bool): if True, exceptions won't be
                recordered (default: False).
        """
        self.executor = ProcessPoolExecutor(max_workers=max_workers)
        self.max_workers = max_workers
        self.verbose = verbose
        self.suppress_exceptions = suppress_exceptions
        self._pbar = None

    @property
    def pbar(self):
        """a tqdm progress bar."""
        if self._pbar is None:
            self._pbar = tqdm(desc=f'Running in parallel with {self.max_workers} workers',
                              unit='task',
                              disable=not self.verbose)
        return self._pbar

    def __enter__(self):
        """Context manager enter.

        Notes:
            Progress bar is initialized via the context
            manager mechanism in order to correctly close
            pbar in case of exiting with an error.
        """
        pbar = self.pbar
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.pbar.close()
        self.shutdown(wait=True)
        self._pbar = None
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
                Can't accept lambdas.
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
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        futures = deque()
        exec_func = partial(applier, func)
        try:
            collection = (item for item in iterable)
        except TypeError:
            raise ValueError("Provided iterable is not iterable.")
        for _ in range(self.max_workers):
            try:
                futures.append(self.executor.submit(exec_func,
                                                    next(collection)))
            except StopIteration:
                break
        results = list()
        exceptions = list()
        while futures:
            future = futures.popleft()
            if future.done():
                self.pbar.update()
                try:
                    futures.append(self.executor.submit(exec_func,
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
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are returned in the same order as
                corresponding arguments were passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=False, applier=apply)

    def starmap(self, func, iterable):
        """Maps many-argument function to an iterable.

        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are returned in the same order as
                corresponding arguments were passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=False, applier=starapply)

    def map_async(self, func, iterable):
        """Maps one-argument function to an iterable in asynchronous way.

        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are not guaranteed to be returned
                in the same order as corresponding arguments were
                passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=True, applier=apply)

    def starmap_async(self, func, iterable):
        """Maps multi-argument function to an iterable in asynchronous way.

        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are not guaranteed to be returned
                in the same order as corresponding arguments were
                passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map_async(func, iterable,
                               asynchronous=True, applier=starapply)


class TimeoutProcessPool:
    """Provides multiprocessing with no exceptions invoked
    and allows for cancelling tasks after timeout.
    Based on pebble.ProcessPool mimics multiprocessing.Pool API
    for mapping function to an iterable containing arguments.
    It also mimics pebble API in terms of timeout supply.

    Attributes:
        executor (ProcessPoolExecutor): multiprocessing instance.
        max_workers (int): number of max_workers for executor
            (default: None, i.e. os.cpu_count()).
        verbose (bool): if True will show progress via tqdm
            (default: True).
        suppress_exceptions: if True, exceptions won't be
            recordered (default: False).

    Notes:
        There are 2 methods for multiprocessing: `map`, `starmap`.
        `map` allows multiprocessing of one-argument function execution.
        `starmap` allows multiprocessing of many-arguments function execution.
        Each method takes function and iterable containing arguments one argument
        (or group of arguments) per item. Iterable can be a collection,
        generator and coroutine as well. `map` and `starmap`
        return results and exceptions in the same order as arguments
        were passed by an iterable.

    Usage:
        def foo(x):
            for _ in range(2000000):
                continue
            if x % 10 == 0:
                raise ValueError(f"{x} is divided by ten.")
            return x
        with TimeoutProcessPool(max_workers=5) as p:
            results, exceptions = p.map(foo, range(100), timeout=None)

    Explanation of usage:
        All results returned by the functions are stored in
        the results list. All the exceptions are stored in
        the exceptions list unless `suppress_exceptions` is
        True.
    """

    def __init__(self,
                 max_workers=None,
                 verbose=True,
                 suppress_exceptions=False):
        """Initializes TimeoutProcessPool instance.
        Args:
            max_workers (int): number of max_workers for executor
                (default: None, i.e. os.cpu_count()).
            verbose (bool): if True will show progress via tqdm
                (default: True).
            suppress_exceptions: if True, exceptions won't be
                recordered (default: False).
        """
        self.executor = ProcessPool(max_workers=max_workers)
        self.max_workers = max_workers
        self.verbose = verbose
        self.suppress_exceptions = suppress_exceptions
        self._pbar = None

    @property
    def pbar(self):
        if self._pbar is None:
            self._pbar = tqdm(desc=f'Running in parallel with {self.max_workers} workers',
                              unit='task',
                              disable=not self.verbose)
        return self._pbar

    def __enter__(self):
        """Context manager enter.

        Notes:
            Progress bar is initialized via the context
            manager mechanism in order to correctly close
            pbar in case of exiting with an error.
        """
        pbar = self.pbar
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.pbar.close()
        self.executor.close()
        self.executor.join()
        self._pbar = None
        return False

    def _map(self, func, iterable, timeout=None, applier=apply):
        """Helper function to sustain API.
        If `self.verbose`is True, than the progress will be
        reported with tqdm.tqdm progress bar.
        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.
            timeout (int): time in seconds to allow for process to
                execute, then kill if not completed. If None,
                waits until complete (default: None).
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
        exec_func = partial(applier, func)
        results, exceptions = list(), list()

        try:
            future = self.executor.map(exec_func,
                                       iterable,
                                       timeout=timeout)
        except TypeError:
            raise ValueError("Provided iterable is not iterable.")

        iterator = future.result()

        while True:
            try:
                result = next(iterator)
                results.append(result)
                self.pbar.update()
            except StopIteration:
                self.pbar.close()
                break
            except Exception as error:
                if not self.suppress_exceptions:
                    exceptions.append(error)
                self.pbar.update()
        return results, exceptions

    def map(self, func, iterable, timeout=None):
        """Maps one-argument function to an iterable.

        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.
            timeout (int): time in seconds to allow for process to
                execute, then kill if not completed. If None,
                waits until complete (default: None).

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are returned in the same order as
                corresponding arguments were passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map(func, iterable, timeout=timeout, applier=apply)

    def starmap(self, func, iterable, timeout=None):
        """Maps many-argument function to an iterable.

        Args:
            func (function): function object to multiprocess.
                Can't accept lambdas.
            iterable (iterable): collection, range, map, zip,
                generator, coroutine object containing arguments
                to pass to func. The iterable is guaranteed to
                be not transformed to list so it is safe to pass
                generators and similar objects.
            timeout (int): time in seconds to allow for process to
                execute, then kill if not completed. If None,
                waits until complete (default: None).

        Returns:
            tuple: two lists. The first contains results
                of successful call, the second contains exceptions
                if `self.suppress_exceptions` is False. Results
                and exceptions are returned in the same order as
                corresponding arguments were passed via iterable.

        Raises:
            ValueError: in case provided iterable is not iterable.
        """
        return self._map(func, iterable, timeout=timeout, applier=starapply)


class ExceptionLogger(Exception):
    """
    Logs exception in a human-readable way.

    Attributes:
        exception (Exception): raised exception.
        variable (any type): any important variable
            to log with the exception.
        message (str): some explanatory message (default: "").
    """

    def __init__(self, exception, variable, message=""):
        """Initializes an ExceptionLogger isntance.

        Args:
            exception (Exception): raised exception.
            variable (any type): any important variable
                to log with the exception.
            message (str): some explanatory message (default: "").

        Returns:
            None
        """
        self.exception = exception
        self.variable = variable
        self.message = message

    def __str__(self):
        return f'The exception "{self.exception}" has been raised. ' \
               f'Important variable is {self.variable}. {self.message}'
