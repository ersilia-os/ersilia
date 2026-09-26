import asyncio


def run_coroutine(coro):
    """
    Run a coroutine to completion from synchronous code.

    Uses ``asyncio.run`` when no event loop is running. Inside a running loop
    (e.g. a Jupyter notebook), ``asyncio.run`` is not allowed, so the loop is
    made re-entrant with ``nest_asyncio`` first. This keeps the patch out of
    processes that never need it, instead of applying it at import time.

    Parameters
    ----------
    coro : coroutine
        The coroutine to run.

    Returns
    -------
    Any
        The coroutine's result.
    """
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)
    import nest_asyncio

    nest_asyncio.apply(loop)
    return loop.run_until_complete(coro)
