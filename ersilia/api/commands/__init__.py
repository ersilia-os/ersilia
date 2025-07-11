from ... import __version__ as __version__
from ... import logger


def ersilia_ali(verbose):
    """
    🦠 Welcome to Ersilia! 💊
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
