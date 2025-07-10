from ... import __version__, logger


def ersilia_ali(verbose):
    """
    🦠 Welcome to Ersilia! 💊
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
