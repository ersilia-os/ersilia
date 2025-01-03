import importlib


def import_extra(mod):
    """Try to import a module, if not found return None"""
    try:
        return importlib.import_module(mod)
    except ImportError:
        return None
