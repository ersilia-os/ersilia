# External imports
# ruff: noqa
import os
from ._version import __version__
import warnings

# Filter out some warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Disable GPU
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# Default variables
from .default import EOS, CONFIG_JSON, INSTALL_STATUS_FILE

# Config
if not os.path.exists(os.path.join(EOS, CONFIG_JSON)):
    from .utils.config import Checker

    Checker().config()

# Environmental variables
os.environ["EOS_HOME"] = EOS

# User profile
from .default import bashrc_cli_snippet

bashrc_cli_snippet(overwrite=False)


# Lazy re-exports of heavy library classes/utilities. Importing them eagerly
# pulled ~300 ms of transitive dependencies into every CLI startup. PEP 562
# __getattr__ defers each import until first access while preserving the
# `from ersilia import X` ergonomics.
_LAZY = {
    "logger": ("ersilia.utils.logging", "logger"),
    "throw_ersilia_exception": (
        "ersilia.utils.exceptions_utils.throw_ersilia_exception",
        "throw_ersilia_exception",
    ),
    "Config": ("ersilia.utils.config", "Config"),
    "ErsiliaBase": ("ersilia.core.base", "ErsiliaBase"),
    "ModelBase": ("ersilia.core.modelbase", "ModelBase"),
    "ErsiliaModel": ("ersilia.core.model", "ErsiliaModel"),
}


def __getattr__(name):
    if name in _LAZY:
        import importlib

        mod_name, attr = _LAZY[name]
        obj = getattr(importlib.import_module(mod_name), attr)
        globals()[name] = obj
        return obj
    raise AttributeError(f"module 'ersilia' has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY.keys()))


# Check status of installs
def check_install_status():
    fn = os.path.join(EOS, INSTALL_STATUS_FILE)
    if not os.path.exists(fn):
        status = None
    else:
        with open(fn, "r") as f:
            status = f.read().strip()
    results = {"install_status_file": fn, "status": status}
    return results


INSTALL_STATUS = check_install_status()["status"]
