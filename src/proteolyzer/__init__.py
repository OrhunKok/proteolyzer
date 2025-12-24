from . import config
from . import utils
from . import transformers
from .utils import Data
import importlib

__all__ = ["config", "utils", "plots", "transformers", 'Data']

__lazy__ = ["aas", "cellenone", "plots", "unimod"]

def __getattr__(name):

    # If user accesses aas or unimod, load both at once
    if name in ("aas", "unimod"):
        for module_name in ("aas", "unimod"):
            if module_name not in globals():
                globals()[module_name] = importlib.import_module(f"{__name__}.{module_name}")
        return globals()[name]

    # Normal lazy-loading for other modules
    if name in __lazy__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    raise AttributeError(f"module {__name__} has no attribute {name}")