import sys
if sys.version_info < (3, 0):
    import patch_vasp
else:
    from . import patch_vasp
from .other_vasp import VaspRelax
