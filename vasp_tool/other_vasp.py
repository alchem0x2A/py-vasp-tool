################################################
# Wrapper for several Vasp-derived classes     #
# including: Relaxation, Dielectric, G0W0, BSE #
################################################

from ase.calculators.vasp import Vasp
from ase.units import GPa
from pymatgen.io.vasp import Vasprun, Outcar

# Some fix for importing
import sys
if sys.version_info < (3, 0):
    import patch_vasp
else:
    from . import patch_vasp


class VaspRelax(Vasp):
    def __init__(self, restart=None, output_template="vasp",
                 track_output=False, **kwargs):
        # First generate using normal Vasp class taking default params
        # Overwrite the parameters by user default
        default_params = {"istart": 0,  # Start from scratch
                          "icharg": 2,  # should be 2 for relaxation
                          "amix": 0.1,   # mixing parameters
                          "bmix": 0.01,
                          "ismear": 0,  # Gaussian smear
                          "sigma": 0.01,  # smearing strength
                          "ediff": 1e-8,  # energy difference
                          "ediffg": 1e-7,
                          "prec": "Accurate",  # precision
                          "lwave": False,      # Do not store wave function
                          "lcharg": False,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 200,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": 2,  # do relaxation
                          "isif": 2,    # default setting for stress tensor (2)
                          "nsw": 50,   # numbers for relaxation
                          "xc": "pbe",  # use PBE
        }
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)




if __name__ == "__main__":
    # Test
    pass
