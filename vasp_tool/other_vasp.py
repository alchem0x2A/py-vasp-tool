################################################
# Wrapper for several Vasp-derived classes     #
# including: Relaxation, Dielectric, G0W0, BSE #
################################################

from ase.calculators.vasp import Vasp
from ase.units import GPa
from pymatgen.io.vasp import Vasprun, Outcar
from ase.dft.kpoints import ibz_points, get_bandpath

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
                          "isif": 3,    # ionic + position
                          "nsw": 200,   # numbers for relaxation
                          "xc": "pbe",  # use PBE
        }
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)

class VaspGround(Vasp):         # Calculate the ground state, always restart
    def __init__(self, restart=False, output_template="vasp",
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
                          "lwave": True,      # Do not store wave function
                          "lcharg": True,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 500,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": -1,  # do relaxation
                          "xc": "pbe",  # use PBE
        }
        # You should provide the kpoints line for yourself!
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)

# Generate Kpoint path using symmetry
def gen_line_path(kpath, lattice_type):
    valid_pts = ibz_points[lattice_type]
    kpoints = []
    for i, p in enumerate(kpath):
        if p == "G":
            p = "Gamma"
        cood = valid_pts[p]
        # Double the kpoints for intermediate points
        if (i == 0) or (i == len(kpath) - 1):
            kpoints.append(cood)
        else:
            [kpoints.append(cood) for i in range(2)]
    # print(kpoints)
    return kpoints

    
class VaspBandStructure(Vasp):         # Calculate the ground state, always restart
    def __init__(self, restart=False, output_template="vasp",
                 track_output=False, kpath=None,
                 lattice_type=None, **kwargs):
        # First generate using normal Vasp class taking default params
        # Overwrite the parameters by user default
        default_params = {"istart": 1,  # Start from scratch
                          "icharg": 11,  # Read the charge density
                          "amix": 0.1,   # mixing parameters
                          "bmix": 0.01,
                          "ismear": 0,  # Gaussian smear
                          "sigma": 0.01,  # smearing strength
                          "ediff": 1e-8,  # energy difference
                          "ediffg": 1e-7,
                          "prec": "Accurate",  # precision
                          "lwave": True,      # Do not store wave function
                          "lcharg": True,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 500,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": -1,  # do relaxation
                          "xc": "pbe0",
                          "lorbit": 11,   # write DOSCAR and PROCAR
                          "algo": "exact",  # Diagonization
                          "kpts_nintersections": 10,
                          "reciprocal": True,
        }
        for key in kwargs:
            default_params[key] = kwargs[key]

        if (kpath is not None) and (lattice_type is not None):
            kpts = gen_line_path(kpath, lattice_type)
            default_params["kpts"] = kpts
            

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)
        self.kpath = kpath
        self.lattice_type = lattice_type

class VaspHybridBandgap(Vasp):         # DFT+HF Bandgap
    def __init__(self, restart=False, output_template="vasp",
                 track_output=False, **kwargs):
        # First generate using normal Vasp class taking default params
        # Overwrite the parameters by user default
        default_params = {"istart": 1,  # Start from previous
                          "icharg": 2,  # should be 2 for relaxation
                          "amix": 0.1,   # mixing parameters
                          "bmix": 0.01,
                          "ismear": 0,  # Gaussian smear
                          "sigma": 0.01,  # smearing strength
                          "ediff": 1e-8,  # energy difference
                          "ediffg": 1e-7,
                          "prec": "Accurate",  # precision
                          "lwave": True,      # Do not store wave function
                          "lcharg": True,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 500,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": -1,  # do relaxation
                          "xc": "pbe0",  # use PBE
                          "algo": "exact",  # or algo = D
                          "precfock": "fast",  # use fast fft grid for fock matrx
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
