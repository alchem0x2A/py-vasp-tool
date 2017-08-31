################################################
# Wrapper for several Vasp-derived classes     #
# including: Relaxation, Dielectric, G0W0, BSE #
################################################

from ase.calculators.vasp import Vasp
from ase.units import GPa
from pymatgen.io.vasp import Vasprun, Outcar
from ase.dft.kpoints import special_points, bandpath
from .paramters import default_parameters
import numpy

# Some fix for importing
import sys
if sys.version_info < (3, 0):
    import patch_vasp
else:
    from . import patch_vasp


# Common class for Vasp from paramters
class VaspGeneral(Vasp):
    def __init__(self, restart=None, output_template="vasp",
                 track_output=True, profile=None, **kwargs):
        if profile is not None:
            if type(profile) == dict:
                dp = profile    # profile as dict
            elif profile in default_parameters:
                dp = default_parameters[profile]
            else:
                raise ValueError("Profile not recognized !")
        dp.update(**kwargs)
        self.profile = profile
        # update the paramter
        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **dp)
        
    def write_bs_kpoints(self, path_string,
                      intersections=40, lattice_type=None):
        if lattice_type is None:
            # TODO use build-in search method
            raise TypeError("The type of lattice cannot be None")
        else:
            self.input_params["kpath"] = path_string
            if self.profile == "bs_DFT":
                line_kpts = gen_line_path(path_string,
                                          lattice_type=lattice_type,
                                          n_int=None)
                self.input_params["kpts"] = line_kpts  # no weights
                self.input_params["kpts_nintersections"] = intersections
                self.input_params["reciprocal"] = True
                print(self.input_params["kpts_nintersections"])
                self.write_kpoints()
            elif self.profile == "bs_hybrid":
                line_kpts = gen_line_path(path_string,
                                          lattice_type=lattice_type,
                                          n_int=intersections)
                # read the ibzk
                kpoints = []
                with open("IBZKPT", "r") as f_ibz:
                    f_ibz.readline()
                    n = int(f_ibz.readline().strip())
                    f_ibz.readline()
                    for i in range(n):
                        l = f_ibz.readline().strip().split()
                        pts = []
                        for j in range(3):
                            pts.append(float(l[j]))
                        pts.append(int(l[-1]))
                        kpoints.append(pts)
                for l in line_kpts:
                    kpoints.append(l + [0])
                self.input_params["kpts"] = kpoints
                self.input_params["reciprocal"] = True
                self.write_kpoints()
            
# Generate Kpoint path using symmetry
def gen_line_path(kpath, lattice_type, n_int=None):
    # n_interp = None for generation of line mode,
    # Otherwise for explicit kpoints
    valid_pts = special_points[lattice_type]
    kpoints = []
    prev = None
    for i, p in enumerate(kpath):
        cood = valid_pts[p]
        # Double the kpoints for intermediate points
        if i == 0:
            kpoints.append(cood)
            prev = cood
        else:
            if n_int is None:  # line mode
                kpoints.append(cood)
                if i != len(kpath) - 1:
                    kpoints.append(cood)
            else:
                n_int = int(n_int)
                assert n_int >= 0
                start = numpy.array(prev)
                interv = (numpy.array(cood) - start) / (n_int + 1)
                for i in range(n_int + 1):
                    pts = start + (i + 1) * interv
                    kpoints.append(pts.tolist())
                prev = cood
    # print(kpoints)
    return kpoints           
        
    
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
                          # "algo": "e",  # Diagonization
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
                          "ediff": 1e-6,  # energy difference
                          "ediffg": 1e-5,
                          "prec": "Accurate",  # precision
                          "lwave": True,      # Do not store wave function
                          "lcharg": False,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 500,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": -1,  # do relaxation
                          "xc": "pbe0",  # use PBE
                          "algo": "Damped",  # or algo = D, not Diag!
                          "precfock": "fast",  # use fast fft grid for fock matrx
        }
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)

class VaspGW(Vasp):         # GW mode
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
                          "ediff": 1e-6,  # energy difference
                          "ediffg": 1e-5,
                          "prec": "Accurate",  # precision
                          "lwave": True,      # Do not store wave function
                          "lcharg": False,    # Do not store the charge density
                          "lvtot": False, # Do not store the local potential
                          "encut": 800, # energy cutoff
                          "nelm": 500,  # max SC steps
                          "nelmin": 4,  # min SC steps
                          "ibrion": -1,  # do relaxation
                          "xc": "pbe0",  # use PBE
                          "algo": "Damped",  # or algo = D, not Diag!
                          "precfock": "fast",  # use fast fft grid for fock matrx
        }
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params)
        
if __name__ == "__main__":
    pass
