################################################
# Wrapper for several Vasp-derived classes     #
# including: Relaxation, Dielectric, G0W0, BSE #
################################################

from ase.calculators.vasp import Vasp
import patch_vasp
from ase.units import GPa
from pymatgen.io.vasp import Vasprun, Outcar


class VaspRelax(Vasp):
    default_params = {"istart": 0,  # Start from scratch
                      "icharge": 2,  # should be 2 for relaxation
                      "amix": 0.1,   # mixing parameters
                      "bmix": 0.01,
                      "ismear": 0,  # Gaussian smear
                      "sigma": 0.01,  # smearing strength
                      "ediff": 1e-8,  # energy difference
                      "ediffg": 1e-7,
                      "prec": "Accurate",  # precision
                      "lwave": False,      # Do not store wave function
                      "lcharge": False,    # Do not store the charge density
                      "lvtot": False, # Do not store the local potential
                      "encut": 800, # energy cutoff
                      "nelm": 200,  # max SC steps
                      "nelmin": 4,  # min SC steps
                      "ibrion": 2,  # do relaxation
                      "isif": 2,    # default setting for stress tensor (2)
                      "nsw": 50,   # numbers for relaxation
                      "xc": "pbe",  # use PBE
    }
    def __init__(self, restart=None, output_template="vasp",
                 track_output=False, **kwargs):
        # First generate using normal Vasp class taking default params
        # Overwrite the parameters by user default
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params,
        )



class VaspHSE(Vasp):
    default_params = {"istart": 0,  # Start from scratch
                      "icharge": 2,  # should be 2 for relaxation
                      "amix": 0.1,   # mixing parameters
                      "bmix": 0.01,
                      "ismear": 0,  # Gaussian smear
                      "sigma": 0.01,  # smearing strength
                      "ediff": 1e-8,  # energy difference
                      "ediffg": 1e-7,
                      "prec": "Accurate",  # precision
                      "lwave": False,      # Do not store wave function
                      "lcharge": False,    # Do not store the charge density
                      "lvtot": False, # Do not store the local potential
                      "encut": 800, # energy cutoff
                      "nelm": 200,  # max SC steps
                      "nelmin": 4,  # min SC steps
                      "ibrion": 2,  # do relaxation
                      "isif": 2,    # default setting for stress tensor (2)
                      "nsw": 50,   # numbers for relaxation
                      "xc": "pbe",  # use PBE
    }
    def __init__(self, restart=None, output_template="vasp",
                 track_output=False, **kwargs):
        # First generate using normal Vasp class taking default params
        # Overwrite the parameters by user default
        for key in kwargs:
            default_params[key] = kwargs[key]

        Vasp.__init__(self, restart=restart,
                      output_template=output_template,
                      track_output=track_output,
                      **default_params,
        )

    # def load_vasprun(self, filename="vasprun.xml"):
    #     self.vasprun = Vasprun(filename)

    # # read the bandgap from vasprun.xml
    # def read_bandgap(self):
    #     if not hasattr(self, "vasprun"):
    #         self.load_vasprun()
    #     # From DOS
    #     dos = self.vasprun.complete_dos
    #     bg_dos = dos.get_gap()
    #     # From Band structure
    #     bs = self.vasprun.get_band_structure()
    #     bg_bs = bs.get_band_gap()
    #     # Return the bandgaps calculated by DOS or band structure
    #     return (bg_dos, bg_bs)

    # def read_extern_stress(self, form="kB",
    #                        filename="OUTCAR"):
    #     stress = None
    #     for line in open(filename):
    #         if line.find('external pressure') != -1:
    #             stress = line.split()[3]
    #             if form != "kB":
    #                 # in GPa
    #                 stress = stress * 0.1 * GPa
    #     return stress

if __name__ == "__main__":
    # Test
    pass
