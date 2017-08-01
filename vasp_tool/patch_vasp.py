#####################################################################
# The patcher for factory ase.calculator.vasp.Vasp class            #
# will change the behavior of the POSCAR writer to use vasp5 format #
#####################################################################

from ase.calculators.vasp.create_input import GenerateVaspInput
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp import Vasprun
import os, os.path, shutil

# Tell vasp calculator to write the POSCAR using vasp5 style
def _new_write_input(self, atoms, directory='./', vasp5=True):
        from ase.io.vasp import write_vasp
        from os.path import join
        write_vasp(join(directory, 'POSCAR'),
                   self.atoms_sorted,
                   symbol_count=self.symbol_count, vasp5=vasp5)
        self.write_incar(atoms, directory=directory)
        self.write_potcar(directory=directory)
        self.write_kpoints(directory=directory)
        self.write_sort_file(directory=directory)

# Hot patch for the GenerateVaspInput class
GenerateVaspInput.write_input = _new_write_input


def _load_vasprun(self, filename="vasprun.xml"):
        self.vasprun = Vasprun(filename)

    # read the bandgap from vasprun.xml
def _read_bandgap(self):
        if not hasattr(self, "vasprun"):
            self.load_vasprun()
        # From DOS
        dos = self.vasprun.complete_dos
        bg_dos = dos.get_gap()
        # From Band structure
        bs = self.vasprun.get_band_structure()
        bg_bs = bs.get_band_gap()
        # Return the bandgaps calculated by DOS or band structure
        return (bg_dos, bg_bs)

def _read_extern_stress(self, form="kB", filename="OUTCAR"):
    stress = None
    for line in open(filename):
        if line.find('external pressure') != -1:
            stress = line.split()[3]
            if form != "kB":
                # in GPa
                stress = stress * 0.1 * GPa
    return stress

def _copy_files(self, tag="tag"):
    # copy_file is supposed to be used only after the calculation!
    if hasattr(self, "tag"):
        tag = self.tag
    
    for fname in ["INCAR", "OUTCAR", "WAVECAR", "CONTCAR"
                  "WAVEDER", "DOSCAR", "vasprun.xml"]:
        if os.path.exists(fname):
            f_new = ".".join((fname, tag))
            shutil.copy(fname, f_new)

# Hot patch to the Vasp class
Vasp.read_bandgap = _read_bandgap
Vasp.load_vasprun = _load_vasprun
Vasp.read_extern_stress = _read_extern_stress
Vasp.copy_files = _copy_files

