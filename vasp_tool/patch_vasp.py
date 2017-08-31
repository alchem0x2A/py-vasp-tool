#####################################################################
# The patcher for factory ase.calculator.vasp.Vasp class            #
# will change the behavior of the POSCAR writer to use vasp5 format #
#####################################################################

from ase.calculators.vasp.create_input import GenerateVaspInput
from ase.calculators.vasp.create_input import bool_keys, int_keys, float_keys
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp import Vasprun
import os
import os.path
import shutil
from ase.io import read

# Tell vasp calculator to write the POSCAR using vasp5 style


def _new_write_input(self, atoms, directory='./', direct=True, vasp5=True):
    from ase.io.vasp import write_vasp
    from os.path import join
    write_vasp(join(directory, 'POSCAR'),
               self.atoms_sorted,
               direct=direct,
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


def _copy_files(self, select_names=None,
                exclude_names=None,
                tag="tag"):
    # copy_file is supposed to be used only after the calculation!
    if hasattr(self, "tag"):
        tag = self.tag
    default_names = ["INCAR", "OUTCAR", "WAVECAR", "CONTCAR",
                     "WAVEDER", "DOSCAR", "vasprun.xml"]
    if exclude_names != None:
        tmp = [p for p in default_names if p not in exclude_names]
        default_names = tmp
    elif select_names != None:
        default_names = select_names

    for fname in default_names:
        if os.path.exists(fname):
            f_new = ".".join((fname, tag))
            shutil.copy(fname, f_new)

# Get the final potential from vasprun.xml


def _get_final_E(self, filename="vasprun.xml"):
    v = Vasprun(filename)
    fe = v.final_energy.real
    return fe

oldrun = Vasp.run
def _run(self):
    # Handle the incomplete BSE vasprun problem
    oldrun(self)
    if os.path.exists("vasprun.xml"):
        with open("vasprun.xml", "rb+") as f:
            f.seek(-12, 2)                        # to the -12
            s = f.read()
            s.decode("utf8")
            if s.strip() != "</modeling>":  # Not the last line
                f.seek(0, 2)                # To the last
                f.write(b"</modeling>\n")
                print("Warning! The vasprun.xml seems incomplete.")


# path for writing kpoints
# taken from jasp
def _write_kpoints(self, directory="", fname=None):
    """Write out the KPOINTS file.
    The KPOINTS file format is as follows:
    line 1: a comment
    line 2: number of kpoints
        n <= 0   Automatic kpoint generation
        n > 0    explicit number of kpoints
    line 3: kpt format
        if n > 0:
            C,c,K,k = cartesian coordinates
            anything else = reciprocal coordinates
        if n <= 0
            M,m,G,g for Monkhorst-Pack or Gamma grid
            anything else is a special case
    line 4: if n <= 0, the Monkhorst-Pack grid
        if n > 0, then a line per kpoint
    line 5: if n <=0 it is the gamma shift
    After the kpts may be tetrahedra, but we do now support that for
    now.
    """
    import numpy as np

    if fname is None:
        fname = os.path.join(directory, 'KPOINTS')

    p = self.input_params

    kpts = p.get('kpts', None)  # this is a list, or None
    # kpts_weight = p.get("kpts_weight", None)  # weights of the kpoints for BS

    if kpts is None:
        NKPTS = None
    elif len(np.array(kpts).shape) == 1:
        NKPTS = 0  # automatic
    else:
        NKPTS = len(p['kpts'])

    # figure out the mode
    if NKPTS == 0 and not p.get('gamma', None):
        MODE = 'm'  # automatic monkhorst-pack
    elif NKPTS == 0 and p.get('gamma', None):
        MODE = 'g'  # automatic gamma monkhorst pack
    # we did not trigger automatic kpoints
    elif p.get('kpts_nintersections', None) is not None:
        MODE = 'l'
    elif p.get('reciprocal', None) is True:
        MODE = 'r'
    else:
        MODE = 'c'

    with open(fname, 'w') as f:
        # line 1 - comment
        comm = 'KPOINTS created by Atomic Simulation Environment\n'
        if p.get("kpath", None) is not None:
            comm = "KPATH: {} \n".format(p.get("kpath", None))
        f.write(comm)
        # line 2 - number of kpts
        if MODE in ['c', 'k', 'm', 'g', 'r']:
            f.write('{}\n'.format(NKPTS))
        elif MODE in ['l']:  # line mode, default intersections is 10
            f.write('{}\n'.format(p.get('kpts_nintersections')))

        # line 3
        if MODE in ['m', 'g', 'l']:
            if MODE == 'm':
                f.write('Monkhorst-Pack\n')  # line 3
            elif MODE == 'g':
                f.write('Gamma\n')
            else:
                f.write("Line mode\n")
        elif MODE in ['c', 'k']:
            f.write('Cartesian\n')
        else:
            f.write('Reciprocal\n')

        # line 4
        if MODE in ['m', 'g']:
            f.write('{0:9f} {1:0.9f} {2:0.9f}\n'.format(*p.get('kpts', (1, 1, 1))))
        elif MODE in ['c', 'k', 'r']:
            for n in range(NKPTS):
                # I assume you know to provide the weights
                f.write('{0:9f} {1:9f} {2:9f} {3:4d}\n'.format(*p['kpts'][n]))
        elif MODE in ['l']:
            if p.get('reciprocal', None) is False:
                f.write('Cartesian\n')
            else:
                f.write('Reciprocal\n')
            for n in range(NKPTS):
                f.write('{0:9f} {1:9f} {2:9f}\n'.format(*p['kpts'][n]))

        # line 5 - only if we are in automatic mode
        if MODE in ['m', 'g']:
            if p.get('gamma', None):
                f.write('{0:9f} {1:9f} {2:9f}\n'.format(*p['gamma']))
            else:
                f.write('0.0 0.0 0.0\n')

# Patch method for get the atoms from previous calculation


def read_atoms_sorted(path=""):
    f_sort = os.path.join(path, 'ase-sort.dat')
    f_contcar = os.path.join(path, "CONTCAR")
    if os.path.isfile(f_sort):
        sort = []
        resort = []
        line = None
        with open(f_sort, 'r') as dat_sort:
            lines = dat_sort.readlines()
        for line in lines:
            data = line.split()
            sort.append(int(data[0]))
            resort.append(int(data[1]))
        atoms = read(f_contcar, format='vasp')[resort]
    else:
        atoms = read(f_contcar, format='vasp')
    return atoms


# Hot patch to the Vasp class
Vasp.read_bandgap = _read_bandgap
Vasp.load_vasprun = _load_vasprun
Vasp.read_extern_stress = _read_extern_stress
Vasp.copy_files = _copy_files
Vasp.read_final_E = _get_final_E
Vasp.write_kpoints = _write_kpoints
Vasp.run = _run

# Add missing keys
bool_keys += ["lusew",
              "ladder",
              "lhartree",
              "lpead",
              "lvdwexpansion",
              "lorbitalreal"]

int_keys += ["antires",
             "omegamax",
]

# Patching the vasprun.xml for BSE calculations
@property
def _converged_electronic(self):
    """
    Checks that electronic step convergence has been reached in the final
    ionic step
    """
    try:
        final_esteps = self.ionic_steps[-1]["electronic_steps"]
    except IndexError:
        return False            # no actual ionic steps
    if 'LEPSILON' in self.incar and self.incar['LEPSILON']:
        i = 1
        to_check = set(['e_wo_entrp', 'e_fr_energy', 'e_0_energy'])
        while set(final_esteps[i].keys()) == to_check:
            i += 1
        return i + 1 != self.parameters["NELM"]
    return len(final_esteps) < self.parameters["NELM"]

@property
def optical_transitions(self):
    # Get optical transitions of BSE calculation
    from xml.etree import ElementTree as ET
    import numpy
    ep = None
    for event, elem in ET.iterparse(self.filename):
        if ("name" in elem.attrib) and (elem.attrib["name"] == "opticaltransitions"):
            ep = elem
            break
    ot_array = []
    for v in ep:
        # print(v)
        ot_array.append(list(map(float, v.text.strip().split())))
    ot_array = numpy.array(ot_array)
    return ot_array

Vasprun.converged_electronic = _converged_electronic
Vasprun.optical_transitions = optical_transitions
