import numpy as np

from ase.calculators.calculator import Calculator
#from ase import Atoms
#from ase.data import atomic_numbers
#from ase.utils import IOContext
#from ase.geometry import get_distances
#from ase.cell import Cell
from ase.calculators.qmmm import SimpleQMMM
from ase.io import read
import sys
from mm_lammps import reset_positions, unit_cell_coords

np.set_printoptions(threshold=sys.maxsize)


def get_qm_gau_calc_cm5_charges(file, atoms):
    qm_output = open(file, 'r')
    for l_no, line in enumerate(qm_output):
        if "Hirshfeld charges, spin densities, dipoles, and CM5 charges" in line:
            break
    ch = np.loadtxt(fname = file, skiprows = l_no+2, max_rows = len(atoms), usecols = 7)
    return ch


def modify_forces(forces, nfrozen):
    for i in range(len(forces)):
        if i not in nfrozen:
            forces[i] = [0]*len(forces[i])
    return forces


def modify_link_atoms_positions(atoms):
    for i in range(len(atoms.nn_lnk_atms_idx)):
        if atoms.arrays["atomtypes"][atoms.nn_lnk_atms_idx[i]].startswith("C") and atoms.arrays["atomtypes"][atoms.nn_lnk_atms_idx_2[i]].startswith("C"):
            scale = 0.74
            atoms.positions[atoms.lnk_atms_idx[i]][0] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][0] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][0]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][0])
            atoms.positions[atoms.lnk_atms_idx[i]][1] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][1] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][1]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][1])
            atoms.positions[atoms.lnk_atms_idx[i]][2] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][2] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][2]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][2])
        else:
            continue
    return atoms


class SmpQMMM(SimpleQMMM):
 # want to make a derived class of SimpleQMMM to run QMMM calcs using Gaussian for QM and LAMMPS for MM
 # changes which are made:
 # 1) The conventional SimpleQMMM runs the QM calculations 2 times, one to compute QM energy while the other is to compute forces,
 #    so basically QM calculations are run twice, which can be very time consuming. While computing QM forces, Gaussian also computes
 #    energy, so the first modification made is to extract both energy and force information from a sigle Gaussian (QM) run.
 # 2) We need partial charges because they are used to compute MM non bonded interactions. Gaussian calculations of the QM subsystem
 #    yields that information (in the form of Mulliken charges, not accurate) but the conventional SimpleQMMM calculator does not output 
 #    them, so the second modification made is to read those data from the QM output file and incorporate them in LAMMPS (MM) calculations.

    """Simple QMMM calculator derived from parent SimpleQMMM specifically for Gaussian (QM) and LAMMPS (MM)"""

    implemented_properties = ['energy', 'forces']

    def __init__(self, selection, qmcalc, mmcalc1, mmcalc2, vacuum=None, **kwargs):

        SimpleQMMM.__init__(self, selection, qmcalc, mmcalc1, mmcalc2, vacuum=None)
        self.__dict__.update(kwargs)


    def calculate(self, atoms, properties, system_changes):
        frozen = []
        atoms.wrap(eps=1e-15)
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.qmatoms is None:
            self.initialize_qm(atoms)

        if hasattr(atoms, "nn_lnk_atms_idx") and hasattr(atoms, "nn_lnk_atms_idx_2"):
            atoms = modify_link_atoms_positions(atoms) # changes link atoms coordinates
        self.qmatoms.positions = atoms.positions[self.selection]
        self.qmatoms.cell = atoms.cell
        self.qmatoms.pbc = False # setting pbc false for qm atoms to performm qm calculations
        #mm1_energy = self.mmcalc1.get_potential_energy(self.qmatoms)
        if len(self.qmatoms) > 0: 
            reset_positions(self.qmatoms, self.qmatoms.get_scaled_positions()[0]) # fold the qm atoms so that they form a complete structure    
        if self.vacuum:
            self.qmatoms.positions += (self.center -
                                       self.qmatoms.positions.mean(axis=0))

        #energy = self.qmcalc.get_potential_energy(self.qmatoms)
        self.mmcalc2.label = atoms.get_chemical_formula()
        self.mmcalc1.label = self.qmatoms.get_chemical_formula()
        qmforces = self.qmcalc.get_forces(self.qmatoms)
        qm_output = read(self.qmcalc.label + '.log', format='gaussian-out')
        qm_charges = get_qm_gau_calc_cm5_charges(self.qmcalc.label + '.log', self.qmatoms)
        atoms.arrays["charges"][self.selection] = qm_charges
        if hasattr(atoms, "lnk_atms_idx"):
            #offset = np.sum(atoms.charges[atoms.lnk_atms_idx])/len(atoms.lnk_atms_idx)
            #atoms.charges[atoms.nn_lnk_atms_idx]+=offset
            atoms.arrays["charges"][atoms.nn_lnk_atms_idx]+=atoms.arrays["charges"][atoms.lnk_atms_idx]
            atoms.arrays["charges"] = np.array([0 if i in atoms.lnk_atms_idx else atoms.arrays["charges"][i] for i in range(len(atoms))])      
        self.qmatoms.charges = qm_charges
        self.qmcalc.calc = qm_output.calc
        self.qmcalc.results = qm_output.calc.results 
        open_file = open('energy-and-forces-data.xyz', "a+")
        energy = self.qmcalc.results["energy"] # gaussian QM energy 
        print("QM Energy from Gaussian of QM subsystem = {}". format(self.qmcalc.results["energy"]), file=open_file)
        mm2_energy = self.mmcalc2.get_potential_energy(atoms)
        print("MM Energy from LAMMPS of entire system = {}". format(mm2_energy), file=open_file)
        energy += mm2_energy
        mm2_forces = self.mmcalc2.get_forces(atoms)
        forces = mm2_forces

        if self.vacuum:
            qmforces -= qmforces.mean(axis=0)
        forces[self.selection] += qmforces

        if hasattr(self,"rt"): # to store the optimization calculation type data (if exists) in the qmatoms object
            self.qmatoms.rt = self.rt

        self.qmatoms.pbc = True # resetting pbc true for qm atoms to performm mm calculations
        self.qmatoms.positions = atoms.positions[self.selection] # resetting the positions back after changes made to the positions during qm calculation
        mm1_energy = self.mmcalc1.get_potential_energy(self.qmatoms)
        mm2_forces = self.mmcalc2.get_forces(atoms) # need to reset mm2_forces, do not understand why it gets altered (need to check)
        print("MM Energy from LAMMPS of QM subsystem = {}". format(mm1_energy), file=open_file)
        energy -= mm1_energy
        print("Total Energy = {}". format(energy), file=open_file)
        mm1_forces = self.mmcalc1.get_forces(self.qmatoms)
        forces[self.selection] -= mm1_forces
        print("QM Forces from Gaussian of QM subsystem = ", file=open_file)
        open_file.writelines(str(qmforces))
        print("", file=open_file)
        print("MM Forces from LAMMPS of entire system = ", file=open_file)
        open_file.writelines(str(mm2_forces))
        print("", file=open_file)
        print("MM Forces from LAMMPS of QM subsystem = ", file=open_file)
        open_file.writelines(str(mm1_forces))
        print("", file=open_file)
        print("Total Forces = ", file=open_file)
        open_file.writelines(str(forces))
        print("", file=open_file)
        print("############################################################################################################################################", file=open_file)
        print("############################################################################################################################################", file=open_file)
        print("", file=open_file)
        open_file.close()
        forces = modify_forces(forces, self.nfrz_idx)
        self.results['energy'] = energy
        self.results['forces'] = forces




