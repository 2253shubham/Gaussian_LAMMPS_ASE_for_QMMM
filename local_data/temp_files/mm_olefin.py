from ase import Atoms
from ase.io import read, write
from ase.calculators.qmmm import SimpleQMMM
from ase.calculators.gaussian import Gaussian
from ase.calculators.gromacs import Gromacs
from mm_lammps import (init_atoms, set_calc_params, reset_positions) 
from mm_lammps import MM_LAMMPS 
import itertools as it
#####################################################################################
# Set-up molecular system, read input files for lammps header, cmds, topology params
# energy unit output is in eV 
# note - lammps output file has units of kcal/mol and gaussian output file has Hartrees

# file = 'structures/MFIt23.pdb'
sys_label = 'MFIt23_with_ethane_QM_part'; #class_label = 'ethane' 
atoms_file = f'structures/{sys_label}.pdb'

#sys_label = 'butane_gas_2'; class_label = 'butane' ; molec = 'C4H10'
#atoms_file = f'structures/{sys_label}.pdb'

# tester = load_system(file)
# get topology and other params for lammps calculation
atoms = read(atoms_file)
reset_positions(atoms, atoms.get_scaled_positions()[0]) # resetting positions of atoms to bring closer to the first atoms
atoms = init_atoms(atoms)
class_label = atoms.get_chemical_formula()
calc_params = set_calc_params(atoms, class_label, class_label)
calc_lmp = MM_LAMMPS(calc_params)
atoms.calc = calc_lmp
atoms.calc.label = atoms.get_chemical_formula()
print(atoms.get_potential_energy())
# print(atoms)
# qm 1butane = -4305.377384680251 eV
# mm 2butane = 17.5610127 eV
# mm 1butane = 0.39886227 eV
# total = -4288.214869611141 eV