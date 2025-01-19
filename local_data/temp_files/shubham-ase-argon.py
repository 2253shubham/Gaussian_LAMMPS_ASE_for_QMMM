from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, Trajectory
from ase import Atoms
import numpy as np
from asemc.calculators.lammps_calculator import LammpsCalculator

def get_lammps_lj():
    return dict(
        lammps_header='''
        units real
        atom_style atomic
        atom_modify map array
        ''',
        atom_types={
            18: 1,
        },
        command='''
        pair_style lj/cut 14.0
        pair_modify mix arithmetic
        pair_coeff * * 0.3844 3.4
        special_bonds coul 0.0 0.0 0.5
        ''',
    )

#check lammpsdata_read_lammps_data
# can specify bond_coeff, pair_coeff and angle_coeff in main.lmp
#atoms = read("benzalgs.pdb")
atoms = read("argon.pdb")
#atoms = read("Ar_gas.pdb")
#atoms.cell = 3e1*np.eye(3, dtype=np.float_)
calc = LammpsCalculator(label='label_str', **(get_lammps_lj()))
atoms.calc = calc
print(atoms.get_potential_energy())

#from asemc.calculators.lammps_calculator import LammpsCalculator
#atoms.calc = LammpsCalculator

# problems with loading topmon
# energy units
# atom_type pair_style error