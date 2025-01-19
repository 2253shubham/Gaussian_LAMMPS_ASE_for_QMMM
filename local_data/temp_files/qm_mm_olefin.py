# sample submission line
# %run qm_mm_olefin.py -sl MFIt23_with_ethane -cl MFIt23_with_ethane -qi 22 23 53 76 54 88 62 36 35 37 77 52 20 21 86 61 89 90 91 92 93 94 95 96 -gl calc/gaussianQMMM_whole -mult 2

# QM/MM formula

 # E = E  (R  ) - E  (R  ) + E  (R   )
 #     QM  QM     MM  QM     MM  all
 #                mmcalc1    mmcalc2


from ase import Atoms
from ase.io import read, write
from ase.calculators.qmmm import SimpleQMMM
from ase.calculators.gaussian import Gaussian
from ase.calculators.gromacs import Gromacs
from mm_lammps import (init_atoms, set_calc_params, reset_positions) 
from mm_lammps import MM_LAMMPS 
import itertools as it
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='QM/MM calculation')
    parser.add_argument('-sl', '--sys_label', action='store',
                        default='MFIt23')
    #parser.add_argument('-cl', '--class_label', action='store',
    #                    default='MFIt23')
    parser.add_argument('-qi', '--qm_idx', action='store', nargs='*', type=int, default=[0,1,2,3,4,5,6,7])
    parser.add_argument('-gl', '--gaussian_label', action='store', default='calc/gaussianQMMM')
    parser.add_argument('-mult', '--multiplicity', action='store', default=1)
    args, unknown = parser.parse_known_args()
    return args


args = parse_args()
print(args)

#####################################################################################
# Set-up molecular system, read input files for lammps header, cmds, topology params
# energy unit output is in eV 
# note - lammps output file has units of kcal/mol and gaussian output file has Hartrees

# file = 'structures/MFIt23.pdb'
sys_label = args.sys_label #class_label = args.class_label ; #molec = 'C4H52O26Si23'
atoms_file = f'structures/{sys_label}.pdb'

#sys_label = 'MFIt23_with_ethane'; class_label = 'MFIt23' ; #molec = 'C4H52O26Si23'
#atoms_file = f'structures/{sys_label}.pdb'

# get topology and other params for lammps calculation
atoms = read(atoms_file)
reset_positions(atoms, atoms.get_scaled_positions()[0]) # resetting positions of atoms to bring closer to the first atoms
atoms = init_atoms(atoms)
class_label2 = atoms.get_chemical_formula()
calc_params2 = set_calc_params(atoms, class_label2, class_label2) # parameters for MM calc2 (whole system)
calc_lmp2 = MM_LAMMPS(calc_params2)#, coeff_cmds = coeff_comds, run_cmds = run_comds)


######################################################################################
# Make QM atoms selection with help of ase gui
#qm_idx = [20, 21, 22, 23, 32, 35, 36, 37, 52, 53, 54, 61, 62, 76, 77, 86, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100] # atom index starts from zero
# qm_idx = [0,1,2,3,4,5,6,7]
# qm_idx = [20, 21, 22, 23, 32, 35, 36, 37, 52, 53, 54, 61, 62, 76, 77, 86, 88, 89, 90, 91, 92] # atom index starts from zero
# qm_idx = [0,1,2,3,4,5,6,7,8,9,10,11,12,13] #,14,15,16,17,18,19,20,21,22,23,24,25,26,27]

qm_idx = args.qm_idx # qm subsystem atom indices
atoms1 = atoms[qm_idx] # qm subsystem atoms object
atoms1 = init_atoms(atoms1)
class_label1 = atoms1.get_chemical_formula()
calc_params1 = set_calc_params(atoms1, class_label1, class_label1) # parameters for MM calc1 (qm subsystem)
calc_lmp1 = MM_LAMMPS(calc_params1)

# Set up (separate) QM and MM calculators
# UPDATE SUB-TOPMON module add gaussian/g16.a03
#prog='/share/apps/gaussian/g16.a03/g16'
CALC_QM = Gaussian(label=args.gaussian_label, 
                     xc='wB97X', 
                     basis='6-31+G*', 
                     scf='maxcycle=200',
                     charge=0,
                     mult=int(args.multiplicity))
#                     command=f'{prog} > qm_out.log')

# to compute spin multiplicity - formula 2S+1 where S is the number of unpaired electrons. Each unpaired electron contributes 1/2 

CALC_MM1 = calc_lmp1
CALC_MM2 = calc_lmp2
#CALC_MM.set_own_params_runs('extra_mdrun_parameters', ' -nt 1 ') # for serial runs

# Set up hybrid calculator: simple subtractive
atoms.calc = SimpleQMMM(qm_idx, CALC_QM, CALC_MM1, CALC_MM2)
# caveats
# 1> pbc is set to false, cannot handle pbc
# 2> when making qm region and runnign qm_calc only the coordinates are passed and no info on the cell parameters are used
# 3> .com file between the ase gaussian run and ase gaussian run in qmmm is different. The latter has a force keyword appended to the header and no cell parameters are present.

print(atoms.get_potential_energy())



