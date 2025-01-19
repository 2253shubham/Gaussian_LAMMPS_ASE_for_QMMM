#!/usr/bin/env python

# command in ipython (depreciated)
#       %run optimize_dd_qm_gaus_mm_lammps.py -c 0 -m 1 -fl1 stdout1.log -ftr1 stdtraj1.traj -fl2 stdout2.log -ftr2 stdtraj2.traj -qi1 15 41
#       266 306 397 398 403 458 575 590 591 592 593 594 598 -qi2 625-646 -fz1 15 41 266 306 397 398 403 458 575 -fz2 576-646 -gl1 gaussian
#       -gl2 gaussian2 -s structure_num_1_of_nc18_TS_whole_MM_ps_dd_actual.pdb -tol1 0.08 -tol2 0.05 -bs1 gen1.basis -bs2 gen2.basis -metd
#       PBEPBE -xyz_fl traj1.xyz -rxn_type ps


# for nC18 pc in LTA at O1 site rs (old)  : 15,41,211,266,306,361,397,398,403,458,529,552,575,576-649 (648)  (atom index starts from zero and for range, the higher limit should be incremented with 1)
# for nC18 pc in LTA at O1 site ps : 15,41,211,266,306,361,397,398,403,458,529,552,575,590,591,592,593,594,598,625-650 (649)
# for nC18 pc in LTA at O1 site ps act : 15,41,266,306,397,398,403,458,575,590,591,592,593,594,598,625-646 (645)

'''Use ASE (git: fea26c) to perform geometry optimization using damped dynamics.'''
# Need two step optimization
# STEP 1 - pVDZ, fmax < 0.08
# STEP 2 - AUG-cc-pVDZ, fmax < 0.05

from __future__ import (print_function, unicode_literals)
from builtins import (str, range)
import argparse
#import mpi4py # including this gives mpi error
from ase.calculators.gaussian import Gaussian
from ase.io import read, Trajectory, write, trajectory
from ase.optimize.optimize import Optimizer
from ase.optimize import FIRE as DampedDynamics
from ase.optimize import BFGS as QuasiNewton
import ase.parallel as mpi
import numpy as np
#from asemc.util import sys
from ase import Atoms
from ase.calculators.qmmm import SimpleQMMM
from ase.calculators.gaussian import Gaussian
from ase.calculators.gromacs import Gromacs
from QM_MM_implemetation.bin_metd.modules.mm_lammps import (init_atoms, set_calc_params, reset_positions, get_charges, get_atom_types_from_pdb_file) 
from QM_MM_implemetation.bin_metd.modules.mm_lammps import MM_LAMMPS 
from QM_MM_implemetation.bin_metd.modules.simpleqmmm_mod import SmpQMMM
from ase import neighborlist
import itertools as it
import os
import numpy as np
import re


def get_atoms_info(atoms):
    atoms = init_atoms(atoms)
    get_charges(atoms)
    #atoms.calc = get_calculator(qm_idx, set_QM_calc(calc_label, method, basis_file), set_MM_calc1(atoms, qm_idx), set_MM_calc2(atoms))
    return atoms


def get_calculator(qm_idx, CALC_QM, CALC_MM1, CALC_MM2, **kwargs):
    return SmpQMMM(qm_idx, CALC_QM, CALC_MM1, CALC_MM2, **kwargs)


def get_gaussian_params(**kwargs):
    params = dict(
        chk='gaussian.chk',
        #xc='PBEPBE', 
        basis='Gen/Auto', 
        #basisfile='gen.basis',
        # SCRF='SMD,Solvent=Water',
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        Integral='Grid=UltraFine', SCF='Tight,IntRep,XQC', Symmetry="None",
        pop='Hirshfeld',
        EmpiricalDispersion="GD3BJ",
        **kwargs
    )

    #if 'nprocshared' not in params:
    #    params['nprocshared'] = sys.core_count

    #if 'mem' not in params and sys.mem_bytes is not None:
    #    params['mem'] = '{}GB'.format(int(sys.mem_bytes*0.6/(1024.**3)/(sys.core_count/params['nprocshared'])))

    print(params)
    return params


def get_link_atoms_nearest_neighbors(atoms):
    lnk_atms = []
    nn_lnk_atms = []
    nn_lnk_atms2 = []
    cutoff = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(
        cutoff, self_interaction=False, bothways=True)
    nl.update(atoms)
    for i in atoms.qm_idx:
        if atoms.arrays["atomtypes"][i]=="HL":
            lnk_atms.append(i)
            ind, offs = nl.get_neighbors(i)
            n1, n2 = nearest_neighbors(atoms, [i, *ind])
            nn_lnk_atms.append(n1)
            nn_lnk_atms2.append(n2)
    atoms.lnk_atms_idx = lnk_atms
    atoms.nn_lnk_atms_idx = nn_lnk_atms
    atoms.nn_lnk_atms_idx_2 = nn_lnk_atms2


def nearest_neighbors(atoms, array):
    ref = array[0]
    leng = []
    for j in array[1:]:
        leng.append(np.linalg.norm(atoms[ref].position-atoms[j].position))
    sort = np.argsort(leng)
    return array[1:][sort[0]], array[1:][sort[1]]


def parse_args():
    parser = argparse.ArgumentParser(description='Damped Dynamics Optimization')
    parser.add_argument('-c', '--charge', action='store', default=0, type=int)
    parser.add_argument('-m', '--multiplicity', action='store', default=1, type=int)
    parser.add_argument('-fl1', '--log_filename1', action='store', nargs='?', type=str, default="stdout") # log file for STEP1
    parser.add_argument('-ftr1', '--traj_filename1', action='store', nargs='?', type=str, default="traj") # traj file for STEP1
    parser.add_argument('-fl2', '--log_filename2',  action='store', nargs='?', type=str, default="") # log file for STEP2
    parser.add_argument('-ftr2', '--traj_filename2', action='store', nargs='?', type=str, default="") # traj file for STEP2
    parser.add_argument('-nproc','--proc-per-image', action='store', default=1, type=int)
    parser.add_argument('-mem', '--memory', action='store', default='1GB', type=str)
    parser.add_argument('-qi1', '--qm_idx_1', action='store', nargs='*', type=int, default=[]) # provide QM atoms
    parser.add_argument('-qi2', '--qm_idx_2', action='store', type=str, help="two numbers separated by a hyphen", default='0-0') # provide range of QM atoms if in sequence
    parser.add_argument('-nfz1', '--nfrz_idx_1', action='store', nargs='*', type=int, default=[]) # provide QM atoms
    parser.add_argument('-nfz2', '--nfrz_idx_2', action='store', type=str, help="two numbers separated by a hyphen", default='0-0') # provide range of QM atoms if in sequence
    parser.add_argument('-gl1', '--gaus_label1', action='store', nargs='?', type=str, default = "gaussian1") # STEP 1 gaussian output file name
    parser.add_argument('-gl2', '--gaus_label2', action='store', nargs='?', type=str, default="") # STEP 2 gaussian output file name
    parser.add_argument('-s', '--struct', action='store', required=True) # input structure to perform optimization on
    parser.add_argument('-tol1', '--tolerance1', action='store', default=0.08, type=float) # force convergence criteria tolerance for STEP 1
    parser.add_argument('-tol2', '--tolerance2', action='store', default=0.05, type=float) # force convergence criteria tolerance for STEP 2
    parser.add_argument('-bs1', '--basis_set1', action='store', nargs='?', type=str, default='gen1.basis') # basis set file for STEP 1
    parser.add_argument('-bs2', '--basis_set2', action='store', nargs='?', type=str, default="") # basis set file for STEP 2
    parser.add_argument('-qm_metd', '--qm_method', action='store', default='PBEPBE') # method to be used for optimization
    parser.add_argument('-xyz_fl1', '--xyz_file_name1', action='store', nargs='?', type=str, default='traj1.xyz') # file to store trajectory in xyz format
    parser.add_argument('-xyz_fl2', '--xyz_file_name2',  action='store', nargs='?', type=str, default="") # file to store trajectory in xyz format
    parser.add_argument('-inp_tr_fl1', '--inp_traj_filename1', action='store', nargs='?', type=str, default="") # input traj file from previous eqlb. run
    parser.add_argument('-inp_tr_fl2', '--inp_traj_filename2', action='store', nargs='?', type=str, default="") # input traj file from previous prd. run  
    #parser.add_argument('-rxn_type', '--reaction_type', action='store', default='rs') # rs for reactant state and ps for product state
    parser.add_argument('-r_st', '--run_style', action='store', type=int, default=1) # 0 for only MM calculations on the entire system, 1 for QM/MM calculation
    parser.add_argument('-coul_calc', '--coul_calc', action='store', type=int, default=0) # flag for coulomb calculations in lammps: 0 is No and 1 is Yes 
    args = parser.parse_args()
    return args


def run_optimizer(dd_optimizer, tol):
#    dd = DampedDynamics(dd_optimizer)
    dd_optimizer.run(fmax=tol, steps=1000)
    return dd_optimizer


def set_calculator(atoms, qm_idx, calc_label, method, basis_file, **kwargs):
    return get_calculator(qm_idx, set_QM_calc(calc_label, method, basis_file), set_MM_calc1(atoms, qm_idx), set_MM_calc2(atoms), **kwargs)


def set_MM_calc1(atoms, qm_idx): # MM calc of QM subsystem
    #qm_idx = args.qm_idx # qm subsystem atom indices
    atoms1 = atoms[qm_idx] # qm subsystem atoms object
    #atoms1 = init_atoms(atoms1)
    coeff_class_label1 = atoms.get_chemical_formula()
    class_label1 = atoms1.get_chemical_formula()
    calc_params1 = set_calc_params(atoms1, class_label1, coeff_class_label1) # parameters for MM calc1 (qm subsystem)
    calc_lmp1 = MM_LAMMPS(calc_params1)
    return calc_lmp1


def set_MM_calc2(atoms): # MM calc of the whole system
    coeff_class_label2 = atoms.get_chemical_formula()
    class_label2 = atoms.get_chemical_formula()
    calc_params2 = set_calc_params(atoms, class_label2, coeff_class_label2) # parameters for MM calc2 (whole system)
    calc_lmp2 = MM_LAMMPS(calc_params2)
    return calc_lmp2


def set_QM_calc(calc_label, method, basis_file): # QM calc of QM subsystem
    if any(i in args.memory for i in ("G", "g")): # if memory specified in GB in slurm
        memr = '{}MB'.format(int(float(*re.findall(r'[0-9]+', args.memory))*1024))
    elif any(i in args.memory for i in ("M", "m")):  # if memory specified in MB in slurm
        memr = '{}MB'.format(int(float(*re.findall(r'[0-9]+', args.memory))))
    else:
        memr = '10000MB'
    calc_params = get_gaussian_params(
    charge=args.charge,
    mult=args.multiplicity,
    basisfile=basis_file,
    xc=method,
    mem=memr,
    nprocshared=args.proc_per_image)
    calc_gau = Gaussian(label=calc_label, **calc_params)
    return calc_gau


def set_optimizer(atoms, filename, traj_file_name):
    dd_optimizer = DampedDynamics(atoms, restart=None, logfile=filename, trajectory=traj_file_name)
    return dd_optimizer


if __name__ == "__main__":
    args = parse_args()
    #if mpi.rank == 0:
    print(args)

    # Calculator parameters
    #calc_params = get_gaussian_params(
    #    charge=args.charge,
    #    multiplicity=args.multiplicity,
    #    nprocshared=args.proc_per_image,
    #)

    #calc_label1 = args.gaus_label1
 
    #if mpi.size > 1:
    #    calc_label = str(mpi.rank) + '-' + calc_label
    #if len(args.tmp) > 0:
    #    calc_label = args.tmp.rstrip('/') + '/' + calc_label

    if args.qm_idx_2:
        qm_idx1 = args.qm_idx_1
        before, after = args.qm_idx_2.split('-')
        qm_idx2 = np.arange(int(before), int(after)).tolist()
        qm_idx = qm_idx1 + qm_idx2
    else:
        qm_idx = args.qm_idx_1

    if args.nfrz_idx_2:
        nfrz_idx1 = args.nfrz_idx_1
        before, after = args.nfrz_idx_2.split('-')
        nfrz_idx2 = np.arange(int(before), int(after)).tolist()
        nfrz_idx = nfrz_idx1 + nfrz_idx2
    else:
        nfrz_idx = args.nfrz_idx_1

    # STEP 1
    atoms = read(args.struct)
    if ("atomtypes" not in atoms.arrays):
        atomtypes = get_atom_types_from_pdb_file(open(args.struct))
        atoms.set_array("atomtypes", np.array(atomtypes))
    #atoms.rt = args.reaction_type
    atoms = get_atoms_info(atoms) # Create the Atoms object
    atoms.qm_idx = qm_idx
    atoms.nfrz_idx = nfrz_idx
    get_link_atoms_nearest_neighbors(atoms)

    if args.run_style==0:
        atoms.qm_idx = []
        qm_idx = []

    qm_atoms = atoms[atoms.qm_idx]
    qm_atoms = get_atoms_info(qm_atoms)

    atoms.coul_calc = args.coul_calc
    qm_atoms.coul_calc = args.coul_calc

    if (args.traj_filename1 != None):

        if (args.inp_traj_filename1 != None):
            atoms_temp = trajectory.TrajectoryReader(args.inp_traj_filename1)
            atoms.set_positions(atoms_temp[-1].positions)
            qm_atoms.set_positions(atoms.positions[atoms.qm_idx])

        atoms.calc = set_calculator(atoms, atoms.qm_idx, args.gaus_label1, args.qm_method, args.basis_set1, nfrz_idx=atoms.nfrz_idx, lnk_atms_idx=atoms.lnk_atms_idx, nn_lnk_atms_idx=atoms.nn_lnk_atms_idx, qm_atoms_info = qm_atoms)
        dd_optimizer1 = set_optimizer(atoms, args.log_filename1, args.traj_filename1)
        dd_optimizer1 = run_optimizer(dd_optimizer1, tol = args.tolerance1)
   
        atoms.write("STEP1_opt.pdb")
    
        k = trajectory.TrajectoryReader(args.traj_filename1)
        write(args.xyz_file_name1, k)

    # STEP2

    if (args.traj_filename2 != None):

        if (args.inp_traj_filename2 != None):
            atoms_temp = trajectory.TrajectoryReader(args.inp_traj_filename2)
            atoms.set_positions(atoms_temp[-1].positions)
            qm_atoms.set_positions(atoms.positions[atoms.qm_idx])

        atoms.calc = set_calculator(atoms, atoms.qm_idx, args.gaus_label2, args.qm_method, args.basis_set2, nfrz_idx=atoms.nfrz_idx, lnk_atms_idx=atoms.lnk_atms_idx, nn_lnk_atms_idx=atoms.nn_lnk_atms_idx, qm_atoms_info = qm_atoms)
        dd_optimizer2 = set_optimizer(atoms, args.log_filename2, args.traj_filename2)
        dd_optimizer2 = run_optimizer(dd_optimizer2, tol = args.tolerance2)   

        atoms.write("STEP2_opt.pdb")
    
        k = trajectory.TrajectoryReader(args.traj_filename2)
        write(args.xyz_file_name2, k)