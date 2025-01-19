from __future__ import (print_function, unicode_literals)
from builtins import (str, range)
import argparse
import copy
#import mpi4py # including this gives mpi error
from ase.calculators.gaussian import Gaussian
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io import read, Trajectory, write, trajectory
from nvtberendsen_mod import NVTBerendsen_mod 
import ase.parallel as mpi
import numpy as np
from plumed import Plumed
from ase.calculators.gaussian import Gaussian
from ase.calculators.plumed import Plumed
from mm_lammps import (init_atoms, set_calc_params, reset_positions, get_charges, get_atom_types_from_pdb_file) 
from mm_lammps import MM_LAMMPS 
from simpleqmmm_mod import SmpQMMM
from atoms_mod import Atoms_mod
from ase import neighborlist
from ase import units
import itertools as it
import os
from ase import units
import re


# for all calculations - energy in eV, time in ase units and length in angs.

def get_atoms_info(atoms):
    atoms = init_atoms(atoms)
    get_charges(atoms)
    #atoms.calc = get_calculator(qm_idx, set_QM_calc(calc_label, method, basis_file), set_MM_calc1(atoms, qm_idx), set_MM_calc2(atoms), qm_idx=qm_idx)
    return atoms


def get_calculator(qm_idx, CALC_QM, CALC_MM1, CALC_MM2, **kwargs):
    return SmpQMMM(qm_idx, CALC_QM, CALC_MM1, CALC_MM2, **kwargs)


def get_gaussian_params(**kwargs):
    params = dict(
        chk='gaussian.chk',
        basis='Gen/Auto', 
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
            nn_lnk_atms.append(ind[0])
            nn_lnk_atms2.append(ind[1])
    atoms.lnk_atms_idx = lnk_atms
    atoms.nn_lnk_atms_idx = nn_lnk_atms
    atoms.nn_lnk_atms_idx_2 = nn_lnk_atms2


def parse_args():
    parser = argparse.ArgumentParser(description='QM/MM with Metadynamics, with Berendsen')
    parser.add_argument('-c', '--charge', action='store', default=0, type=int)
    parser.add_argument('-m', '--multiplicity', action='store', default=1, type=int)
    #parser.add_argument('-fl', '--log_filename', action='store', default="stdout") # log file 
    parser.add_argument('-e_tr_fl', '--eqlb_traj_filename', action='store', default="eqlb.traj") # traj file for eqlb. run
    parser.add_argument('-p_tr_fl', '--prd_traj_filename', action='store', default="prd.traj") # traj file for prd. run
    parser.add_argument('-inp_e_tr_fl', '--inp_eqlb_traj_filename', action='store', nargs='?', default="") # input traj file from previous eqlb. run
    parser.add_argument('-inp_p_tr_fl', '--inp_prd_traj_filename', action='store', nargs='?', default="") # input traj file from previous prd. run
    #parser.add_argument('--proc-per-image', action='store',default=1, type=int)
    parser.add_argument('-mem', '--memory', action='store', default='1GB', type=str)   
    parser.add_argument('-nproc', '--proc-per-image', action='store', default=1, type=int)  
    parser.add_argument('-qi1', '--qm_idx_1', action='store', nargs='*', type=int, default=[]) # provide QM atoms
    parser.add_argument('-qi2', '--qm_idx_2', action='store', type=str, help="two numbers separated by a hyphen", default='0-0') # provide range of QM atoms if in sequence
    parser.add_argument('-nfz1', '--nfrz_idx_1', action='store', nargs='*', type=int, default=[]) # provide QM atoms
    parser.add_argument('-nfz2', '--nfrz_idx_2', action='store', type=str, help="two numbers separated by a hyphen", default='0-0') # provide range of QM atoms if in sequence
    parser.add_argument('-gl', '--gaus_label', action='store', type=str, default = "gaussian") # gaussian output file name
    parser.add_argument('-s1', '--struct1', action='store', required=True) # input structure to perform optimization on
    parser.add_argument('-bs', '--basis_set', action='store', default='gen1.basis') # basis set file
    parser.add_argument('-qm_metd', '--qm_method', action='store', default='PBEPBE') # method to be used for QM optimization
    parser.add_argument('-metdy_setup_1', '--metdy_setup_file_1', action='store', nargs='?', default='plumed_1.dat') # method to be used for optimization
    parser.add_argument('-metdy_setup_2', '--metdy_setup_file_2', action='store', nargs='?', default='plumed_2.dat') # method to be used for optimization
    parser.add_argument('-xyz_fl', '--xyz_file_name', action='store', default='traj.xyz') # file to store production trajectory in xyz format
    parser.add_argument('-esteps', '--eqlb_md_steps', action='store', type=int, default=1000) # total time steps to equilibriate the system
    parser.add_argument('-psteps', '--prd_md_steps', action='store', type=int, default=50000) # total time steps to run production runs
    parser.add_argument('-temp', '--temperature', action='store', type=float, default=600) # temperature in K
    parser.add_argument('-ts', '--time_step', action='store', type=float, default=5) # time step in femtoseconds
    parser.add_argument('-r_st', '--run_style', action='store', type=int, default=1) # 0 for only MM calculations on the entire system, 1 for QM/MM calculation
    parser.add_argument('-coul_calc', '--coul_calc', action='store', type=int, default=0) # flag for coulomb calculations in lammps: 0 is No and 1 is Yes 
    args = parser.parse_args()
    return args


def run_md(md, steps):
    md.run(steps)
    return md


def set_calculator(atoms, qm_idx, calc_label, method, basis_file, **kwargs):
    return get_calculator(qm_idx, set_QM_calc(calc_label, method, basis_file), set_MM_calc1(atoms, qm_idx), set_MM_calc2(atoms), **kwargs)
 

def set_MM_calc1(atoms, qm_idx): # MM calc of QM subsystem
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


def set_md(atoms, timestep, temp, traj_file_name): # run eqlb runs
    md = NVTBerendsen_mod(atoms, timestep, temp, taut=10, trajectory=traj_file_name, fixcm=False)  
    return md


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

    atoms = read(args.struct1)
    atoms = Atoms_mod(atoms)
    MaxwellBoltzmannDistribution(atoms, temperature_K=args.temperature)
    atoms.arrays["momenta"][[i for i in range(len(atoms.arrays["momenta"])) if i not in nfrz_idx]] = 0
    set_temp = atoms.get_temperature()
    if ("atomtypes" not in atoms.arrays):
        atomtypes = get_atom_types_from_pdb_file(open(args.struct1))
        atoms.set_array("atomtypes", np.array(atomtypes))
    atoms = get_atoms_info(atoms) # Create the Atoms object
    #ini_temp = 1000
    #MaxwellBoltzmannDistribution(atoms, temperature_K=ini_temp)
    #scale_temp = args.temperature/atoms.get_temperature()
    #while not (scale_temp < 1.1 and scale_temp > 0.9):
    #    ini_temp *= scale_temp
    #    MaxwellBoltzmannDistribution(atoms, temperature_K=ini_temp)
    #    atoms.arrays["momenta"][[i for i in range(len(atoms.arrays["momenta"])) if i not in nfrz_idx]] = 0
    #    scale_temp = args.temperature/atoms.get_temperature()
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

    if (args.inp_eqlb_traj_filename != None):
        atoms_temp = trajectory.TrajectoryReader(args.inp_eqlb_traj_filename)
        atoms.set_positions(atoms_temp[-1].positions)
        atoms.set_momenta(atoms_temp[-1].get_momenta())
        qm_atoms.set_positions(atoms.positions[atoms.qm_idx])
        qm_atoms.set_momenta(atoms.get_momenta()[atoms.qm_idx])

    if args.eqlb_md_steps > 0:
        atoms.calc = set_calculator(atoms, atoms.qm_idx, args.gaus_label, args.qm_method, args.basis_set, nfrz_idx=atoms.nfrz_idx, lnk_atms_idx=atoms.lnk_atms_idx, nn_lnk_atms_idx=atoms.nn_lnk_atms_idx, qm_atoms_info = qm_atoms) # Create the Atoms object)
        if (args.metdy_setup_file_1 != None):
            setup = open(args.metdy_setup_file_1, "r").read().splitlines() # time in ase units, length in angs., energy in eV
            energy_calc = atoms.calc
            if (args.inp_eqlb_traj_filename == None): 
                atoms.calc = Plumed(calc=energy_calc, input=setup, timestep=args.time_step, atoms=atoms, kT=args.temperature*units.kB) # updating calculator for metadynamics (plumed has ps time units)
            else:
                atoms.calc = Plumed(calc=energy_calc, input=setup, timestep=args.time_step, atoms=atoms, kT=args.temperature*units.kB, restart=True) # updating calculator for metadynamics (plumed has ps time units) 
                atoms.calc.istep = len(trajectory.TrajectoryReader(args.inp_eqlb_traj_filename)) - 1
       
        # equilibriation runs
        eqlb_md_run = set_md(atoms, args.time_step, set_temp, args.eqlb_traj_filename)
        eqlb_md_run = run_md(eqlb_md_run, args.eqlb_md_steps) # behaving wierdly with "wrap" in nvtberendsen_mod probably with plumed calc
        ########
        k = trajectory.TrajectoryReader(args.eqlb_traj_filename)
        write("eqlb_traj.xyz", k)

    #p = read(args.eqlb_traj_filename) # last eqlb. configuration
    #p.write(args.struct2)
    atom_types_ret = atoms.arrays["atomtypes"]

    #atoms2 = copy.deepcopy(atoms) # gives issues with atoms object has plumed calc
    #qm_atoms2 = copy.deepcopy(qm_atoms)
    if "k" in locals():
        #atoms2 = read(args.struct1)
        #if ("atomtypes" not in atoms2.arrays):
        #    atomtypes = get_atom_types_from_pdb_file(open(args.struct))
        #    atoms2.set_array("atomtypes", np.array(atomtypes))
        #atoms2.arrays["atomtypes"] = atom_types_ret
        #atoms2 = get_atoms_info(atoms2)# Create the Atoms object
        atoms.set_positions(k[-1].positions)
        atoms.set_momenta(k[-1].get_momenta())
        atoms.write("eqlb_config.pdb")
        qm_atoms.set_positions(atoms.positions[atoms.qm_idx])
        qm_atoms.set_momenta(atoms.get_momenta()[atoms.qm_idx])
        
    #MaxwellBoltzmannDistribution(atoms2, temperature_K=args.ini_temp)
    #atoms2.arrays["momenta"][[i for i in range(len(atoms2.arrays["momenta"])) if i not in nfrz_idx]] = 0
    #atoms2.qm_idx = qm_idx
    #atoms2.nfrz_idx = nfrz_idx
    #get_link_atoms_nearest_neighbors(atoms2)

    #qm_atoms2 = atoms2[atoms2.qm_idx]
    #qm_atoms2 = get_atoms_info(qm_atoms2)

    if (args.inp_prd_traj_filename != None):
        atoms_temp = trajectory.TrajectoryReader(args.inp_prd_traj_filename)
        atoms.set_positions(atoms_temp[-1].positions)
        atoms.set_momenta(atoms_temp[-1].get_momenta())
        qm_atoms.set_positions(atoms.positions[atoms.qm_idx])
        qm_atoms.set_momenta(atoms.get_momenta()[atoms.qm_idx])

    if args.prd_md_steps > 0 :
        atoms.calc = set_calculator(atoms, atoms.qm_idx, args.gaus_label, args.qm_method, args.basis_set, nfrz_idx=atoms.nfrz_idx, lnk_atms_idx=atoms.lnk_atms_idx, nn_lnk_atms_idx=atoms.nn_lnk_atms_idx, qm_atoms_info = qm_atoms) # Create the Atoms object)
        if (args.metdy_setup_file_2 != None):
            setup = open(args.metdy_setup_file_2, "r").read().splitlines() # time in ase units, length in angs., energy in eV
            energy_calc = atoms.calc
            if (args.inp_prd_traj_filename == None):
                atoms.calc = Plumed(calc=energy_calc, input=setup, timestep=args.time_step, atoms=atoms, kT=args.temperature*units.kB) # updating calculator for metadynamics (plumed has ps time units)
            else:
                atoms.calc = Plumed(calc=energy_calc, input=setup, timestep=args.time_step, atoms=atoms, kT=args.temperature*units.kB, restart=True) # updating calculator for metadynamics (plumed has ps time units) 
                atoms.calc.istep = len(trajectory.TrajectoryReader(args.inp_prd_traj_filename)) - 1
        
        # production runs    
        prd_md_run = set_md(atoms, args.time_step, set_temp, args.prd_traj_filename)
        prd_md_run = run_md(prd_md_run, args.prd_md_steps)
        # to convert .traj file to .xyz file
        k = trajectory.TrajectoryReader(args.prd_traj_filename)
        write(args.xyz_file_name, k)

        atoms.set_positions(k[-1].positions)
        atoms.set_momenta(k[-1].get_momenta())
        atoms.write("prod_config.pdb")

