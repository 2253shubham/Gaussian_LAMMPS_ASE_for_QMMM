#!/usr/bin/env python

'''Use ASE (git: fea26c) to perform a transition-state search.'''


from __future__ import (print_function, unicode_literals)
from builtins import (str, range)
import argparse
#import mpi4py
from ase.calculators.gaussian import Gaussian
from ase.io import read, Trajectory, write, trajectory
from ase.neb import NEB
from ase.optimize import FIRE as DampedDynamics
from ase.optimize import BFGS as QuasiNewton
import ase.parallel as mpi
import numpy as np
from asemc.util import sys # needs asemc!!
from ase import Atoms
from simpleqmmm_mod import SmpQMMM
from ase.calculators.qmmm import SimpleQMMM
from ase.calculators.gaussian import Gaussian
from ase.calculators.gromacs import Gromacs
from Gaussian_LAMMPS_ASE_for_QMMM.modules.mm_lammps import (init_atoms, set_calc_params, reset_positions, get_atom_types_from_pdb_file) 
from Gaussian_LAMMPS_ASE_for_QMMM.modules.mm_lammps import MM_LAMMPS 
import itertools as it
import os


def get_atoms_from_file(qm_idx, frz_idx, filename, calc_label, id):
    atoms = read(filename)
    if ("atomtypes" not in atoms.arrays):
        atomtypes = get_atom_types_from_pdb_file(open(filename))
        atoms.set_array("atomtypes", np.array(atomtypes))
    #reset_positions(atoms, atoms.get_scaled_positions()[0]) # resetting positions of atoms to bring closer to the first atoms
    atoms = init_atoms(atoms)
    atoms.frz_idx = frz_idx
    atoms.calc = get_calculator(qm_idx, set_QM_calc(calc_label), set_MM_calc1(atoms, qm_idx, id), set_MM_calc2(atoms, id))
    #atoms.cell = 3e1*np.eye(3, dtype=np.float_)
    # atoms.pbc = (True, True, True)
    return atoms


def get_calculator(qm_idx, CALC_QM, CALC_MM1, CALC_MM2):
    return SmpQMMM(qm_idx, CALC_QM, CALC_MM1, CALC_MM2)


def get_gaussian_params(**kwargs):
    params = dict(
        chk='gaussian.chk',
        xc='PBEPBE', basis='Gen/Auto', basisfile='gen.basis',
        # SCRF='SMD,Solvent=Water',
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        Integral='Grid=UltraFine', SCF='Tight,IntRep,XQC', #Symmetry=None,
        **kwargs
    )

    if 'nprocshared' not in params:
        params['nprocshared'] = sys.core_count

    if 'mem' not in params and sys.mem_bytes is not None:
        params['mem'] = '{}GB'.format(int(sys.mem_bytes*0.6/(1024.**3)/(sys.core_count/params['nprocshared'])))

    return params


def neb_data(path, distance, energies, fname):
    '''Write out reaction coordinates and energies.'''

    data = []
    for i, apath in enumerate(path):
        data.append([np.sqrt(distance(apath, path[0])[0]),
                     energies[i]-energies[0]])
    np.savetxt(fname, np.asarray(data), delimiter=' ')


def parse_args():
    parser = argparse.ArgumentParser(description='Automatic transition state'
                                     ' search (ReaxFF/Gaussian).')
    parser.add_argument('-c', '--charge', action='store',
                        default=0, type=int)
    parser.add_argument('-m', '--multiplicity', action='store',
                        default=1, type=int)
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-ni', '--neb_l', action='store', default=-1, type=int)
    parser.add_argument('-ngc', '--neb-guide-count', action='store',
                        default=0, type=int)
    parser.add_argument('ngp', '--neb-guide-prefix', action='store',
                        default='neb')
    parser.add_argument('ngs', '--neb-guide-suffix', action='store',
                        default='.pdb')
    parser.add_argument('-ps', '--prod-struct', action='store')
    parser.add_argument('-ppi', '--proc-per-image', action='store',
                        default=1, type=int)
    #parser.add_argument('-mem', '--memory', action='store', default=1, type=int) 
    parser.add_argument('-qi1', '--qm_idx_1', action='store', nargs='*', type=int, default=[0,1,2,3,4,5,6,7])
    parser.add_argument('-qi2', '--qm_idx_2', action='store', type=str, help="two numbers separated by a hyphen")
    parser.add_argument('-fz1', '--frz_idx_1', action='store', required=True, nargs='*', type=int, default=[0,1,2,3,4,5,6,7]) # provide QM atoms
    parser.add_argument('-fz2', '--frz_idx_2', action='store', type=str, help="two numbers separated by a hyphen") # provide range of QM atoms if in sequence
    parser.add_argument('-rs', '--react_struct', action='store', required=True)
    parser.add_argument('-gl1', '--gaus_label1', action='store', type=str, default = "gaussian1") # STEP 1 gaussian output file name
    parser.add_argument('-gl2', '--gaus_label2', action='store', type=str, default = "gaussian2") # STEP 2 gaussian output file name
    parser.add_argument('-bs1', '--basis_set1', action='store', required=True, default='gen1.basis') # basis set file for STEP 1
    parser.add_argument('-bs2', '--basis_set2', action='store', required=True, default='gen2.basis') # basis set file for STEP 2
    parser.add_argument('-metd', '--method', action='store', required=True, default='PBEPBE') # method to be used for optimization
    parser.add_argument('--tmp', action='store', default='')
    args = parser.parse_args()
    return args


def prepare_neb(atoms1, atoms2, basis_set, method, nimages, qm_idx, frz_idx, calc_label,
                nguide=0, prefix='neb', suffix='.pdb'):
    images = [atoms1]
    guide_images = [0]

    for i in range(nimages):
        if i < nguide:
            guide_images.append(i+1)
            images.append(get_atoms_from_file(qm_idx, frz_idx, prefix+str(i+1)+suffix,
                                              calc_label.format(id=i+1), i+1))
        else:
            at = atoms1.copy()
            at.set_calculator(get_calculator(qm_idx, set_QM_calc(calc_label.format(id=i+1), method, basis_set), set_MM_calc1(at, qm_idx, i+1), set_MM_calc2(at, i+1)))
            images.append(at)

    images.append(atoms2)
    guide_images.append(-1)

    neb = NEB(images, k=5, climb=False, method='improvedtangent',
              remove_rotation_and_translation=False, parallel=True)
    neb.interpolate(guide_images=guide_images)

    if mpi.rank == 0:
        write_neb(images)

    return neb


def run_neb(neb, tol=0.08):
    tol_stage1 = tol

    dd = DampedDynamics(neb)
    qn = QuasiNewton(neb)

    for i, im in enumerate(neb.images[1:-1]):
        if i == mpi.rank:
            traj = Trajectory('neb{}.traj'.format(i+1), 'w', im, master=True)
            dd.attach(traj)
            qn.attach(traj)

    dd.run(fmax=tol_stage1, steps=1)

    if mpi.rank == 0:
        print('*Add "Guess=Read"')

    for im in neb.images[1:-1]:
        im.calc.set(Guess='Read')

    dd.run(fmax=tol_stage1, steps=500)

    # if mpi.rank == 0:
    #     print('== Switch to QuasiNewton (BFGS) ==')

    # qn.run(fmax=tol, steps=500)

    return neb


def set_MM_calc1(atoms, qm_idx, idx): # MM calc of QM subsystem
    #qm_idx = args.qm_idx # qm subsystem atom indices
    atoms1 = atoms[qm_idx] # qm subsystem atoms object
    atoms1 = init_atoms(atoms1)
    coeff_class_label1 = '{id}/atoms.get_chemical_formula()'
    class_label1 = atoms1.get_chemical_formula()
    calc_params1 = set_calc_params(atoms1, class_label1, coeff_class_label1.format(id=idx)) # parameters for MM calc1 (qm subsystem)
    calc_lmp1 = MM_LAMMPS(calc_params1)
    return calc_lmp1


def set_MM_calc2(atoms, idx): # MM calc of the whole system
    coeff_class_label2 = '{id}/atoms.get_chemical_formula()'
    class_label2 = atoms.get_chemical_formula()
    calc_params2 = set_calc_params(atoms, class_label2, coeff_class_label2.format(id=idx)) # parameters for MM calc2 (whole system)
    calc_lmp2 = MM_LAMMPS(calc_params2)#, 
    return calc_lmp2


def set_QM_calc(calc_label, method, basis_file): # QM calc of QM subsystem
    calc_params = get_gaussian_params(
    charge=args.charge,
    mult=args.multiplicity,
    basisfile=basis_file,
    xc=method,
    nprocshared=args.proc_per_image)
    calc_gau = Gaussian(label=calc_label, **calc_params)
    return calc_gau


def write_neb(images=None): # needs asemc!!
    '''Convert an ASE trajectory to an ARC movie file.'''

    if images is None:
        # images = read('neb.traj@-18:')
        images = read('neb.traj@:18')

    from asemc.io.helper import write_movie
    # write_movie('neb.arc', path, xyz.atomtypes,
    #             SimulationSystem.calc.get_boxvec(), shift=[0, 0, 0.5])
    write_movie('neb.res', [x.positions for x in images],
                [x.get_chemical_symbols() for x in images],
                [x.cell for x in images], filetype='res',
                shift=[0.5, 0.5, 0.5])


def write_movie_alt(filename, nimages, index):
    #open('latest_update.pdb', 'w').close()
    for i in range(nimages):
        tr = trajectory.TrajectoryReader("neb"+str(i+1)+".traj")
        write(filename, tr[index], append=True)


if __name__ == "__main__":
    args = parse_args()
    if mpi.rank == 0:
        print(args)

    # Calculator parameters
    #calc_params = get_gaussian_params(
    #    charge=args.charge,
    #    multiplicity=args.multiplicity,
    #    nprocshared=args.proc_per_image,
    #)

    calc_label = '{id}/gaussian'
    if mpi.size > 1:
        calc_label = str(mpi.rank) + '-' + calc_label
    if len(args.tmp) > 0:
        calc_label = args.tmp.rstrip('/') + '/' + calc_label

    if args.qm_idx_2:
        qm_idx1 = args.qm_idx_1
        before, after = args.qm_idx_2.split('-')
        qm_idx2 = np.arange(int(before), int(after)).tolist()
        qm_idx = qm_idx1 + qm_idx2
    else:
        qm_idx = args.qm_idx_1

    if args.frz_idx_2:
        frz_idx1 = args.frz_idx_1
        before, after = args.frz_idx_2.split('-')
        frz_idx2 = np.arange(int(before), int(after)).tolist()
        frz_idx = frz_idx1 + frz_idx2
    else:
        frz_idx = args.frz_idx_1

    # Create the Atoms object
    min1 = get_atoms_from_file(qm_idx, frz_idx, args.react_struct,
                               calc_label.format(id=0), 0)

    if args.prod_struct is not None:
        min2 = get_atoms_from_file(qm_idx, frz_idx, args.prod_struct,
                                   calc_label.format(id=args.neb_l+1), args.neb_l+1)

        if args.neb_l > 2:
            neb = prepare_neb(min1, min2, args.basis_set1,
                              args.method,
                              nimages=args.neb_l,
                              qm_idx=qm_idx,
                              frz_idx=frz_idx,
                              calc_label=calc_label,
                              nguide=args.neb_guide_count,
                              prefix=args.neb_guide_prefix,
                              suffix=args.neb_guide_suffix)
            if not args.dry_run:
                if mpi.size != args.neb_l:
                    raise RuntimeError('# of MPI ranks != # of NEB images')
                neb = run_neb(neb)

    k = trajectory.TrajectoryReader("neb"+str(3)+".traj")
    for i in range(len(k)):
        write_movie_alt("neb_step_"+str(i+1)+"_movie.pdb",args.neb_l,i)
