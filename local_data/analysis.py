from ase import neighborlist
from ase.geometry.analysis import Analysis
import numpy as np
from ase.io import read, Trajectory, write, trajectory
import QM_MM_implemetation.bin_metd.modules.mm_lammps as mm_lammps

def init_atoms(atoms):
    cutoff = neighborlist.natural_cutoffs(atoms)
    neigh_list = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=False)
    neigh_list.update(atoms)
    topology_basic = Analysis(atoms)
    atoms.ase_bonds = topology_basic.unique_bonds
    atoms.ase_angles = topology_basic.unique_angles
    atoms.ase_dihedrals = topology_basic.unique_dihedrals
    return atoms, topology_basic


def get_mags(ana, atoms):
    bonds = []
    angles = []
    dihedrals = []
    for i in range(len(atoms.bonds)):
        bonds.append(ana.get_bond_value(0,[atoms.bonds[i][1]-1, atoms.bonds[i][2]-1]))
    for i in range(len(atoms.angles)):
        angles.append(ana.get_angle_value(0,[atoms.angles[i][1]-1, atoms.angles[i][2]-1, atoms.angles[i][3]-1]))
    for i in range(len(atoms.dihedrals)):
        dihedrals.append(ana.get_dihedral_value(0,[atoms.dihedrals[i][1]-1, atoms.dihedrals[i][2]-1, atoms.dihedrals[i][3]-1, atoms.dihedrals[i][4]-1]))
    return bonds, angles, dihedrals


traj_n = trajectory.TrajectoryReader("eqlb.traj")

atom0 = read("../../md_rs.pdb")
atom0 = mm_lammps.init_atoms(atom0)
atom0, ana0 = init_atoms(atom0)
bonds0, angles0, dihedrals0 = get_mags(ana0,atom0)

coor3 = traj_n[3].positions
atom3 = read("../../md_rs.pdb")
atom3.set_positions(coor3)
atom3 = mm_lammps.init_atoms(atom3)
atom3, ana3 = init_atoms(atom3)
bonds3, angles3, dihedrals3 = get_mags(ana3,atom3)

coor8 = traj_n[8].positions
atom8 = read("../../md_rs.pdb")
atom8.set_positions(coor8)
atom8 = mm_lammps.init_atoms(atom8)
atom8, ana8 = init_atoms(atom8)
bonds8, angles8, dihedrals8 = get_mags(ana8,atom8)

coor10 = traj_n[10].positions
atom10 = read("../../md_rs.pdb")
atom10.set_positions(coor10)
atom10 = mm_lammps.init_atoms(atom10)
atom10, ana10 = init_atoms(atom10)
bonds10, angles10, dihedrals10 = get_mags(ana10,atom10)


bd3C1H = [i for i in range(len(atom3.bonds)) if atom3.bonds[i][0]==9]

bd8C1H = [i for i in range(len(atom8.bonds)) if atom8.bonds[i][0]==8]

bd0C1H = [i for i in range(len(atom0.bonds)) if atom0.bonds[i][0]==8]

bd3C2H = [i for i in range(len(atom3.bonds)) if atom3.bonds[i][0]==11]

bd0C2H = [i for i in range(len(atom0.bonds)) if atom0.bonds[i][0]==10]

bd8C2H = [i for i in range(len(atom8.bonds)) if atom8.bonds[i][0]==10]

bd3C1C2 = [i for i in range(len(atom3.bonds)) if atom3.bonds[i][0]==10]

bd8C1C2 = [i for i in range(len(atom8.bonds)) if atom8.bonds[i][0]==9]

bd0C1C2 = [i for i in range(len(atom0.bonds)) if atom0.bonds[i][0]==9]

bdl0C1H = np.array(bonds0)[bd0C1H]

bdl3C1H = np.array(bonds3)[bd3C1H]

bdl8C1H = np.array(bonds8)[bd8C1H]

bdl0C2H = np.array(bonds0)[bd0C2H]

bdl3C2H = np.array(bonds3)[bd3C2H]

bdl8C2H = np.array(bonds8)[bd8C2H]

bdl0C1C2 = np.array(bonds0)[bd0C1C2]

bdl3C1C2 = np.array(bonds3)[bd3C1C2]

bdl8C1C2 = np.array(bonds8)[bd8C1C2]

db0_8C1H = abs(bdl0C1H-bdl8C1H)

db0_3C1H = abs(bdl0C1H-bdl3C1H)

db3_8C1H = abs(bdl3C1H-bdl8C1H)

db0_8C2H = abs(bdl0C2H-bdl8C2H)

db0_3C2H = abs(bdl0C2H-bdl3C2H)

db3_8C2H = abs(bdl3C2H-bdl8C2H)

e1 = np.loadtxt(fname = "Tot_energy.txt")
e2 = np.loadtxt(fname = "MM_entire_energy.txt")
e3 = np.loadtxt(fname = "MM_QM_sub_energy.txt")
e4 = np.loadtxt(fname = "QM_QM_sub_energy.txt")

er1 = np.loadtxt(fname = "../02_19_2024_qmmm_run_expanse/Tot_energy.txt")
er2 = np.loadtxt(fname = "../02_19_2024_qmmm_run_expanse/MM_entire_energy.txt")
er3 = np.loadtxt(fname = "../02_19_2024_qmmm_run_expanse/MM_QM_sub_energy.txt")
er4 = np.loadtxt(fname = "../02_19_2024_qmmm_run_expanse/QM_QM_sub_energy.txt")

e1 = e1 - er1[0]
e2 = e2 - er2[0]
e3 = e3 - er3[0]
e4 = e4 - er4[0]

plt.plot(e1, "blue")
plt.plot(e2, "green")
plt.plot(e3, "red")
plt.plot(e4, "orange")