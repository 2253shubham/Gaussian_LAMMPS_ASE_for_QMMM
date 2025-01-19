import os
import glob
import numpy as np
from ase import Atoms
import lammps as lammps
from cProfile import run
from ase import neighborlist
from ase.io import read, write
from ase.geometry.analysis import Analysis
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.lammps import convert
from ase import units
from ase.geometry import wrap_positions
from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_masses
import warnings
import sys

######################################################################################


def read_atom_line_pdb(line_full):
    """
    Read atom line from pdb format
    HETATM    1  H14 ORTE    0       6.301   0.693   1.919  1.00  0.00        H
    """

    line = line_full.rstrip("\n")
    type_atm = line[0:6]
    if type_atm == "ATOM  " or type_atm == "HETATM":

        name = line[12:16].strip()

        altloc = line[16]
        resname = line[17:21]
        # chainid = line[21]        # Not used

        resseq = int(line[22:26].split()[0])  # sequence identifier
        # icode = line[26]          # insertion code, not used

        # atomic coordinates
        try:
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=np.float64,
            )
        except ValueError:
            raise ValueError("Invalid or missing coordinate(s)")

        # occupancy & B factor
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = None  # Rather than arbitrary zero or one

        if occupancy is not None and occupancy < 0:
            warnings.warn("Negative occupancy in one or more atoms")

        try:
            bfactor = float(line[60:66])
        except ValueError:
            # The PDB use a default of zero if the data is missing
            bfactor = 0.0

        # segid = line[72:76] # not used
        symbol = line[76:78].strip().upper()

    else:
        raise ValueError("Only ATOM and HETATM supported")

    return symbol, name, altloc, resname, coord, occupancy, bfactor, resseq


def get_atom_types_from_pdb_file(filename):
    atomtypes = []
    for line in filename.readlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            line_info = read_atom_line_pdb(line)
            atomtypes.append(line_info[1])
    return atomtypes


def init_atoms(atoms):
    # atoms = read(f'structures/{file}')
    # atoms = read(file)
    cutoff = neighborlist.natural_cutoffs(atoms)
    neigh_list = neighborlist.NeighborList(
        cutoff, self_interaction=False, bothways=False
    )
    neigh_list.update(atoms)
    # bonds = neigh_list.get_connectivity_matrix()
    topology_basic = Analysis(atoms)
    atoms.ase_bonds = topology_basic.unique_bonds
    atoms.ase_angles = topology_basic.unique_angles
    atoms.ase_dihedrals = topology_basic.unique_dihedrals
    # atoms.charges = np.zeros(len(atoms)) # initially all atoms are assigned zero charges
    # atoms.center(vacuum=4.0)
    get_atom_types(atoms)
    get_bond_types(atoms)
    get_angle_types(atoms)
    get_dihedral_types(atoms)
    return atoms


def create_atoms_in_lmp(atoms):
    # get positions of each atom to use lammps command "create atoms" WITH atom type accd to molecular info
    atom_cmds = []
    for i, pos in enumerate(atoms.positions):
        atom_cmds.append(
            "create_atoms {} single {} {} {}".format(
                atoms.atom_type_list[i], pos[0], pos[1], pos[2]
            )
        )
    for i in range(len(atoms)):
        atom_cmds.append(
            "set atom {} charge {}".format(i + 1, atoms.arrays["charges"][i])
        )
    return atom_cmds


def get_atom_types(atoms):
    # loop over atoms to get types, attach to the atoms object so we can use to create bonds
    atom_types_ase = {}
    atom_types_pdb = atoms.arrays["atomtypes"]
    uniq_atom_type = np.unique(atom_types_pdb)  # unique atom label from PDB file
    for i in range(
        len(uniq_atom_type)
    ):  # build dictionary to attach to ASE-atoms object
        atom_types_ase[uniq_atom_type[i]] = i + 1
    atoms.atom_types = atom_types_ase
    atoms.atom_type_list = [
        atoms.atom_types[
            atom_types_pdb[i]
        ]  # assign atom types to each atom in ASE-atoms object
        for i in range(len(atom_types_pdb))
    ]


def get_atomtype_masses(atoms):
    molec = atoms.get_chemical_formula()
    mass_cmds = []
    type_masses = {}
    # mass_lbl = 'all_masses.{}'.format(molec)
    mass_dict = getattr(all_masses, molec)()
    for i in atoms.atom_types:
        type_masses[i] = mass_dict[i]
    # for i, sym in enumerate(atoms.calc.parameters.atom_types):
    for i, sym in enumerate(atoms.atom_types):
        mass = type_masses[sym]
        mass_cmds.append("mass {} {}".format(i + 1, mass))
    # atoms.calc.parameters.atom_type_masses = type_masses
    return type_masses, mass_cmds


def get_bond_types(atoms):
    # loop over bond matrix to identify and map bond types
    # want to reverse index the atom TYPES in each connected bond, return bond types based on the atom TYPES bonded together (i.e., ATOM 1 is carbon, but atom TYPE 1 is ch2 because it's a carbon on ch2)
    # atom TYPES come from atoms.atom_type_list
    unq_bonded_atom_types = []
    bond_list_with_type = []
    for i in range(len(atoms.ase_bonds[0])):
        for j in range(len(atoms.ase_bonds[0][i])):
            bonded_atoms = [i, atoms.ase_bonds[0][i][j]]
            bonded_atom_types = [atoms.atom_type_list[k] for k in bonded_atoms]
            if (
                bonded_atom_types in unq_bonded_atom_types
                or bonded_atom_types[::-1] in unq_bonded_atom_types
            ):  # look for existing bond type
                bond_type = (
                    unq_bonded_atom_types.index(bonded_atom_types) + 1
                    if bonded_atom_types in unq_bonded_atom_types
                    else unq_bonded_atom_types.index(bonded_atom_types[::-1]) + 1
                )
            else:
                bond_type = (
                    len(unq_bonded_atom_types) + 1
                )  # new bond type is always +1 existing bond types
                unq_bonded_atom_types.append(bonded_atom_types)
            # print(i, bond_type, bonded_atoms, bonded_atom_types, unq_bonded_atom_types)
            bond_list_with_type.append([bond_type, *[k + 1 for k in bonded_atoms]])
    atoms.bonds = np.array(bond_list_with_type)  # desired format for LAMMPS
    atoms.unq_bonds = np.array(
        unq_bonded_atom_types
    )  # list of atom TYPES in unique bonds; len(unq_bonds) = # of types of bonds


def get_angle_types(atoms):
    # loop over angle matrix to identify and map angle types
    unq_angular_atoms = []
    angle_list_with_type = []
    for i in range(len(atoms.ase_angles[0])):
        for j in range(len(atoms.ase_angles[0][i])):
            angular_atoms = [i, *atoms.ase_angles[0][i][j]]
            angular_atom_types = [atoms.atom_type_list[k] for k in angular_atoms]
            if (
                angular_atom_types in unq_angular_atoms
                or angular_atom_types[::-1] in unq_angular_atoms
            ):
                ang_type = (
                    unq_angular_atoms.index(angular_atom_types) + 1
                    if angular_atom_types in unq_angular_atoms
                    else unq_angular_atoms.index(angular_atom_types[::-1]) + 1
                )
            else:
                ang_type = len(unq_angular_atoms) + 1
                unq_angular_atoms.append(angular_atom_types)
            angle_list_with_type.append([ang_type, *[k + 1 for k in angular_atoms]])
    atoms.angles = np.array(angle_list_with_type)
    atoms.unq_angles = np.array(unq_angular_atoms)


def get_dihedral_types(atoms):
    # loop over dihedral matrix to identify and map dihedral types
    unq_dihedral_atoms = []
    dihedral_list_with_type = []
    for i in range(len(atoms.ase_dihedrals[0])):
        for j in range(len(atoms.ase_dihedrals[0][i])):
            dihedral_atoms = [i, *atoms.ase_dihedrals[0][i][j]]
            dihedral_atom_type = [atoms.atom_type_list[k] for k in dihedral_atoms]
            if (
                dihedral_atom_type in unq_dihedral_atoms
                or dihedral_atom_type[::-1] in unq_dihedral_atoms
            ):
                dih_type = (
                    unq_dihedral_atoms.index(dihedral_atom_type) + 1
                    if dihedral_atom_type in unq_dihedral_atoms
                    else unq_dihedral_atoms.index(dihedral_atom_type[::-1]) + 1
                )
            else:
                dih_type = len(unq_dihedral_atoms) + 1
                unq_dihedral_atoms.append(dihedral_atom_type)
            dihedral_list_with_type.append([dih_type, *[k + 1 for k in dihedral_atoms]])
    atoms.dihedrals = np.array(dihedral_list_with_type)
    atoms.unq_dihedrals = np.array(unq_dihedral_atoms)


def get_impropers(atoms):
    # get impropers
    return impropers


def get_improper_types(impropers):
    # loop over improper matrix to identify and map improper types
    return impropers_with_type


def create_coeff_file(atoms, coeff_class_label):
    filename1 = "coeff/" + coeff_class_label + ".coeff"
    filename2 = "coeff/" + coeff_class_label + ".del"
    out = open(filename1, "w")
    out.truncate(0)
    out2 = open(filename2, "w")
    out2.truncate(0)
    infs_bonds = []
    from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_atoms_LJ
    from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_bonds
    from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_angles
    from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_dihedrals

    for i, sym in enumerate(atoms.atom_types):
        LJ_dict = getattr(all_atoms_LJ, sym)()
        eps = LJ_dict["epsln"]
        sig = LJ_dict["sigma"]
        print("pair_coeff {} {} {} {}".format(i + 1, i + 1, eps, sig), file=out)
    for i in range(len(atoms.unq_bonds)):
        a1 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_bonds[i][0])
        ]
        a2 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_bonds[i][1])
        ]
        bond = a1 + "_" + a2
        bond_dict = getattr(all_bonds, bond)()
        k = bond_dict["k_value"]
        eq_b = bond_dict["eqlb_dist"]
        if eq_b < 0:
            infs_bonds.append(bond)
            print("bond_coeff {} {} {}".format(i + 1, k, eq_b), file=out)
            print("delete_bonds all bond {} remove special".format(i + 1), file=out2)
        else:
            print("bond_coeff {} {} {}".format(i + 1, k, eq_b), file=out)
    for i in range(len(atoms.unq_angles)):
        a1 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_angles[i][0])
        ]
        a2 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_angles[i][1])
        ]
        a3 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_angles[i][2])
        ]
        angle = a1 + "_" + a2 + "_" + a3
        if any(pattern in angle for pattern in infs_bonds):
            print("angle_coeff {} {} {}".format(i + 1, 0, 0), file=out)
            print("delete_bonds all angle {} remove special".format(i + 1), file=out2)
        else:
            angle_dict = getattr(all_angles, angle)()
            k = angle_dict["k_value"]
            eq_a = angle_dict["eqlb_ang"]
            print("angle_coeff {} {} {}".format(i + 1, k, eq_a), file=out)
    for i in range(len(atoms.unq_dihedrals)):
        a1 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_dihedrals[i][0])
        ]
        a2 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_dihedrals[i][1])
        ]
        a3 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_dihedrals[i][2])
        ]
        a4 = list(atoms.atom_types.keys())[
            list(atoms.atom_types.values()).index(atoms.unq_dihedrals[i][3])
        ]
        dihedral = a1 + "_" + a2 + "_" + a3 + "_" + a4
        if any(pattern in dihedral for pattern in infs_bonds):
            print("dihedral_coeff {} {} {} {} {}".format(i + 1, 0, 0, 0, 0), file=out)
            print(
                "delete_bonds all dihedral {} remove special".format(i + 1), file=out2
            )
        else:
            dihedral_dict = getattr(all_dihedrals, dihedral)()
            p1 = dihedral_dict["param_1"]
            p2 = dihedral_dict["param_2"]
            p3 = dihedral_dict["param_3"]
            p4 = dihedral_dict["param_4"]
            print(
                "dihedral_coeff {} {} {} {} {}".format(i + 1, p1, p2, p3, p4), file=out
            )
    out.close()
    out2.close()
    return filename1, filename2


def get_charges(atoms):
    atoms.arrays["charges"] = np.zeros(
        len(atoms)
    )  # initially all atoms are assigned zero charges
    from Gaussian_LAMMPS_ASE_for_QMMM.local_data import all_atoms_LJ

    for i in range(len(atoms.arrays["atomtypes"])):
        LJ_dict = getattr(all_atoms_LJ, atoms.arrays["atomtypes"][i])()
        atoms.arrays["charges"][i] = LJ_dict["q"]


def reset_positions(atoms_obj, ref):
    for j in range(3):
        atoms_obj.positions[:, j][
            atoms_obj.get_scaled_positions()[:, j] - ref[j] > 0.5
        ] = (
            atoms_obj.get_scaled_positions()[:, j][
                atoms_obj.get_scaled_positions()[:, j] - ref[j] > 0.5
            ]
            - 1
        ) * atoms_obj.cell.array.diagonal()[
            j
        ]
        atoms_obj.positions[:, j][
            atoms_obj.get_scaled_positions()[:, j] - ref[j] < -0.5
        ] = (
            atoms_obj.get_scaled_positions()[:, j][
                atoms_obj.get_scaled_positions()[:, j] - ref[j] < -0.5
            ]
            + 1
        ) * atoms_obj.cell.array.diagonal()[
            j
        ]


def reset_charges(atoms):
    reset_charges_cmds = []
    for i in atoms.qm_idx:
        reset_charges_cmds.append(
            "set atom {} charge {}".format(i + 1, atoms.arrays["charges"][i])
        )
    return reset_charges_cmds


def set_calc_params(atoms, class_label, coeff_class_label):
    header = []
    lcmds = []
    coeff_comds = []
    run_comds = []
    lmp_files = glob.glob(f"./*/{class_label}.*")
    for lf in lmp_files:
        ext = os.path.splitext(lf)[1]
        lines = open(lf)
        print(sys.version.split()[0])
        #    if (sys.version.split()[0]!="3.9.16"):
        #        match ext: # supported in python 3.10 >
        #            case ".lhdr":
        #                header = set_header(lines)
        #            case ".lcmd":
        #                lcmds = set_cmds(lines)
        #            #case ".coeff":
        #            #    coeff_comds = set_coeff(lines)
        #            case ".lrun":
        #                run_comds = set_run_cmds(lines)
        #    else: # probably using lower python version (need to correct this) # currently set to "3.9.16"
        if ext == ".lhdr":
            header = set_header(lines)
        elif ext == ".lcmd":
            lcmds = set_cmds(lines)
        # elif ext == ".coeff":
        #    coeff_comds = set_coeff(lines)
        elif ext == ".lrun":
            run_comds = set_run_cmds(lines)
    # coeff_comds = set_coeff(open(create_coeff_file(atoms, coeff_class_label)))
    return dict(
        lammps_header=header,
        lmpcmds=lcmds,
        create_box=False,
        create_atoms=False,
        boundary=False,
        run_cmds=run_comds,
    )  # , atom_types=atoms.atom_types, run_cmds=run_comds)#, coeff_cmds=coeff_comds, coeff_cmds, run_cmds # log_file=f'{sys_label}_aselammps.log'


def set_cmds(commands):
    c = []
    for i in commands:
        c.append(i)
    return c


def set_coeff(coeff):
    coeff_cmds = []
    for i in coeff:
        coeff_cmds.append(i)
    return coeff_cmds
    # add error flag: if number of bonds/angles/dihedrals doesn't match len(unq_bonds/angles/dihedrals), throw an error


def set_header(hdr):
    header_cmds = []
    for i in hdr:
        header_cmds.append(i)
    return header_cmds
    # add error flag: if header doesn't match expected format, throw error


def set_run_cmds(rcmds):
    run_cmds = []
    for i in rcmds:
        run_cmds.append(i)
    return run_cmds


def set_topology(atoms):
    bond_cmds = []
    angle_cmds = []
    dihed_cmds = []
    for t, i1, i2 in atoms.bonds:
        bond_cmds.append("create_bonds single/bond {} {} {} ".format(t, i1, i2))
    for t, i1, i2, i3 in atoms.angles:
        angle_cmds.append("create_bonds single/angle {} {} {} {}".format(t, i1, i2, i3))
    for t, i1, i2, i3, i4 in atoms.dihedrals:
        dihed_cmds.append(
            "create_bonds single/dihedral {} {} {} {} {}".format(t, i1, i2, i3, i4)
        )
    return bond_cmds, angle_cmds, dihed_cmds


def unit_cell_coords(atoms_obj):
    for j in range(3):
        atoms_obj.positions[:, j][
            atoms_obj.positions[:, j] > atoms_obj.cell.array.diagonal()[j]
        ] = (
            atoms_obj.positions[:, j][
                atoms_obj.positions[:, j] > atoms_obj.cell.array.diagonal()[j]
            ]
            - atoms_obj.cell.array.diagonal()[j]
        )
        atoms_obj.positions[:, j][atoms_obj.positions[:, j] < 0] = (
            atoms_obj.positions[:, j][atoms_obj.positions[:, j] < 0]
            + atoms_obj.cell.array.diagonal()[j]
        )


class MM_LAMMPS(LAMMPSlib):
    #     # WANT TO MAKE MM_LAMMPS A SUBCLASS OF THE PARENT CLASS LAMMPSLIB CALCULATOR TO ACCESS ALL THESE FUNCTIONS
    #     # object-based programming: this means we can use our lammps object BUILT ON TOP OF the lammpslib calculator
    #     # python tutorial for object based programming
    # def __init__(self, *args, **kwargs):
    #     LAMMPSlib.__init__(self, *args, **kwargs)

    _energy_unit = units.kcal / units.mol
    _time_unit = units.fs
    _length_unit = units.Angstrom
    _force_unit = _energy_unit / _length_unit
    _pressure_unit = 1.01325e5 * units.Pascal
    _velocity_unit = _length_unit / _time_unit
    _mass_unit = 1e-3 * units.kg / units.mol

    default_parameters = dict(
        atom_types=None,
        atom_type_masses=None,
        log_file=None,
        lammps_name="",
        keep_alive=False,
        new_system=True,  # new lammps_system to which the same calculator is attached
        lammps_header=[
            "units metal",
            "atom_style atomic",
            "atom_modify map array sort 0 0",
        ],
        amendments=None,
        post_changebox_cmds=None,
        boundary=True,
        create_box=False,
        create_atoms=False,
        read_molecular_info=False,
        coeff_cmds=None,
        run_cmds=None,
        comm=None,
    )

    def __init__(self, calc_params):
        LAMMPSlib.__init__(self, **calc_params)

    #    startup_options  = ['-echo', 'screen', '-log', f'{sys_label}_mmlammps.log']
    #    self.lmp = lammps.lammps('', startup_options)
    #    self.start_lammps()

    # def init_calc(self, sys_label):
    #    calc = LAMMPSlib(**(calc_params))
    #    startup_options  = ['-echo', 'screen', '-log', f'{sys_label}_mmlammps.log']
    #    self.lmp = lammps.lammps('', startup_options)
    # self.initialized = False
    #    self.start_lammps()
    #    atoms.calc = calc
    #    return atoms.calc

    def attach_all_cmds(self, cmds):
        for i in cmds:
            self.lmp.command(i)

    def hackaround_change_box(self, atoms):
        # do all the commands that WOULD happen if create box was TRUE because we pass create_box in as FALSE so the system can be initialized correctly --> 2 functions: box tilt large and "initialize box" code section
        box_cmds = []
        box_dim = np.empty(3)
        n_types = len(atoms.atom_types)
        nbond_types = len(atoms.unq_bonds)
        nangle_types = len(atoms.unq_angles)
        ndihed_types = len(atoms.unq_dihedrals)
        for i in range(len(atoms.cell)):
            box_dim[i] = np.linalg.norm(atoms.cell[i])
        xh, yh, zh = box_dim
        xl, yl, zl = np.zeros(3)
        region = "region mybox prism {} {} {} {} {} {} {} {} {}".format(
            xl, xh, yl, yh, zl, zh, 0, 0, 0
        )
        create_box = "create_box {} mybox bond/types {} angle/types {} dihedral/types {} extra/bond/per/atom {} extra/angle/per/atom {} extra/dihedral/per/atom {}".format(
            n_types, nbond_types, nangle_types, ndihed_types, 50, 50, 50
        )  # change extra/bond/angle/dihedral as per case to case basis
        box_cmds.append(region)
        box_cmds.append(create_box)
        if hasattr(atoms, "coul_calc"):
            if atoms.coul_calc == 1:
                pairstyle_cmd = "pair_style lj/cut/coul/long 12.0"  # change cutoff accordingly to system size
                kspace_cmd = "kspace_style pppm 1.0e-4"
                box_cmds.append(pairstyle_cmd)
                box_cmds.append(kspace_cmd)
        return box_cmds
        # atoms.calc.lmp.command('box tilt large')

    def initialise_lammps(self, atoms):
        # Initialising commands
        if not hasattr(atoms, "bonds"):
            atoms = init_atoms(atoms)
            # get_charges(atoms)
        self.parameters.atom_types = atoms.atom_types

        if self.parameters.boundary:
            # if the boundary command is in the supplied commands use that
            # otherwise use atoms pbc
            for cmd in self.parameters.lmpcmds:
                if "boundary" in cmd:
                    break
            else:
                self.lmp.command("boundary " + self.lammpsbc(atoms))

        # Initialize cell
        self.set_cell(
            atoms, change=self.parameters.create_box
        )  # different from original lammpslib function

        if self.parameters.atom_types is None:
            # if None is given, create from atoms object in order of appearance
            s = atoms.get_chemical_symbols()
            _, idx = np.unique(s, return_index=True)
            s_red = np.array(s)[np.sort(idx)].tolist()
            self.parameters.atom_types = {j: i + 1 for i, j in enumerate(s_red)}

        # Initialize box
        if self.parameters.create_box:
            # count number of known types
            n_types = len(self.parameters.atom_types)
            create_box_command = "create_box {} cell".format(n_types)
            self.lmp.command(create_box_command)

        # Initialize the atoms with their types
        # positions do not matter here
        if self.parameters.create_atoms:
            self.lmp.command("echo none")  # don't echo the atom positions
            self.rebuild(atoms)
            self.lmp.command("echo log")  # turn back on
        else:
            self.previous_atoms_numbers = atoms.numbers.copy()

        # execute the user commands
        box_cmds = self.hackaround_change_box(atoms)
        atom_cmds = create_atoms_in_lmp(atoms)
        type_masses, mass_cmds = get_atomtype_masses(atoms)
        reset_charges_cmds = []
        #    if hasattr(atoms,"qm_idx"):
        #        reset_charges_cmds = reset_charges(atoms)
        bond_cmds, angle_cmds, dihed_cmds = set_topology(atoms)
        f1, f2 = create_coeff_file(atoms, self.label)
        self.parameters.coeff_cmds = set_coeff(open(f1))  # manual edit made in SMPQMMM
        self.parameters.del_cmds = set_coeff(open(f2))  # manual edit made in SMPQMMM

        all_cmds = (
            box_cmds
            + atom_cmds
            + mass_cmds
            + self.parameters.coeff_cmds
            + bond_cmds
            + angle_cmds
            + dihed_cmds
            + self.parameters.del_cmds
            + self.parameters.run_cmds
        )  # +reset_charges_cmds+self.parameters.run_cmds
        # print(all_cmds)
        for i in all_cmds:
            # print("[FROM ASE] : "+ i)
            self.lmp.command(i)
        # for cmd in self.parameters.lmpcmds:
        # print("[FROM ASE] : "+ cmd)
        # self.lmp.command(cmd)
        self.parameters.atom_type_masses = type_masses

        # Set masses after user commands, e.g. to override
        # EAM-provided masses
        # for sym in self.parameters.atom_types:
        #    mass = self.parameters.atom_type_masses[sym]
        #   self.lmp.command('mass %d %.30f' % (
        #        self.parameters.atom_types[sym],
        #        convert(mass, "mass", "ASE", self.units)))

        # Define force & energy variables for extraction
        self.lmp.command("variable pxx equal pxx")
        self.lmp.command("variable pyy equal pyy")
        self.lmp.command("variable pzz equal pzz")
        self.lmp.command("variable pxy equal pxy")
        self.lmp.command("variable pxz equal pxz")
        self.lmp.command("variable pyz equal pyz")

        # I am not sure why we need this next line but LAMMPS will
        # raise an error if it is not there. Perhaps it is needed to
        # ensure the cell stresses are calculated
        self.lmp.command(
            "thermo_style custom step etotal ke temp pe epair emol ebond eangle edihed eimp evdwl ecoul elong etail press vol pxx pyy pzz pxy pxz pyz"
        )
        self.lmp.command("thermo_modify norm no flush no line multi")
        self.lmp.command("variable fx atom fx")
        self.lmp.command("variable fy atom fy")
        self.lmp.command("variable fz atom fz")

        # do we need this if we extract from a global ?
        self.lmp.command("variable pe equal pe")

        self.lmp.command("neigh_modify delay 0 every 1 check yes")

        self.initialized = True

    def propagate(
        self,
        atoms,
        properties,
        system_changes,
        n_steps,
        dt=None,
        dt_not_real_time=False,
        velocity_field=None,
    ):
        """ "atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these five: 'positions', 'numbers', 'cell',
            'pbc', 'charges' and 'magmoms'.
        """
        if len(system_changes) == 0:
            print("enters here")
            return

        self.coord_transform = None

        if not self.started:
            # self.parameters.log_file = f'{atoms.sys_label}_aselammps.log'
            if hasattr(
                atoms, "rt"
            ):  # check if the atoms object has some information on the optimization calculation type
                self.parameters.log_file = f"{atoms.get_chemical_formula()}_aselammps_{atoms.rt}.log"  # hash out this line if want to print MM results on the screen, manual edit
            else:
                self.parameters.log_file = f"{atoms.get_chemical_formula()}_aselammps.log"  # hash out this line if want to print MM results on the screen, manual edit
            self.start_lammps()
        if not self.initialized:
            self.initialise_lammps(atoms)
        else:  # still need to reset cell
            # NOTE: The whole point of ``post_changebox_cmds`` is that they're
            # executed after any call to LAMMPS' change_box command.  Here, we
            # rely on the fact that self.set_cell(), where we have currently
            # placed the execution of ``post_changebox_cmds``, gets called
            # after this initial change_box call.

            # Apply only requested boundary condition changes.  Note this needs
            # to happen before the call to set_cell since 'change_box' will
            # apply any shrink-wrapping *after* it's updated the cell
            # dimensions
            if "pbc" in system_changes:
                change_box_str = "change_box all boundary {}"
                change_box_cmd = change_box_str.format(self.lammpsbc(atoms))
                self.lmp.command(change_box_cmd)

            # Reset positions so that if they are crazy from last
            # propagation, change_box (in set_cell()) won't hang.
            # Could do this only after testing for crazy positions?
            # Could also use scatter_atoms() to set values (requires
            # MPI comm), or extra_atoms() to get pointers to local
            # data structures to zero, but then we would have to be
            # careful with parallelism.
            self.lmp.command("set atom * x 0.0 y 0.0 z 0.0")
            self.set_cell(atoms, change=True)
        """    

        if self.parameters.atom_types is None:
            raise NameError("atom_types are mandatory.")

        do_rebuild = (not np.array_equal(atoms.numbers,
                                         self.previous_atoms_numbers)
                      or ("numbers" in system_changes))
        if not do_rebuild:
            do_redo_atom_types = not np.array_equal(
                atoms.numbers, self.previous_atoms_numbers)
        else:
            do_redo_atom_types = False

        self.lmp.command('echo none')  # don't echo the atom positions
        if do_rebuild:
            self.rebuild(atoms)
        elif do_redo_atom_types:
            self.redo_atom_types(atoms)
        self.lmp.command('echo log')  # switch back log

        self.set_lammps_pos(atoms)
      
        if self.parameters.amendments is not None:
            for cmd in self.parameters.amendments:
                self.lmp.command(cmd)

        if n_steps > 0:
            if velocity_field is None:
                vel = convert(
                    atoms.get_velocities(),
                    "velocity",
                    "ASE",
                    self.units)
            else:
                # FIXME: Do we need to worry about converting to lammps units
                # here?
                vel = atoms.arrays[velocity_field]

            # If necessary, transform the velocities to new coordinate system
            if self.coord_transform is not None:
                vel = np.dot(self.coord_transform, vel.T).T

            # Convert ase velocities matrix to lammps-style velocities array
            lmp_velocities = list(vel.ravel())

            # Convert that lammps-style array into a C object
            c_double_array = (ctypes.c_double * len(lmp_velocities))
            lmp_c_velocities = c_double_array(*lmp_velocities)
            self.lmp.scatter_atoms('v', 1, 3, lmp_c_velocities)

        # Run for 0 time to calculate
        if dt is not None:
            if dt_not_real_time:
                self.lmp.command('timestep %.30f' % dt)
            else:
                self.lmp.command('timestep %.30f' %
                                 convert(dt, "time", "ASE", self.units))
        """
        self.lmp.command("run %d" % n_steps)

        if n_steps > 0:
            # TODO this must be slower than native copy, but why is it broken?
            pos = np.array([x for x in self.lmp.gather_atoms("x", 1, 3)]).reshape(-1, 3)
            if self.coord_transform is not None:
                pos = np.dot(pos, self.coord_transform)

            # Convert from LAMMPS units to ASE units
            pos = convert(pos, "distance", self.units, "ASE")

            atoms.set_positions(pos)

            vel = np.array([v for v in self.lmp.gather_atoms("v", 1, 3)]).reshape(-1, 3)
            if self.coord_transform is not None:
                vel = np.dot(vel, self.coord_transform)
            if velocity_field is None:
                atoms.set_velocities(convert(vel, "velocity", self.units, "ASE"))

        # Extract the forces and energy
        self.results["energy"] = convert(
            self.lmp.extract_variable("pe", None, 0), "energy", self.units, "ASE"
        )
        # self.results['energy'] = self.lmp.extract_variable('pe', None, 0)
        self.results["free_energy"] = self.results["energy"]

        stress = np.empty(6)
        stress_vars = ["pxx", "pyy", "pzz", "pyz", "pxz", "pxy"]

        for i, var in enumerate(stress_vars):
            stress[i] = self.lmp.extract_variable(var, None, 0)

        stress_mat = np.zeros((3, 3))
        stress_mat[0, 0] = stress[0]
        stress_mat[1, 1] = stress[1]
        stress_mat[2, 2] = stress[2]
        stress_mat[1, 2] = stress[3]
        stress_mat[2, 1] = stress[3]
        stress_mat[0, 2] = stress[4]
        stress_mat[2, 0] = stress[4]
        stress_mat[0, 1] = stress[5]
        stress_mat[1, 0] = stress[5]
        if self.coord_transform is not None:
            stress_mat = np.dot(
                self.coord_transform.T, np.dot(stress_mat, self.coord_transform)
            )
        stress[0] = stress_mat[0, 0]
        stress[1] = stress_mat[1, 1]
        stress[2] = stress_mat[2, 2]
        stress[3] = stress_mat[1, 2]
        stress[4] = stress_mat[0, 2]
        stress[5] = stress_mat[0, 1]

        # self.results['stress'] = -stress
        self.results["stress"] = convert(-stress, "pressure", self.units, "ASE")

        # definitely yields atom-id ordered force array
        # f = np.array(self.lmp.gather_atoms("f", 1, 3)).reshape(-1, 3)
        f = convert(
            np.array(self.lmp.gather_atoms("f", 1, 3)).reshape(-1, 3),
            "force",
            self.units,
            "ASE",
        )

        if self.coord_transform is not None:
            self.results["forces"] = np.dot(f, self.coord_transform)
        else:
            self.results["forces"] = f.copy()

        # otherwise check_state will always trigger a new calculation
        self.atoms = atoms.copy()

        if not self.parameters.keep_alive:
            self.lmp.close()
            if self.parameters.new_system:
                self.started = False
                self.initialized = False
                self.lmp = None
                # print("############################################################################################################################################")
                # print("############################################################################################################################################")

    def rebuild(self, atoms):
        try:
            n_diff = len(atoms.numbers) - len(self.previous_atoms_numbers)
        except Exception:  # XXX Which kind of exception?
            n_diff = len(atoms.numbers)

        if n_diff > 0:
            if any([("reax/c" in cmd) for cmd in self.parameters.lmpcmds]):
                self.lmp.command("pair_style lj/cut 2.5")
                self.lmp.command("pair_coeff * * 1 1")

                for cmd in self.parameters.lmpcmds:
                    if (
                        ("pair_style" in cmd)
                        or ("pair_coeff" in cmd)
                        or ("qeq/reax" in cmd)
                    ):
                        self.lmp.command(cmd)

            cmd = "create_atoms 1 random {} 1 NULL".format(n_diff)
            self.lmp.command(cmd)
        elif n_diff < 0:
            cmd = "group delatoms id {}:{}".format(
                len(atoms.numbers) + 1, len(self.previous_atoms_numbers)
            )
            self.lmp.command(cmd)
            cmd = "delete_atoms group delatoms"
            self.lmp.command(cmd)

        self.redo_atom_types(atoms)

    def redo_atom_types(self, atoms):
        current_types = set(
            (i + 1, self.parameters.atom_types[sym])
            for i, sym in enumerate(atoms.arrays["atomtypes"].tolist())
        )

        try:
            previous_types = set(
                (i + 1, self.parameters.atom_types[sym])
                for i, sym in enumerate(self.previous_atoms_types)
            )
        except Exception:  # XXX which kind of exception?
            previous_types = set()

        for i, i_type in current_types - previous_types:
            cmd = "set atom {} type {}".format(i, i_type)
            self.lmp.command(cmd)

        self.previous_atoms_numbers = atoms.numbers.copy()
        self.previous_atoms_types = atoms.arrays["atomtypes"].tolist()

    #    def start_lammps(self):
    #       LAMMPSlib.start_lammps(self)

    def start_lammps(self):
        # Only import lammps when running a calculation
        # so it is not required to use other parts of the
        # module

        if self.parameters.log_file is None:
            cmd_args = ["-echo", "screen", "-log", "none"]
        else:
            cmd_args = [
                "-echo",
                "screen",
                "-log",
                "none",
                "-screen",
                self.parameters.log_file,
            ]

        self.cmd_args = cmd_args

        if self.lmp is None:
            self.lmp = lammps.lammps(
                self.parameters.lammps_name, self.cmd_args, comm=self.parameters.comm
            )

        # Run header commands to set up lammps (units, etc.)
        for cmd in self.parameters.lammps_header:
            # print("[FROM ASE] : "+ cmd)
            self.lmp.command(cmd)

        for cmd in self.parameters.lammps_header:
            if "units" in cmd:
                self.units = cmd.split()[1]

        if "lammps_header_extra" in self.parameters:
            if self.parameters.lammps_header_extra is not None:
                for cmd in self.parameters.lammps_header_extra:
                    self.lmp.command(cmd)

        self.started = True
