from __future__ import unicode_literals


def get_params_CH4():
    return dict(
        molecule_atom_mapper=[0,  1,  2,  3,  4],
        molecule_atomic_numbers=[6,  1,  1,  1,  1],
        molecule_type_smarts='[#6](~[#1])(~[#1])(~[#1])~[#1]',
    )


def get_params_C2H6():
    return dict(
        molecule_atom_mapper=[0,  1,  2,  3,  4,  5,  6,  7],
        molecule_atomic_numbers=[6,  1,  1,  1,  6,  1,  1,  1],
        molecule_type_smarts='[#6](~[#1])(~[#1])(~[#1])~[#6](~[#1])(~[#1])~[#1]',
    )


def get_params_H2O():
    return dict(
        molecule_atom_mapper=[0, 1, 2],
        molecule_atomic_numbers=[8, 1, 1],
        molecule_type_smarts='[#8](~[#1])~[#1]',
    )


def get_params_H3O():
    return dict(
        molecule_atom_mapper=[0, 1, 2, 3],
        molecule_atomic_numbers=[8, 1, 1, 1],
        molecule_type_smarts='[#8](~[#1])(~[#1])~[#1]',
    )


def get_params_H5O2():
    return dict(
        molecule_atom_mapper=[0,  1,  2,  3,  4,  5,  6],
        molecule_atomic_numbers=[8,  1,  1,  1,  8,  1,  1],
        molecule_type_smarts='[#8](~[#1])(~[#1])~[#1]~[#8](~[#1])~[#1]',
    )
