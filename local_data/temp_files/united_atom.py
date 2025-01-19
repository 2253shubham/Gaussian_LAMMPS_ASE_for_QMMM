from __future__ import unicode_literals


def get_params_CH4():
    return dict(
        molecule_atom_mapper=[0],
        molecule_atom_valency=[4],
        molecule_atomic_numbers=[6],
        molecule_type_smarts='[#6](~[#1])(~[#1])(~[#1])~[#1]',
    )


def get_params_C4H10():
    return dict(
        molecule_atom_mapper=[0, 4, 7, 10],
        molecule_atom_valency=[4, 4, 4, 4],
        molecule_atomic_numbers=[6, 6, 6, 6],
        molecule_type_smarts='[#6](~[#1])(~[#1])(~[#1])~[#6](~[#1])(~[#1])~[#6](~[#1])(~[#1])~[#6](~[#1])(~[#1])(~[#1])',
    )
