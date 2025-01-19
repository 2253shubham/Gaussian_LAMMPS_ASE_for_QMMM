from __future__ import unicode_literals
from ase.data import atomic_numbers
import numpy as np
def get_params_element(elem):
    number = atomic_numbers[elem]
    return dict(molecule_atom_mapper=[0], molecule_atomic_numbers=[number], molecule_type_smarts='[#{}]'.format(number))
class ParameterBuilder(dict):
    molecule_fields = ['molecule_atom_mapper', 'molecule_atomic_numbers']
    def __init__(self, **kwargs):
        self.params = dict(
            checkpoint_frequency=60,
            hmc_adjust_frequency=100,
            hmc_steps=50,
            hmc_target_acceptance=0.95,
            hmc_timestep=10.0,
            md_ensemble='nvt',
            md_steps=100,
            md_timestep=0.5,
            molecule_atom_mapper=[],
            molecule_atomic_numbers=[],
            molecule_num_atoms=[],
            molecule_type_smarts=[],
            msmc_adjust_frequency=20,
            msmc_steps=10,
            msmc_target_acceptance=0.35,
            mspmc_pool_max_size=2000,
            volume_adjust_frequency=10,
            volume_target_acceptance=0.4)
        self.params.update(kwargs)
    def add_atoms(self, elem):
        try:
            for e in elem:
                p = get_params_element(e)
                self.add_molecule(p)
        except TypeError:
            p = get_params_element(elem)
            self.add_molecule(p)
        return self
    def add_molecule(self, params):
        for key in self.molecule_fields:
            val = params.pop(key)
            self.params[key].append(list(val))
            num_atoms = len(val)
        self.params['molecule_num_atoms'].append(num_atoms)
        self.params['molecule_type_smarts'].append(params.pop('molecule_type_smarts'))
        self.params.update(params)
        return self
    def build(self):
        p = self.params.copy()
        max_len = None
        for key in self.molecule_fields:
            if isinstance(p[key], np.ndarray):
                continue
            if max_len is None:
                max_len = max([len(x) for x in p[key]])
            p[key] = np.array([x+[-1]*(max_len-len(x)) for x in p[key]], dtype=np.int32)
        p['molecule_num_atoms'] = np.array(p['molecule_num_atoms'], dtype=np.int32)
        return p
    def update(self, d):
        self.params.update(d)
        return self