from __future__ import unicode_literals
from builtins import zip
from future.utils import viewitems
import logging
import os
import numpy as np
from ase import units
from ase.calculators.calculator import Calculator
from ase.constraints import FixAtoms
from ase.data import atomic_masses
import lammps


_logger = logging.getLogger('mc.calculator')


class LammpsCalculator(Calculator):
    '''Potentials provided by LAMMPS.

    This calculator does not convert properties to ASE
    standard units (Ang for length, amu for mass, and eV for
    energy), as different potentials may have different
    conventions (e.g., ReaxFF potential files use "real"
    units).

    '''

    implemented_properties = ['energy', 'forces', 'stress',
                              'pressure', 'temperature']

    default_parameters = {
        'ftol': 0.003,
        'lammps_header': '''
        units real
        atom_style charge
        atom_modify map array
        ''',
        'pbc': (True, True, True),
    }

    # these conversion factors correspond to the 'real'
    # units in LAMMPS
    _energy_unit = units.kcal / units.mol
    _time_unit = units.fs
    _length_unit = units.Angstrom
    _force_unit = _energy_unit / _length_unit
    _pressure_unit = 1.01325e5 * units.Pascal
    _velocity_unit = _length_unit / _time_unit
    _mass_unit = 1e-3 * units.kg / units.mol

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

        self.cell = np.array([30.0, 30.0, 30.0], dtype=np.float_)
        self.setup_lammps()

    def __del__(self):
        if self.lmp is not None:
            _logger.debug('Deleting LAMMPS pointer.')
            self.lmp.close()
            self.lmp = None

    def __getstate__(self):
        state = self.__dict__.copy()

        for key in ['lmp', '_positions_ctypes',
                    '_atomtypes_ctypes', '_velocities_ctypes']:
            state[key] = None

        return state

    def __setstate__(self, state):
        self.__dict__ = state
        natoms = self.natoms
        positions = self.positions
        atomtypes = self.atomtypes
        velocities = self.velocities
        self.setup_lammps()
        if natoms > 0:
            self.set_number_atoms(natoms)
            self.set_lammps_arrays(positions=positions, atomtypes=atomtypes,
                                   velocities=velocities)

    def calculate(self, atoms, properties, system_changes):
        '''Calculate properties using LAMMPS.'''

        _logger.debug('In calculate, system_changes = %s', system_changes)

        if len(system_changes) != 0:
            Calculator.calculate(self, atoms, properties, system_changes)
            self.set_configuration(positions=atoms.positions,
                                   numbers=atoms.numbers,
                                   cell=atoms.cell.diagonal(),
                                   velocities=atoms.get_velocities())
            if self.natoms == 0:
                self.zero_properties()
                return
            self.lmp.command('run 0 post no')
            self.get_lammps_properties()

    def clean(self):
        pass

    def get_boxvec(self):
        '''Get cell dimensions from LAMMPS.'''

        boxxlo = self.lmp.extract_global('boxxlo', 1)
        boxxhi = self.lmp.extract_global('boxxhi', 1)
        boxylo = self.lmp.extract_global('boxylo', 1)
        boxyhi = self.lmp.extract_global('boxyhi', 1)
        boxzlo = self.lmp.extract_global('boxzlo', 1)
        boxzhi = self.lmp.extract_global('boxzhi', 1)
        return np.array([boxxhi-boxxlo, boxyhi-boxylo, boxzhi-boxzlo],
                        dtype=np.float_)

    def get_lammps_arrays(self, data=['positions', 'types', 'velocities']):
        '''Get numpy arrays from LAMMPS.

        gather_atoms() returns ctypes arrays, which are cast
        to numpy arrays that share the same memory.

        '''

        # _logger.debug('Getting LAMMPS arrays.')

        if 'positions' in data:
            self._positions_ctypes = self.lmp.gather_atoms('x', 1, 3)
            self.natoms = len(self._positions_ctypes) // 3
            self.positions = np.ctypeslib.as_array(
                self._positions_ctypes).reshape((self.natoms, 3))

        if 'types' in data:
            self._atomtypes_ctypes = self.lmp.gather_atoms('type', 0, 1)
            self.natoms = len(self._atomtypes_ctypes)
            self.atomtypes = np.ctypeslib.as_array(
                self._atomtypes_ctypes, shape=self.natoms)

        if 'velocities' in data:
            self._velocities_ctypes = self.lmp.gather_atoms('v', 1, 3)
            self.natoms = len(self._velocities_ctypes) // 3
            self.velocities = np.ctypeslib.as_array(
                self._velocities_ctypes).reshape((self.natoms, 3))

    def get_lammps_energies(self):
        '''Get the energies of the system using variables.'''

        self.resultsv['energy'] = (self.lmp.extract_variable('pec', 'all', 0) *
                                   self._energy_unit)
        self.resultsv['ke'] = (self.lmp.extract_variable('kec', 'all', 0) *
                               self._energy_unit)
        self.resultsv['etotal'] = (self.lmp.extract_variable(
            'etotalc', 'all', 0) * self._energy_unit)

    def get_lammps_properties(self, properties=implemented_properties):
        '''Get properties from LAMMPS.'''

        # _logger.debug('Getting LAMMPS properties.')

        if 'energy' in properties:
            self.results['energy'] = (self.lmp.extract_compute(
                'thermo_pe', 0, 0) * self._energy_unit)

        if 'temperature' in properties:
            self.results['temperature'] = self.lmp.extract_compute(
                'thermo_temp', 0, 0)

        if 'pressure' in properties:
            self.results['pressure'] = (self.lmp.extract_compute(
                'thermo_press', 0, 0) * self._pressure_unit)

        if 'stress' in properties:
            s = self.lmp.extract_compute('thermo_press', 0, 1)
            stress = np.empty(6, dtype=np.float_)
            stress[0] = s[0]  # xx
            stress[1] = s[1]  # yy
            stress[2] = s[2]  # zz
            stress[3] = s[5]  # yz
            stress[4] = s[4]  # xz
            stress[5] = s[3]  # xy
            self.results['stress'] = (stress * self._pressure_unit)

        if 'forces' in properties:
            self.results['forces'] = (np.ctypeslib.as_array(
                self.lmp.gather_atoms('f', 1, 3)).reshape((self.natoms, 3)) *
                                      self._force_unit)

    def get_lammps_stresses(self):
        '''Get the stresses of the system using variables.'''

        stress = np.empty(6, dtype=np.float_)
        for i, var in enumerate(['pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy']):
            stress[i] = self.lmp.extract_variable(var, 'all', 0)
        self.resultsv['stress'] = (stress * self._pressure_unit)
        self.resultsv['pressure'] = (self.lmp.extract_variable(
            'press', 'all', 0) * self._pressure_unit)

    def get_pressure(self, atoms=None):
        '''Get the isotropic pressure of the system.'''

        return (self.get_property('pressure', atoms),)

    def get_temperature(self, atoms=None):
        '''Get the kinetic temperature of the system.'''

        return (self.get_property('temperature', atoms),)

    def init_variables(self):
        '''Instruct LAMMPS to assign properties to variables.

        This is needed to use get_lammps_energies() or
        get_lammps_stresses(), but it is preferable to use
        get_lammps_properties(), which calls the
        extract_compute() API instead of extract_variable().
        Variables are only assigned after a time step, which
        means a 'run 0' command is necessary.

        '''

        # energy variables
        self.lmp.command('variable pec equal pe')
        self.lmp.command('variable kec equal ke')
        self.lmp.command('variable etotalc equal etotal')

        # pressure variables; these require 'run post yes'
        self.lmp.command('variable pxx equal pxx')
        self.lmp.command('variable pyy equal pyy')
        self.lmp.command('variable pzz equal pzz')
        self.lmp.command('variable pxy equal pxy')
        self.lmp.command('variable pxz equal pxz')
        self.lmp.command('variable pyz equal pyz')
        self.lmp.command('variable press equal press')

        self.resultsv = {}

    def initialize_velocity(self, atoms, temperature, velocity_seed=-1,
                            group_selector=None):
        '''Initialize velocities and send back to atoms.'''

        if velocity_seed < 0:
            velocity_seed = np.random.random_integers(10000)

        if group_selector is None:
            lammps_group = 'all'
        else:
            self.lmp.command('group initv id ' + group_selector)
            lammps_group = 'initv'

        self.lmp.command('velocity {} create {} {}'.
                         format(lammps_group, temperature, velocity_seed))
        self.get_lammps_arrays()
        atoms.set_velocities(self.velocities*self._velocity_unit)

        if group_selector is not None:
            self.lmp.command('group initv clear')

    def minimize(self, atoms, algorithm='hftn', dmax=0.2,
                 etol=1e-6, ftol=0.0, maxiter=100, maxeval=5000):
        '''Perform a geometry optimization.

        Updated coordinates are passed back to the atoms
        object. 'atom_modify sort 0 0' is used to make sure
        that LAMMPS does not spatially sort the particles.

        '''

        self.set_configuration(positions=atoms.positions,
                               numbers=atoms.numbers,
                               cell=atoms.cell.diagonal())
        if self.natoms == 0:
            self.zero_properties()
            return

        # freeze atoms in the specified range by zeroing
        # forces on them
        nfrozen = 0
        if atoms.constraints:
            frozen = np.zeros(self.natoms, dtype=np.bool_)
            for c in atoms.constraints:
                if isinstance(c, FixAtoms):
                    frozen[c.index] = True

            # find continuous frozen indices so that we can
            # use as few LAMMPS commands as possible
            diffs = np.diff(np.hstack(([0], frozen, [0])))
            frozen_starts, = np.where(diffs > 0)
            frozen_ends, = np.where(diffs < 0)

            for s, e in zip(frozen_starts, frozen_ends):
                self.lmp.command('group tmp id <> {} {}'.format(s+1, e))
                self.lmp.command('group constrained union tmp')
                self.lmp.command('group tmp clear')

            self.lmp.command(
                'fix constraining constrained setforce 0.0 0.0 0.0')

        # optimization will stop if the 2-norm (length) of
        # the global force vector is less than the ftol
        ftol *= self.natoms - nfrozen

        # set neighbor list parameters, optimization algorithm
        self.lmp.command('atom_modify sort 0 0')
        self.lmp.command('neighbor 2.0 bin')
        self.lmp.command('neigh_modify delay 0 every 1 check yes once no')
        self.lmp.command('min_style {}'.format(algorithm))
        self.lmp.command('min_modify dmax {}'.format(dmax))

        # optimize to within specified criteria
        self.lmp.command('minimize {} {} {} {}'.
                         format(etol, ftol, maxiter, maxeval))

        if atoms.constraints:
            self.lmp.command('unfix constraining')
            self.lmp.command('group constrained clear')

        atoms.positions = np.ctypeslib.as_array(self.lmp.gather_atoms(
            'x', 1, 3)).reshape((self.natoms, 3))

    def run_md(self, atoms, nstep=1, temperature=1.0, pressure=1.0,
               ensemble='nve', timestep=0.5, integrator='verlet',
               T_damp=500, T_chain=1, T_loop=1, velocity_seed=None,
               P_damp=1000, P_chain=3, P_loop=1, isotropic=True, mtk=True):
        '''Perform a molecular dynamics run.

        Updated coordinates and velocities are passed back
        to the atoms object. 'atom_modify sort 0 0' is used
        to make sure that LAMMPS does not spatially sort the
        particles.

        '''

        self.set_configuration(positions=atoms.positions,
                               numbers=atoms.numbers,
                               cell=atoms.cell.diagonal(),
                               velocities=atoms.get_velocities())
        if self.natoms == 0:
            self.zero_properties()
            return

        # set ensemble and related parameters
        ensemble = ensemble.lower()
        if ensemble == 'nve':
            self.lmp.command('fix ensemble all nve')
        else:
            temp_str = (
                ' temp {} {} {} tchain {} tloop {} '.
                format(temperature, temperature, T_damp, T_chain, T_loop))

            if ensemble == 'nvt':
                self.lmp.command('fix ensemble all nvt' + temp_str)
            elif ensemble == 'npt':
                # whether the x, y, and z dimensions are
                # driven independently
                if isotropic:
                    isotropic = 'iso'
                else:
                    isotropic = 'aniso'

                # whether to add the correction terms due to
                # Martyna, Tuckerman, and Klein
                if mtk:
                    mtk = 'yes'
                else:
                    mtk = 'no'

                self.lmp.command(
                    'fix ensemble all npt' + temp_str +
                    '{} {} {} {} mtk {} pchain {} ploop {}'.
                    format(isotropic, pressure, pressure, P_damp, mtk, P_chain,
                           P_loop))
            else:
                raise RuntimeError('Ensemble not supported!')

        # set neighbor list parameters, time step, and integrator
        self.lmp.command('atom_modify sort 0 0')
        self.lmp.command('neighbor 2.0 bin')
        self.lmp.command('neigh_modify delay 10 every 1 check yes once no')
        self.lmp.command('timestep {}'.format(timestep))
        self.lmp.command('run_style {}'.format(integrator))

        # reinitialize velocity
        if velocity_seed is not None:
            self.initialize_velocity(temperature, velocity_seed)

        # run specified number of steps
        self.lmp.command('run {}'.format(nstep))
        self.lmp.command('unfix ensemble')

        self.get_lammps_properties()
        self.get_lammps_arrays()
        atoms.positions = self.positions
        atoms.set_velocities(self.velocities*self._velocity_unit)
        self.atoms = atoms.copy()

    def restore_state(self):
        '''Restore previous results and the atoms object.'''

        try:
            backup = self.backup
            self.results = backup['results']
            self.atoms = backup['atoms']
        except (AttributeError, TypeError, KeyError):
            raise RuntimeError('No state saved!')

    def save_state(self):
        '''Save current results and the associated atoms object.'''

        self.backup = {
            'results': self.results.copy(),
            'atoms': self.atoms,
        }

    def set_cell(self, cell=None, pbc=None, change_box=True):
        '''Set simulation cell parameters.

        Use the boundary or change_box command to set
        periodic boundary conditions and cell dimensions.

        '''

        box_cmd = ''

        if pbc is not None:
            box_cmd += 'boundary ' + ' '.join(
                ['p' if x else 's' for x in pbc])

        if cell is not None:
            box_cmd += ''.join([
                ' {} final 0 {}'.format(axis, cell[i])
                for i, axis in enumerate(['x', 'y', 'z'])])
            change_box = True
            self.cell = cell

        if len(box_cmd) > 0:
            if change_box:
                box_cmd = 'change_box all ' + box_cmd
            self.lmp.command(box_cmd)

    def set_configuration(self, positions=None, numbers=None, cell=None,
                          velocities=None, restart=True):
        '''Set system configuration in LAMMPS.'''

        if cell is not None:
            self.set_cell(cell=cell)

        if positions is not None and len(positions) != self.natoms:
            if restart:
                self.lmp.close()
                self.setup_lammps()
            self.set_number_atoms(len(positions))

        if numbers is None:
            atomtypes = None
        else:
            atom_types = self.parameters['atom_types']
            atomtypes = [atom_types[atomno] for atomno in numbers]

        velocities = (None if velocities is None
                      else (velocities/self._velocity_unit))

        self.set_lammps_arrays(positions=positions, atomtypes=atomtypes,
                               velocities=velocities)

    def set_lammps_arrays(self, positions=None, atomtypes=None,
                          velocities=None):
        '''Set LAMMPS arrays.'''

        if positions is not None:
            self.positions[:] = positions
            self.lmp.scatter_atoms('x', 1, 3, self._positions_ctypes)

        if atomtypes is not None:
            self.atomtypes[:] = atomtypes
            self.lmp.scatter_atoms('type', 0, 1, self._atomtypes_ctypes)

        if velocities is not None:
            self.velocities[:] = velocities
            self.lmp.scatter_atoms('v', 1, 3, self._velocities_ctypes)

    def set_number_atoms(self, num_new_atoms, remove=None):
        '''Change the number of atoms in the system.'''

        _logger.debug('Calculator N = %d -> %d.',
                      self.natoms, num_new_atoms)

        if num_new_atoms > self.natoms:
            self.lmp.command('create_atoms 1 random {} 1 NULL'.
                             format(num_new_atoms-self.natoms))
        elif num_new_atoms < self.natoms:
            if remove is None:
                remove = [num_new_atoms+1, self.natoms]
            self.lmp.command('group extra id <> {} {}'.
                             format(remove[0], remove[1]))
            self.lmp.command('delete_atoms group extra')
            self.lmp.command('group extra clear')
        else:
            return

        self.get_lammps_arrays()
        if self.natoms != num_new_atoms:
            raise RuntimeError('Number of atoms mismatch!')

    def setup_lammps(self):
        '''Create a LAMMPS instance.'''

        self.lmp = None

        try:
            atom_types = self.parameters['atom_types']
        except KeyError:
            raise KeyError('atom_types must be provided to map ASE atom'
                           '  symbols to LAMMPS atom types in pair_style.')

        try:
            command = self.parameters['command']
        except KeyError:
            raise KeyError('command must be provided to specify at least'
                           ' atom_style and pair_style.')

        if self.label is None:
            startup_options = ['-echo', 'screen', '-log', 'none']
        else:
            startup_options = ['-echo', 'screen', '-log', 'none',
                               '-screen', os.path.join(self.directory,
                                                       self.prefix)]

        self.lmp = lammps.lammps('', startup_options)
        # self.lmp.file(inputfile)

        # parameters that must be set before simulation box
        # is defined, including units and atom_styles
        for cmd in self.parameters['lammps_header'].splitlines():
            self.lmp.command(cmd.strip())

        # set up periodic boundary conditions
        self.set_cell(pbc=self.parameters['pbc'], change_box=False)

        # set up simulation cell and atom types
        self.lmp.command('region cell block {} {} {} {} {} {} units box'.
                         format(0.0, self.cell[0], 0.0, self.cell[1],
                                0.0, self.cell[2]))
        self.lmp.command('create_box {} cell'.format(len(atom_types)))
        for anum, atype in viewitems(atom_types):
            self.lmp.command('mass {} {}'.format(atype, atomic_masses[anum]))
        self.natoms = 0
        self.get_lammps_arrays()

        # pair_style and other output commands: thermo,
        # dump, restart, etc.
        for cmd in command.splitlines():
            self.lmp.command(cmd.strip())

        self.lmp.command('compute_modify thermo_temp dynamic yes')

    def zero_properties(self):
        self.results = {
            'energy': 0.0,
            'temperature': 0.0,
            'pressure': 0.0,
            'stress': np.zeros(6, dtype=np.float_),
            'forces': np.zeros((self.natoms, 3), dtype=np.float_)
        }
