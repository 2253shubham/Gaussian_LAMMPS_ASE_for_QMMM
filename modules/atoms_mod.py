import copy

from ase.atoms import Atoms


class Atoms_mod(Atoms):
    # want to make a derived class of SimpleQMMM to run QMMM calcs using Gaussian for QM and LAMMPS for MM
    # changes which are made:

    ase_objtype = "atoms"

    def __init__(
        self,
        symbols=None,
        positions=None,
        numbers=None,
        tags=None,
        momenta=None,
        masses=None,
        magmoms=None,
        charges=None,
        scaled_positions=None,
        cell=None,
        pbc=None,
        celldisp=None,
        constraint=None,
        calculator=None,
        info=None,
        velocities=None,
    ):

        Atoms.__init__(
            self,
            symbols=symbols,
            positions=positions,
            numbers=numbers,
            tags=tags,
            momenta=momenta,
            masses=masses,
            magmoms=magmoms,
            charges=charges,
            scaled_positions=scaled_positions,
            cell=cell,
            pbc=pbc,
            celldisp=celldisp,
            constraint=constraint,
            calculator=calculator,
            info=info,
            velocities=velocities,
        )

    def copy(self):
        """Return a copy."""
        atoms = Atoms.copy(self)
        s = set(dir(atoms))
        diff = [x for x in dir(self) if x not in s]
        for i in diff:
            setattr(atoms, i, getattr(self, i))
        return atoms
