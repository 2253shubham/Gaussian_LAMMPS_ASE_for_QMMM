from ase.md.nvtberendsen import NVTBerendsen
from ase.parallel import world
import numpy as np
from QM_MM_implemetation.bin_metd.modules.mm_lammps import init_atoms


def modify_link_atoms_positions(atoms):
    for i in range(len(atoms.nn_lnk_atms_idx)):
        if atoms.arrays["atomtypes"][atoms.nn_lnk_atms_idx[i]].startswith("C") and atoms.arrays["atomtypes"][atoms.nn_lnk_atms_idx_2[i]].startswith("C"):
            scale = 0.74
            atoms.positions[atoms.lnk_atms_idx[i]][0] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][0] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][0]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][0])
            atoms.positions[atoms.lnk_atms_idx[i]][1] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][1] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][1]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][1])
            atoms.positions[atoms.lnk_atms_idx[i]][2] = atoms.positions[atoms.nn_lnk_atms_idx_2[i]][2] + scale*(atoms.positions[atoms.nn_lnk_atms_idx[i]][2]-atoms.positions[atoms.nn_lnk_atms_idx_2[i]][2])
        else:
            continue
    return atoms


class NVTBerendsen_mod(NVTBerendsen):
 # want to make a derived class of SimpleQMMM to run QMMM calcs using Gaussian for QM and LAMMPS for MM
 # changes which are made:

    def __init__(self, atoms, timestep, temperature=None, taut=None,
                 fixcm=True, *, temperature_K=None,
                 trajectory=None, logfile=None, loginterval=1,
                 communicator=world, append_trajectory=False):
        
        NVTBerendsen.__init__(self, atoms, timestep, temperature=temperature, taut=taut,
                 fixcm=fixcm, temperature_K=temperature_K,
                 trajectory=trajectory, logfile=logfile, loginterval=loginterval,
                 communicator=communicator, append_trajectory=append_trajectory)

    
    def step(self, forces=None):
        """Move one timestep forward using Berenden NVT molecular dynamics."""
        self.scale_velocities()

        if not hasattr(self.atoms, "bonds"):
            self.atoms = init_atoms(self.atoms)

        # one step velocity verlet
        atoms = self.atoms

        if forces is None:
            forces = atoms.get_forces(md=True)

        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * forces

        if self.fix_com:
            # calculate the center of mass
            # momentum and subtract it
            psum = p.sum(axis=0) / float(len(p))
            p = p - psum

        self.atoms.set_positions(
            self.atoms.get_positions() +
            self.dt * p / self.atoms.get_masses()[:, np.newaxis])
        
        if hasattr(self.atoms, "nn_lnk_atms_idx") and hasattr(self.atoms, "nn_lnk_atms_idx_2"):
            self.atoms = modify_link_atoms_positions(self.atoms) # changes link atoms coordinates  
            atoms = modify_link_atoms_positions(atoms) # changes link atoms coordinates

        self.atoms.wrap(eps=1e-15)
        atoms.wrap(eps=1e-15)
  
        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.  For the same reason, we
        # cannot use self.masses in the line above.

        self.atoms.set_momenta(p)
        forces = self.atoms.get_forces(md=True)
        atoms.set_momenta(self.atoms.get_momenta() + 0.5 * self.dt * forces)
        return forces
