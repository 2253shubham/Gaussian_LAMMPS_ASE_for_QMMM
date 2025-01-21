# Description of Scripts

## 1. `mm_lammps.py` 
Located in [`modules`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/modules)

### Overview
Derived class of ASE’s `LAMMPSlib`, which is ASE’s LAMMPS calculator ([reference](https://wiki.fysik.dtu.dk/ase/ase/calculators/lammpslib.html)). While the default `LAMMPSlib` is suitable for atomic systems, it has certain limitations:
- It cannot handle molecular systems or import bonded force-field parameters.

### Modifications
- Functions are added to compute bonds, angles, and dihedrals in molecular systems.
- Users must provide force-field parameter information in separate files.

---

## 2. `simpleqmmm_mod.py`
Located in [`modules`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/modules)

### Overview
Derived class of `SimpleQMMM`, the default class for performing QM/MM calculations using the subtractive scheme ([reference](https://wiki.fysik.dtu.dk/ase/ase/calculators/qmmm.html#simple-subtractive-qmmm-calculations)).

### Modifications
1. **Optimized QM Calculations**: Extract both energy and force information from a single Gaussian (QM) run instead of running QM calculations twice.
2. **Partial Charge Handling**: Reads Mulliken charges from the Gaussian QM output file for MM non-bonded interaction calculations.
3. **Framework Atom Freezing**: Adds a function to set the forces of selected atoms to zero, e.g., for freezing zeolite framework atoms.

---

## 3. `geo_opt.py`

### Overview
Script to perform geometry optimization of a given structure using QM/MM. Implements the damped dynamics algorithm for optimization.

### Inputs
- Charge, multiplicity, memory, and `nproc` for Gaussian calculations
- Atom indices for the QM region
- Atom indices fixed during simulation
- Structure input for optimization
- Basis set and functionals for Gaussian calculations
- Force convergence criteria (tolerance)

### HPC Support
Scripts available in [`bash_scripts`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/bash_scripts):
- `sub-qmmm-opt.script`
- `sub-qmmm-opt-perlmutter.script` (supports NERSC Perlmutter; can be modified for other HPCs)

### Running Calculations
1. Create a folder and copy the structure file (`.pdb` format) for optimization.
2. Link the required codes:
   - `geo_opt.py`, `mm_lammps.py`, `simpleqmmm_mod.py`, force-field parameter files
3. Include the input file `ref_data_opt.txt` from [`parameters`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/parameters).
4. Run locally or submit as a batch job on HPC. Ensure resource allocation aligns with Gaussian parameters (`memory` and `nproc`).

---

## 4. `nvtberendsen_mod.py`
Located in [`modules`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/modules)

### Overview
Derived class of `NVT_berendsen`, an ASE class for MD simulations ([reference](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nvtberendsen)).

### Modification
- Resets link atom positions along bonds broken when creating QM/MM boundaries. This prevents issues with dangling bonds during MD simulations.

---

## 5. `qmmm_metd.py`
### Overview
Script for QM/MM MD simulations assisted by metadynamics. Includes two stages:
1. Equilibration run
2. Production run with metadynamics (implemented via the Python `plumed` package)

### Inputs
- Charge, multiplicity, memory, and `nproc` for Gaussian calculations
- Atom indices for QM region
- Atom indices fixed during simulation
- Structure input for optimization
- Basis set and functionals for Gaussian calculations
- Force convergence criteria (tolerance)
- Metadynamics input files
- Number of steps for equilibration and production
- Temperature

### HPC Support
Scripts available in [`bash_scripts`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/bash_scripts):
- `sub-qmmm-metd.script`
- `sub-qmmm-metd-perlmutter.script` (supports NERSC Perlmutter; can be modified for other HPCs)

### Running Calculations
1. Create a folder and copy the optimized reactant state structure (output of `geo_opt.py`).
2. Link the required codes:
   - `qmmm_metd.py`, `nvtberendsen_mod.py`, `atoms_mod.py`, `mm_lammps.py`, `simpleqmmm_mod.py`, force-field parameter files
3. Include the input file `ref_data.txt` from [`parameters`](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/tree/main/parameters).
4. If using `plumed`, create the required input file (`plumed.dat`).
5. Submit the folder as a batch job on HPC, ensuring appropriate resource allocation for Gaussian parameters.

---
