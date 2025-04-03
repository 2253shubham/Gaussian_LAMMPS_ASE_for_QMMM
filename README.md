# QMMM Calculations implementing Gaussian and LAMMPS with ASE

<p class="center-content"> 
  <img src="https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/blob/main/docs/Gaussian_LAMMPS_ASE.png" alt=""/>
</p>

[ASE + Gaussian ONIOM](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM) methodology, although is better than standalone Gaussian ONIOM approach, it suffers from insufficient options for available force-fields to model MM region and inability to implement periodic boundary conditions.

To address these, another algorithm was developed, which combines Gaussian with LAMMPS, integrated with ASE. 

The reason to include LAMMPS is twofold— it allows for implementation of periodic boundary conditions and that it supports more robust force-fields like OPLS, TraPPE, Amber. Thus, not only can we perform transition state search via NEB, but it also opens the avenue to perform enhanced sampling calculations like metadynamics/umbrella sampling.

This repository contains the scripts to perform QMMM calculations and implement them in geometry optimization and transition state searches. Gaussian is used to compute QM energy and forces, whereas the MM counterparts are computed using LAMMPS. The total energy and forces are then computed, which is utilized by ASE’s optimization algorithms to modify the structure. 


## Documentation

[Documentation](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/blob/main/docs/documentation.md)



## Authors

- [@shubham](https://github.com/2253shubham)


## Environment Variables

You may need to have an installed Gaussian software or if not, install/compile Gaussian ([version - g16](https://gaussian.com/gaussian16)) in your local system (or where you plan to run this code). Following which, the required environment variables can be set by the followiong lines:
```bash
module spider gaussian      # list installed Gaussian versions
module load gaussian        # load default version
```

You will also need [LAMMPS](https://docs.lammps.org/Manual.html) installed/compiled in your running unit. The environment variables can be set using the following:
```bash
module spider LAMMPS      # list installed LAMMPS versions
module load LAMMPS        # load default version
```

Required python version  - 3.9 or higher \
Download and install [requirements.txt](https://github.com/2253shubham/Gaussian_LAMMPS_ASE_for_QMMM/blob/main/requirements.txt) to install python library dependencies.
```bash
pip install -r /path/to/requirements.txt
```
