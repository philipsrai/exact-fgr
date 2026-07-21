# exact-fgr

**exact-fgr** is an open-source Fortran package for calculating nonadiabatic electron-transfer rate constants using the **Exact Fermi Golden Rule (FGR)** method.

The package computes rate constants using the following methods:

1. Exact Fermi Golden Rule (Exact FGR)
2. High-Temperature (HT) approximation
3. Short-Time (ST) approximation
4. Marcus theory

Two programs are included:

- **exact_fgr_sgpt.f90** – Performs single-point calculations.
- **exact_fgr_grid.f90** – Performs large-scale parameter-grid calculations for constructing Marcus validity maps.

Both programs support OpenMP shared-memory parallelization. The parameter-grid program can also be divided into independent jobs for efficient execution on HPC clusters.

## Repository Structure

```text
exact-fgr/
├── src/
├── docs/
└── examples/
```

## Getting Started

Compile the programs from the `src` directory using either **gfortran** or **ifort**.

### GNU Fortran (gfortran)

```bash
gfortran -O3 -fopenmp exact_fgr_sgpt.f90 -o exact_fgr
gfortran -O3 -fopenmp exact_fgr_grid.f90 -o exact_fgr_grid
```

### Intel Fortran (ifort)

```bash
ifort -O3 -qopenmp exact_fgr_sgpt.f90 -o exact_fgr
ifort -O3 -qopenmp exact_fgr_grid.f90 -o exact_fgr_grid
```

Run the executables using

```bash
./exact_fgr
```

or

```bash
./exact_fgr_grid
```

## Documentation

Complete documentation, including compilation instructions, input parameters, numerical convergence strategy, and descriptions of the output files, is available in **docs/User_Manual.md**.

## References

Citation information will be added later.
