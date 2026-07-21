# exact-fgr

**exact-fgr** is an open-source Fortran package for calculating nonadiabatic electron-transfer rate constants using the **Fermi Golden Rule (FGR)** method.

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
./exact_fgr > out.log &
```

or

```bash
./exact_fgr_grid > out.log &
```

where `out.log` is a user-defined file that stores the terminal output.

## Documentation

Complete documentation, including compilation instructions, input parameters, numerical convergence strategy, and descriptions of the output files, is available in **docs/User_Manual.md**.

## References

Citation information will be added later.
