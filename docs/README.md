# User Manual

## 1. Introduction

The Exact FGR package is a Fortran program developed to calculate nonadiabatic electron-transfer rate constants using the Exact Fermi Golden Rule (FGR) method. The software is intended for studying electron-transfer reactions in condensed-phase systems, where interactions between the electronic states and the surrounding environment play an important role.

In addition to the Exact FGR method, the package also calculates rate constants using three approximations in the Exact FGR rate expression:

- High-Temperature (HT) approximation
- Short-Time (ST) approximation
- Both approximations at a time, i.e., Marcus theory

This allows users to compare the approximate methods directly with the numerically converged Exact FGR results and identify the conditions under which these approximations are valid.

The package contains two independent programs.

- **exact_fgr.f90** performs calculations for a single set of physical parameters. It is useful for studying individual systems and comparing the different theoretical methods.
- **exact_fgr_grid.f90** performs calculations over a three-dimensional parameter space consisting of the reorganization energy (λ), bath relaxation parameter (η), and reaction free energy (ΔG). This program is designed for large-scale calculations, such as generating Marcus validity maps and analyzing the performance of approximate theories over a wide range of parameters.

Both programs support OpenMP shared-memory parallelization to speed up the numerical calculations. The grid program also allows the parameter space to be divided into independent jobs, making it suitable for execution on high-performance computing (HPC) clusters.

This manual explains how to compile and run the programs, describes the input parameters and output files, and discusses the numerical methodology and convergence procedure used in the calculations.

---

## 2. Features

The Exact FGR package provides the following features:

- Calculates nonadiabatic electron-transfer rate constants using the Exact Fermi Golden Rule (FGR) method.
- Computes rate constants using the High-Temperature (HT), Short-Time (ST), and Marcus approximations.
- Performs both single-point and parameter-grid calculations.
- Automatically converges the numerical parameters for reliable results.
- Supports OpenMP parallel execution for faster calculations.
- Allows parameter-grid partitioning for large calculations on HPC clusters.
- Generates detailed output and convergence log files for further analysis.

---

## 3. Directory Structure

The repository is organized as follows:

```text
exact-fgr/
│
├── README.md
│
├── src/
│   ├── exact_fgr_sgpt.f90
│   ├── exact_fgr_grid.f90
│   └── README.md
│
├── docs/
│   ├── README.md
│   └── User_Manual.md
│
└── examples/
    └── README.md
```

- **src/** contains the Fortran source code.
- **docs/** contains the software documentation and user manual.
- **examples/** contains example calculations and descriptions of the output files.

---

## 4. Compilation

The source code is written in standard Fortran and supports OpenMP parallelization. It can be compiled using either the GNU Fortran compiler (**gfortran**) or the Intel Fortran compiler (**ifort**).

### Using GNU Fortran (gfortran)

Compile the single-point program:

```bash
gfortran -O3 -fopenmp exact_fgr_sgpt.f90 -o exact_fgr
```

Compile the parameter-grid program:

```bash
gfortran -O3 -fopenmp exact_fgr_grid.f90 -o exact_fgr_grid
```

### Using Intel Fortran (ifort)

Compile the single-point program:

```bash
ifort -O3 -qopenmp exact_fgr_sgpt.f90 -o exact_fgr
```

Compile the parameter-grid program:

```bash
ifort -O3 -qopenmp exact_fgr_grid.f90 -o exact_fgr_grid
```

Here, `-O3` enables compiler optimizations for improved performance, while `-fopenmp` (gfortran) and `-qopenmp` (ifort) enable OpenMP parallelization.

---

## 5. Running the Programs

Before running the programs, set the number of OpenMP threads according to the available CPU cores. For example, to use 8 threads:

```bash
export OMP_NUM_THREADS=8
```

### Running the Single-Point Program

Execute the single-point program using

```bash
./exact_fgr > out.log &
```

### Running the Parameter-Grid Program

Execute the parameter-grid program using

```bash
./exact_fgr_grid > out.log &
```
where `out.log` is a user-defined file that stores the terminal output.

The program performs calculations over the specified parameter grid and generates the corresponding output and convergence log files.

### Running on Multiple Jobs

For large parameter-grid calculations, the grid can be divided into independent jobs by modifying the `DG_start` and `DG_end` variables in `exact_fgr_grid.f90`. Each job computes a different portion of the parameter space, allowing multiple jobs to run simultaneously on an HPC cluster.

---

## 6. Input Parameters

The physical parameters and numerical settings are specified directly in the source code before compilation.

### Physical Parameters

| Parameter | Description | Unit |
|-----------|-------------|------|
| `lambda_cm` | Reorganization energy (λ) | cm⁻¹ |
| `DG_cm` | Reaction free energy (ΔG) | cm⁻¹ |
| `Vc_cm` | Electronic coupling (V) | cm⁻¹ |
| `Temp` | Temperature | K |
| `eta_cm` (single-point) | Bath relaxation parameter (η) | cm⁻¹ |

### Grid Parameters

The parameter-grid program scans over a range of physical parameters.

| Parameter | Description |
|-----------|-------------|
| `Nlambda` | Number of λ values |
| `Neta` | Number of η values |
| `NDG` | Number of ΔG values |
| `lambda_min`, `lambda_max` | Range of λ |
| `eta_min`, `eta_max` | Range of η |
| `DG_min`, `DG_max` | Range of ΔG |

### Convergence Criteria

The numerical calculations use the following convergence criteria:

- The recovered reorganization energy must agree with the input value within 5%.
- The Exact FGR rate constant is considered converged when successive calculations differ by less than 10%.

> **Note:** The reorganization energy recovery criterion is used only to obtain a good initial estimate of the frequency integration parameters (`wmax` and `dw`). The final numerical convergence is determined from the convergence of the Exact FGR rate constant.

---

## 7. Output Files

The programs generate output files containing the calculated rate constants and numerical convergence information.

### Single-Point Program

The single-point program displays the calculated quantities directly on the terminal. These include the calculated rate constants, recovered reorganization energy, convergence information, and other numerical details. No output files are generated.

### Parameter-Grid Program

The parameter-grid program generates the following files.

| File | Description |
|------|-------------|
| `validity_map.dat` | Calculated rate constants and other quantities for every point in the parameter grid. |
| `convergence.log` | Numerical convergence details for each grid point. |

The `convergence.log` file records the convergence history, including the values of `tmax`, `dt`, `wmax`, `dw`, and the recovered reorganization energy.

The data in `validity_map.dat` can be used for further analysis, visualization, and comparison of the Exact FGR, HT, ST, and Marcus rate constants.

---

## 8. Numerical Convergence Strategy

The Exact FGR rate constant is obtained by numerically evaluating the frequency and time integrals. To ensure reliable results, the program automatically converges the numerical parameters instead of requiring the user to choose them manually.

The convergence procedure consists of the following steps.

### 1. Frequency Grid Initialization

The program first estimates suitable values of `wmax` and `dw` by recovering the input reorganization energy. The frequency grid is accepted when the recovered reorganization energy agrees with the input value within 5%.

### 2. Automatic Convergence

Starting from the initial numerical parameters, the program performs successive convergence steps by doubling `tmax` and `wmax`, while halving `dt` and `dw`. After each step, the Exact FGR rate constant is recalculated and compared with the value from the previous step. The convergence process continues until the change in the Exact FGR rate constant is less than 10%.

The convergence history is written to `convergence.log`, allowing the user to monitor the evolution of `tmax`, `dt`, `wmax`, `dw`, the recovered reorganization energy, and the Exact FGR rate constant throughout the calculation.

---

## 9. Troubleshooting

Under construction.
