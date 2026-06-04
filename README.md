# Stokes Periodization — Numerical Tests

This repository contains the boundary integral solver and the numerical
experiments for periodic Stokes flow on the unit cell `[0,1]^3`, where
singly-, doubly-, and triply-periodic problems are handled by a precomputed
periodization operator built on top of the free-space fast multipole method.

The solver and all reported experiments accompany the paper:

> **A scalable Ewald-free BIE framework for periodic Stokes flow via hierarchical proxy sums**
> [arXiv:2605.30805](http://arxiv.org/abs/2605.30805)

The drivers in `test/` reproduce the figures and tables of the paper's
numerical-results section (see [Reproducing the paper results](#reproducing-the-paper-results)).

## Requirements

- C++17 compiler with OpenMP
- MPI (e.g. OpenMPI or MPICH, providing `mpicxx` and `mpirun`)
- Autotools (to build PVFMM)
- BLAS/LAPACK (e.g. Intel MKL) and, optionally, FFTW

On the Flatiron Institute clusters the required modules are listed in
`pvfmm_modules`; `source pvfmm_modules` before building or running.

## Building

```bash
# 1. Clone the repository together with its submodules (CSBQ and PVFMM)
git clone --recurse-submodules https://github.com/TianyueLi75/stokes-periodize-numtest.git
cd stokes-periodize-numtest

# If you already cloned without --recurse-submodules:
git submodule update --init --recursive

# 2. Build PVFMM (one time)
cd extern/pvfmm
./autogen.sh
./configure CXXFLAGS="-march=native -O3"
make -j
cd ../..

# 3. Build the test binaries into bin/
make test            # all test binaries, or build one with: make <target> (see below)
```

`make` reads PVFMM's compiler and library settings from
`extern/pvfmm/MakeVariables`. Edit the `CXXFLAGS` in the top-level `Makefile`
to point at your BLAS/LAPACK (and FFTW) installation if you are not using the
default MKL configuration.

## Running

Each binary is invoked through `mpirun`. Bind one OpenMP thread per core and
pass the per-test command-line arguments documented in the header comment of
the corresponding source file:

```bash
export OMP_NUM_THREADS=16
mpirun -n <Nproc> --map-by slot:pe=$OMP_NUM_THREADS ./bin/<binary> <args...>

# Example: triply-periodic self-convergence on a 25-sphere suspension
make test_selfconv
mpirun -n 1 --map-by numa:pe=$OMP_NUM_THREADS ./bin/test_convergence 2 3 25 0
```

Output `.vtu` visualization files are written to `vis/`, reference solutions
and logs to `out/`, and cached preconditioner matrices to `data/`. The
particle configurations used by the suspension tests are also provided under
`data/`.

## Reproducing the paper results

| Make target              | Binary                | Paper artifact |
|--------------------------|-----------------------|----------------|
| `make test_selfconv`     | `test_convergence`    | Self-convergence figure (channel+sphere, wall-bound particles, suspension) |
| `make test_manufactured_soln` | `test_manufactured_soln` | Manufactured-solution digits-of-accuracy table and the `(N_p, N_f)` parameter grid |
| `make timing_precomp`    | `precompute_time`     | Precomputation-time vs. multipole-order figure |
| `make timing_periodization` | `periodization_time` | Periodization-overhead table |
| `make timing`            | `timing`               | Weak/strong scaling table and the disturbance-flow streamline figures |
| `make examples`          | `examples`            | Converging-diverging channel and wall-bound loop visualizations |
| `make test_peri`         | `test_periodicity`    | Periodicity sanity check (velocity match across opposing cell faces) |

The exact command lines (and argument values) used to generate each result are
collected in `scripts/`:

- `scripts/test_runs.sh` — self-convergence, manufactured solution, precomputation,
  periodization overhead, streamlines, and the visualization examples;
- `scripts/weak_scaling.sh`, `scripts/strong_scaling.sh` — the scaling study
  (Slurm batch scripts; the per-process particle count and discretization are set inside).

Every driver's header comment documents its command-line arguments, the
quantity it computes, and the underlying formulation.

## Notes

- **Figure 7b is not reproduced by these drivers.** This figure reports the
  accuracy of the singly-periodic solver on the spherical-suspension geometry
  of the manufactured-solution test (see paper Section 4.3) as a function of the number of
  image levels $N_{\text{levels}}$. Because PVFMM must be recompiled for each
  value of $N_{\text{levels}}$, the data points were collected manually rather
  than from a single automated run.

- **Two precomputations occur on the first run.** First, the periodization
  operator is precomputed once per multipole order $m$ and cached to disk
  (see paper Sections 3 and 4.4), so subsequent runs at the same $m$ reuse it. Second, the
  block preconditioner is built and written to file for each discretization pair
  $(N_p, N_f)$. The preconditioner step may be disabled by commenting out the
  corresponding code; improving the efficiency of its construction is left to
  future work.