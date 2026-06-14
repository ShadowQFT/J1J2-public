# J1J2 Simulation Data Release

This repository contains the simulation source code and the raw data needed for
figures shown in the paper.

## Contents

- `sim/`: Fortran simulation sources.
- `vendor/`: bundled FFT sources required by the simulation build.
- `data/raw/`: raw simulation outputs needed for figures shown in the paper.

The private development repository also contains analysis scripts, CLI tooling,
cluster automation, generated plots, and paper drafting files. Those files are
not part of this public release repository.

## Build Notes

The main simulation source is `sim/J1J2_6tv.f90`. It depends on
`sim/corr_rand_3d.f90` and bundled FFT sources under `vendor/fft2d/`.

Build commands are environment-specific because the production runs used an HPC
cluster and MPI Fortran compiler. A typical local compiler invocation starts
from:

```bash
mpif90 -cpp sim/corr_rand_3d.f90 sim/J1J2_6tv.f90 \
  vendor/fft2d/fftsg3d.f vendor/fft2d/fftsg.f -o J1J2
```

Simulation parameters are compile-time constants in `sim/J1J2_6tv.f90`.

## Data Layout

Raw data are organized under `data/raw/`. File naming conventions and run
metadata are encoded in the raw output headers and directory names.

Generated figures and private analysis tooling are not included in this public
repository.
