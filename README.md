# J1J2 Simulation Code and Data

This repo holds the simulation code and data files used for the J1-J2 random
field simulations.

## Contents

- `sim/`: Fortran simulation source code.
- `vendor/`: bundled FFT source files used by the simulation.
- `data/raw/`: simulation output files.

## Build Notes

The main simulation source is `sim/J1J2_6tv.f90`. It depends on
`sim/corr_rand_3d.f90` and bundled FFT sources under `vendor/fft2d/`.

Build commands are environment-specific. A typical local MPI Fortran build
starts from:

```bash
mpif90 -cpp sim/corr_rand_3d.f90 sim/J1J2_6tv.f90 \
  vendor/fft2d/fftsg3d.f vendor/fft2d/fftsg.f -o J1J2
```

Simulation parameters are compile-time constants in `sim/J1J2_6tv.f90`.

## Data Layout

Data files are organized under `data/raw/`. File names, directory names, and
headers carry the run settings and metadata. Each run folder also has a
`run_meta.txt` file with the parameters for that run.
