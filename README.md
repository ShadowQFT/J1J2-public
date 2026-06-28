# J1-J2 simulation code and data

This repo has the Fortran code and data files used for the J1-J2 random-field
simulation figures.

## Contents

- `sim/`: Fortran simulation source code.
- `vendor/`: bundled FFT source files for generating correlated random fields.
- `data/raw/`: simulation output files.

## Build Notes

The main simulation source is `sim/J1J2_6tv.f90`. It depends on
`sim/corr_rand_3d.f90` and FFT sources under `vendor/fft2d/`.

Build commands are environment-specific. A typical local MPI Fortran build
starts from:

```bash
mpif90 -cpp sim/corr_rand_3d.f90 sim/J1J2_6tv.f90 \\
  vendor/fft2d/fftsg3d.f vendor/fft2d/fftsg.f -o J1J2
```

Simulation parameters are compile-time constants in `sim/J1J2_6tv.f90`.

## Run Parameters

Simulation parameters are compile-time constants in `sim/J1J2_6tv.f90`.  The values used in
a data set are listed in 'run_meta.txt.'

A few parameters need extra explanation:

- `PHI`: Random-field amplitude used by the code.  In the paper, the disorder
  strength is written as \(W\).  Because the paper defines the local nematic
  variable with a factor of \(1/4\), \(W = 4\,\texttt{PHI}\). 
  Thus runs with `PHI=2` correspond to \(W=8\).

- `XI`: Disorder correlation length.  In the paper this is written as
  \(\xi_d\).  All runs use `XI=3`.

- `START`: Initial spin configuration.  `HOT` means a random spin
  configuration.  `STR` means a stripe-ordered initial state.

- `STRAXIS`: Stripe direction for `START=STR`.  The listed value is the axis that
the stripes run parallel to. So `STRAXIS=X` means the stripes run along the \(x\)-axis, 
meaning there is alternation along the \(y\)-axis. 

- `NEQ`: Number of equilibration sweeps. In parallel tempering mode, this is the number of
"parallel-tempering cycles," with each cycle having a fixed number of sweeps then one swap attempt.

- `NMESS`: Number of measurement sweeps. Same parallel-tempering definition as `NEQ`.

- `PTSWAPEVER`: Number of sweeps between swap attempts.  For the paper data, swaps were attempted every five sweeps.

- `NREPLICA`: Number of temperatures/replicas in the parallel-tempering ladder.

- `NCONF`: Number of disorder configurations.

The main paper data set uses \(J_1=-J_2=J_\perp=1\), no site dilution,
`PHI=2` (\(W=8\)), `XI=3`, and parallel tempering over the temperature range
\(T=8.0\) to \(11.5\).

## Data Layout

Data files are organized under `data/raw/`. File names, directory names, and
headers indicate run settings and metadata. Each run folder also has a
`run_meta.txt` file with the parameters for that run.
