# Surge 3D Brady Engine

This is a modular C++17 3D hydrodynamics engine scaffold for building toward
Brady-style Stokesian Dynamics. Brady's monolayer is treated as a constrained
projection of a full 3D six-degree-of-freedom particle system.

## Modules

- `core.hpp`: particles, dense block conventions, strain flattening, periodic box helpers.
- `coefficients.*`: high-order coefficient table loader hook.
- `near_field.*`: two-body Jeffrey/Onishi pair-resistance bridge.
- `ewald.*`: open-space RPY and triply periodic Stokeslet Ewald translational mobility.
- `solver.*`: dense hydrodynamic assembly and force-balance solve.
- `time_integrator.*`: Euler and RK4 time integration.
- `projection.*`: Brady monolayer projection to `[Ux, Uy, Omega_z]`.
- `surge_pair_resistance.hpp`: local pair-resistance formulas used by the near-field module.
- `examples/two_sphere_trajectory.cpp`: trajectory driver.
- `examples/smoke_checks.cpp`: finite-value, symmetry, Ewald, projection, solve, and RK4 checks.
- `scripts/plot_trajectory.py`: CSV trajectory plotter.

## Requirements

On macOS:

```bash
brew install cmake eigen python
python3 -m pip install matplotlib
```

The `requirements.txt` file records the same project dependencies.

## Build

```bash
cmake -S . -B build
cmake --build build
```

Without CMake:

```bash
make EIGEN3_INCLUDE="$(brew --prefix eigen)/include/eigen3"
```

## Verify

```bash
./build/smoke_checks
```

Expected output:

```text
smoke checks passed
```

## Run

```bash
./build/two_sphere_trajectory two_sphere_trajectory.csv
python3 scripts/plot_trajectory.py two_sphere_trajectory.csv two_sphere_trajectory.png
```

## Scientific Status

Implemented now:

- clean 3D architecture;
- dense matrix assembly and solve with Eigen;
- Euler and RK4 time stepping;
- pairwise near-field resistance bridge using the Jeffrey/Onishi-derived header;
- open-space RPY translational mobility;
- triply periodic Stokeslet Ewald translational mobility backend;
- Brady-style 2D monolayer projection from the 3D generalized coordinates;
- smoke tests for finite matrix values, resistance symmetry, Ewald execution,
  projection size, force-balance solution, and RK4 stepping.

This is not yet a full paper-grade Brady/Stokesian Dynamics solver. Passing the
included checks proves that the software plumbing is coherent and runnable; it
does not prove numerical agreement with Brady, Jeffrey, or Onishi.

Still required for research-grade agreement:

- replace the printed `f_0...f_11` scalar truncation with recurrence/table-backed
  high-order Jeffrey/Onishi coefficients;
- implement the full far-field grand mobility beyond translation-translation:
  translation-rotation, rotation-rotation, and force/torque/stresslet couplings;
- implement the exact Stokesian Dynamics correction
  `R = R_far + R_2B_exact - R_2B_far`;
- validate every scalar family and tensor block against the source tables;
- add Brady's background vorticity convention for simple shear trajectories;
- add proper periodic or shear-periodic boundary updates for production runs.

