# TOML writing instructions for calculation parameters

This document shows the structure of the TOML file for the calculation parameters, units, and simulation settings.

## Header

The file contains "commented" info about the version of ADSIM used.
```toml
# ADSIM calculation parameters file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"
```

## Units and Dimensions

The following section defines the unit system used in the simulation.
```toml
[units]
geometry_unit = "m"
mass_unit = "kg"
pressure_unit = "kPa"
temperature_unit = "K"
time_unit = "s"
```

## Gravity Vector

Define gravity acceleration magnitude and components.
```toml
[gravity]
gravity_magnitude = 9.81
gravity_x_component = 0.0
gravity_y_component = -1.0
```

For Richards water flow, `gravity_x_component` and `gravity_y_component` define
the gravity direction used by the solver. Set both components to `0.0` to run a
water-flow calculation without gravity.

`gravity_magnitude` is a physical unit-conversion constant used to compute
saturated hydraulic conductivity, $K_{sat}=k\rho_w g/\mu_w$, and water pressure,
$P_w=\rho_w g h$. It must remain at the physical gravitational acceleration
(for example, `9.81` in SI units), even when the directional components are zero.
Setting `gravity_magnitude = 0.0` drives the computed $K_{sat}$ to zero.

The `[solver] gravity` flag does not control gravity in the Richards solver. It is
currently used only by `log_analysis_type` for the analysis banner. In particular,
do not use `gravity = 0` under `[solver]` to infer that Richards gravity is off;
the direction components above are authoritative.

For Richards water flow, an unspecified boundary has the natural zero-flux
condition everywhere in the domain, with no positional exceptions. Unit-gradient
free drainage must be requested explicitly with the `free_drainage_bc` mesh block;
an unconstrained bottom boundary is not treated as free drainage.

## Solver Settings

Define the solver type for the problem.
```toml
[solver]
solver_type = "2D-Plane"
```

## Time Stepping Data

Parameters for time-dependent calculations.
```toml
[time_stepping]
total_simulation_time = 1.0
time_per_step = 0.01
courant_number = 0.98
```

## Data Saving and Probing

Settings for data output and monitoring specific nodes/elements.
```toml
[data_saving]
data_saving_interval = 1.0

[probing]
# Nodes to follow
number_of_nodes = 3
nodes_to_probe = [1, 2, 3]

# Elements to follow
number_of_elements = 2
elements_to_probe = [1, 2]
```

## Complete Example

```toml
# ADSIM calculation parameters file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"

[units]
geometry_unit = "m"
mass_unit = "kg"
pressure_unit = "kPa"
temperature_unit = "K"
time_unit = "s"

[gravity]
gravity_magnitude = 9.81
gravity_x_component = 0.0
gravity_y_component = -1.0

[solver]
solver_type = "2D-Plane"

[time_stepping]
total_simulation_time = 100.0
time_per_step = 0.1
courant_number = 0.98

[data_saving]
data_saving_interval = 1.0

[probing]
# Nodes to follow during simulation
number_of_nodes = 3
nodes_to_probe = [10, 25, 50]

# Elements to follow during simulation (data at Gaussian points)
number_of_elements = 2
elements_to_probe = [5, 15]
```
