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
