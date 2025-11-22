# Time Step Module Documentation

## Overview
The `time_step.jl` module calculates the critical time step and manages time stepping for ADSIM simulations based on stability criteria from diffusion, advection, and chemical reaction processes.

## Module Structure

### Data Structure: `TimeStepData`
Stores all time stepping information:
- `critical_dt`: Critical time step from stability criteria [s]
- `actual_dt`: Actual time step used (critical_dt × Courant number) [s]
- `total_time`: Total simulation time requested by user [s]
- `num_steps`: Number of time steps required
- `courant_number`: Courant number (0 < C_N ≤ 1)
- `h_min`: Minimum characteristic element size [m]

## Critical Time Step Calculation

The critical time step is calculated as the minimum of three time scales:

### 1. Diffusive Time Scale
```
Δt_diffusion = (h_min² × τ) / (θ_g × D_max)
```
Where:
- `h_min`: Minimum characteristic element length [m]
- `τ`: Tortuosity factor [-]
- `θ_g = n - θ_w`: Volumetric gas content [-]
- `D_max`: Maximum diffusion coefficient [m²/s]

### 2. Advective Time Scale
```
Δt_advection = h_min² × (μ_g / (C_g^i × K × T × R))_min
```
Where:
- `μ_g`: Minimum gas dynamic viscosity [Pa·s]
- `C_g^i`: Maximum initial gas concentration [mol/m³]
- `K`: Intrinsic permeability [m²]
- `T`: Temperature [K]
- `R`: Universal gas constant = 8.314 J/(mol·K)

### 3. Reactive Time Scale
```
Δt_reaction = 1 / (κ × θ_w × C_CO2_max)
```
Where:
- `κ`: Reaction rate [1/s]
- `θ_w`: Volumetric water content [-]
- `C_CO2_max`: Maximum CO2 concentration [mol/m³]

### Final Time Step
```
Δt = C_N × min{Δt_diffusion, Δt_advection, Δt_reaction}
```
Where `C_N` is the Courant number (0 < C_N ≤ 1).

## Key Functions

### `calculate_element_characteristic_length(mesh, elem_id)`
Calculates the characteristic length of a quadrilateral element as the minimum distance between opposite side midpoints.

### `find_minimum_characteristic_length(mesh)`
Finds the minimum characteristic length across all elements in the mesh.

### `calculate_critical_time_step(mesh, materials, T_ref)`
Calculates the critical time step based on the three stability criteria.

### `calculate_time_step_info(mesh, materials, calc_params)`
Main function that:
1. Calculates the critical time step
2. Applies the Courant number to get the actual time step
3. Determines the number of time steps needed for the total simulation time
4. Returns a `TimeStepData` structure with all information

### `print_time_step_info(time_data, time_unit)`
Prints formatted time step information to the console.

## Usage in Kernel

The module is integrated into `kernel.jl` at Step 6:

```julia
include("time_step.jl")

# ...

# Step 6: Calculate time step information
time_data = calculate_time_step_info(mesh, materials, calc_params)
log_print("   ✓ Minimum characteristic length: $(time_data.h_min) m")
log_print("   ✓ Critical time step: $(time_data.critical_dt) s")
log_print("   ✓ Courant number: $(time_data.courant_number)")
log_print("   ✓ Actual time step: $(time_data.actual_dt) s")
log_print("   ✓ Number of time steps: $(time_data.num_steps)")
```

## Helper Functions

The module includes several helper functions to extract maximum/minimum values:
- `get_maximum_diffusion_coefficient(materials)`
- `get_minimum_gas_viscosity(materials)`
- `get_maximum_initial_concentration(mesh, NGases)`
- `get_maximum_co2_concentration(mesh, materials)`
- `get_minimum_permeability_ratio(mesh, materials, T_ref)`
- `get_maximum_reaction_parameters(mesh, materials, C_co2_max)`
- `get_minimum_gas_volume_fraction(materials)`

## Example Output

```
================================================================
TIME STEP INFORMATION
================================================================
Minimum characteristic length: 0.012919358897127609 m
Critical time step: 1.7120719658157846e-9 s
Courant number: 0.98
Actual time step: 1.677830526499469e-9 s
Total simulation time: 10 s
Number of time steps: 5960077518
================================================================
```

## Notes

1. The module automatically detects CO2 in the gas dictionary (case-insensitive)
2. Reference temperature is calculated as the mean of initial temperatures
3. All stability criteria are evaluated to ensure stable explicit time integration
4. The number of time steps is rounded up to ensure complete coverage of the simulation time
