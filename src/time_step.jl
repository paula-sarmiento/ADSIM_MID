#------------------------------------------------------------------------------
# ADSIM Time Step Module
# This module contains functions to calculate the critical time step and
# manage time stepping for ADSIM FEM calculations
#------------------------------------------------------------------------------

"""
TimeStepData

Structure to store time stepping information for the simulation.

# Fields
- `critical_dt::Float64`: Critical time step calculated from stability criteria [s]
- `actual_dt::Float64`: Actual time step used (critical_dt × Courant number) [s]
- `total_time::Float64`: Total simulation time requested by user [s]
- `time_per_step::Float64`: Time per load step (when data is saved) [s]
- `num_steps::Int`: Number of time steps required for the simulation
- `num_steps_per_load::Int`: Number of time steps per load step
- `num_load_steps::Int`: Number of load steps in the simulation
- `courant_number::Float64`: Courant number (0 < C_N ≤ 1)
- `h_min::Float64`: Minimum characteristic element size [m]
"""
mutable struct TimeStepData
    critical_dt::Float64
    actual_dt::Float64
    total_time::Float64
    time_per_step::Float64
    num_steps::Int
    num_steps_per_load::Int
    num_load_steps::Int
    courant_number::Float64
    h_min::Float64
    
    function TimeStepData()
        new(0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0)
    end
end


"""
    calculate_element_characteristic_length(mesh::MeshData, elem_id::Int)::Float64

Calculate the characteristic length of a quadrilateral element.
Uses the minimum distance between opposite sides of the element.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `elem_id::Int`: Element ID (1-based)

# Returns
- `Float64`: Characteristic length h [m]

# Formula
For a quadrilateral element, h is computed as the minimum of:
- Distance between midpoints of sides 1-2 and 3-4
- Distance between midpoints of sides 2-3 and 4-1
"""
function calculate_element_characteristic_length(mesh::MeshData, elem_id::Int)::Float64
    # Get the four nodes of the element
    nodes = mesh.elements[elem_id, :]
    
    # Get coordinates of all four nodes
    x1, y1 = mesh.coordinates[nodes[1], 1], mesh.coordinates[nodes[1], 2]
    x2, y2 = mesh.coordinates[nodes[2], 1], mesh.coordinates[nodes[2], 2]
    x3, y3 = mesh.coordinates[nodes[3], 1], mesh.coordinates[nodes[3], 2]
    x4, y4 = mesh.coordinates[nodes[4], 1], mesh.coordinates[nodes[4], 2]
    
    # Calculate midpoints of opposite sides
    # Side 1-2 midpoint
    mx12 = (x1 + x2) / 2.0
    my12 = (y1 + y2) / 2.0
    
    # Side 3-4 midpoint
    mx34 = (x3 + x4) / 2.0
    my34 = (y3 + y4) / 2.0
    
    # Side 2-3 midpoint
    mx23 = (x2 + x3) / 2.0
    my23 = (y2 + y3) / 2.0
    
    # Side 4-1 midpoint
    mx41 = (x4 + x1) / 2.0
    my41 = (y4 + y1) / 2.0
    
    # Distance between opposite side midpoints
    d1 = sqrt((mx34 - mx12)^2 + (my34 - my12)^2)
    d2 = sqrt((mx41 - mx23)^2 + (my41 - my23)^2)
    
    # Return minimum distance as characteristic length
    return min(d1, d2)
end


"""
    find_minimum_characteristic_length(mesh::MeshData)::Float64

Find the minimum characteristic length across all elements in the mesh.

# Arguments
- `mesh::MeshData`: Mesh data structure

# Returns
- `Float64`: Minimum characteristic length h_min [m]
"""
function find_minimum_characteristic_length(mesh::MeshData)::Float64
    h_min = Inf
    
    for elem_id in 1:mesh.num_elements
        h = calculate_element_characteristic_length(mesh, elem_id)
        h_min = min(h_min, h)
    end
    
    return h_min
end


"""
    suggest_implicit_water_dt(soil::SWRCModel, mesh::MeshData, F::Float64=50.0) -> Float64

Suggest an implicit time step for water-only (Richards equation) solver.

For water-only simulations without gas coupling, compute an aggressive time step
based on the explicit stability criterion scaled by amplification factor F.

# Arguments
- `soil::SWRCModel`: Any SWRC model instance (VanGenuchten, ConstantSoil, etc.)
- `mesh::MeshData`: Mesh data structure
- `F::Float64`: Time step amplification factor (default: 50.0)

# Returns
- `Float64`: Suggested time step Δt_impl [s]

# Formula
```
    dt_suggested = F × 0.5 × l_min² / (4 × D_max)
    where: D_max = K_h(soil, 0.0) / C_min
           K_h(soil, 0.0) = K_sat (saturation conductivity)
           C_min = 1e-4 m⁻¹ (regularization constant)
           l_min = minimum element characteristic length [m]
```

# Notes
- Backward Euler (implicit) is unconditionally stable in theory
- Δt is limited by precision and Picard convergence, not stability
- Uses polymorphic dispatch: works for any SWRC model without conditionals
- If Picard converges in ≤3 iterations: can increase F
- If Picard needs >7 iterations or diverges: reduce F to F/2

# Example
```julia
dt_suggested = suggest_implicit_water_dt(soil_model, mesh, 50.0)
```
"""
function suggest_implicit_water_dt(soil::SWRCModel, mesh::MeshData, F::Float64=50.0)
    # Regularization constant for C_moist (same as in swrc_models.jl)
    C_min = 1.0e-4  # [m⁻¹]
    CFL = 0.5
    
    # Get minimum element characteristic length
    l_min = find_minimum_characteristic_length(mesh)
    
    # Use polymorphic dispatch to compute K_sat and C_min from soil model
    # K_h at saturation (h=0) gives K_sat for any SWRC model
    K_sat = K_h(soil, 0.0)
    
    # C_moist at saturation (h=0) - for regularized version
    C_actual = C_moist(soil, 0.0)
    
    # Use maximum of regularization constant and actual C
    C_eff = max(C_min, C_actual)
    
    # Compute D_max = K_sat / C_eff
    D_max = K_sat / C_eff
    
    # Compute suggested time step: dt = F * CFL * l_min² / (4 * D_max)
    dt_suggested = F * CFL * l_min^2 / (4.0 * D_max)
    
    return dt_suggested
end


"""
    get_maximum_diffusion_coefficient(materials::MaterialData)::Float64

Get the maximum diffusion coefficient among all gases.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Maximum diffusion coefficient D_max [m²/s]
"""
function get_maximum_diffusion_coefficient(materials::MaterialData)::Float64
    D_max = 0.0
    
    for gas_name in materials.gas_dictionary
        gas = materials.gases[gas_name]
        D_max = max(D_max, gas.diff_coefficient)
    end
    
    return D_max
end


"""
    get_minimum_gas_viscosity(materials::MaterialData)::Float64

Get the minimum dynamic viscosity among all gases.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Minimum gas dynamic viscosity μ_g [Pa·s]
"""
function get_minimum_gas_viscosity(materials::MaterialData)::Float64
    μ_min = Inf
    
    for gas_name in materials.gas_dictionary
        gas = materials.gases[gas_name]
        μ_min = min(μ_min, gas.dynamic_viscosity)
    end
    
    return μ_min
end


"""
    get_maximum_initial_concentration(mesh::MeshData, NGases::Int)::Float64

Get the maximum initial gas concentration across all elements and gases.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `NGases::Int`: Number of gas species

# Returns
- `Float64`: Maximum initial concentration C_g^i [mol/m³]
"""
function get_maximum_initial_concentration(mesh::MeshData, NGases::Int)::Float64
    C_max = 0.0
    
    for (elem_id, concentrations) in mesh.initial_concentrations
        for gas_idx in 1:NGases
            C_max = max(C_max, concentrations[gas_idx])
        end
    end
    
    return C_max
end


"""
    get_maximum_co2_concentration(mesh::MeshData, materials::MaterialData)::Float64

Get the maximum CO2 concentration from initial conditions.
Assumes CO2 is one of the gases in the gas dictionary.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure

# Returns
- `Float64`: Maximum CO2 concentration C_CO2_max [mol/m³]
"""
function get_maximum_co2_concentration(mesh::MeshData, materials::MaterialData)::Float64
    # Find CO2 index in gas dictionary
    co2_idx = 0
    for (idx, gas_name) in enumerate(materials.gas_dictionary)
        if uppercase(gas_name) == "CO2" || uppercase(gas_name) == "CO₂"
            co2_idx = idx
            break
        end
    end
    
    # If CO2 not found, return 0
    if co2_idx == 0
        return 0.0
    end
    
    # Find maximum CO2 concentration
    C_co2_max = 0.0
    for (elem_id, concentrations) in mesh.initial_concentrations
        C_co2_max = max(C_co2_max, concentrations[co2_idx])
    end
    
    return C_co2_max
end


"""
    get_minimum_permeability_ratio(mesh::MeshData, materials::MaterialData, T_ref::Float64)::Float64

Get the minimum value of (μ_g / (C_g^i × K × T × R)) across all elements.
This is needed for the advective time scale calculation.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `T_ref::Float64`: Reference temperature [K]

# Returns
- `Float64`: Minimum permeability ratio [s/m²]

# Note
Uses universal gas constant R = 8.314 J/(mol·K)
"""
function get_minimum_permeability_ratio(mesh::MeshData, materials::MaterialData, T_ref::Float64)::Float64
    ratio_min = Inf
    
    # Get minimum gas viscosity
    μ_min = get_minimum_gas_viscosity(materials)
    
    # Get maximum initial concentration
    C_max = get_maximum_initial_concentration(mesh, length(materials.gas_dictionary))
    
    # Loop through all elements to find minimum ratio
    for elem_id in 1:mesh.num_elements
        # Get material for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get soil properties
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            
            # Calculate permeability ratio for this element
            # Using maximum concentration and minimum viscosity to get minimum ratio
            K = soil.intrinsic_permeability
            
            if K > 0.0 && C_max > 0.0
                ratio = μ_min / (4* C_max * K * T_ref * R)
                ratio_min = min(ratio_min, ratio)
            end
        end
    end
    
    return ratio_min
end


"""
    get_maximum_reaction_parameters(mesh::MeshData, materials::MaterialData, C_co2_max::Float64)::Float64

Get the maximum value of (κ × θ_w × C_CO2_max) for the reactive time scale.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `C_co2_max::Float64`: Maximum CO2 concentration [mol/m³]

# Returns
- `Float64`: Maximum reaction parameter [1/s]
"""
function get_maximum_reaction_parameters(mesh::MeshData, materials::MaterialData, C_co2_max::Float64)::Float64
    param_max = 0.0
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get soil properties
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            
            # Calculate reaction parameter
            κ = soil.reaction_rate  # Reaction rate [1/s]
            θ_w = soil.porosity * soil.saturation  # Volumetric water content [-]
            
            if C_co2_max > 0.0
                param = κ * θ_w * C_co2_max
                param_max = max(param_max, param)
            end
        end
    end
    
    return param_max
end


"""
    get_minimum_gas_volume_fraction(materials::MaterialData)::Float64

Get the minimum gas volume fraction θ_g = n - θ_w across all soils.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Minimum gas volume fraction θ_g [-]
"""
function get_minimum_gas_volume_fraction(materials::MaterialData)::Float64
    θ_g_min = Inf
    
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        n = soil.porosity
        θ_w = n * soil.saturation
        θ_g = n - θ_w
        
        θ_g_min = min(θ_g_min, θ_g)
    end
    
    return θ_g_min
end


"""
    calculate_critical_time_step(mesh::MeshData, materials::MaterialData, calc_params::Dict, T_ref::Float64)::Tuple{Float64, String}

Calculate the critical time step based on simulation type and stability criteria.

## Water-Only (Richards) Case
For water-flow-only simulations (no gas transport/reaction), use implicit Backward Euler
time stepping with polymorphic SWRC dispatch:
```
    Δt_crit = F × 0.5 × h_min² / (4 × D_max)
    where: D_max = K_h(soil, 0.0) / C_min
           K_h(soil, 0.0) = K_sat (saturation conductivity, any SWRC model)
           C_min = 1e-4 m⁻¹ (regularization constant)
           F = implicit_water_dt_factor (default 50.0)
```

## Multi-Physics Case
For coupled simulations, compute minimum of three stability criteria:
1. Diffusive time scale: h_min² × τ / (θ_g × D_max)
2. Advective time scale: h_min² × (μ_g / (C_g^i × K × T × R))_min
3. Reactive time scale: 1 / (κ × θ_w × C_CO2_max)

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `calc_params::Dict`: Calculation parameters (solver_settings, time_stepping)
- `T_ref::Float64`: Reference temperature [K]

# Returns
- `Tuple{Float64, String}`: (Critical time step [s], limiting scale description)
"""
function calculate_critical_time_step(mesh::MeshData, materials::MaterialData, calc_params::Dict, T_ref::Float64)::Tuple{Float64, String}
    # Check if this is water-flow-only case (Richards equation, no transport/reaction)
    solver_settings = get(calc_params, "solver_settings", Dict())
    water_flow = get(solver_settings, "water_flow", 0)
    diffusion = get(solver_settings, "diffusion", 0)
    advection = get(solver_settings, "advection", 0)
    reaction = get(solver_settings, "reaction_kinetics", 0)
    
    # WATER-ONLY CASE: Use implicit water dt calculation with polymorphic dispatch
    if water_flow == 1 && diffusion == 0 && advection == 0 && reaction == 0
        F = get(get(calc_params, "time_stepping", Dict()), "implicit_water_dt_factor", 50.0)
        
        # Get first soil with water properties for dt calculation
        dt_critical = Inf
        limiting_scale = "Water-only (implicit, no soil properties found)"
        
        for soil_name in materials.soil_dictionary
            soil = materials.soils[soil_name]
            if hasproperty(soil, :water) && soil.water !== nothing
                swrc_model = soil.water.swrc_model_instance
                if swrc_model !== nothing
                    dt_critical = suggest_implicit_water_dt(swrc_model, mesh, F)
                    F_int = isfinite(F) ? Int(round(F)) : 50
                    limiting_scale = "Implicit-water (F=$(F_int)× explicit, unconditionally stable)"
                    return dt_critical, limiting_scale
                end
            end
        end
        
        return dt_critical, limiting_scale
    end
    
    # MULTI-PHYSICS CASE: Calculate traditional critical time steps
    # Find minimum characteristic length
    h_min = find_minimum_characteristic_length(mesh)
    
    # Get maximum diffusion coefficient
    D_max = get_maximum_diffusion_coefficient(materials)
    
    # Get minimum gas volume fraction
    θ_g_min = get_minimum_gas_volume_fraction(materials)
    
    # Calculate diffusive time scale
    # Note: τ (tortuosity) is stored as granular_tortuosity
    # For simplicity, use the minimum tortuosity from all soils
    τ_min = Inf
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        τ_min = min(τ_min, soil.granular_tortuosity)
    end
    
    dt_diffusion = Inf
    if D_max > 0.0 && θ_g_min > 0.0
        dt_diffusion = (h_min^2 * τ_min) / (4 * θ_g_min * D_max)
    end
    
    # Calculate advective time scale
    dt_advection = Inf
    permeability_ratio_min = get_minimum_permeability_ratio(mesh, materials, T_ref)
    if permeability_ratio_min < Inf
        dt_advection = h_min^2 * permeability_ratio_min
    end
    
    # Calculate reactive time scale
    dt_reaction = Inf
    C_co2_max = get_maximum_co2_concentration(mesh, materials)
    reaction_param_max = get_maximum_reaction_parameters(mesh, materials, C_co2_max)
    if reaction_param_max > 0.0
        dt_reaction = 1.0 / (2 * reaction_param_max)
    end
    
    # Determine limiting time scale based on which is smallest
    limiting_scale = "Unknown"
    if dt_diffusion <= dt_advection && dt_diffusion <= dt_reaction
        limiting_scale = "Diffusive"
    elseif dt_advection <= dt_diffusion && dt_advection <= dt_reaction
        limiting_scale = "Advective"
    else
        limiting_scale = "Reactive"
    end
    
    # Return minimum of all three time scales and the limiting scale
    return (min(dt_diffusion, dt_advection, dt_reaction), limiting_scale)
end


"""
    calculate_time_step_info(mesh::MeshData, materials::MaterialData, calc_params::Dict)::TimeStepData

Calculate all time stepping information for the simulation.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `calc_params::Dict`: Calculation parameters including Courant number and total time

# Returns
- `TimeStepData`: Structure containing all time step information

# Example
```julia
time_info = calculate_time_step_info(mesh, materials, calc_params)
println("Critical time step: ", time_info.critical_dt, " s")
println("Actual time step: ", time_info.actual_dt, " s")
println("Number of steps: ", time_info.num_steps)
```
"""
function calculate_time_step_info(mesh::MeshData, materials::MaterialData, calc_params::Dict)::TimeStepData
    time_data = TimeStepData()
    
    # Get reference temperature from initial conditions (use maximum temperature)
    T_ref = 298.15  # Default to 298.15 K (25°C) if not specified

    T_max = 0.0
    for (elem_id, temp) in mesh.initial_temperature
        T_max = max(T_max, temp)
    end
    if T_max > 0.0
        T_ref = T_max
    end

    
    # Calculate critical time step (detects water-only case internally)
    time_data.critical_dt, limiting_scale = calculate_critical_time_step(mesh, materials, calc_params, T_ref)
    
    # Get Courant number from calculation parameters
    time_data.courant_number = calc_params["time_stepping"]["courant_number"]
    
    # Calculate actual time step (apply Courant number)
    time_data.actual_dt = time_data.critical_dt * time_data.courant_number
    
    # Get time per load step (when data is saved)
    time_data.time_per_step = calc_params["time_stepping"]["time_per_step"]
    
    # Store minimum characteristic length
    time_data.h_min = find_minimum_characteristic_length(mesh)
    
    # Calculate number of time steps per load step
    time_data.num_steps_per_load = ceil(Int, time_data.time_per_step / time_data.actual_dt)
    
    # Get total simulation time
    time_data.total_time = calc_params["time_stepping"]["total_simulation_time"]
    
    # Calculate number of load steps
    time_data.num_load_steps = ceil(Int, time_data.total_time / time_data.time_per_step)
    
    # Calculate total number of time steps for entire simulation
    time_data.num_steps = time_data.num_steps_per_load * time_data.num_load_steps
    
    return time_data, limiting_scale
end


# Export all public functions and types
export TimeStepData
export calculate_element_characteristic_length, find_minimum_characteristic_length
export suggest_implicit_water_dt
export get_maximum_diffusion_coefficient, get_minimum_gas_viscosity
export get_maximum_initial_concentration, get_maximum_co2_concentration
export calculate_critical_time_step, calculate_time_step_info
