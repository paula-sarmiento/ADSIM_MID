#------------------------------------------------------------------------------
# ADSIM Materials Reader Module
# This module contains functions to read and parse .toml material files for 
# ADSIM FEM calculations
#------------------------------------------------------------------------------

using TOML
include("swrc_models.jl")  # Import SWRC model factory

"""
GasProperties

Structure to store gas-specific properties.

# Fields
- `name::String`: Gas name
- `dynamic_viscosity::Float64`: Dynamic viscosity [Pa·s]
- `molar_mass::Float64`: Molar mass [g/mol]
- `diff_coefficient::Float64`: Diffusion coefficient [m²/s]
"""
mutable struct GasProperties
    name::String
    dynamic_viscosity::Float64
    molar_mass::Float64
    diff_coefficient::Float64
    
    function GasProperties(name::String)
        new(name, 0.0, 0.0, 0.0)
    end
end


"""
LiquidProperties

Structure to store liquid phase properties.

# Fields
- `dynamic_viscosity::Float64`: Dynamic viscosity [Pa·s]
- `density::Float64`: Density [kg/m³]
- `specific_heat::Float64`: Specific heat [J/(kg·K)]
"""
mutable struct LiquidProperties
    dynamic_viscosity::Float64
    density::Float64
    specific_heat::Float64
    
    function LiquidProperties()
        new(0.0, 0.0, 0.0)
    end
end


"""
WaterSoilProperties

Structure to store water-specific soil properties for Richards equation (vadose zone water flow).
Includes soil-water retention characteristics (SWRC), hydraulic conductivity, and water state variables.

# Fields
- `residual_water_content::Float64`: Residual water content [-]

- `swrc_model::String`: SWRC model name ("None", "Van_Genuchten", "Cavalcante")
- `swrc_max_anw::Float64`: Maximum air-water interfacial area [L⁻¹]
- `swrc_saturation_max_anw::Float64`: Saturation at maximum anw [-]
- `swrc_vg_alpha::Float64`: Van Genuchten fitting parameter α [1/m]
- `swrc_vg_n::Float64`: Van Genuchten fitting parameter n [-]
- `swrc_cav_delta::Float64`: Cavalcante fitting parameter δ [1/m]
- `theta_s::Float64`: Saturated water content [-]
- `theta_r::Float64`: Residual water content for SWRC [-]
- `K_sat::Float64`: Saturated hydraulic conductivity [m/s]
- `K_sat_x::Float64`: Saturated hydraulic conductivity x-direction (anisotropic) [m/s]
- `K_sat_y::Float64`: Saturated hydraulic conductivity y-direction (anisotropic) [m/s]

- `K_h::Function`: Closure K(h) - hydraulic conductivity as function of matric head [m/s]
- `theta_h::Function`: Closure θ(h) - volumetric water content as function of matric head [-]
- `h_theta::Function`: Closure h(θ) - matric head as function of volumetric water content [m]
- `c_s::Function`: Closure c_s(h) - water capacity as function of matric head [1/m]
- `D_w::Function`: Closure D_w(h) - water diffusivity as function of matric head [m²/s]
"""
mutable struct WaterSoilProperties
    residual_water_content::Float64
    swrc_model::String
    swrc_max_anw::Float64
    swrc_saturation_max_anw::Float64
    swrc_vg_alpha::Float64
    swrc_vg_n::Float64
    swrc_cav_delta::Float64
    theta_s::Float64
    theta_r::Float64
    K_sat::Float64
    K_sat_x::Float64
    K_sat_y::Float64
    K_h::Function
    theta_h::Function
    h_theta::Function
    c_s::Function
    D_w::Function
    
    function WaterSoilProperties()
        # Dummy closures (will be replaced if hydraulic properties are provided)
        dummy_func = (x) -> 0.0
        new(0.0, "None", 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            dummy_func, dummy_func, dummy_func, dummy_func, dummy_func)
    end
end


"""
SoilProperties

Structure to store soil material properties including mechanical, thermal, chemical, and hydraulic characteristics.

# Fields
- `name::String`: Soil name
- `specific_gravity::Float64`: Specific gravity of soil solids [-]
- `porosity::Float64`: Porosity [-]
- `saturation::Float64`: Degree of saturation [-]
- `granular_tortuosity::Float64`: Granular tortuosity [-]
- `intrinsic_permeability::Float64`: Intrinsic permeability [m²]
- `intrinsic_permeability_x::Float64`: Intrinsic permeability x-direction (anisotropic) [m²]
- `intrinsic_permeability_y::Float64`: Intrinsic permeability y-direction (anisotropic) [m²]
- `lime_content::Float64`: Lime content [-]
- `residual_lime::Float64`: Residual lime concentration [mol/m³] (reactive chemistry property)
- `reaction_rate::Float64`: Chemical reaction rate [1/s]
- `specific_heat_solids::Float64`: Specific heat of solids [J/(kg·K)]
- `water::WaterSoilProperties`: Water-related properties (SWRC, conductivity, closures)
"""
mutable struct SoilProperties
    name::String
    specific_gravity::Float64
    porosity::Float64
    saturation::Float64
    granular_tortuosity::Float64
    intrinsic_permeability::Float64
    intrinsic_permeability_x::Float64      # Anisotropic x-direction [m²]
    intrinsic_permeability_y::Float64      # Anisotropic y-direction [m²]
    lime_content::Float64
    residual_lime::Float64                 # Residual lime (reactive chemistry property)
    reaction_rate::Float64
    specific_heat_solids::Float64
    water::WaterSoilProperties             # Water-specific properties (SWRC, K_h, closures)
    
    function SoilProperties(name::String)
        new(name, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, WaterSoilProperties())
    end
end


"""
MaterialData

Structure to store all material properties including gases, liquids, and soils.

# Fields
- `gas_dictionary::Vector{String}`: List of gas names
- `gases::Dict{String, GasProperties}`: Gas properties indexed by name
- `liquid::LiquidProperties`: Liquid phase properties
- `soil_dictionary::Vector{String}`: List of soil names
- `soils::Dict{String, SoilProperties}`: Soil properties indexed by name
"""
mutable struct MaterialData
    gas_dictionary::Vector{String}
    gases::Dict{String, GasProperties}
    liquid::LiquidProperties
    soil_dictionary::Vector{String}
    soils::Dict{String, SoilProperties}
    
    function MaterialData()
        new(String[], 
            Dict{String, GasProperties}(),
            LiquidProperties(),
            String[],
            Dict{String, SoilProperties}())
    end
end


"""
read_materials_file(filename::String) -> MaterialData

Read an ADSIM materials .toml file and return a MaterialData structure containing
all material properties.

# Arguments
- `filename::String`: Path to the .toml materials file

# Returns
- `MaterialData`: Structure containing all material properties

# Example
```julia
materials = read_materials_file("materials.toml")
println("Number of gases: ", length(materials.gas_dictionary))
println("Number of soils: ", length(materials.soil_dictionary))
```
"""
function read_materials_file(filename::String)
    materials = MaterialData()
    
    # Read TOML file
    data = TOML.parsefile(filename)
    
    # Parse gas dictionary
    materials.gas_dictionary = data["gas_dictionary_"]
    
    # Parse gas properties
    parse_gas_properties!(materials, data["gas"])
    
    # Parse liquid properties
    parse_liquid_properties!(materials, data["liquid"])
    
    # Parse soil dictionary
    materials.soil_dictionary = data["soil_dictionary_"]
    
    # Parse soil properties
    parse_soil_properties!(materials, data["soil"])
    
    return materials
end


"""
parse_gas_properties!(materials::MaterialData, gas_data::Dict)

Parse gas properties from TOML data and store in MaterialData structure.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `gas_data::Dict`: Dictionary containing gas property data from TOML
"""
function parse_gas_properties!(materials::MaterialData, gas_data::Dict)
    for gas_name in materials.gas_dictionary
        gas_props = GasProperties(gas_name)
        gas_info = gas_data[gas_name]
        
        # Read gas properties
        gas_props.dynamic_viscosity = Float64(gas_info["dynamic_viscosity"])
        gas_props.molar_mass = Float64(gas_info["molar_mass"])
        gas_props.diff_coefficient = Float64(gas_info["diff_coefficient"])
        
        materials.gases[gas_name] = gas_props
    end
end


"""
parse_liquid_properties!(materials::MaterialData, liquid_data::Dict)

Parse liquid properties from TOML data and store in MaterialData structure.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `liquid_data::Dict`: Dictionary containing liquid property data from TOML
"""
function parse_liquid_properties!(materials::MaterialData, liquid_data::Dict)
    # Read liquid properties with reasonable defaults if not specified (0.0)
    # Default values match properties of water at room temperature (25°C, ~1 atm)
    materials.liquid.dynamic_viscosity = Float64(get(liquid_data, "dynamic_viscosity", 0.0))
    if materials.liquid.dynamic_viscosity == 0.0
        materials.liquid.dynamic_viscosity = 1.0e-3  # Pa·s (default: water viscosity)
    end
    
    materials.liquid.density = Float64(get(liquid_data, "density", 0.0))
    if materials.liquid.density == 0.0
        materials.liquid.density = 1000.0  # kg/m³ (default: water density)
    end
    
    materials.liquid.specific_heat = Float64(get(liquid_data, "specific_heat", 0.0))
    if materials.liquid.specific_heat == 0.0
        materials.liquid.specific_heat = 4186.0  # J/(kg·K) (default: water specific heat)
    end
end


"""
parse_soil_properties!(materials::MaterialData, soil_data::Dict)

Parse soil properties from TOML data and store in MaterialData structure.
Includes support for hydraulic properties that generate SWRC model closures.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `soil_data::Dict`: Dictionary containing soil property data from TOML
"""
function parse_soil_properties!(materials::MaterialData, soil_data::Dict)
    for soil_name in materials.soil_dictionary
        soil_props = SoilProperties(soil_name)
        soil_info = soil_data[soil_name]
        
        # Read physical properties
        soil_props.specific_gravity = Float64(soil_info["specific_gravity"])
        soil_props.porosity = Float64(soil_info["porosity"])
        soil_props.saturation = Float64(soil_info["saturation"])
        soil_props.granular_tortuosity = Float64(soil_info["granular_tortuosity"])
        soil_props.intrinsic_permeability = Float64(soil_info["intrinsic_permeability"])
        
        # Anisotropic permeability with fallback logic:
        # If x,y are provided (non-zero), use them; otherwise use isotropic value
        k_x = Float64(get(soil_info, "intrinsic_permeability_x", 0.0))
        k_y = Float64(get(soil_info, "intrinsic_permeability_y", 0.0))
        
        if k_x > 0.0
            soil_props.intrinsic_permeability_x = k_x
        else
            soil_props.intrinsic_permeability_x = soil_props.intrinsic_permeability
        end
        
        if k_y > 0.0
            soil_props.intrinsic_permeability_y = k_y
        else
            soil_props.intrinsic_permeability_y = soil_props.intrinsic_permeability
        end
        
        soil_props.lime_content = Float64(soil_info["lime_content"])
        soil_props.residual_lime = Float64(soil_info["residual_lime"])
        soil_props.reaction_rate = Float64(get(soil_info, "lime_reaction_rate", 0.0))
        
        # Read thermal properties
        soil_props.specific_heat_solids = Float64(soil_info["specific_heat_solids"])
        
        # Read SWRC properties (optional, backward compatible)
        soil_props.water.swrc_model = String(get(soil_info, "swrc_model", "None"))
        soil_props.water.swrc_max_anw = Float64(get(soil_info, "swrc_max_anw", 0.0))
        soil_props.water.swrc_saturation_max_anw = Float64(get(soil_info, "swrc_saturation_max_anw", 0.0))
        soil_props.water.swrc_vg_alpha = Float64(get(soil_info, "swrc_vg_alpha", 0.0))
        soil_props.water.swrc_vg_n = Float64(get(soil_info, "swrc_vg_n", 0.0))
        soil_props.water.swrc_cav_delta = Float64(get(soil_info, "swrc_cav_delta", 0.0))
        
        # Store residual_water_content in water struct
        soil_props.water.residual_water_content = Float64(get(soil_info, "residual_water_content", 0.0))
        

        # Initialize water flow parameters 
        # Initialize SWRC parameters - closures will be created at runtime
        if soil_props.water.swrc_model != "None"
            soil_props.water.theta_s = soil_props.porosity          
            soil_props.water.theta_r = soil_props.water.residual_water_content  
            soil_props.water.K_sat = 0.0  # Placeholder - computed at runtime
            
            # Set dummy closures (will be replaced in compute_K_sat_runtime!())
            dummy_func = (x) -> 0.0
            soil_props.water.K_h = dummy_func
            soil_props.water.theta_h = dummy_func
            soil_props.water.h_theta = dummy_func
            soil_props.water.c_s = dummy_func
            soil_props.water.D_w = dummy_func
        else
            soil_props.water.swrc_model = "None"
        end
        
        materials.soils[soil_name] = soil_props
    end
end


"""
get_gas_properties(materials::MaterialData, gas_name::String) -> Union{GasProperties, Nothing}

Get the properties for a specific gas by name.

# Arguments
- `materials::MaterialData`: Material data structure
- `gas_name::String`: Name of the gas

# Returns
- `Union{GasProperties, Nothing}`: Gas properties or nothing if not found

# Example
```julia
co2 = get_gas_properties(materials, "CO2")
if co2 !== nothing
    println("CO2 viscosity: ", co2.dynamic_viscosity)
end
```
"""
function get_gas_properties(materials::MaterialData, gas_name::String)
    return get(materials.gases, gas_name, nothing)
end


"""
get_soil_properties(materials::MaterialData, soil_name::String) -> Union{SoilProperties, Nothing}

Get the properties for a specific soil by name.

# Arguments
- `materials::MaterialData`: Material data structure
- `soil_name::String`: Name of the soil

# Returns
- `Union{SoilProperties, Nothing}`: Soil properties or nothing if not found

# Example
```julia
soil = get_soil_properties(materials, "Soil 1")
if soil !== nothing
    println("Porosity: ", soil.porosity)
end
```
"""
function get_soil_properties(materials::MaterialData, soil_name::String)
    return get(materials.soils, soil_name, nothing)
end


"""
get_num_gases(materials::MaterialData) -> Int

Get the total number of gases defined in the material data.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `Int`: Number of gases
"""
function get_num_gases(materials::MaterialData)
    return length(materials.gas_dictionary)
end


"""
get_num_soils(materials::MaterialData) -> Int

Get the total number of soils defined in the material data.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `Int`: Number of soils
"""
function get_num_soils(materials::MaterialData)
    return length(materials.soil_dictionary)
end


"""
get_liquid_properties(materials::MaterialData) -> LiquidProperties

Get the liquid phase properties.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `LiquidProperties`: Liquid properties

# Example
```julia
liquid = get_liquid_properties(materials)
println("Liquid density: ", liquid.density)
```
"""
function get_liquid_properties(materials::MaterialData)
    return materials.liquid
end


"""
validate_swrc_parameters(materials::MaterialData) -> Bool

Validate that SWRC model parameters are properly specified when a model is selected.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `Bool`: true if validation passes

# Throws
- `ErrorException`: If SWRC model is selected but required parameters are missing or invalid
"""
function validate_swrc_parameters(materials::MaterialData)
    for (soil_name, soil) in materials.soils
        if soil.water.swrc_model != "None"
            # Validate Van_Genuchten model
            if soil.water.swrc_model == "Van_Genuchten"
                if soil.water.swrc_vg_alpha <= 0.0 || soil.water.swrc_vg_n <= 0.0
                    error("""
                    SWRC Validation Error: '$soil_name' uses Van_Genuchten model but parameters are invalid.
                    
                    Required parameters for Van_Genuchten model:
                    - swrc_vg_alpha > 0 (current: $(soil.water.swrc_vg_alpha))
                    - swrc_vg_n > 0 (current: $(soil.water.swrc_vg_n))
                    
                    Please verify the SWRC parameters in your materials definition.
                    """)
                end
            
            # Validate Cavalcante model
            elseif soil.water.swrc_model == "Cavalcante"
                if soil.water.swrc_cav_delta <= 0.0
                    error("""
                    SWRC Validation Error: '$soil_name' uses Cavalcante model but parameters are invalid.
                    
                    Required parameters for Cavalcante model:
                    - swrc_cav_delta > 0 (current: $(soil.water.swrc_cav_delta))
                    
                    Please verify the SWRC parameters in your materials definition.
                    """)
                end
            
            else
                error("""
                SWRC Validation Error: '$soil_name' has unknown SWRC model: '$(soil.water.swrc_model)'
                
                Valid SWRC models are:
                - "None"
                - "Van_Genuchten"
                - "Cavalcante"
                """)
            end
        end
    end
    return true
end


"""
compute_K_sat_runtime!(materials::MaterialData, calc_params::Dict)

Compute saturated hydraulic conductivity K_sat for all soils at runtime,
after gravity and liquid properties are known.

K_sat = (k_intrinsic * ρ_water * g) / μ_water

# Arguments
- `materials::MaterialData`: Material data structure
- `calc_params::Dict`: Calculation parameters (contains gravity magnitude)

# Notes
This must be called after read_materials_file and after calc_params are loaded,
since K_sat depends on both intrinsic permeability (from materials) and gravity 
(from calculation parameters).
"""
function compute_K_sat_runtime!(materials::MaterialData, calc_params::Dict)
    # Get gravitational acceleration
    g = calc_params["gravity"]["magnitude"]
    
    # Get liquid properties
    rho_w = materials.liquid.density  # [kg/m³]
    mu_w = materials.liquid.dynamic_viscosity  # [Pa·s]
    
    # Loop through all soils and compute K_sat
    for (soil_name, soil) in materials.soils
        # Always compute K_sat - it's a fundamental property of the soil
        # Isotropic case (default)
        k_intrinsic = soil.intrinsic_permeability  # [m²]
        
        if k_intrinsic > 0.0
            soil.water.K_sat = (k_intrinsic * rho_w * g) / mu_w  # [m/s]
        else
            # Default fallback for zero or missing intrinsic permeability
            soil.water.K_sat = 1.0e-7
        end
        
        # Anisotropic case: compute directional K_sat values
        k_intrinsic_x = soil.intrinsic_permeability_x  # [m²]
        k_intrinsic_y = soil.intrinsic_permeability_y  # [m²]
        
        if k_intrinsic_x > 0.0
            soil.water.K_sat_x = (k_intrinsic_x * rho_w * g) / mu_w  # [m/s]
        else
            soil.water.K_sat_x = soil.water.K_sat  # Use isotropic value as fallback
        end
        
        if k_intrinsic_y > 0.0
            soil.water.K_sat_y = (k_intrinsic_y * rho_w * g) / mu_w  # [m/s]
        else
            soil.water.K_sat_y = soil.water.K_sat  # Use isotropic value as fallback
        end
        
        # CRITICAL: Recreate SWRC closures with the newly-computed K_sat
        # The closures were created with placeholder K_sat; now update them
        if soil.water.swrc_model != "None"
            try
                swrc_params = Dict{String, Float64}(
                    "K_sat" => soil.water.K_sat,
                    "theta_s" => soil.water.theta_s,
                    "theta_r" => soil.water.theta_r
                )
                
                if soil.water.swrc_model == "Van_Genuchten"
                    swrc_params["alpha"] = soil.water.swrc_vg_alpha
                    swrc_params["n_param"] = soil.water.swrc_vg_n
                elseif soil.water.swrc_model == "Cavalcante"
                    swrc_params["delta"] = soil.water.swrc_cav_delta
                end
                
                # Recreate SWRC closures with updated K_sat
                swrc_model = create_swrc_model(soil.water.swrc_model, swrc_params)
                soil.water.K_h = swrc_model.K_h
                soil.water.theta_h = swrc_model.theta_h
                soil.water.h_theta = swrc_model.h_theta
                soil.water.c_s = swrc_model.c_s
                soil.water.D_w = swrc_model.D_w
            catch error_obj
                # Log but don't fail - use previous closures
                @warn "Failed to update SWRC closures for soil '$soil_name' after K_sat computation: $(string(error_obj))"
            end
        end
    end
end


# Export all public functions and types
export MaterialData, GasProperties, LiquidProperties, SoilProperties
export read_materials_file, get_gas_properties, get_soil_properties
export get_num_gases, get_num_soils, get_liquid_properties
export validate_swrc_parameters, compute_K_sat_runtime!
