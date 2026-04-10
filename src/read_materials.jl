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
- `swrc_model_instance::Union{SWRCModel, Nothing}`: SWRC model instance for dispatch-based evaluation
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
    swrc_model_instance::Union{SWRCModel, Nothing}
    
    function WaterSoilProperties()
        # Dummy closures (will be replaced if hydraulic properties are provided)
        dummy_func = (x) -> 0.0
        new(0.0, "None", 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            dummy_func, dummy_func, dummy_func, dummy_func, dummy_func, nothing)
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
        raw_swrc_model = String(get(soil_info, "swrc_model", "None"))
        soil_props.water.swrc_model = replace(raw_swrc_model, " " => "_")
        soil_props.water.swrc_max_anw = Float64(get(soil_info, "swrc_max_anw", 0.0))
        soil_props.water.swrc_saturation_max_anw = Float64(get(soil_info, "swrc_saturation_max_anw", 0.0))
        soil_props.water.swrc_vg_alpha = Float64(get(soil_info, "swrc_vg_alpha", 0.0))
        soil_props.water.swrc_vg_n = Float64(get(soil_info, "swrc_vg_n", 0.0))
        soil_props.water.swrc_cav_delta = Float64(get(soil_info, "swrc_cav_delta", 0.0))
        
        # Store residual_water_content in water struct
        soil_props.water.residual_water_content = Float64(get(soil_info, "residual_water_content", 0.0))
        

        # Initialize water flow parameters 
        # Initialize SWRC parameters - model instance will be created at runtime
        if soil_props.water.swrc_model != "None"
            soil_props.water.theta_s = soil_props.porosity          
            soil_props.water.theta_r = soil_props.water.residual_water_content  
            soil_props.water.K_sat = 0.0  # Placeholder - computed at runtime
            soil_props.water.swrc_model_instance = nothing  # Will be assigned in compute_K_sat_runtime!()
        else
            soil_props.water.swrc_model = "None"
            soil_props.water.swrc_model_instance = nothing
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
        
        # CRITICAL: Create SWRC model instance with the newly-computed K_sat
        # Create model struct for method dispatch (struct+dispatch pattern)
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
                
                # Create SWRC model struct instance with updated K_sat (option: add directional K_s_x, K_s_y)
                swrc_params["K_sat_x"] = soil.water.K_sat_x
                swrc_params["K_sat_y"] = soil.water.K_sat_y
                
                soil.water.swrc_model_instance = create_swrc_model(soil.water.swrc_model, swrc_params)
            catch error_obj
                # Log but don't fail - leave as nothing
                @warn "Failed to create SWRC model instance for soil '$soil_name' after K_sat computation: $(string(error_obj))"
            end
             # Post-creation validation: fail hard if model was requested but not created
            if soil.water.swrc_model_instance === nothing 
                error("""SWRC Model Initialization Error: Soil '$soil_name' has swrc_model='$(soil.water.swrc_model)'
                        but the model instance could not be created.

                        Check that:
                        - The SWRC model name is valid: "Van_Genuchten" or "Cavalcante"
                        - Required parameters are non-zero (α, n for VG; δ for Cavalcante)
                        - Intrinsic permeability is positive (needed for K_sat computation)
                        - Liquid density and viscosity are positive

                        Current values:
                        K_sat = $(soil.water.K_sat)
                        theta_s = $(soil.water.theta_s), theta_r = $(soil.water.theta_r)
                        vg_alpha = $(soil.water.swrc_vg_alpha), vg_n = $(soil.water.swrc_vg_n)
                        cav_delta = $(soil.water.swrc_cav_delta)
                """)
            end
        end
    end
end


"""
    get_water_model(water_props::WaterSoilProperties) -> SWRCModel

Create a Richards equation water SWRC model from WaterSoilProperties.

Returns appropriate SWRC model instance (Van Genuchten, Cavalcante, etc.)
for use with the implicit Richards solver.

# Arguments
- `water_props::WaterSoilProperties`: Water material properties with SWRC model instance

# Returns
- `Union{SWRCModel, Nothing}`: Model struct instance or nothing if no model is assigned

# Example
\`\`\`julia
water_props = soil.water
model = get_water_model(water_props)
if model !== nothing
    K_val = K_h(model, h)      # Method dispatch
    theta_val = theta(model, h)
    C_val = C_moist(model, h)
end
\`\`\`

# Supported models
- `VanGenuchten`: Van Genuchten-Mualem model (ADSIM swrc_models.jl)
- `CavalcanteZornberg`: Cavalcante-Zornberg exponential model (ADSIM swrc_models.jl)
- `LinearSoil`: Linear retention curve (testing/verification only)
- `ConstantSoil`: Constant K and C coefficients (testing/verification only)
- `nothing`: No model assigned
"""
function get_water_model(water_props::WaterSoilProperties)::Union{SWRCModel, Nothing}
    return water_props.swrc_model_instance
end



# ══════════════════════════════════════════════════════════════════════════════
# Element-level water property helpers (used by Richards solvers)
# ══════════════════════════════════════════════════════════════════════════════
 
"""
    ElementWaterProps
 
Water properties for a single element, extracted from ADSIM materials.
Packages the SWRC model instance and anisotropic K_sat values for efficient
access during element-level assembly (avoids repeated dictionary lookups).
 
# Fields
- `model::SWRCModel` — SWRC model instance (for theta, C_moist, K_h, Se dispatch)
- `K_sat_x::Float64` — Saturated hydraulic conductivity in x [L/T]
- `K_sat_y::Float64` — Saturated hydraulic conductivity in y [L/T]
- `K_sat::Float64` — Isotropic K_sat (used to extract k_r = K_h/K_sat)
"""
struct ElementWaterProps
    model   :: SWRCModel
    K_sat_x :: Float64
    K_sat_y :: Float64
    K_sat   :: Float64
end
 
"""
    get_element_water_props(mesh, materials, elem_id::Int) → ElementWaterProps
 
Extract water properties for element `elem_id` from ADSIM materials.
Falls back to first soil if no material assignment exists.
 
# Arguments
- `mesh`: Mesh data structure
- `materials`: MaterialData structure
- `elem_id::Int`: Element ID
 
# Returns
- `ElementWaterProps`: Packaged water properties for this element
"""
function get_element_water_props(mesh, materials, elem_id::Int) :: ElementWaterProps
    mat_idx = get_element_material(mesh, elem_id)
    if mat_idx === nothing
        mat_idx = 1
    end
 
    soil_name = materials.soil_dictionary[mat_idx]
    soil = materials.soils[soil_name]
 
    model = soil.water.swrc_model_instance
    if model === nothing
        error("Element $elem_id: soil '$soil_name' has no SWRC model instance. " *
              "Ensure compute_K_sat_runtime! was called and swrc_model ≠ 'None'.")
    end
 
    return ElementWaterProps(model, soil.water.K_sat_x, soil.water.K_sat_y, soil.water.K_sat)
end
 
"""
    precompute_element_water_props(mesh, materials) → Vector{ElementWaterProps}
 
Precompute water properties for ALL elements at once.
Called once before the time loop to avoid repeated dictionary lookups during assembly.
 
# Arguments
- `mesh`: Mesh data structure
- `materials`: MaterialData structure
 
# Returns
- `Vector{ElementWaterProps}`: One entry per element
"""
function precompute_element_water_props(mesh, materials) :: Vector{ElementWaterProps}
    return [get_element_water_props(mesh, materials, e) for e in 1:mesh.num_elements]
end
 
 
# Export all public functions and types
export MaterialData, GasProperties, LiquidProperties, SoilProperties
export read_materials_file, get_gas_properties, get_soil_properties
export get_num_gases, get_num_soils, get_liquid_properties
export validate_swrc_parameters, compute_K_sat_runtime!, get_water_model
export ElementWaterProps, get_element_water_props, precompute_element_water_props
 