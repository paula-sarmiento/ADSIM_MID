#------------------------------------------------------------------------------
# ADSIM variable initialization Module
# This module contains functions to initialize variables for ADSIM
# FEM calculations
#------------------------------------------------------------------------------

using Base.Threads

#------------------------------------------------------------------------------
# Global simulation variables - shared across all modules
#------------------------------------------------------------------------------

# Dimensional parameters
global NDim::Int = 0
global Nnodes::Int = 0
global Nelements::Int = 0
global NSoils::Int = 0
global NGases::Int = 0

# Boundary node influence lengths
global boundary_node_influences::Dict{Int, Float64} = Dict{Int, Float64}()

# State variables
# Gas transport state
global C_g::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global P::Vector{Float64} = Float64[]
global T::Vector{Float64} = Float64[]
global v::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global P_boundary::Matrix{Int} = Matrix{Int}(undef, 0, 0)
global λ_bc::Vector{Float64} = Float64[]

# Reactive species state
global C_lime::Vector{Float64} = Float64[]
global C_caco3::Vector{Float64} = Float64[]
global C_lime_residual::Vector{Float64} = Float64[]
global Caco3_max::Vector{Float64} = Float64[]

# Water flow state
global h::Vector{Float64} = Float64[]
global theta_w::Vector{Float64} = Float64[]
global S_r::Vector{Float64} = Float64[]
global P_water::Vector{Float64} = Float64[]
global v_water::Matrix{Float64} = zeros(Float64, 0, 2)

# Time derivatives
global dC_g_dt::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global dT_dt::Vector{Float64} = Float64[]
global dC_lime_dt::Vector{Float64} = Float64[]
global dtheta_dt::Vector{Float64} = Float64[]

#Analysis variables for soil carbonation
global binder_content::Vector{Float64} = Float64[]
global degree_of_carbonation::Vector{Float64} = Float64[]

#------------------------------------------------------------------------------
# Initialize all variables
#------------------------------------------------------------------------------
"""
    zero_variables!(mesh, materials)

Zero all simulation variables based on mesh and material data.
This function should be called once at the start of the simulation.
The exclamation mark indicates it modifies global variables.

# Arguments
- `mesh`: Mesh data structure containing node and element information
- `materials`: Material data structure containing soil and gas dictionaries
"""
function zero_variables!(mesh, materials)
    global NDim, Nnodes, Nelements, NSoils, NGases
    global C_g, P, T, v, P_boundary, λ_bc, boundary_node_influences
    global C_lime, C_caco3, C_lime_residual, binder_content, degree_of_carbonation, Caco3_max
    global dC_g_dt, dT_dt, dC_lime_dt, dtheta_dt
    global h, theta_w, S_r, P_water, v_water
  
    # Set dimensions
    NDim = 2  # Number of spatial dimensions - TODO: generalize for 3D
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    NSoils = length(materials.soil_dictionary)
    NGases = length(materials.gas_dictionary)
    
    # Allocate and initialize state variables
    C_g = zeros(Float64, Nnodes, NGases)
    P = zeros(Float64, Nnodes)
    T = zeros(Float64, Nnodes)
    v = zeros(Float64, Nnodes, NDim)
    P_boundary = ones(Int, Nnodes, NGases)  # 1 = free node, 0 = concentration BC node
    λ_bc = zeros(Float64, Nnodes)  # Lagrange multipliers for pressure BCs
    
    # Calculate and store boundary node influence lengths
    boundary_influences = get_boundary_node_influences(mesh)
    boundary_node_influences = boundary_influences.node_influences
    
    # Allocate and initialize reactive species
    C_lime = zeros(Float64, Nnodes)
    C_caco3 = zeros(Float64, Nnodes)
    C_lime_residual = zeros(Float64, NSoils)
    
    # Allocate and initialize time derivatives
    dC_g_dt = zeros(Float64, Nnodes, NGases)
    dT_dt = zeros(Float64, Nnodes)
    dC_lime_dt = zeros(Float64, Nnodes)    
    dtheta_dt = zeros(Float64, Nnodes)        

    # Allocate analysis variables
    binder_content = zeros(Float64, Nnodes)
    degree_of_carbonation = zeros(Float64, Nnodes)
    Caco3_max = zeros(Float64, Nnodes)

    # Allocate and initialize water state variables (Richards equation)
    h = zeros(Float64, Nnodes)           # Matric head [m]
    theta_w = zeros(Float64, Nnodes)     # Volumetric water content [-]
    S_r = zeros(Float64, Nnodes)         # Water saturation [-]
    P_water = zeros(Float64, Nnodes)     # Water pressure [Pa]
    v_water = zeros(Float64, Nnodes, NDim)  # Water velocity [m/s]
    
end


#------------------------------------------------------------------------------
# Apply initial conditions
#------------------------------------------------------------------------------
"""
    apply_initial_concentrations!(mesh::MeshData)

Apply initial gas concentrations from mesh data to the global C_g array.
This function reads element-based initial concentrations from the mesh and 
assigns them to all nodes within those elements.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial concentration data

# Note
- Modifies global variable `C_g`
- If a node belongs to multiple elements with different initial conditions,
  the last element's value will be used
"""
function apply_initial_concentrations!(mesh)
    global C_g, NGases
    
    # Apply element-based initial concentrations to nodes
    for (elem_id, concentrations) in mesh.initial_concentrations
        # Get nodes of this element
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Apply concentrations to each node of the element
        for node_id in element_nodes
            @threads for gas_idx in 1:NGases
                C_g[node_id, gas_idx] = concentrations[gas_idx]
            end
        end
    end    
end


"""
    apply_initial_temperature!(mesh::MeshData)

Apply initial temperatures from mesh data to the global T array.
This function reads element-based initial temperatures from the mesh and 
assigns them to all nodes within those elements.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial temperature data

# Note
- Modifies global variable `T`
- If a node belongs to multiple elements with different initial conditions,
  the last element's value will be used
"""
function apply_initial_temperature!(mesh)
    global T
    
    # Apply element-based initial temperatures to nodes
    for (elem_id, temperature) in mesh.initial_temperature
        # Get nodes of this element
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Apply temperature to each node of the element
        for node_id in element_nodes
            T[node_id] = temperature
        end
    end
end


"""
    apply_concentration_bc!(mesh::MeshData)

Apply concentration boundary conditions from mesh data to the global C_g array.
This function sets fixed gas concentrations at nodes where concentration 
boundary conditions are specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing concentration BC data

# Note
- Modifies global variable `C_g`
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_concentration_bc!(mesh)
    global C_g, NGases, P_boundary
    
    # Apply nodal concentration boundary conditions
    for (node_id, concentrations) in mesh.concentration_bc
        @threads for gas_idx in 1:NGases
            C_g[node_id, gas_idx] = concentrations[gas_idx]
            P_boundary[node_id, gas_idx] = 0  # Mark node as having concentration BC
        end
    end
end


"""
    apply_pressure_bc!(mesh::MeshData)

Apply absolute pressure boundary conditions from mesh data to the global P array.
This function sets fixed pressures at nodes where pressure boundary conditions 
are specified and restricts all gases at those nodes.

# Arguments
- `mesh::MeshData`: Mesh data structure containing pressure BC data

# Note
- Modifies global variables `P` and `P_boundary`
- P_boundary is set to 0 for all gases at pressure BC nodes
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_pressure_bc!(mesh)
    global P, P_boundary
    
    # Apply nodal pressure boundary conditions
    for (node_id, pressure) in mesh.absolute_pressure_bc
        P[node_id] = pressure         
        # Restrict all gases at this node
        #P_boundary[node_id, :] .= 0  # Mark node as having pressure BC
    end
end


"""
    apply_initial_lime_concentration!(mesh::MeshData, materials)

Apply initial lime concentrations from material properties to the global C_lime array.
This function loops through all elements, gets their assigned material, retrieves
the lime content from that material, and assigns it to all nodes in the element.

# Arguments
- `mesh::MeshData`: Mesh data structure containing element-material assignments
- `materials`: Material data structure containing soil properties with lime content

# Note
- Modifies global variable `C_lime`
- If a node belongs to multiple elements with different materials,
  the last element's value will be used
"""
function apply_initial_lime_concentration!(mesh, materials)
    global C_lime, C_lime_residual, Caco3_max
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material index for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get the soil name from the soil dictionary
            soil_name = materials.soil_dictionary[material_idx]
            
            # Get the soil properties for this material
            soil_props = get_soil_properties(materials, soil_name)
            
            if soil_props !== nothing
                # Get lime content from material
                β_l = soil_props.lime_content
                G_s = soil_props.specific_gravity
                n=soil_props.porosity                
                M_lime=74.093   # Molar mass of Ca(OH)2 in g/mol
                #Calculate lime concentration in mol/m^3 
                lime_concentration= (β_l * G_s * (1 - n) * 1e6 ) / M_lime #Asumes ρ_w= 1000 kg/m^3  

                #Calculatte reidual lime 
                residual_percent= soil_props.residual_lime
                C_lime_residual[material_idx] = residual_percent * lime_concentration

                #Calculate Caco3 max for degree of carbonation
                Caco3_max_concentration= lime_concentration  
                
                # Get nodes of this element
                element_nodes = get_element_nodes(mesh, elem_id)
                
                # Assign lime content and Caco3_max to each node of the element
                for node_id in element_nodes
                    C_lime[node_id] = lime_concentration
                    Caco3_max[node_id] = Caco3_max_concentration
                end
            end
        end
    end
end


"""
    apply_partial_pressure_bc!(mesh::MeshData)

Apply partial pressure boundary conditions from mesh data to the global C_g array.
This function sets gas concentrations at nodes where partial pressure boundary 
conditions are specified, using the ideal gas law: C_g[i] = P_partial[i] / (R * T).
It also marks these nodes in P_boundary to prevent the solver from updating them.

# Arguments
- `mesh::MeshData`: Mesh data structure containing partial pressure BC data

# Note
- Modifies global variables `C_g` and `P_boundary`
- Uses ideal gas law: P_partial = C_g * R * T, where R = 8.314 J/(mol·K)
- P_boundary is set to 0 for all gases at partial pressure BC nodes
- Concentrations will be dynamically updated in solver to maintain partial pressure
"""
function apply_partial_pressure_bc!(mesh)
    global C_g, NGases, P_boundary, T
    
    R = 8.314  # Universal gas constant [J/(mol·K)]
    
    # Apply nodal partial pressure boundary conditions
    for (node_id, partial_pressures) in mesh.partial_pressure_bc
        @threads for gas_idx in 1:NGases
            # Calculate concentration from partial pressure using ideal gas law
            # P_partial = C_g * R * T  =>  C_g = P_partial / (R * T)
            C_g[node_id, gas_idx] = partial_pressures[gas_idx] / (R * T[node_id])
            P_boundary[node_id, gas_idx] = 0  # Mark node as having partial pressure BC
        end
    end
end


"""
    apply_initial_water_volumetric_content!(mesh::MeshData, materials)

Apply initial volumetric water content from mesh data to global theta_w array.
Element-based specification - assigns values to all nodes of the element.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial water content data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- Element-level specification - prioritizes volumetric content over pressure head
- Recovers h from θ via SWRC inversion
"""
function apply_initial_water_volumetric_content!(mesh, materials)
    global theta_w, h
    
    for (elem_id, theta_ic) in mesh.initial_volumetric_content
        element_nodes = get_element_nodes(mesh, elem_id)
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            soil_name = materials.soil_dictionary[material_idx]
            soil_props = get_soil_properties(materials, soil_name)
            
            if soil_props !== nothing
                for node_id in element_nodes
                    theta_w[node_id] = theta_ic
                    h[node_id] = soil_props.water.h_theta(theta_ic)
                end
            end
        end
    end
end


"""
    apply_initial_water_pressure_head!(mesh::MeshData, materials)

Apply initial pressure head (matric head) from mesh data to global h array.
Element-based specification - converts h → θ via SWRC for element interior.
Only applied if volumetric content not specified for same element.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial pressure head data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- SECONDARY priority: Only used if initial_volumetric_content not specified
- Converts h → θ via SWRC inversion
"""
function apply_initial_water_pressure_head!(mesh, materials)
    global theta_w, h
    
    for (elem_id, h_ic) in mesh.initial_pressure_head
        # Skip if this element already has volumetric content IC (higher priority)
        if !haskey(mesh.initial_volumetric_content, elem_id)
            element_nodes = get_element_nodes(mesh, elem_id)
            material_idx = get_element_material(mesh, elem_id)
            
            if material_idx !== nothing
                soil_name = materials.soil_dictionary[material_idx]
                soil_props = get_soil_properties(materials, soil_name)
                
                if soil_props !== nothing && hasproperty(soil_props, :water)
                    wmodel = soil_props.water.swrc_model_instance
                    if wmodel !== nothing
                        for node_id in element_nodes
                            theta_w[node_id] = theta(wmodel, h_ic)
                            h[node_id] = h_inv(wmodel, theta_w[node_id])
                        end
                    end
                end
            end
        end
    end
end


"""
    apply_water_initial_conditions!(mesh::MeshData, materials)

Apply all water initial conditions from mesh and materials (t=0).
Calls individual functions in priority order.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure

# Note
- Modifies global variables: `theta_w`, `h`
- Priority: volumetric_content IC > pressure_head IC
"""
function apply_water_initial_conditions!(mesh, materials)
    apply_initial_water_volumetric_content!(mesh, materials)  # Priority 1
    apply_initial_water_pressure_head!(mesh, materials)       # Priority 2
end


"""
    apply_water_boundary_conditions!(mesh::MeshData, materials)

Apply all water boundary conditions from mesh (at t=0).
Calls individual BC functions in priority order.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure

# Note
- Modifies global variables: `theta_w`, `h`, `q_boundary_water`
- Priority: volumetric_content BC > pressure_head BC > flux BC
"""
function apply_water_boundary_conditions!(mesh, materials)
    apply_water_volumetric_content_bc!(mesh, materials)  # Priority 1
    apply_water_pressure_head_bc!(mesh, materials)       # Priority 2
    apply_water_flux_bc!(mesh)                           # Priority 3
end


"""
    apply_water_volumetric_content_bc!(mesh::MeshData, materials)

Apply volumetric water content boundary conditions at t=0.
This is the PRIMARY water BC - if specified, it takes precedence over pressure head BCs.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- HIGHEST priority: Directly sets θ values
- Recovers h from θ via SWRC inversion for consistency
"""
function apply_water_volumetric_content_bc!(mesh, materials)
    global theta_w, h
    
    # Apply volumetric content BCs - convert θ to h (primary state for BCs)
    for (node_id, theta_bc) in mesh.volumetric_content_bc
        h[node_id] = theta_bc  # First convert and store as h (primary)
        
        # Recover θ from h via SWRC for consistency
        material_idx = get_node_material(mesh, node_id)
        if material_idx !== nothing
            soil_name = materials.soil_dictionary[material_idx]
            soil_props = get_soil_properties(materials, soil_name)
            if soil_props !== nothing && hasproperty(soil_props, :water)
                wmodel = soil_props.water.swrc_model_instance
                if wmodel !== nothing
                    h[node_id] = h_inv(wmodel, theta_bc)  # Convert θ to h
                    theta_w[node_id] = theta(wmodel, h[node_id])  # Derive θ from h
                end
            end
        end
    end
end


"""
    apply_water_pressure_head_bc!(mesh::MeshData, materials)

Apply pressure head boundary conditions at t=0.
This is the SECONDARY water BC - applied only if volumetric content is not specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- SECONDARY priority: Only applied if node doesn't have volumetric_content_bc
- Converts prescribed h → θ via SWRC
"""
function apply_water_pressure_head_bc!(mesh, materials)
    global theta_w, h
    
    # Apply pressure head BCs (only if not already set by volumetric content BC)
    for (node_id, h_bc) in mesh.pressure_head_bc
        # Skip if this node already has a volumetric content BC (higher priority)
        if !haskey(mesh.volumetric_content_bc, node_id)
            material_idx = get_node_material(mesh, node_id)
            
            if material_idx !== nothing
                soil_name = materials.soil_dictionary[material_idx]
                soil_props = get_soil_properties(materials, soil_name)
                
                if soil_props !== nothing && hasproperty(soil_props, :water)
                    wmodel = soil_props.water.swrc_model_instance
                    if wmodel !== nothing
                        # Convert prescribed h → θ via SWRC
                        theta_w[node_id] = theta(wmodel, h_bc)
                        h[node_id] = h_bc
                    end
                end
            end
        end
    end
end


"""
    apply_water_flux_bc!(mesh)

Apply water flux (Neumann) boundary conditions at t=0.
This is the TERTIARY water BC - least restrictive, allows flexibility.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data

# Note
- LOWEST priority: Applied if neither volumetric_content_bc nor pressure_head_bc exist
- Flux BCs initialize target volumetric flux q [m/s] for boundary nodes
"""
function apply_water_flux_bc!(mesh)
    global q_flux_water
    
    # Apply volumetric water flux boundary conditions
    for (node_id, flux_bc) in mesh.liquid_discharge_bc
        # Skip if this node has a Dirichlet BC (higher priority)
        if !haskey(mesh.volumetric_content_bc, node_id) && !haskey(mesh.pressure_head_bc, node_id)
            q_flux_water[node_id] = flux_bc
        end
    end
end


"""
    enforce_water_volumetric_content_bc!(mesh::MeshData, materials)

Enforce volumetric water content Dirichlet boundary conditions during time-stepping.
This is the PRIMARY water BC - if specified, it takes precedence over pressure head BCs.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- HIGHEST priority: Directly enforces θ values
- Volumetric content BCs: Directly enforces theta_w[node_id] = θ_bc
- Recovers h from θ via SWRC inversion for consistency
"""
function enforce_water_volumetric_content_bc!(mesh, materials)
    global theta_w, h
    
    # Helper function: find node's material (use first adjacent element's material)
    function get_node_material_enforce(mesh_data, node_id)
        for elem_id in 1:mesh_data.num_elements
            element_nodes = get_element_nodes(mesh_data, elem_id)
            if node_id in element_nodes
                return get_element_material(mesh_data, elem_id)
            end
        end
        return nothing
    end
    
    # Enforce volumetric content BCs directly
    for (node_id, theta_bc) in mesh.volumetric_content_bc
        theta_w[node_id] = theta_bc
        
        # Recover h from θ via SWRC inversion for consistency
        material_idx = get_node_material_enforce(mesh, node_id)
        if material_idx !== nothing
            soil_name = materials.soil_dictionary[material_idx]
            soil_props = get_soil_properties(materials, soil_name)
            if soil_props !== nothing && hasproperty(soil_props, :h_theta)
                h[node_id] = soil_props.h_theta(theta_bc)
            end
        end
    end
end


"""
    enforce_water_pressure_head_bc!(mesh::MeshData, materials)

Enforce pressure head (matric head) Dirichlet boundary conditions during time-stepping.
This is the SECONDARY water BC - applied only if volumetric content is not specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`
- SECONDARY priority: Only applied if node doesn't have volumetric_content_bc
- Pressure head BCs: Convert prescribed h → θ via SWRC
- h[node_id] = h_bc (maintains prescribed pressure head)
"""
function enforce_water_pressure_head_bc!(mesh, materials)
    global theta_w, h
    
    # Helper function: find node's material (use first adjacent element's material)
    function get_node_material_enforce(mesh_data, node_id)
        for elem_id in 1:mesh_data.num_elements
            element_nodes = get_element_nodes(mesh_data, elem_id)
            if node_id in element_nodes
                return get_element_material(mesh_data, elem_id)
            end
        end
        return nothing
    end
    
    # Enforce pressure head BCs (only if not already set by volumetric content BC)
    for (node_id, h_bc) in mesh.pressure_head_bc
        # Skip if this node already has a volumetric content BC (higher priority)
        if !haskey(mesh.volumetric_content_bc, node_id)
            material_idx = get_node_material_enforce(mesh, node_id)
            
            if material_idx !== nothing
                soil_name = materials.soil_dictionary[material_idx]
                soil_props = get_soil_properties(materials, soil_name)
                
                if soil_props !== nothing && hasproperty(soil_props, :theta_h)
                    # Convert prescribed h → θ via SWRC
                    theta_w[node_id] = soil_props.theta_h(h_bc)
                    h[node_id] = h_bc
                end
            end
        end
    end
end


"""
    enforce_water_flux_bc!(mesh)

Enforce water flux (Neumann) boundary conditions during time-stepping.
This is the TERTIARY water BC - least restrictive, allows flexibility.

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data

# Note
- Modifies global variables: `q_flux_water` (via interaction with solver)
- LOWEST priority: Applied if neither volumetric_content_bc nor pressure_head_bc exist
- Flux BCs: Maintain specified volumetric flux q [m/s]
"""
function enforce_water_flux_bc!(mesh)
    global q_flux_water
    
    # Enforce volumetric water flux boundary conditions
    for (node_id, flux_bc) in mesh.liquid_discharge_bc
        # Skip if this node has a Dirichlet BC (higher priority)
        if !haskey(mesh.volumetric_content_bc, node_id) && !haskey(mesh.pressure_head_bc, node_id)
            q_flux_water[node_id] = flux_bc
        end
    end
end


"""
    enforce_water_dirichlet_bc!(mesh::MeshData, materials)

Apply all water Dirichlet boundary conditions with proper priority hierarchy.
This function calls the individual BC enforcement functions in priority order:
1. Volumetric content BC (highest priority)
2. Pressure head BC (secondary)
3. Flux BC (lowest priority)

# Arguments
- `mesh::MeshData`: Mesh data structure containing water BC data
- `materials`: Material data structure containing SWRC model closures

# Note
- Modifies global variables: `theta_w`, `h`, `q_flux_water`
- Called after each time step to re-enforce BC values
- Priority hierarchy ensures no conflicting BCs are applied to same node
"""
function enforce_water_dirichlet_bc!(mesh, materials)
    # Apply BCs in priority order
    enforce_water_volumetric_content_bc!(mesh, materials)  # Priority 1: θ
    enforce_water_pressure_head_bc!(mesh, materials)       # Priority 2: h
    enforce_water_flux_bc!(mesh)                           # Priority 3: q (Neumann)
end


"""
    apply_all_initial_conditions!(mesh::MeshData, materials)

Apply all initial conditions and boundary conditions from mesh and material data.
This is a convenience function that calls all individual application functions.

# Arguments
- `mesh::MeshData`: Mesh data structure containing all initial and boundary condition data
- `materials`: Material data structure containing soil and gas properties

# Note
- Modifies global variables: `C_g`, `T`, `P`, `C_lime`, `theta_w`, `h`, `q_boundary_water`
- Call this after `zero_variables!()` to set up the initial state
"""
function apply_all_initial_conditions!(mesh, materials)
    apply_initial_concentrations!(mesh)
    apply_initial_temperature!(mesh)
    apply_concentration_bc!(mesh)
    apply_partial_pressure_bc!(mesh)
    apply_pressure_bc!(mesh)
    apply_initial_lime_concentration!(mesh, materials)
    apply_water_initial_conditions!(mesh, materials)
    apply_water_boundary_conditions!(mesh, materials)
    
    println("\nAll initial conditions and BCs applied successfully")
end


