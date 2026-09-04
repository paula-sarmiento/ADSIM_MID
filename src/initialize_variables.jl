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

# Water BC masking (1 = free node, 0 = Dirichlet BC node)
global P_boundary_water::Vector{Int} = Int[]

# Water flux boundary conditions
global q_flux_water::Vector{Float64} = Float64[]

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
    global h, theta_w, S_r, P_water, v_water, P_boundary_water, q_flux_water
  
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
    
    # Initialize water BC mask (1 = free, 0 = BC constrained)
    P_boundary_water = ones(Int, Nnodes)
    
    # Initialize water flux boundary conditions
    q_flux_water = zeros(Float64, Nnodes)
    
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
            
            if soil_props !== nothing && soil_props.water.swrc_model_instance !== nothing
                wmodel = soil_props.water.swrc_model_instance
                for node_id in element_nodes
                    theta_w[node_id] = theta_ic
                    h[node_id] = h_inv(wmodel, theta_ic)
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
    apply_water_initial_conditions!(mesh, materials)

Apply all water initial conditions (t=0).

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties

# Priority
1. volumetric_content IC
2. pressure_head IC
"""
function apply_water_initial_conditions!(mesh, materials)
    apply_initial_water_volumetric_content!(mesh, materials)  # Priority 1
    apply_initial_water_pressure_head!(mesh, materials)       # Priority 2
end


"""
    apply_water_boundary_conditions!(mesh, materials)

Apply all water boundary conditions (t=0).

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties

# Priority
1. volumetric_content BC
2. pressure_head BC
3. flux BC
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
function get_node_material(mesh_data, node_id)
    """Find the material index for a node by locating its first adjacent element."""
    for elem_id in 1:mesh_data.num_elements
        element_nodes = get_element_nodes(mesh_data, elem_id)
        if node_id in element_nodes
            return get_element_material(mesh_data, elem_id)
        end
    end
    return nothing
end


function apply_water_volumetric_content_bc!(mesh, materials)
    global theta_w, h
    
    # Apply volumetric content BCs - convert θ to h (primary state for BCs)
    for (node_id, theta_bc) in mesh.volumetric_content_bc
        # Recover h from θ via SWRC inversion
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
    apply_water_pressure_head_bc!(mesh, materials)

Apply pressure head boundary conditions (t=0).

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties with SWRC

# Notes
- Modifies global theta_w and h
- Converts prescribed h → θ via SWRC
- Priority 2 (after volumetric content)
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

Apply water flux (Neumann) boundary conditions (t=0).

# Arguments
- `mesh`: Mesh structure with liquid_discharge_bc data

# Notes
- Modifies global q_flux_water
- Priority 3 (lowest - after head/content BCs)
- `q_flux_water` is no longer consumed by the Richards residual.
- `liquid_discharge_bc` reaches Richards through `assemble_water_neumann!`.
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
    apply_neumann_edge_richards!(R, node_i, node_j, q_bar, coords)

Add Neumann flux BC to residual for edge (node_i, node_j).

Edge-based assembly via 2-point Gauss quadrature.
**Note:** Assumes one constant `q_bar` per edge and is superseded in the Richards
solver by `assemble_water_neumann!`, which integrates nodally varying values.

# Arguments
- `R`: Residual vector (modified in-place)
- `node_i`, `node_j`: Edge node IDs
- `q_bar`: Prescribed volumetric flux [m/s]
- `coords`: Node coordinates [n_nodes, 2]

# Sign Convention
- q > 0: inflow (increases h)
- q < 0: outflow (decreases h)
"""
function apply_neumann_edge_richards!(
    R      :: Vector{Float64},
    node_i :: Int,
    node_j :: Int,
    q_bar  :: Float64,
    coords :: Matrix{Float64}
)
    # ── Compute edge length ──────────────────────────────────────────
    xi, yi = coords[node_i, 1], coords[node_i, 2]
    xj, yj = coords[node_j, 1], coords[node_j, 2]
    l_e = sqrt((xj - xi)^2 + (yj - yi)^2)

    # ── 2-point Gauss quadrature: ξ = ±1/√3, w = 1 ────────────────
    # Linear shape functions: N_i(ξ) = (1-ξ)/2, N_j(ξ) = (1+ξ)/2
    inv_sqrt3 = 1.0 / sqrt(3.0)
    for s in (-inv_sqrt3, inv_sqrt3)
        # Shape function values at Gauss point s
        N_i_s = 0.5 * (1.0 - s)  # Weight at node i
        N_j_s = 0.5 * (1.0 + s)  # Weight at node j
        
        # Integrate: ∫ q_bar * N dξ * (l_e / 2)
        # [l_e / 2 converts from isoparametric [-1,1] to physical]
        flux_contribution = q_bar * (l_e / 2.0)
        
        R[node_i] += N_i_s * flux_contribution
        R[node_j] += N_j_s * flux_contribution
    end
    
    return nothing
end


"""
    apply_water_dirichlet_bc!(mesh, materials)

Mark water Dirichlet BC nodes in P_boundary_water.

Sets P_boundary_water[node_id] = 0 for prescribed nodes.

# Priority
1. volumetric_content_bc
2. pressure_head_bc

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties
"""
function apply_water_dirichlet_bc!(mesh::MeshData, materials)
    global P_boundary_water
    
    # Priority 1: Volumetric content BC (highest priority - used before normalization)
    for (node_id, theta_bc) in mesh.volumetric_content_bc
        P_boundary_water[node_id] = 0  # Node is constrained
    end
    
    # Priority 2: Pressure head BC (includes converted volumetric_content_bc after normalization)
    for (node_id, h_bc) in mesh.pressure_head_bc
        # Only mark if not already marked by Priority 1
        if !haskey(mesh.volumetric_content_bc, node_id)
            P_boundary_water[node_id] = 0  # Node is constrained
        end
    end
    
    return nothing
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
    enforce_water_pressure_head_bc!(mesh, materials)

Enforce pressure head Dirichlet BCs during time-stepping.

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties

# Notes
- Modifies global theta_w and h
- Priority 2 (after volumetric content)
- Converts h → θ via SWRC
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
    SEEPAGE_MAX_FLIPS

Switch budget per node per time step for the potential seepage face. A node that has
changed state this many times is frozen for the remainder of the step; see
`resolve_seepage_face!`.
"""
const SEEPAGE_MAX_FLIPS = 10


"""
    init_seepage_face!(mesh) -> BitVector

Initialise the seepage-face partition with every candidate node dry.

Candidates start unsaturated and are opened by `resolve_seepage_face!` as the solution
requires. Starting closed is the conservative choice: a node that should be
discharging is opened on the first iteration that pushes its pressure head positive,
whereas starting open would inject water at nodes that are genuinely dry.

Throws if a candidate also carries a prescribed head, because the two conditions would
contend for the same node.
"""
function init_seepage_face!(mesh)
    for node in mesh.seepage_face_bc
        haskey(mesh.pressure_head_bc, node) &&
            error("node $node is both a seepage-face candidate and a prescribed-head node")
        haskey(mesh.volumetric_content_bc, node) &&
            error("node $node is both a seepage-face candidate and a volumetric-content node")
    end
    return falses(mesh.num_nodes)
end


"""
    resolve_seepage_face!(active, h_new, reaction, candidates, flips, max_flips) -> Int

Update the potential-seepage-face partition and return how many nodes switched state.

A seepage face is a boundary whose wet extent is part of the solution rather than an
input. Each candidate node is held in one of two states, chosen so that the two
physical conditions of a face open to the atmosphere — `h <= 0` everywhere, and no
inflow anywhere — both end up satisfied:

  * **discharging**: constrained to `h = 0`, water free to leave. If its reaction shows
    flow *entering* the domain, the face would be sucking water out of the air, so the
    node is released to zero flux.
  * **unsaturated**: left to the natural zero-flux condition. If the solve pushes its
    pressure head above zero, the node cannot stay dry, so it begins discharging.

Switching mutates the global `P_boundary_water` mask that the Richards assembly reads,
so the next assembly sees the revised partition. Iterating to a fixed point is what
locates the exit point of the free surface.

`flips` counts state changes within the current time step. A node that reaches
`max_flips` is frozen in its present state, which prevents the partition oscillating
between two configurations that each appear wrong from the other. The caller is
expected to report freezing rather than treat it as convergence.

# Arguments
- `active::BitVector`: Per-node discharging flag, mutated in place
- `h_new::Vector{Float64}`: Current pressure head, mutated at nodes that start discharging
- `reaction::Vector{Float64}`: Reactions from the most recent assembly. Physical inward
  flow at node `i` is `-reaction[i]`, and is only meaningful where the node was masked
  as Dirichlet, which is exactly the discharging set.
- `candidates`: Node ids forming the potential seepage face
- `flips::Vector{Int}`: Per-node switch counter for this time step, mutated in place
- `max_flips::Int`: Switch budget per node before freezing

# Returns
- `Int`: Number of nodes whose state changed. Zero means the partition is consistent.
"""
function resolve_seepage_face!(active::BitVector, h_new::Vector{Float64},
                               reaction::Vector{Float64}, candidates,
                               flips::Vector{Int}, max_flips::Int)::Int
    global P_boundary_water
    switched = 0

    for node in candidates
        flips[node] >= max_flips && continue

        if active[node]
            -reaction[node] > 0.0 || continue      # still discharging, leave it
            active[node] = false
            P_boundary_water[node] = 1             # release to natural zero flux
        else
            h_new[node] > 0.0 || continue          # still dry, leave it
            active[node] = true
            P_boundary_water[node] = 0             # hold at atmospheric pressure
            h_new[node] = 0.0
        end

        flips[node] += 1
        switched += 1
    end

    return switched
end


"""
    enforce_water_flux_bc!(mesh)

Enforce water flux (Neumann) BCs during time-stepping.

# Arguments
- `mesh`: Mesh structure with liquid_discharge_bc data

# Notes
- Modifies global q_flux_water
- Priority 3 (lowest - after head/content BCs)
- `q_flux_water` is no longer consumed by the Richards residual.
- `liquid_discharge_bc` reaches Richards through `assemble_water_neumann!`.
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
    enforce_water_dirichlet_bc!(mesh, materials)

Apply all water Dirichlet BCs with priority hierarchy.

# Arguments
- `mesh`: Mesh structure
- `materials`: Material properties

# Priority
1. volumetric_content_bc
2. pressure_head_bc
3. flux_bc

# Notes
- Modifies global theta_w, h, q_flux_water
- Called after each time step
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


"""
    apply_initial_conditions_water!(mesh, materials)

Apply water-only initial conditions (both IC and BC).
Preferred for water-only simulations (e.g., Richards solver verification tests).

This function applies water initial conditions and boundary conditions
without triggering multi-physics initializations (gas, temperature, lime, etc.).

# Arguments
- `mesh`: Mesh structure containing water IC/BC data
- `materials`: Material properties containing SWRC models

# Behavior
Applies in priority order:
1. Initial conditions (volumetric_content > pressure_head)
2. Boundary conditions (volumetric_content_bc > pressure_head_bc > flux_bc)

# Usage
For water-only tests:
```julia
apply_initial_conditions_water!(mesh, materials)
initialize_all_flows!(mesh, materials, mesh.num_nodes, 0)  # 0 gases
```
"""
function apply_initial_conditions_water!(mesh, materials)
    # Apply water initial conditions (priority: volumetric_content > pressure_head)
    apply_initial_water_volumetric_content!(mesh, materials)
    apply_initial_water_pressure_head!(mesh, materials)
    
    # Apply water boundary conditions (priority: volumetric_content_bc > pressure_head_bc > flux_bc)
    apply_water_boundary_conditions!(mesh, materials)
end


# ══════════════════════════════════════════════════════════════════════════════
# Implicit Richards solver helpers
# These functions support the implicit Richards solver but operate on
# ADSIM's global water state variables, so they belong here.
# ══════════════════════════════════════════════════════════════════════════════

"""
    build_dirichlet_lists(mesh, materials) → (dbc_nodes, dbc_vals)

Build Dirichlet BC node/value lists from ADSIM mesh for the implicit Richards solver.
Respects ADSIM's BC priority hierarchy:
  1. volumetric_content_bc (highest — convert θ→h via SWRC h_inv)
  2. pressure_head_bc (direct h value)

# Arguments
- `mesh`: Mesh data structure containing water BC data
- `materials`: MaterialData structure containing SWRC models

# Returns
- `dbc_nodes::Vector{Int}`: Node IDs with Dirichlet BCs
- `dbc_vals::Vector{Float64}`: Prescribed h values at those nodes
"""
function build_dirichlet_lists(mesh, materials)
    dbc_dict = Dict{Int, Float64}()

    # Priority 2: pressure_head_bc → direct h value
    for (node_id, h_val) in mesh.pressure_head_bc
        dbc_dict[node_id] = h_val
    end

    # Priority 1: volumetric_content_bc → convert θ to h (overrides pressure_head)
    for (node_id, theta_bc) in mesh.volumetric_content_bc
        # Find material for this node (use first adjacent element)
        mat_idx = nothing
        for e in 1:mesh.num_elements
            for a in 1:4
                if mesh.elements[e, a] == node_id
                    mat_idx = get_element_material(mesh, e)
                    break
                end
            end
            mat_idx !== nothing && break
        end

        if mat_idx !== nothing
            soil_name = materials.soil_dictionary[mat_idx]
            soil = materials.soils[soil_name]
            wmodel = soil.water.swrc_model_instance
            if wmodel !== nothing
                dbc_dict[node_id] = h_inv(wmodel, theta_bc)
            end
        end
    end

    dbc_nodes = collect(keys(dbc_dict))
    dbc_vals  = [dbc_dict[n] for n in dbc_nodes]

    return dbc_nodes, dbc_vals
end


"""
    update_water_globals!(elem_props, mesh, e_g, cache, rho_w, g_mag)

Update global water state variables (theta_w, S_r, P_water, v_water)
from the global h vector after a Richards solver time step.

Uses element water properties for multi-material support and
precomputed cache for efficient gradient computation.

# Arguments
- `elem_props::Vector{ElementWaterProps}`: Precomputed element water properties
- `mesh`: Mesh data structure
- `e_g::Vector{Float64}`: Gravity unit vector [gx, gy]
- `cache`: RichardsCache with precomputed Bp, detJ, Np
- `rho_w::Float64`: Liquid water density [kg/m³]
- `g_mag::Float64`: Gravity magnitude [m/s²]

# Note
- Modifies global variables: `theta_w`, `S_r`, `P_water`, `v_water`
- Reads global variable: `h`
- Uses anisotropic Darcy velocity: v = -K(h)(∇h - e_g)
"""
function update_water_globals!(elem_props, mesh, e_g::Vector{Float64},
                                cache, rho_w::Float64, g_mag::Float64)
    global h, theta_w, S_r, P_water, v_water

    Nn = mesh.num_nodes

    # Reset velocity accumulator
    fill!(v_water, 0.0)
    node_counts = zeros(Int, Nn)

    for e in 1:mesh.num_elements
        epr = elem_props[e]
        nodes = mesh.elements[e, :]
        h_e = h[nodes]

        # Update scalar fields (first element to touch a node wins)
        for a in 1:4
            nid = nodes[a]
            if node_counts[nid] == 0
                theta_w[nid] = theta(epr.model, h[nid])
                S_r[nid]     = Se(epr.model, h[nid])
                P_water[nid] = rho_w * g_mag * h[nid]
            end
            node_counts[nid] += 1
        end

        # Velocity: average from Gauss points via shape functions
        for p in 1:4
            Np = cache.Np[p]
            hp = dot(Np, h_e)

            # Anisotropic conductivity directly from SWRC model dispatch
            Kx = K_h_x(epr.model, hp)
            Ky = K_h_y(epr.model, hp)

            grad_h_x = cache.Bp[e,p,1,1]*h_e[1] + cache.Bp[e,p,1,2]*h_e[2] +
                       cache.Bp[e,p,1,3]*h_e[3] + cache.Bp[e,p,1,4]*h_e[4]
            grad_h_y = cache.Bp[e,p,2,1]*h_e[1] + cache.Bp[e,p,2,2]*h_e[2] +
                       cache.Bp[e,p,2,3]*h_e[3] + cache.Bp[e,p,2,4]*h_e[4]

            # Darcy velocity: v = -K (∇h - e_g)
            vx = -Kx * (grad_h_x - e_g[1])
            vy = -Ky * (grad_h_y - e_g[2])

            # For 2x2 Q4 quadrature, sum_p N_a(xi_p, eta_p) = 1 for each local node,
            # so each element contributes one velocity and node_counts tracks elements.
            for a in 1:4
                nid = nodes[a]
                v_water[nid, 1] += Np[a] * vx
                v_water[nid, 2] += Np[a] * vy
            end
        end
    end

    # Average velocity
    for i in 1:Nn
        if node_counts[i] > 0
            v_water[i, 1] /= node_counts[i]
            v_water[i, 2] /= node_counts[i]
        end
    end
end


"""
    compute_darcy_nodes_extrapolation(h::Vector{Float64}, mesh, elem_props::Vector{ElementWaterProps}, 
                                       cache, e_g::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}}

Compute Darcy velocity at every mesh node using Gauss-point extrapolation and nodal averaging.

# Algorithm
1. Build 4×4 extrapolation matrix E = P_C * inv(P_G) where:
   - P_G has rows [1, ξ_p, η_p, ξ_p*η_p] for each Gauss point p
   - P_C has rows [1, ξ_a, η_a, ξ_a*η_a] for each corner node a
2. For each element:
   - Compute velocity at 4 Gauss points
   - Extrapolate to 4 corner nodes via E
   - Accumulate contributions to global nodes
3. Average nodal contributions by node count

# Arguments
- `h::Vector{Float64}` : Nodal pressure heads, length num_nodes
- `mesh` : Mesh data with fields: coordinates (num_nodes×2), elements (num_elements×4), num_nodes, num_elements
- `elem_props::Vector{ElementWaterProps}` : Element SWRC properties (one per element)
- `cache` : RichardsCache with fields: Bp (num_elements×4×2×4), Np (4 vectors), detJ, A_e, weights
- `e_g::Vector{Float64}` : Gravity unit vector, e.g. [0.0, -1.0]

# Returns
- `(q_x, q_y)::Tuple{Vector{Float64}, Vector{Float64}}` : Darcy velocity components at each node

# Notes
- Uses anisotropic conductivity: K_x and K_y from model dispatch
- Velocity formula: v = -K(∇h - e_g) at Gauss points, then extrapolated to corners
- Follows zero-allocation pattern: all buffers pre-allocated outside element loop
"""
function compute_darcy_nodes_extrapolation(h::Vector{Float64}, mesh, 
                                            elem_props::Vector{ElementWaterProps}, 
                                            cache, e_g::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}}
    
    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 1: Build extrapolation matrix E = P_C * inv(P_G)
    # ─────────────────────────────────────────────────────────────────────────
    
    E = ShapeFunctions.build_q4_extrapolation_matrix()
    
    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 2: Allocate buffers (zero-allocation pattern)
    # ─────────────────────────────────────────────────────────────────────────
    
    Nn = mesh.num_nodes
    q_x = zeros(Float64, Nn)
    q_y = zeros(Float64, Nn)
    cnt = zeros(Int, Nn)
    
    # Pre-allocate element-local buffers
    h_e = zeros(Float64, 4)
    qG_x = zeros(Float64, 4)        # Velocity at Gauss points
    qG_y = zeros(Float64, 4)
    q_corners_x = zeros(Float64, 4) # Extrapolated velocity at corners
    q_corners_y = zeros(Float64, 4)
    grad_h = zeros(Float64, 2)      # Gradient [grad_h_x, grad_h_y]
    
    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 3: Element loop - compute velocity at Gauss points, extrapolate
    # ─────────────────────────────────────────────────────────────────────────
    
    for e in 1:mesh.num_elements
        # Gather nodal pressure heads for this element
        for a in 1:4
            h_e[a] = h[mesh.elements[e, a]]
        end
        
        # Get element SWRC model
        epr = elem_props[e]
        
        # Compute velocity at each Gauss point
        for p in 1:4
            # Interpolate pressure head at Gauss point p
            h_p = dot(cache.Np[p], h_e)
            
            # Get anisotropic hydraulic conductivity
            K_p_x = K_h_x(epr.model, h_p)
            K_p_y = K_h_y(epr.model, h_p)
            
            # Compute gradient via shape function derivatives
            # grad_h_x = sum(∂N_a/∂x * h_e[a])
            # grad_h_y = sum(∂N_a/∂y * h_e[a])
            grad_h[1] = 0.0
            grad_h[2] = 0.0
            for a in 1:4
                grad_h[1] += cache.Bp[e, p, 1, a] * h_e[a]
                grad_h[2] += cache.Bp[e, p, 2, a] * h_e[a]
            end
            
            # Compute Darcy velocity at Gauss point: v = -K(∇h - e_g)
            qG_x[p] = -K_p_x * (grad_h[1] - e_g[1])
            qG_y[p] = -K_p_y * (grad_h[2] - e_g[2])
        end
        
        # Extrapolate from Gauss points to corner nodes
        mul!(q_corners_x, E, qG_x)
        mul!(q_corners_y, E, qG_y)
        
        # Accumulate contributions to global nodes
        for a in 1:4
            nid = mesh.elements[e, a]
            q_x[nid] += q_corners_x[a]
            q_y[nid] += q_corners_y[a]
            cnt[nid] += 1
        end
    end
    
    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 4: Nodal averaging
    # ─────────────────────────────────────────────────────────────────────────
    
    for i in 1:Nn
        if cnt[i] > 0
            q_x[i] /= cnt[i]
            q_y[i] /= cnt[i]
        end
    end
    
    return q_x, q_y
end