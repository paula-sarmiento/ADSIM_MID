#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.x.x
# Author: Paula Sarmiento
#______________________________________________________

#______________________________________________________
# Fully explicit solver functions for ADSIM (Water Flow)
# Solves Richards equation using forward Euler time integration scheme
#______________________________________________________

using Statistics
using Base.Threads
using Printf
using LinearAlgebra

"""
    assemble_lumped_mass_vector_water!(M::Vector{Float64}, mesh, materials, h)

Assemble lumped mass vector for water flow.
Mass lumping accounts for nonlinear water capacity c_s(h).

# Arguments
- `M::Vector{Float64}`: Lumped mass vector to fill [Nnodes]
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `h::Vector{Float64}`: Current matric head [m]

# Formula
For each element: M_e = ∫ c_s(h) N dΩ ≈ c_s(h) × A_e / 4
where c_s(h) = dθ_w/dh is nonlinear capacity [1/m]
"""
function assemble_lumped_mass_vector_water!(M::Vector{Float64}, mesh, materials, h)
    fill!(M, 0.0)
    
    #______________________________________________________
    # Loop over all elements
    #______________________________________________________
    @threads for e in 1:mesh.num_elements
        # Get element nodes
        nodes = mesh.elements[e, :]
        
        # Get material properties for this element
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing
            error("Element $e has no material assigned. Check mesh material definitions.")
        end
        
        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        
        #______________________________________________________
        # Evaluate water capacity at element nodes (nonlinear!)
        #______________________________________________________
        h_e = [h[nodes[i]] for i in 1:4]
        c_s_e = [soil.c_s(h_e[i]) for i in 1:4]
        c_s_mean = mean(c_s_e)
        
        # Calculate element area using Gaussian quadrature
        A_e = 0.0
        for p in 1:4  # 4 Gauss points
            detJ = ShapeFunctions.get_detJ(e, p)
            w = ShapeFunctions.shape_funcs.gauss_weights[p]
            A_e += w * detJ
        end
        
        # Distribute mass equally to all 4 nodes (lumped mass)
        M_node = c_s_mean * A_e / 4.0
        
        # Add contribution to each node
        for i in 1:4
            node_id = nodes[i]
            M[node_id] += M_node
        end
    end
end

"""
    assemble_stiffness_matrices_water(mesh, materials, h, K_geometric)

Assemble element stiffness matrices for water flow including nonlinear K(h).

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `h::Vector{Float64}`: Current matric head [m]
- `K_geometric::Vector{Matrix{Float64}}`: Precomputed geometric stiffness matrices

# Returns
- `Vector{Matrix{Float64}}`: Element stiffness matrices with K(h) applied [Nelements][4×4]

# Formula
K_e(h) = ∫ K(h) ∇N · ∇N dΩ ≈ K(h) × K_geometric_e
"""
function assemble_stiffness_matrices_water(mesh, materials, h, K_geometric)
    Nelements = mesh.num_elements
    K_elements_water = Vector{Matrix{Float64}}(undef, Nelements)
    
    #______________________________________________________
    # Loop over all elements
    #______________________________________________________
    @threads for e in 1:Nelements
        # Get element nodes
        nodes = mesh.elements[e, :]
        
        # Get material properties for this element
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing
            error("Element $e has no material assigned. Check mesh material definitions.")
        end
        
        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        
        #______________________________________________________
        # Evaluate nonlinear K(h) at element nodes
        #______________________________________________________
        h_e = [h[nodes[i]] for i in 1:4]
        K_h_e = [soil.K_h(h_e[i]) for i in 1:4]
        K_h_mean = mean(K_h_e)
        
        # Apply K(h) to geometric stiffness
        K_elements_water[e] = K_h_mean * K_geometric[e]
    end
    
    return K_elements_water
end

"""
    check_courant_stability_water(mesh, materials, h, dt)

Verify Courant number for water flow stability.
Fo = D_max × Δt / (Δx_min)² ≤ 0.25

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `h::Vector{Float64}`: Current matric head [m]
- `dt::Float64`: Current time step [s]

# Returns
- `Tuple`: (is_stable::Bool, Fo::Float64, dt_max::Float64)

# Formula
Fo = max(D_w(h)) × Δt / (min(Δx))²
D_w(h) = K(h) / c_s(h) = water diffusivity [m²/s]
"""
function check_courant_stability_water(mesh, materials, h, dt)
    D_max = 0.0
    dx_min = Inf
    
    #______________________________________________________
    # Find maximum diffusivity and minimum element size
    #______________________________________________________
    @threads for e in 1:mesh.num_elements
        nodes = mesh.elements[e, :]
        
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing
            continue
        end
        
        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        
        # Maximum diffusivity in element
        h_e = [h[nodes[i]] for i in 1:4]
        D_e_all = [soil.D_w(h_e[i]) for i in 1:4]
        D_e_max = maximum(D_e_all)
        
        if isfinite(D_e_max) && D_e_max > 0.0
            D_max = max(D_max, D_e_max)
        end
        
        # Minimum element dimension
        for i in 1:3
            for j in i+1:4
                dist = norm(mesh.coords[nodes[i]] - mesh.coords[nodes[j]])
                dx_min = min(dx_min, dist)
            end
        end
    end
    
    # Compute Courant number
    Fo = D_max * dt / (dx_min^2)
    dt_max = 0.25 * (dx_min^2) / D_max
    is_stable = Fo <= 0.25
    
    return (is_stable, Fo, dt_max)
end

"""
    fully_explicit_richards_solver(mesh, materials, calc_params, time_data, project_name, log_print, initial_state=nothing)

Main fully explicit solver for water flow in porous media.
Solves the Richards equation using forward Euler time integration.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `calc_params`: Calculation parameters dictionary
- `time_data`: Time stepping data structure
- `project_name::String`: Project name for output files
- `log_print`: Function for logging output
- `initial_state`: Optional checkpoint state to resume from

# Governing Equation
∂θ_w/∂t = ∇·[K(h)∇h] + ∂K(h)/∂z

where:
- h = matric head [m]
- θ_w = volumetric water content [-]
- K(h) = nonlinear hydraulic conductivity [m/s]

# Time Integration (Forward Euler)
h^(n+1) = h^n + Δt × (1/M) × q_net^n

where:
- M = lumped mass matrix (depends on c_s(h))
- q_net = capillary + gravitational fluxes
"""
function fully_explicit_richards_solver(mesh, materials, calc_params, time_data, project_name, log_print, initial_state=nothing)
    log_print("\n[7/8] Starting fully explicit Richards solver for water flow")
    log_print("   Using $(Threads.nthreads()) threads for parallel execution")
    
    #______________________________________________________
    # Access global variables for water
    #______________________________________________________
    global h, theta_w, S_r, P_water, v_water
    global dh_dt, q_flux_water
    
    # Track warnings during simulation
    negative_head_warned = false
    courant_warning_issued = false
    
    #______________________________________________________
    # Physical constants
    #______________________________________________________
    ρ_water = 1000.0  # Water density [kg/m³]
    g_magnitude = 9.81  # Gravity magnitude [m/s²]
    γ_water = ρ_water * g_magnitude  # Specific weight [N/m³]
    
    #______________________________________________________
    # Get dimensions and settings
    #______________________________________________________
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    
    # Get solver settings
    solver_settings = calc_params["solver_settings"]
    calculate_capillary = solver_settings["capillary_flow"] == 1
    calculate_gravity = solver_settings["gravity_flow"] == 1
    
    # Get gravity vector
    gravity_params = calc_params["gravity"]
    g_magnitude = gravity_params["magnitude"]
    g_x = gravity_params["x_component"]
    g_y = gravity_params["y_component"]
    g_vector = [g_x, g_y]  # Unit vector
    
    # Time stepping parameters
    dt = time_data.actual_dt
    total_time = time_data.total_time
    load_step_time = time_data.time_per_step
    num_steps = time_data.num_steps
    
    #______________________________________________________
    # Initialize storage arrays
    #______________________________________________________
    M = zeros(Float64, Nnodes)  # Lumped mass vector (changes with h!)
    K_geometric = assemble_element_stiffness_matrices(mesh)  # Geometric part only
    
    # Initialize water variables from TOML or defaults
    if initial_state !== nothing
        # Load from checkpoint
        h = copy(initial_state.h)
        current_time = initial_state.current_time
        output_counter = initial_state.output_counter
        next_output_time = initial_state.next_output_time
        log_print("      Continuing from time: $(current_time) $(calc_params["units"]["time_unit"])")
    else
        # Start from initial conditions
        h_init = get(calc_params["initial_conditions"], "h_initial", -1.0)
        h = fill(h_init, Nnodes)
        
        current_time = 0.0
        next_output_time = load_step_time
        output_counter = 1
        
        log_print("      Load step 0 (0.0%)")
        # write_output_vtk_water(mesh, materials, 0, 0.0, project_name)
    end
    
    # Initialize secondary water variables
    theta_w = zeros(Float64, Nnodes)
    S_r = zeros(Float64, Nnodes)
    P_water = zeros(Float64, Nnodes)
    v_water = zeros(Float64, Nnodes, 2)
    dh_dt = zeros(Float64, Nnodes)
    q_flux_water = zeros(Float64, Nnodes)
    
    #______________________________________________________
    # Main time stepping loop
    #______________________________________________________
    save_data = false
    for step in 1:num_steps
        
        #______________________________________________________
        # Update nonlinear mass matrix based on current h
        #______________________________________________________
        assemble_lumped_mass_vector_water!(M, mesh, materials, h)
        
        # Check for zero or very small masses
        if any(M .<= 1.0e-15)
            error("Lumped mass vector contains zero or negative values at step $step!")
        end
        
        #______________________________________________________
        # Update nonlinear stiffness matrices based on current K(h)
        #______________________________________________________
        K_elements_water = assemble_stiffness_matrices_water(mesh, materials, h, K_geometric)
        
        #______________________________________________________
        # Check Courant stability
        #______________________________________________________
        (is_stable, Fo, dt_max) = check_courant_stability_water(mesh, materials, h, dt)
        if !is_stable && !courant_warning_issued
            log_print("⚠ WARNING: Courant Fo = $(round(Fo, digits=3)) > 0.25 (potentially unstable)")
            log_print("  Recommendation: reduce Δt to $(round(dt_max, sigdigits=3)) seconds")
            courant_warning_issued = true
        end
        
        #______________________________________________________
        # Reset flow vector
        #______________________________________________________
        q_flux_water .= 0.0
        
        #______________________________________________________
        # Calculate capillary and gravitational fluxes
        #______________________________________________________
        @threads for e in 1:Nelements
            # Get element nodes
            nodes = mesh.elements[e, :]
            
            # Get material properties
            material_idx = get_element_material(mesh, e)
            if material_idx === nothing
                error("Element $e has no material assigned.")
            end
            
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            
            # Get element head values
            h_e = [h[nodes[i]] for i in 1:4]
            
            #__________________________________________
            # Capillary flux: q_cap = -K(h) ∇h
            #__________________________________________
            if calculate_capillary
                q_aux = zeros(4)
                
                # Evaluate K(h) at nodes and average
                K_h_e = [soil.K_h(h_e[i]) for i in 1:4]
                K_h_mean = mean(K_h_e)
                
                # q_cap = -K_h_mean * K_geometric[e] * h_e
                q_flux = -K_h_mean * K_geometric[e] * h_e
                
                for i in 1:4
                    node_id = nodes[i]
                    q_flux_water[node_id] += q_flux[i]
                end
            end
            
            #__________________________________________
            # Gravitational flux: q_grav = K(h) at z-direction
            #__________________________________________
            if calculate_gravity
                # Loop over Gauss points for gravitational term
                for p in 1:4
                    # Get shape functions and derivatives
                    N_p = ShapeFunctions.shape_funcs.N[p]
                    B = ShapeFunctions.get_B(p)
                    invJ = ShapeFunctions.get_invJ(e, p)
                    detJ = ShapeFunctions.get_detJ(e, p)
                    w = ShapeFunctions.shape_funcs.gauss_weights[p]
                    
                    # Evaluate K(h) at Gauss point
                    h_gp = N_p' * h_e
                    K_h_gp = soil.K_h(h_gp)
                    
                    # Gravitational flux: K(h) * ∂z component
                    dN_dx = B * invJ  # [4 nodes, 2 coords]
                    q_grav_contribution = K_h_gp * dN_dx[:, 2] * w * detJ  # z-direction (index 2)
                    
                    for i in 1:4
                        node_id = nodes[i]
                        q_flux_water[node_id] += g_magnitude * q_grav_contribution[i]
                    end
                end
            end
        end
        
        #______________________________________________________
        # Calculate rate of change: dh/dt = (q_net) / c_s(h)
        #______________________________________________________
        @threads for i in 1:Nnodes
            if M[i] > 1.0e-15
                dh_dt[i] = q_flux_water[i] / M[i]
            else
                dh_dt[i] = 0.0
            end
        end
        
        #______________________________________________________
        # Forward Euler: h^(n+1) = h^n + Δt × dh/dt
        #______________________________________________________
        @threads for i in 1:Nnodes
            h[i] += dt * dh_dt[i]
            
            # Enforce minimum head (optional: avoid extreme unsaturation)
            h_min = -10.0  # Minimum head [m]
            if h[i] < h_min
                if !negative_head_warned
                    log_print("Warning: Head below $h_min m at step $step. Clipping.")
                    negative_head_warned = true
                end
                h[i] = h_min
            end
        end
        
        #______________________________________________________
        # Update secondary water variables from SWRC
        #______________________________________________________
        @threads for i in 1:Nnodes
            # Get material at node
            material_idx = get_node_material(mesh, i)
            soil = materials.soils[materials.soil_dictionary[material_idx]]
            
            # Update water content and saturation
            theta_w[i] = soil.theta_h(h[i])
            S_r[i] = theta_w[i] / soil.porosity
            
            # Calculate water pressure: u_w = γ_w × h
            P_water[i] = γ_water * h[i]  # [Pa]
        end
        
        #______________________________________________________
        # Calculate water velocities using Darcy's law
        #______________________________________________________
        v_water .= 0.0
        
        for e in 1:Nelements
            nodes = mesh.elements[e, :]
            material_idx = get_element_material(mesh, e)
            if material_idx === nothing
                continue
            end
            
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            
            h_e = [h[nodes[i]] for i in 1:4]
            
            # Loop over Gauss points
            for p in 1:4
                N_p = ShapeFunctions.shape_funcs.N[p]
                B = ShapeFunctions.get_B(p)
                invJ = ShapeFunctions.get_invJ(e, p)
                detJ = ShapeFunctions.get_detJ(e, p)
                w = ShapeFunctions.shape_funcs.gauss_weights[p]
                
                # Evaluate h at Gauss point
                h_gp = N_p' * h_e
                K_h_gp = soil.K_h(h_gp)
                
                # Pressure gradient: ∇P = γ_w × ∇h
                dN_dx = B * invJ
                grad_h = dN_dx' * h_e
                grad_P = γ_water * grad_h
                
                # Darcy velocity: v = -K/μ × ∇P, but for water typically written as v = -K × ∇h (no viscosity)
                v_gp = -K_h_gp * grad_h
                
                # Mass weight
                mass_weight = w * detJ
                
                # Distribute to nodes
                for i in 1:4
                    node_id = nodes[i]
                    v_water[node_id, 1] += v_gp[1] * N_p[i] * mass_weight  # x-component
                    v_water[node_id, 2] += v_gp[2] * N_p[i] * mass_weight  # z-component
                end
            end
        end
        
        # Normalize velocities by mass
        for i in 1:Nnodes
            if M[i] > 1.0e-15
                v_water[i, :] ./= M[i]
            end
        end
        
        #______________________________________________________
        # Update time
        #______________________________________________________
        current_time += dt
        
        #______________________________________________________
        # Check if save output
        #______________________________________________________
        if save_data || step == num_steps
            progress = 100.0 * step / num_steps
            log_print(@sprintf("      Load Step %d (%.1f%%), Time = %.4e %s",
                              output_counter, progress, current_time, 
                              calc_params["units"]["time_unit"]))
            
            # write_output_vtk_water(mesh, materials, output_counter, current_time, project_name)
            
            next_output_time += load_step_time
            output_counter += 1
            save_data = false
            negative_head_warned = false
        end
        
        #______________________________________________________
        # Adjust dt to hit next output time exactly
        #______________________________________________________
        if current_time + dt > next_output_time
            dt = next_output_time - current_time
            save_data = true
        else
            dt = time_data.actual_dt
            save_data = false
        end
    end
    
    log_print("   ✓ Time integration completed")
    log_print(@sprintf("   ✓ Final time: %.4e %s", current_time, calc_params["units"]["time_unit"]))
    
    # Return checkpoint data
    return (current_time=current_time, output_counter=output_counter, 
            next_output_time=next_output_time, h=h)
end
