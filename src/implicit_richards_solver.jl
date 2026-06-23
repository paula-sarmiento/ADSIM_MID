#______________________________________________________
# ADSIM: Implicit Richards Equation Solver
# Q4 Isoparametric · Backward Euler · Picard Iteration
# Following Celia et al. (1990)
# Authors: Paula Sarmiento, Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Implicit FEM solver for water flow in porous media
# Solves Richards equation via Picard iteration
#______________________________________________________

using LinearAlgebra
using SparseArrays
using Printf

using .ShapeFunctions
using .WriteVTK

#=
IMPLEMENTATION NOTES:
=====================

This module implements an implicit finite element solver for the mixed-form
Richards equation: ∂θ/∂t + ∇·[-K(h)(∇h + e_g)] = 0

Key features:
  • Backward Euler time discretization (unconditionally stable)
  • Picard iteration for nonlinearity (K and θ depend on h)
  • Lumped mass for discrete mass conservation
  • Anisotropic conductivity via SWRC model dispatch
  • Gravity enters residual only (not matrix)

Dependencies:
  swrc_models.jl          → K_h_x, K_h_y, theta, C_moist
  read_materials.jl       → ElementWaterProps
  initialize_variables.jl → BC enforcement, global variables
  shape_functions.jl      → Precomputed shape function data
  write_vtk.jl            → VTK output
=#

#______________________________________________________
# Picard Iteration Defaults (Richards Solver)
#______________________________________________________
const PICARD_TOL_DEFAULT = 1e-8           # Convergence tolerance for ||δh||
const PICARD_MAX_ITER_DEFAULT = 100       # Maximum Picard iterations
const PICARD_RELAXATION_DEFAULT = 1.0     # Relaxation parameter ω ∈ (0, 1]
const MIN_CONDUCTIVITY = 1e-12            # Minimum K(h) to prevent ill-conditioning [m/s]


# Element matrices: lumped capacity + anisotropic stiffness

"""
    element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_new, hᵉ_old, eprops, Δt, e_g, cache, e)

Compute element system matrix and residual for one Q4 element using Backward Euler.

# Arguments
- `Aᵉ::Matrix{Float64}`: Element matrix [4×4] (output, zeroed at start)
- `Rᵉ::Vector{Float64}`: Element residual [4] (output, zeroed at start)
- `hᵉ_new::Vector{Float64}`: Current pressure head at 4 nodes [timestep n+1]
- `hᵉ_old::Vector{Float64}`: Previous pressure head at 4 nodes [timestep n]
- `eprops::ElementWaterProps`: Element SWRC properties
- `Δt::Float64`: Time step
- `e_g::Vector{Float64}`: Gravity vector [gₓ, gᵧ]
- `cache::RichardsCache`: Precomputed shape functions, Jacobians
- `e::Int`: Element index

# Key Feature: SWRC Model Dispatch
Conductivity K(h) is NONLINEAR and depends on current pressure head hᵉ_new.
Uses model dispatch for K_h_x(model, hp) and K_h_y(model, hp) at each Gauss point.
→ Requires recalculation in EACH Picard iteration (called inside picard_richards! loop)

This is different from fully_explicit_solver.jl (which uses constant D_g).

# Notes
- Modified in-place for performance
- 2×2 Gauss quadrature (4 points per element)
- Anisotropic conductivity K = diag(Kₓ, Kᵧ) from SWRC model dispatch
- Capacity C_moist(h) is also nonlinear (recalculated at each Picard iteration)
- Lumped mass for discrete mass conservation
- Minimum conductivity clipping prevents ill-conditioning in dry regions
- References: Celia et al. (1990)
"""
function element_matrices_aniso!(
    Aᵉ      :: Matrix{Float64},
    Rᵉ      :: Vector{Float64},
    hᵉ_new  :: Vector{Float64},
    hᵉ_old  :: Vector{Float64},
    eprops   :: ElementWaterProps,
    Δt       :: Float64,
    e_g      :: Vector{Float64},
    cache    :: RichardsCache,
    e        :: Int
)
    fill!(Aᵉ, 0.0)
    fill!(Rᵉ, 0.0)

    model = eprops.model

    for p in 1:4
        Np = cache.Np[p]
        wp = cache.weights[p]
        dJ = cache.detJ[e, p]
        w  = wp * dJ

        # Interpolate h at Gauss point
        hp = Np[1]*hᵉ_new[1] + Np[2]*hᵉ_new[2] + Np[3]*hᵉ_new[3] + Np[4]*hᵉ_new[4]

        # Anisotropic conductivity directly from SWRC model dispatch
        Kx = K_h_x(model, hp)
        Ky = K_h_y(model, hp)
        
        # Clip conductivity to minimum to prevent ill-conditioning in dry regions
        Kx = max(Kx, MIN_CONDUCTIVITY)
        Ky = max(Ky, MIN_CONDUCTIVITY)

        # LHS: anisotropic stiffness
        @inbounds for a in 1:4
            Bxa = cache.Bp[e, p, 1, a]
            Bya = cache.Bp[e, p, 2, a]
            for b in 1:4
                Bxb = cache.Bp[e, p, 1, b]
                Byb = cache.Bp[e, p, 2, b]
                Aᵉ[a, b] += (Kx * Bxa * Bxb + Ky * Bya * Byb) * w
            end
        end

        # RHS: internal flux + gravity
        grad_h_x = cache.Bp[e,p,1,1]*hᵉ_new[1] + cache.Bp[e,p,1,2]*hᵉ_new[2] +
                   cache.Bp[e,p,1,3]*hᵉ_new[3] + cache.Bp[e,p,1,4]*hᵉ_new[4]
        grad_h_y = cache.Bp[e,p,2,1]*hᵉ_new[1] + cache.Bp[e,p,2,2]*hᵉ_new[2] +
                   cache.Bp[e,p,2,3]*hᵉ_new[3] + cache.Bp[e,p,2,4]*hᵉ_new[4]

        @inbounds for a in 1:4
            Bxa = cache.Bp[e, p, 1, a]
            Bya = cache.Bp[e, p, 2, a]
            Rᵉ[a] += -(Bxa * Kx * grad_h_x + Bya * Ky * grad_h_y) * w
            Rᵉ[a] += +(Bxa * Kx * e_g[1]   + Bya * Ky * e_g[2])   * w
        end
    end

    # Lumped capacity (ensures discrete mass conservation)
    ML_ii = cache.A_e[e] / 4.0
    @inbounds for a in 1:4
        Ca = C_moist(model, hᵉ_new[a])
        Aᵉ[a, a] += Ca / Δt * ML_ii
        θw_new = theta(model, hᵉ_new[a])
        θw_old = theta(model, hᵉ_old[a])
        Rᵉ[a] += -(1.0 / Δt) * ML_ii * (θw_new - θw_old)
    end

    return nothing
end


# Global assembly + sparsity pattern

"""
    build_richards_sparsity(mesh::MeshData)

Build sparsity pattern from mesh connectivity.

# Arguments
- `mesh::MeshData`: Mesh structure

# Returns
- `SparseMatrixCSC{Float64, Int}`: Sparse matrix with precomputed sparsity pattern
"""
function build_richards_sparsity(mesh::MeshData)::SparseMatrixCSC{Float64, Int}
    I_vec = Int[]
    J_vec = Int[]
    for e in 1:mesh.num_elements
        nodes = mesh.elements[e, :]
        for a in 1:4, b in 1:4
            push!(I_vec, nodes[a])
            push!(J_vec, nodes[b])
        end
    end
    return sparse(I_vec, J_vec, zeros(Float64, length(I_vec)),
                  mesh.num_nodes, mesh.num_nodes)
end

"""
    assemble_richards!(A, R, h_new, h_old, mesh, elem_props, Δt, e_g, cache)

Assemble global sparse system matrix A and residual vector R with boundary conditions.

# Arguments
- `A::SparseMatrixCSC`: Global sparse matrix (modified in-place)
- `R::Vector{Float64}`: Global residual vector (modified in-place)
- `h_new::Vector{Float64}`: Current pressure head [timestep n+1]
- `h_old::Vector{Float64}`: Previous pressure head [timestep n]
- `mesh::MeshData`: Mesh data structure
- `elem_props::Vector{ElementWaterProps}`: Element SWRC properties
- `Δt::Float64`: Time step
- `e_g::Vector{Float64}`: Gravity vector [gₓ, gᵧ]
- `cache::RichardsCache`: Precomputed shape function data
- `bot_edges::Vector{Tuple{Int,Int,Int}}`: Bottom boundary edges

# Notes
- Dirichlet BCs: A[i,:] = 0, A[i,i] = 1, R[i] = 0
- Neumann BCs: pre-computed and added to R
- Assumes sparsity pattern pre-built (build_richards_sparsity)
"""
function assemble_richards!(
    A         :: SparseMatrixCSC{Float64, Int},
    R         :: Vector{Float64},
    h_new     :: Vector{Float64},
    h_old     :: Vector{Float64},
    mesh      :: MeshData,
    elem_props :: Vector{ElementWaterProps},
    Δt        :: Float64,
    e_g       :: Vector{Float64},
    cache     :: RichardsCache,
    bot_edges :: Vector{Tuple{Int,Int,Int}}
)
    fill!(A.nzval, 0.0)
    fill!(R, 0.0)

    Aᵉ      = zeros(Float64, 4, 4)
    Rᵉ      = zeros(Float64, 4)
    nodes    = zeros(Int, 4)
    hᵉ_new = zeros(Float64, 4)
    hᵉ_old = zeros(Float64, 4)

    for e in 1:mesh.num_elements
        @inbounds for a in 1:4
            nodes[a] = mesh.elements[e, a]
        end
        @inbounds for a in 1:4
            hᵉ_new[a] = h_new[nodes[a]]
            hᵉ_old[a] = h_old[nodes[a]]
        end

        element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_new, hᵉ_old,
                                 elem_props[e], Δt, e_g, cache, e)

        @inbounds for a in 1:4
            I = nodes[a]
            R[I] += Rᵉ[a]
            for b in 1:4
                A[I, nodes[b]] += Aᵉ[a, b]
            end
        end
    end

    # Apply Dirichlet BC masking: zero matrix rows and residual at BC nodes
    global P_boundary_water
    n = mesh.num_nodes
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1] - 1)
            i = A.rowval[k]
            if P_boundary_water[i] == 0   # BC node: enforce identity row
                A.nzval[k] = (i == j) ? 1.0 : 0.0
            end
        end
    end
    for i in 1:n
        if P_boundary_water[i] == 0
            R[i] = 0.0
        end
    end

    # Apply Neumann boundary flux contributions
    global q_flux_water
    for i in 1:mesh.num_nodes
        if P_boundary_water[i] != 0 && q_flux_water[i] != 0.0
            R[i] += q_flux_water[i]
        end
    end

    # Free-drainage gravity flux at bottom edges
    for (ni, nj, e_idx) in bot_edges
        model_e = elem_props[e_idx].model
        K_ni  = K_h_y(model_e, h_new[ni])
        K_nj  = K_h_y(model_e, h_new[nj])
        q_bot = (K_ni + K_nj) / 2.0 * e_g[2]
        xi = mesh.coordinates[ni, 1];  yi = mesh.coordinates[ni, 2]
        xj = mesh.coordinates[nj, 1];  yj = mesh.coordinates[nj, 2]
        l_e = sqrt((xj - xi)^2 + (yj - yi)^2)
        R[ni] += q_bot * l_e / 2.0
        R[nj] += q_bot * l_e / 2.0
    end

    return nothing
end




# Neumann BC assembly: flow BCs are assembled a priori at solver startup
# via zero_flow_vectors_water!() + apply_boundary_flows_water!().
# Pre-computed q_boundary_water[i] is added to the residual in assemble_richards!().
# Edge-based alternative: apply_neumann_edge_richards!() in initialize_variables.jl is
# called by apply_water_flux_bc!() for liquid-discharge BCs (nodal flux BCs go through q_flux_water).


# Picard iteration solver (one time step)

"""
    picard_richards!(h_new, h_old, mesh, elem_props, Δt, e_g, A, cache; kwargs...)

Picard iteration solver for one time step.

# Arguments
- `h_new::Vector{Float64}`: Current pressure head [modified in-place, timestep n+1]
- `h_old::Vector{Float64}`: Previous time step pressure head [timestep n]
- `mesh::MeshData`: Mesh structure
- `elem_props::Vector{ElementWaterProps}`: Element SWRC properties
- `Δt::Float64`: Time step
- `e_g::Vector{Float64}`: Gravity vector [gₓ, gᵧ]
- `A::SparseMatrixCSC`: Global sparse matrix (reused for efficiency)
- `cache::RichardsCache`: Precomputed shape function data
- `bot_edges::Vector{Tuple{Int,Int,Int}}`: Bottom boundary edges
- `tol::Float64`: Convergence tolerance (default: 1e-8)
- `max_iter::Int`: Maximum iterations (default: 100)
- `ω::Float64`: Relaxation parameter (default: 1.0)

# Returns
- `Tuple{Int, Vector{Float64}}`: Number of iterations, convergence history

# Key Feature: Re-enforce Dirichlet BCs After Each Picard Update
After the Picard update hᵉ_new += ωδh, we EXPLICITLY re-enforce boundary conditions.
This is CRITICAL for convergence because:
  1. The Picard update δh can violate Dirichlet BCs even though row-zeroing is applied
  2. Pressure head at BC nodes must stay at prescribed values through iterations
  3. Without re-enforce, iterative error accumulates at boundaries

See line ~388 ("KEY (from Colab)" comment).

# Notes
- Solves A(h)δh = R(h) via fixed-point iteration  
- K(h) and C(h) recomputed each iteration (nonlinear)
- Backward Euler unconditionally stable (solution accepted even if max_iter reached)
- References: Celia et al. (1990)
"""
function picard_richards!(
    h_new      :: Vector{Float64},
    h_old      :: Vector{Float64},
    mesh       :: MeshData,
    elem_props :: Vector{ElementWaterProps},
    Δt         :: Float64,
    e_g        :: Vector{Float64},
    A          :: SparseMatrixCSC{Float64, Int},
    cache      :: RichardsCache,
    bot_edges  :: Vector{Tuple{Int,Int,Int}};
    tol        :: Float64 = PICARD_TOL_DEFAULT,
    max_iter   :: Int     = PICARD_MAX_ITER_DEFAULT,
    ω          :: Float64 = PICARD_RELAXATION_DEFAULT
) :: Tuple{Int, Vector{Float64}}
    N = mesh.num_nodes
    R = zeros(Float64, N)
    
    # Store initial Dirichlet BC values from mesh
    dbc_nodes = collect(keys(mesh.pressure_head_bc))
    dbc_vals  = [mesh.pressure_head_bc[i] for i in dbc_nodes]
    
    res_history = Float64[]

    # Cold start: reset to previous time step and apply Dirichlet BCs
    h_new .= h_old
    for (i, val) in zip(dbc_nodes, dbc_vals)
        h_new[i] = val
    end

    for m in 1:max_iter
        # Step 1: Assemble nonlinear system (A and R depend on current h)
        assemble_richards!(A, R, h_new, h_old, mesh, elem_props, Δt, e_g, cache, bot_edges)

        # Step 2: Solve for pressure head increments
        delta = A \ R
        delta_norm = maximum(abs.(delta))
        push!(res_history, delta_norm)
        
        # Log convergence every 10 iterations or at convergence/final
        if m % 10 == 1 || delta_norm < tol || m == max_iter
            log_print("      Picard iteration $m: ||δh|| = $(round(delta_norm; sigdigits=3)) m")
        end

        # Step 3: Check convergence
        if delta_norm < tol
            return m, res_history
        end

        # Step 4: Update pressure head
        h_new .+= ω .* delta
        
        # ✅ KEY (from Colab): Re-enforce Dirichlet BCs after update
        # This ensures boundary nodes stay at prescribed values through iterations
        for (i, val) in zip(dbc_nodes, dbc_vals)
            h_new[i] = val
        end
    end

    # Non-convergence: Backward Euler is unconditionally stable so solution is accepted
    @warn "Picard did not converge in $max_iter iterations. ||δ|| = $(res_history[end])"
    return max_iter, res_history
end



# ══════════════════════════════════════════════════════════════════════════════
# MAIN SOLVER — interface matching ADSIM kernel
# ══════════════════════════════════════════════════════════════════════════════

"""
    implicit_richards_solver(mesh, materials, calc_params, time_data,
                              project_name, log_print, cache, elem_props, initial_state)

Implicit FEM solver for Richards equation using Backward Euler.

# Arguments
- `mesh`: MeshData with nodes, elements, BCs
- `materials`: MaterialData with soil SWRC models
- `calc_params`: Dict with gravity, units, solver settings
- `time_data`: TimeData with dt, num_steps, critical_dt
- `project_name`: String for output file naming
- `log_print`: Function for logging (e.g., println or file IO)
- `cache::RichardsCache`: Precomputed shape function data
- `elem_props::Vector{ElementWaterProps}`: Precomputed element SWRC models
- `initial_state`: Optional state from previous stage (checkpoint)

# Returns
- `NamedTuple`: Final state with current_time, output_counter, next_output_time

# Notes
- Backward Euler: unconditionally stable time discretization
- Picard iteration: 5-20 iterations per step typical
- Lumped mass matrix for discrete mass conservation
- Gravity enters residual only (not matrix coefficients)
- References: Celia et al. (1990)
"""
function implicit_richards_solver(mesh::MeshData, materials::MaterialData, calc_params::Dict, time_data::TimeStepData,
                                   project_name::String, log_print::Function, cache::RichardsCache, elem_props::Vector{ElementWaterProps}, initial_state=nothing) :: NamedTuple
    
    global h, theta_w, S_r, P_water, v_water

    log_print("   Starting implicit Richards solver (Picard iteration)")

    # Extract parameters
    dt = time_data.actual_dt
    load_step_time = time_data.time_per_step

    gx = calc_params["gravity"]["x_component"]
    gy = calc_params["gravity"]["y_component"]
    e_g = [gx, gy]

    rho_w = materials.liquid.density
    g_mag = calc_params["gravity"]["magnitude"]

    log_print(@sprintf("   Δt = %.4e %s, total_time = %.4e %s",
                        dt, calc_params["units"]["time_unit"],
                        calc_params["time_stepping"]["total_simulation_time"],
                        calc_params["units"]["time_unit"]))
    log_print(@sprintf("   Gravity = [%.2f, %.2f], |g| = %.2f", gx, gy, g_mag))

    # Build sparsity pattern
    # Note: Cache and elem_props are precomputed in kernel.jl and passed as arguments
    A = build_richards_sparsity(mesh)
    log_print(@sprintf("   ✓ Sparse matrix: %d nonzeros (%.2f%% density)",
                        nnz(A), 100.0 * nnz(A) / mesh.num_nodes^2))

    # Mark Dirichlet BCs in P_boundary_water masking matrix
    apply_water_dirichlet_bc!(mesh, materials)
    n_bc_nodes = count(P_boundary_water .== 0)
    log_print("   ✓ Dirichlet BC nodes marked: $n_bc_nodes nodes")

    # Initialize water flow boundary conditions
    apply_water_flux_bc!(mesh)
    n_neumann_nodes = count(q_flux_water .!= 0.0)
    log_print("   ✓ Neumann BC nodes initialized: $n_neumann_nodes nodes")

    # Precompute free-drainage bottom edges (Colab pattern)
    # Edges at y_min where both nodes are free (not Dirichlet).
    # K·e_g[2] gravity flux is applied here each Picard iteration.
    y_min = minimum(mesh.coordinates[:, 2])
    y_max = maximum(mesh.coordinates[:, 2])
    y_extent = y_max - y_min
    tol_y = y_extent * 1e-10  # Relative tolerance — works at any mesh scale
    bot_edges = Tuple{Int,Int,Int}[]
    seen_edges = Set{Tuple{Int,Int}}()
    for e in 1:mesh.num_elements
        for (a, b) in ((1,2), (2,3), (3,4), (4,1))
            ni = mesh.elements[e, a]
            nj = mesh.elements[e, b]
            if mesh.coordinates[ni, 2] < y_min + tol_y &&
               mesh.coordinates[nj, 2] < y_min + tol_y &&
               P_boundary_water[ni] != 0 &&
               P_boundary_water[nj] != 0
                key = (min(ni,nj), max(ni,nj))
                if key ∉ seen_edges
                    push!(seen_edges, key)
                    push!(bot_edges, (ni, nj, e))
                end
            end
        end
    end
    log_print("   ✓ Free-drainage bottom edges: $(length(bot_edges))")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PHASE 1: Transport Solver Setup (Decoupled One-Way)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Transport solver runs AFTER Richards in each time step.
    # Phase 1: Decoupled — k_g and C_eq are user-specified constants.
    # Initialize theta_w_old for transport solver (frozen water content)
    theta_w_old = similar(theta_w)
    theta_w_old .= theta_w
    
    # Check if transport is enabled (future: read from config)
    enable_transport = false  # TODO: Read from calc_params["transport_settings"]
    
    if enable_transport
        global C_aq_gas  # Dissolved gas concentration [mol/m³]
        
        # Build transport sparsity pattern (identical to Richards mesh)
        A_transport = build_richards_sparsity(mesh)
        b_transport = zeros(Float64, mesh.num_nodes)
        
        # Initialize transport BCs (TODO: read from mesh/config)
        P_boundary_aq = similar(P_boundary_water)          # Dirichlet mask for transport
        mesh_partial_pressure_bc = Dict()                  # Dict of partial pressure BCs
        
        # TODO: apply_transport_dirichlet_bc!(mesh)
        # TODO: apply_transport_flux_bc!(mesh)
        n_transport_bc = count(P_boundary_aq .== 0)
        log_print("   ✓ Transport solver initialized (decoupled, Phase 1)")
        log_print("   ✓ Transport BC nodes marked: $n_transport_bc nodes")
        
        # Transport parameters (Phase 1: constant, user-specified)
        # TODO: Read D_h from config file
        D_h_transport = 1.0e-9                                        # m²/s, placeholder effective diffusivity
        kg_transport = ones(Float64, mesh.num_nodes) .* 1.0e-6        # 1/s, placeholder mass transfer coefficient [1/s]
        
        transport_params = AqueousTransportParams(D_h=D_h_transport)
    else
        A_transport = nothing
        b_transport = nothing
        transport_params = nothing
        log_print("   ⚠ Transport solver DISABLED (set enable_transport=true to activate)")
    end

    # Time tracking
    if initial_state !== nothing
        current_time     = initial_state.current_time
        output_counter   = initial_state.output_counter
        next_output_time = initial_state.next_output_time
        log_print("      Continuing from time: $(current_time) $(calc_params["units"]["time_unit"])")
        log_print("      Next output at: $(next_output_time) $(calc_params["units"]["time_unit"])")
    else
        log_print("      Load step 0 (0.0%)")
        update_water_globals!(elem_props, mesh, e_g, cache, rho_w, g_mag)
        
        # Write initial condition (step 0)
        output_dir = "output"
        filename = joinpath(output_dir, project_name * "_water")
        write_vtk_file_water(filename, 0, 0.0, mesh, h, theta_w, S_r, P_water, v_water)
        
        current_time     = 0.0
        next_output_time = Float64(load_step_time)
        output_counter   = 1
    end

    # Time loop
    h_new = copy(h)
    total_time = calc_params["time_stepping"]["total_simulation_time"]
    nonlinear_solver = get(calc_params, "nonlinear_solver", Dict{String, Any}())
    picard_tol = Float64(get(nonlinear_solver, "picard_tolerance", 1.0e-8))
    picard_max_iter = Int(get(nonlinear_solver, "picard_max_iter", 100))
    picard_relaxation = Float64(get(nonlinear_solver, "picard_relaxation", 1.0))

    log_print(@sprintf("   Picard settings: tol = %.3e, max_iter = %d, ω = %.2f",
                        picard_tol, picard_max_iter, picard_relaxation))
    step_count = 0
    picard_history = Int[]  # n_iter per output step

    while current_time < total_time - 1e-10
        step_count += 1

        # Clamp dt to hit the next output time or total_time
        dt_step = min(dt, next_output_time - current_time, total_time - current_time)
        # Note: at_output checked AFTER the step based on time reached,
        # not whether dt was shortened (which misses exact-hit cases).

        # Picard iteration
        n_iter, res_history = picard_richards!(h_new, h, mesh, elem_props, dt_step, e_g,
                                   A, cache, bot_edges;
                                   tol=picard_tol, max_iter=picard_max_iter, ω=picard_relaxation)

        # Accept step
        h .= h_new
        current_time += dt_step
        enforce_water_dirichlet_bc!(mesh, materials)        
        # Save previous water content for transport solver
        theta_w_old .= theta_w
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # PHASE 1: Transport solver step (after water flow)
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Reads globals set by update_water_globals!:
        #   - theta_w (just updated)
        #   - v_water (just updated)
        # Maintains transport globals via update_transport_globals!:
        #   - C_aq_gas (updates from transport_step!)
        
        if enable_transport
            # Aqueous concentration transport (Phase 1: decoupled, no feedback)
            # Loop over gas species and solve independent concentration transport
            for gas_idx in 1:NGases
                C_aq_gas[:, gas_idx] = aqueous_concentration_solver(
                    A_transport, b_transport, mesh, cache, materials,
                    P_boundary_aq, mesh_partial_pressure_bc,
                    transport_params,
                    theta_w, theta_w_old, v_water,
                    C_aq_gas[:, gas_idx], kg_transport, dt_step, gas_idx)
            end
        end

        # Output at scheduled times and at final time
        at_output = current_time >= next_output_time - 1e-10
        at_end    = current_time >= total_time - 1e-10
        if at_output || at_end
            push!(picard_history, n_iter)
            progress = 100.0 * current_time / total_time
            log_print(@sprintf("      Load Step %d (%.1f%%), Time = %.4e %s, Picard = %d",
                                output_counter, progress, current_time,
                                calc_params["units"]["time_unit"], n_iter))

            update_water_globals!(elem_props, mesh, e_g, cache, rho_w, g_mag)

            output_dir = "output"
            filename = joinpath(output_dir, project_name * "_water")
            write_vtk_file_water(filename, output_counter, current_time,
                                  mesh, h, theta_w, S_r, P_water, v_water)

            next_output_time += load_step_time
            output_counter += 1
        end
    end

    log_print("   ✓ Richards solver completed")
    log_print(@sprintf("   ✓ Final time: %.4e %s", current_time, calc_params["units"]["time_unit"]))

    return (current_time=current_time, output_counter=output_counter, next_output_time=next_output_time, picard_history=picard_history)
end
