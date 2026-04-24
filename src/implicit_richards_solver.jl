#______________________________________________________
# ADSIM: Implicit Richards Equation Solver
# Q4 Isoparametric Elements · Backward Euler · Picard Iteration
# Following Celia et al. (1990), Water Resources Research 26(7):1483–1496
#
# Anisotropic K(h) via SWRC model dispatch:
#   K_h_x(model, h) and K_h_y(model, h) from swrc_models.jl
# Gravity enters residual only — never in the system matrix
# Lumped mass on temporal residual for exact discrete mass conservation
#
# Dependencies (from other ADSIM files):
#   swrc_models.jl          → SWRCModel, theta, C_moist, K_h_x, K_h_y, Se
#   read_materials.jl       → ElementWaterProps, precompute_element_water_props
#   initialize_variables.jl → build_dirichlet_lists, update_water_globals!,
#                              enforce_water_dirichlet_bc!, global h/theta_w/...
#   shape_functions.jl      → ShapeFunctions module (get_N, get_B, get_invJ, get_detJ)
#   write_vtk.jl            → write_vtk_file_water
#
# Authors: Paula Sarmiento, Luis Zambrano-Cruzatty
#______________________________________________________

using LinearAlgebra
using SparseArrays
using Printf

using .ShapeFunctions
using .WriteVTK


# Element matrices: lumped capacity + anisotropic stiffness

"""
    element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_curr, hᵉ_prev, eprops, Δt, e_g, cache, e)

Compute element system matrix Aᵉ (4×4) and residual Rᵉ (4) for one Q4 element.

# Physical Model
Mixed-form Richards equation: ∂θ/∂t + ∇·[-K(h)(∇h + e_g)] = 0

Weak form (Backward Euler, implicit):
- LHS: ∫(∇Nₐ · K∇N_b) dΩ + ∫(Nₐ · C(h)/Δt · N_b) dΩ
- RHS: -∫(∇Nₐ · K∇h) dΩ - ∫(∇Nₐ · K·e_g) dΩ - ∫(Nₐ · (θ_curr-θ_prev)/Δt) dΩ

# Integration
2×2 Gauss quadrature (4 points per element):
- At each point: interpolate h, compute K(h) and C(h)
- K is anisotropic: K = [Kₓ  0]  (Kₓ, Kᵧ from SWRC model)
                         [ 0 Kᵧ]
- Shape derivatives precomputed in cache

# Key Details
- Lumped mass for temporal term (ensures discrete mass conservation)
- Gravity vector e_g = [0, -1] enters residual only (not matrix)
- K and C recomputed each Picard iteration (nonlinear)

# Arguments
- Aᵉ[4,4]: Element matrix (output, zeroed at start)
- Rᵉ[4]: Element residual (output, zeroed at start)
- hᵉ_curr[4]: Current pressure head at 4 nodes
- hᵉ_prev[4]: Previous pressure head at 4 nodes
- eprops: Element SWRC properties
- Δt: Time step
- e_g: Gravity vector
- cache: Precomputed shape functions, Jacobians
- e: Element index

# Notes
- Modified in-place for performance
- Dispatcher to SWRC model for K_h_x, K_h_y, theta, C_moist
- References: Celia et al. (1990)
"""
function element_matrices_aniso!(
    Aᵉ      :: Matrix{Float64},
    Rᵉ      :: Vector{Float64},
    hᵉ_curr :: Vector{Float64},
    hᵉ_prev :: Vector{Float64},
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
        hp = Np[1]*hᵉ_curr[1] + Np[2]*hᵉ_curr[2] + Np[3]*hᵉ_curr[3] + Np[4]*hᵉ_curr[4]

        # Anisotropic conductivity directly from SWRC model dispatch
        Kx = K_h_x(model, hp)
        Ky = K_h_y(model, hp)

        # ── LHS: anisotropic stiffness ────────────────────────────────
        @inbounds for a in 1:4
            Bxa = cache.Bp[e, p, 1, a]
            Bya = cache.Bp[e, p, 2, a]
            for b in 1:4
                Bxb = cache.Bp[e, p, 1, b]
                Byb = cache.Bp[e, p, 2, b]
                Aᵉ[a, b] += (Kx * Bxa * Bxb + Ky * Bya * Byb) * w
            end
        end

        # ── RHS: internal flux + gravity ──────────────────────────────
        grad_h_x = cache.Bp[e,p,1,1]*hᵉ_curr[1] + cache.Bp[e,p,1,2]*hᵉ_curr[2] +
                   cache.Bp[e,p,1,3]*hᵉ_curr[3] + cache.Bp[e,p,1,4]*hᵉ_curr[4]
        grad_h_y = cache.Bp[e,p,2,1]*hᵉ_curr[1] + cache.Bp[e,p,2,2]*hᵉ_curr[2] +
                   cache.Bp[e,p,2,3]*hᵉ_curr[3] + cache.Bp[e,p,2,4]*hᵉ_curr[4]

        @inbounds for a in 1:4
            Bxa = cache.Bp[e, p, 1, a]
            Bya = cache.Bp[e, p, 2, a]
            Rᵉ[a] += -(Bxa * Kx * grad_h_x + Bya * Ky * grad_h_y) * w
            Rᵉ[a] += +(Bxa * Kx * e_g[1]   + Bya * Ky * e_g[2])   * w
        end
    end

    # ── Lumped capacity ───────────────────────────────────────────────
    ML_ii = cache.A_e[e] / 4.0
    @inbounds for a in 1:4
        Ca = C_moist(model, hᵉ_curr[a])
        Aᵉ[a, a] += Ca / Δt * ML_ii
        θ_curr = theta(model, hᵉ_curr[a])
        θ_prev = theta(model, hᵉ_prev[a])
        Rᵉ[a] += -(1.0 / Δt) * ML_ii * (θ_curr - θ_prev)
    end

    return nothing
end


# Global assembly + sparsity pattern

"""
    build_richards_sparsity(mesh) → SparseMatrixCSC

Build sparsity pattern from mesh connectivity. Called once.
"""
function build_richards_sparsity(mesh)
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
    assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache)

Assemble global sparse system matrix A and residual vector R.
Applies Dirichlet and Neumann boundary conditions.

# Algorithm
1. Loop over all elements
2. Compute element matrix and residual via element_matrices_aniso!()
3. Scatter-add contributions to global A and R
4. Apply Dirichlet masking (zero rows at prescribed nodes)
5. Add Neumann flux contributions to residual

# Boundary Conditions
- Dirichlet (prescribed head): A[i,:] = 0, A[i,i] = 1, R[i] = 0
  (Prevents updates at prescribed nodes)
- Neumann (prescribed flux): Add to R at non-Dirichlet nodes

# Arguments
- A: Global sparse matrix (modified in-place)
- R: Global residual vector (modified in-place)
- h_curr: Current pressure head
- h_prev: Previous pressure head
- mesh: Mesh connectivity
- elem_props: Element SWRC properties
- Δt: Time step
- e_g: Gravity vector
- cache: Precomputed shape function data

# Notes
- Assumes sparsity pattern already built (build_richards_sparsity)
- P_boundary_water and q_boundary_water must be pre-populated
- BCs applied once per assembly (not in Picard loop)
"""
function assemble_richards!(
    A         :: SparseMatrixCSC{Float64, Int},
    R         :: Vector{Float64},
    h_curr    :: Vector{Float64},
    h_prev    :: Vector{Float64},
    mesh,
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
    hᵉ_curr = zeros(Float64, 4)
    hᵉ_prev = zeros(Float64, 4)

    for e in 1:mesh.num_elements
        @inbounds for a in 1:4
            nodes[a] = mesh.elements[e, a]
        end
        @inbounds for a in 1:4
            hᵉ_curr[a] = h_curr[nodes[a]]
            hᵉ_prev[a] = h_prev[nodes[a]]
        end

        element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_curr, hᵉ_prev,
                                 elem_props[e], Δt, e_g, cache, e)

        @inbounds for a in 1:4
            I = nodes[a]
            R[I] += Rᵉ[a]
            for b in 1:4
                A[I, nodes[b]] += Aᵉ[a, b]
            end
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Apply P_boundary_water masking: zero matrix rows and residual at BC nodes
    # Single CSC pass — correct even when structurally stored entries are zero.
    # ─────────────────────────────────────────────────────────────────
    global P_boundary_water
    n = mesh.num_nodes
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1] - 1)
            i = A.rowval[k]
            if P_boundary_water[i, 1] == 0   # BC node: enforce identity row
                A.nzval[k] = (i == j) ? 1.0 : 0.0
            end
        end
    end
    for i in 1:n
        if P_boundary_water[i, 1] == 0
            R[i] = 0.0
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Phase 2: Apply Neumann boundary flux contributions
    # Add prescribed flux to residual at non-Dirichlet nodes
    # ─────────────────────────────────────────────────────────────────
    global q_flux_water
    for i in 1:mesh.num_nodes
        if P_boundary_water[i, 1] != 0 && q_flux_water[i] != 0.0
            R[i] += q_flux_water[i]
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Phase 3: Free-drainage gravity flux at bottom edges
    # Arithmetic mean of K at the two edge nodes (Colab approach).
    # q_bot = K_avg * e_g[2]  (negative for downward gravity → outflow)
    # Contribution: R[ni] += q_bot * l_e/2,  R[nj] += q_bot * l_e/2
    # ─────────────────────────────────────────────────────────────────
    for (ni, nj, e_idx) in bot_edges
        model_e = elem_props[e_idx].model
        K_ni  = K_h_y(model_e, h_curr[ni])
        K_nj  = K_h_y(model_e, h_curr[nj])
        q_bot = (K_ni + K_nj) / 2.0 * e_g[2]
        xi = mesh.coordinates[ni, 1];  yi = mesh.coordinates[ni, 2]
        xj = mesh.coordinates[nj, 1];  yj = mesh.coordinates[nj, 2]
        l_e = sqrt((xj - xi)^2 + (yj - yi)^2)
        R[ni] += q_bot * l_e / 2.0
        R[nj] += q_bot * l_e / 2.0
    end

    return nothing
end




#[LZC-RESOLVED] Neumann BC assembly: Flow BCs ARE assembled a priori at solver startup
# via zero_flow_vectors_water!() + apply_boundary_flows_water!()
# Pre-computed q_boundary_water[i] is added to the residual in assemble_richards!()
# Alternative: apply_neumann_edge_richards!() (now in initialize_variables.jl) for edge-based
# assembly (not used — kept for reference)


# Picard iteration solver (one time step)

"""
    picard_richards!(h_curr, h_prev, mesh, elem_props, Δt, e_g, A, cache; kwargs...)

Picard iteration solver for one time step.

Solves the implicit Richards equation A(h) * delta_h = R(h) via fixed-point iteration.

# Algorithm
1. ASSEMBLE: Compute A and R at current h (nonlinear in h)
2. SOLVE: delta_h = A \\ R
3. CHECK: if ||delta_h||_infinity < tol, converged
4. UPDATE: h_new = h_curr + omega * delta_h
5. REPEAT until convergence or max_iter reached

# Arguments
- h_curr: Current pressure head (modified in-place)
- h_prev: Previous time step pressure head
- mesh: Mesh structure
- elem_props: Element SWRC properties
- Δt: Time step
- e_g: Gravity vector [gx, gy]
- A: Global sparse matrix (reused for efficiency)
- cache: Precomputed shape function data
- tol = 1e-8: Convergence tolerance
- max_iter = 100: Maximum iterations (safety limit)
- ω = 1.0: Relaxation parameter (Newton step size)

# Returns
Integer: Number of iterations performed

# Notes
- Convergence criterion: all nodes satisfy |delta_h[i]| <= tol
- K(h) and C(h) are recomputed each iteration (nonlinear)
- Backward Euler is unconditionally stable (solution accepted even if max_iter reached)
- References: Celia et al. (1990), standard Picard iteration theory
"""
function picard_richards!(
    h_curr     :: Vector{Float64},
    h_prev     :: Vector{Float64},
    mesh,
    elem_props :: Vector{ElementWaterProps},
    Δt         :: Float64,
    e_g        :: Vector{Float64},
    A          :: SparseMatrixCSC{Float64, Int},
    cache      :: RichardsCache,
    bot_edges  :: Vector{Tuple{Int,Int,Int}};
    tol        :: Float64 = 1e-8,
    max_iter   :: Int     = 100,
    ω          :: Float64 = 1.0
)
    N = mesh.num_nodes
    R = zeros(Float64, N)

    # ── Cold start: reset to previous time step (Colab pattern) ────────
    h_curr .= h_prev

    for m in 1:max_iter
        # Step 1: Assemble nonlinear system (A and R depend on current h)
        assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache, bot_edges)

        # Step 2: Solve for pressure head increments
        delta = A \ R
        delta_norm = maximum(abs.(delta))

        # Step 3: Update pressure head
        h_curr .+= ω .* delta

        # Step 4: Check convergence
        if delta_norm < tol
            return m
        end
    end

    # Non-convergence: Backward Euler is unconditionally stable so solution is accepted
    @warn "Picard did not converge in $max_iter iterations"
    return max_iter
end



# ══════════════════════════════════════════════════════════════════════════════
# MAIN SOLVER — interface matching ADSIM kernel
# ══════════════════════════════════════════════════════════════════════════════

"""
    implicit_richards_solver(mesh, materials, calc_params, time_data,
                              project_name, log_print, cache, elem_props, initial_state)

Implicit FEM solver for the Richards equation: ∂θ/∂t + ∇·[-K(∇h + e_g)] = 0

**Discretization:**
  • Time: Backward Euler (implicit, unconditionally stable)
  • Space: Q4 isoparametric elements, 2×2 Gauss quadrature
  • Nonlinearity: Picard iteration (5-20 iterations per step typical)
  • Mass: Lumped (ensures discrete mass conservation)

**Boundary Conditions:**
  • Dirichlet (prescribed h): Enforced via row masking in matrix A
  • Neumann (prescribed flux): Pre-computed and added to residual
  • Both applied ONCE per assembly, NOT in Picard loop (correct per Celia 1990)

**Architecture (ADSIM Pattern):**
  • Phase 1 (kernel.jl, before solver):
    - Step 1-3: Read mesh, materials, calc_params
    - Step 3.4: Precompute element SWRC properties → elem_props
    - Step 6: Initialize shape functions
    - Step 6.5: Build Richards cache (shape derivatives, Jacobians) → cache
  • Phase 2 (this function):
    - Initialize boundary conditions (once)
    - Time loop with Picard iteration

**Arguments:**
  - `mesh`: MeshData with nodes, elements, BCs
  - `materials`: MaterialData with soil SWRC models
  - `calc_params`: Dict with gravity, units, solver settings
  - `time_data`: TimeData with dt, num_steps, critical_dt
  - `project_name`: String for output file naming
  - `log_print`: Function for logging (e.g., println or file IO)
  - `cache::RichardsCache`: Precomputed shape function data (from kernel Step 6.5)
  - `elem_props::Vector{ElementWaterProps}`: Precomputed element SWRC models (from kernel Step 3.4)
  - `initial_state`: Optional state from previous stage (checkpoint restoration)

**Returns:**
  NamedTuple with:
    - `current_time`: Final simulation time
    - `output_counter`: Number of output steps written
    - `next_output_time`: Next scheduled output time (for multi-stage)

**Key Properties:**
  ✓ Gravity enters residual only (not matrix coefficients)
  ✓ Anisotropic K via SWRC dispatch: K_h_x(model, h), K_h_y(model, h)
  ✓ No solution reset between Picard iterations (correct formulation)
  ✓ Handles multi-stage calculations via checkpoints

**Implementation References:**
  • Celia et al. (1990) "An Efficient Iterative Scheme for Heterogeneous Porous Media"
  • ADSIM Architecture: kernel.jl (full workflow), swrc_models.jl (soil models)
  • Pattern follows: fully_explicit_solver.jl (multi-solver framework)

**See Also:**
  • IMPLICIT_RICHARDS_SOLVER_ARCHITECTURE.md (detailed documentation)
  • test_gravity_math.jl (verify gravity sign/magnitude)
  • test_part_1_1_check.jl (validate preprocessing)
"""
function implicit_richards_solver(mesh, materials, calc_params, time_data,
                                   project_name, log_print, cache, elem_props, initial_state=nothing)
    
    global h, theta_w, S_r, P_water, v_water

    log_print("   Starting implicit Richards solver (Picard iteration)")

    # ── Extract parameters ────────────────────────────────────────────
    dt = time_data.actual_dt
    load_step_time = calc_params["data_saving_interval"]

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

    # ── Build sparsity pattern ────────────────────────────────────────
    # Note: Cache and elem_props are precomputed in kernel.jl and passed as arguments
    A = build_richards_sparsity(mesh)
    log_print(@sprintf("   ✓ Sparse matrix: %d nonzeros (%.2f%% density)",
                        nnz(A), 100.0 * nnz(A) / mesh.num_nodes^2))

    # ── Mark Dirichlet BCs in P_boundary_water masking matrix ─────────
    apply_water_dirichlet_bc!(mesh, materials)
    n_bc_nodes = count(P_boundary_water .== 0)
    log_print("   ✓ Dirichlet BC nodes marked: $n_bc_nodes nodes")

    # ── Initialize water flow boundary conditions ─────────────────────
    apply_water_flux_bc!(mesh)
    n_neumann_nodes = count(q_flux_water .!= 0.0)
    log_print("   ✓ Neumann BC nodes initialized: $n_neumann_nodes nodes")

    # ── Precompute free-drainage bottom edges (Colab pattern) ─────────
    # Edges at y_min where both nodes are free (not Dirichlet).
    # K·e_g[2] gravity flux is applied here each Picard iteration.
    y_min = minimum(mesh.coordinates[:, 2])
    tol_y = 1e-10
    bot_edges = Tuple{Int,Int,Int}[]
    seen_edges = Set{Tuple{Int,Int}}()
    for e in 1:mesh.num_elements
        for (a, b) in ((1,2), (2,3), (3,4), (4,1))
            ni = mesh.elements[e, a]
            nj = mesh.elements[e, b]
            if mesh.coordinates[ni, 2] < y_min + tol_y &&
               mesh.coordinates[nj, 2] < y_min + tol_y &&
               P_boundary_water[ni, 1] != 0 &&
               P_boundary_water[nj, 1] != 0
                key = (min(ni,nj), max(ni,nj))
                if key ∉ seen_edges
                    push!(seen_edges, key)
                    push!(bot_edges, (ni, nj, e))
                end
            end
        end
    end
    log_print("   ✓ Free-drainage bottom edges: $(length(bot_edges))")


    # ── Time tracking ─────────────────────────────────────────────────
    if initial_state !== nothing
        current_time     = initial_state.current_time
        output_counter   = initial_state.output_counter
        next_output_time = initial_state.next_output_time
        log_print("      Continuing from time: $(current_time) $(calc_params["units"]["time_unit"])")
        log_print("      Next output at: $(next_output_time) $(calc_params["units"]["time_unit"])")
    else
        log_print("      Load step 0 (0.0%)")
        update_water_globals!(elem_props, mesh, e_g, cache, rho_w, g_mag)
        
        # ── Write initial condition (step 0) ──────────────────────────────
        output_dir = "output"
        filename = joinpath(output_dir, project_name * "_water")
        write_vtk_file_water(filename, 0, 0.0, mesh, h, theta_w, S_r, P_water, v_water)
        
        current_time     = 0.0
        next_output_time = Float64(load_step_time)
        output_counter   = 1
    end

    # ── Time loop ─────────────────────────────────────────────────────
    h_new = copy(h)
    total_time = calc_params["time_stepping"]["total_simulation_time"]
    step_count = 0

    while current_time < total_time - 1e-10
        step_count += 1

        # ── Clamp dt to hit the next output time or total_time ────────
        dt_step = min(dt, next_output_time - current_time, total_time - current_time)
        # Note: at_output checked AFTER the step based on time reached,
        # not whether dt was shortened (which misses exact-hit cases).

        # ── Picard iteration ──────────────────────────────────────────
        n_iter = picard_richards!(h_new, h, mesh, elem_props, dt_step, e_g,
                                   A, cache, bot_edges;
                                   tol=1e-4, max_iter=100, ω=1.0)

        # ── Accept step ───────────────────────────────────────────────
        h .= h_new
        current_time += dt_step
        enforce_water_dirichlet_bc!(mesh, materials)

        # ── Output at scheduled times and at final time ───────────────
        at_output = current_time >= next_output_time - 1e-10
        at_end    = current_time >= total_time - 1e-10
        if at_output || at_end
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

    return (current_time=current_time, output_counter=output_counter, next_output_time=next_output_time)
end
