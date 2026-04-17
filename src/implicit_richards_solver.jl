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


# ══════════════════════════════════════════════════════════════════════════════
# RichardsCache: precomputed physical derivatives from ShapeFunctions
# ══════════════════════════════════════════════════════════════════════════════

"""
    RichardsCache

Precomputed element data for the implicit Richards solver.
Built once from ADSIM's ShapeFunctions module before the time loop.

Convention: Bp[e, p, :, :] is [2×4] where row 1 = ∂N/∂x, row 2 = ∂N/∂y.
This is the transpose of ADSIM's dN_dx [4×2].
"""
struct RichardsCache
    Bp      :: Array{Float64, 4}       # [ne, 4, 2, 4]
    detJ    :: Matrix{Float64}         # [ne, 4]
    Np      :: Vector{Vector{Float64}} # [4][4]
    A_e     :: Vector{Float64}         # [ne]
    weights :: Vector{Float64}         # [4]
end

"""
    build_richards_cache(mesh) → RichardsCache

Build precomputed cache from ADSIM's ShapeFunctions module.
Must be called AFTER `initialize_shape_functions!(mesh)` in kernel.jl.
"""
function build_richards_cache(mesh) :: RichardsCache
    ne = mesh.num_elements
    weights = ShapeFunctions.get_gauss_weights()
    Np = [copy(ShapeFunctions.get_N(p)) for p in 1:4]

    Bp_all   = zeros(Float64, ne, 4, 2, 4)
    detJ_all = zeros(Float64, ne, 4)
    A_e_all  = zeros(Float64, ne)

    for e in 1:ne
        area = 0.0
        for p in 1:4
            B_iso = ShapeFunctions.get_B(p)
            invJ  = ShapeFunctions.get_invJ(e, p)
            dJ    = ShapeFunctions.get_detJ(e, p)

            dN_dx = B_iso * invJ   # [4×2]
            for a in 1:4
                Bp_all[e, p, 1, a] = dN_dx[a, 1]   # ∂Na/∂x
                Bp_all[e, p, 2, a] = dN_dx[a, 2]   # ∂Na/∂y
            end

            detJ_all[e, p] = dJ
            area += dJ * weights[p]
        end
        A_e_all[e] = area
    end

    return RichardsCache(Bp_all, detJ_all, Np, A_e_all, weights)
end


# ══════════════════════════════════════════════════════════════════════════════
# Element matrices: lumped capacity + anisotropic stiffness
# ══════════════════════════════════════════════════════════════════════════════

"""
    element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_curr, hᵉ_prev, eprops, Δt, e_g, cache, e)

Compute element system matrix Aᵉ (4×4) and residual Rᵉ (4) in-place.
Anisotropic: uses K_h_x(model, h) and K_h_y(model, h) from swrc_models.jl.
Gravity enters RESIDUAL ONLY. Lumped mass for exact mass conservation.
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


# ══════════════════════════════════════════════════════════════════════════════
# Global assembly + sparsity pattern
# ══════════════════════════════════════════════════════════════════════════════

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
"""
function assemble_richards!(
    A      :: SparseMatrixCSC{Float64, Int},
    R      :: Vector{Float64},
    h_curr :: Vector{Float64},
    h_prev :: Vector{Float64},
    mesh,
    elem_props :: Vector{ElementWaterProps},
    Δt     :: Float64,
    e_g    :: Vector{Float64},
    cache  :: RichardsCache
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
    # This prevents flux contributions from affecting solution at prescribed nodes
    # (since δ_i = 0 when both A row and R[i] are zeroed, and h_i stays at prescribed value)
    # ─────────────────────────────────────────────────────────────────
    global P_boundary_water
    for i in 1:mesh.num_nodes
        if P_boundary_water[i, 1] == 0  # BC node
            # Zero the matrix row i (replace with identity: A[i,i] = 1, A[i,j≠i] = 0)
            # This must be done before zeroing the residual to maintain linear system structure
            rows, cols, vals = findnz(A)
            for k in eachindex(rows)
                if rows[k] == i
                    A.nzval[k] = rows[k] == cols[k] ? 1.0 : 0.0
                end
            end
            # Zero the residual at BC nodes
            R[i] = 0.0
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Phase 2: Apply Neumann boundary flux contributions
    # Add prescribed flux to residual at non-Dirichlet nodes
    # ─────────────────────────────────────────────────────────────────
    global q_boundary_water
    for i in 1:mesh.num_nodes
        # Only add flux at nodes WITHOUT Dirichlet BCs
        if P_boundary_water[i, 1] != 0  # Free node (not Dirichlet)
            if q_boundary_water[i] != 0.0  # Has prescribed flux
                R[i] += q_boundary_water[i]
                # Sign convention:
                #   q > 0 = inflow (increases h)
                #   q < 0 = outflow (decreases h)
            end
        end
    end

    return nothing
end


#[LZC: Assemble flow BC a priori. See how flow is assemble for carbonation equation and then use directly.]
"""
    apply_neumann_edge_richards!(R, node_i, node_j, q_bar, coords) 

Add Neumann flux to residual for edge (node_i, node_j).
q_bar > 0 = flux INTO domain. 2-point Gauss on edge.
"""
function apply_neumann_edge_richards!(
    R      :: Vector{Float64},
    node_i :: Int,
    node_j :: Int,
    q_bar  :: Float64,
    coords :: Matrix{Float64}
)
    xi, yi = coords[node_i, 1], coords[node_i, 2]
    xj, yj = coords[node_j, 1], coords[node_j, 2]
    l_e = sqrt((xj - xi)^2 + (yj - yi)^2)

    inv_sqrt3 = 1.0 / sqrt(3.0)
    for s in (-inv_sqrt3, inv_sqrt3)
        R[node_i] += 0.5 * (1.0 - s) * q_bar * (l_e / 2.0)
        R[node_j] += 0.5 * (1.0 + s) * q_bar * (l_e / 2.0)
    end
    return nothing
end


# ══════════════════════════════════════════════════════════════════════════════
# Picard iteration (one time step)
# ══════════════════════════════════════════════════════════════════════════════

"""
    picard_richards!(h_curr, h_prev, mesh, elem_props, Δt, e_g,
                      A, cache; tol, max_iter, ω)

Picard iteration for one time step. Returns number of iterations.
Dirichlet BCs are enforced via P_boundary_water masking in assemble_richards!().
"""
function picard_richards!(
    h_curr     :: Vector{Float64},
    h_prev     :: Vector{Float64},
    mesh,
    elem_props :: Vector{ElementWaterProps},
    Δt         :: Float64,
    e_g        :: Vector{Float64},
    A          :: SparseMatrixCSC{Float64, Int},
    cache      :: RichardsCache;
    tol        :: Float64 = 1e-8,
    max_iter   :: Int     = 100,
    ω          :: Float64 = 1.0
)
    N = mesh.num_nodes
    R = zeros(Float64, N)

    h_curr .= h_prev

    for m in 1:max_iter
        assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache)
        # Masking is now handled in assemble_richards!() via P_boundary_water

        δ = A \ R
        δ_norm = maximum(abs.(δ))

        if δ_norm < tol
            return m
        end

        h_curr .+= ω .* δ
    end

    @warn "Picard did not converge in $max_iter iterations"
    return max_iter
end


# ══════════════════════════════════════════════════════════════════════════════
# VTK output wrapper
# ══════════════════════════════════════════════════════════════════════════════

"""
    write_richards_output(mesh, output_counter, current_time, project_name)

Write VTK output using ADSIM global water state variables.
"""
function write_richards_output(mesh, output_counter::Int, current_time::Float64,
                                project_name::String)
    global h, theta_w, S_r, P_water, v_water

    output_dir = "output"
    filename = joinpath(output_dir, project_name * "_water")
    write_vtk_file_water(filename, output_counter, current_time,
                          mesh, h, theta_w, S_r, P_water, v_water)
end


# ══════════════════════════════════════════════════════════════════════════════
# MAIN SOLVER — interface matching ADSIM kernel
# ══════════════════════════════════════════════════════════════════════════════

"""
    implicit_richards_solver(mesh, materials, calc_params, time_data,
                              project_name, log_print, cache, elem_props, initial_state)

Implicit FEM solver for the Richards equation using Backward Euler + Picard iteration.
Interface matches ADSIM's fully_explicit_diffusion_solver for drop-in replacement.

# Arguments
- `cache::RichardsCache`: Precomputed shape function and Jacobian data from kernel
- `elem_props::Vector{ElementWaterProps}`: Precomputed element SWRC models from kernel

# Returns
- NamedTuple: (current_time, output_counter, next_output_time)

# Note
These precomputed data structures are built in kernel.jl to avoid redundant
computation if the solver is called multiple times or chained with other solvers.
"""
function implicit_richards_solver(mesh, materials, calc_params, time_data,
                                   project_name, log_print, cache, elem_props, initial_state=nothing)
    
    global h, theta_w, S_r, P_water, v_water

    log_print("   Starting implicit Richards solver (Picard iteration)")

    # ── Extract parameters ────────────────────────────────────────────
    dt = time_data.actual_dt
    n_steps = time_data.num_steps
    load_step_time = calc_params["data_saving_interval"]

    gx = calc_params["gravity"]["x_component"]
    gy = calc_params["gravity"]["y_component"]
    e_g = [gx, gy]

    rho_w = materials.liquid.density
    g_mag = calc_params["gravity"]["magnitude"]

    log_print(@sprintf("   Δt = %.4e %s, n_steps = %d",
                        dt, calc_params["units"]["time_unit"], n_steps))
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
    zero_flow_vectors_water!(mesh.num_nodes)
    boundary_influences = get_boundary_node_influences(mesh)
    apply_boundary_flows_water!(mesh, q_boundary_water, boundary_influences.node_influences)
    n_neumann_nodes = count(q_boundary_water .!= 0.0)
    log_print("   ✓ Neumann BC nodes initialized: $n_neumann_nodes nodes")

    log_print(@sprintf("   ✓ Initial h: min=%.4f, max=%.4f", minimum(h), maximum(h)))

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
        write_richards_output(mesh, 0, 0.0, project_name)
        current_time     = 0.0
        next_output_time = load_step_time
        output_counter   = 1
    end

    # ── Time loop ─────────────────────────────────────────────────────
    h_new = copy(h)
    save_data = false

    for step in 1:n_steps
        n_iter = picard_richards!(h_new, h, mesh, elem_props, dt, e_g,
                                   A, cache;
                                   tol=1e-8, max_iter=100, ω=1.0)

        h .= h_new
        current_time += dt

        # Re-enforce BCs on global water variables (ensures exact BC values after Picard)
        enforce_water_dirichlet_bc!(mesh, materials)

        if save_data || step == n_steps
            progress = 100.0 * step / n_steps
            log_print(@sprintf("      Load Step %d (%.1f%%), Time = %.4e %s, Picard = %d",
                                output_counter, progress, current_time,
                                calc_params["units"]["time_unit"], n_iter))

            update_water_globals!(elem_props, mesh, e_g, cache, rho_w, g_mag)
            write_richards_output(mesh, output_counter, current_time, project_name)

            next_output_time += load_step_time
            output_counter += 1
            save_data = false
        end

        if current_time + dt > next_output_time
            dt = next_output_time - current_time
            save_data = true
        else
            dt = time_data.actual_dt
            save_data = false
        end
    end

    log_print("   ✓ Richards solver completed")
    log_print(@sprintf("   ✓ Final time: %.4e %s", current_time, calc_params["units"]["time_unit"]))

    return (current_time=current_time, output_counter=output_counter, next_output_time=next_output_time)
end
