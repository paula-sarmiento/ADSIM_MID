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

const CAPACITY_FLOOR_DEFAULT = 1.0e-4   # [1/m] lower bound on C(h) for matrix conditioning


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
- c_floor: Lower bound on capacity used only for the matrix diagonal
- floor_bound: Element-local flags indicating where the floor was active

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
    e        :: Int,
    c_floor  :: Float64,
    floor_bound :: Vector{Bool}
)
    fill!(Aᵉ, 0.0)
    fill!(Rᵉ, 0.0)
    fill!(floor_bound, false)

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
        C_exact = C_moist(model, hᵉ_curr[a])
        floor_bound[a] = C_exact < c_floor
        Ca = max(C_exact, c_floor)
        Aᵉ[a, a] += Ca / Δt * ML_ii
        θ_curr = theta(model, hᵉ_curr[a])
        θ_prev = theta(model, hᵉ_prev[a])
        Rᵉ[a] += -(1.0 / Δt) * ML_ii * (θ_curr - θ_prev)
    end

    return nothing
end

"""
    compute_total_water(h, mesh, elem_props, cache) -> Float64

Total stored water W = Σ_i M^L_ii θ(h_i), where M^L_ii is the lumped nodal mass
(A_e/4 summed over the elements touching node i). Uses the same lumping as
`element_matrices_aniso!`, so the discrete balance is exact.

# Arguments
- `h::Vector{Float64}`: Nodal pressure head
- `mesh`: Mesh connectivity
- `elem_props::Vector{ElementWaterProps}`: Element SWRC properties
- `cache::RichardsCache`: Precomputed element areas

# Notes
- At a node shared by different materials, the model of the first touching
  element is used.
"""
function compute_total_water(h::Vector{Float64}, mesh,
                             elem_props::Vector{ElementWaterProps},
                             cache::RichardsCache)::Float64
    M_lump = zeros(Float64, mesh.num_nodes)
    first_element = zeros(Int, mesh.num_nodes)

    for e in 1:mesh.num_elements
        ML_ii = cache.A_e[e] / 4.0
        for a in 1:4
            node = mesh.elements[e, a]
            M_lump[node] += ML_ii
            first_element[node] == 0 && (first_element[node] = e)
        end
    end

    W = 0.0
    for i in 1:mesh.num_nodes
        first_element[i] == 0 && error("Node $i is not connected to any element")
        model = elem_props[first_element[i]].model
        W += M_lump[i] * theta(model, h[i])
    end
    return W
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

function build_boundary_edge_ownership(mesh)
    edge_ownership = Dict{Tuple{Int,Int},Tuple{Int,Int,Int,Int}}()
    for element in 1:mesh.num_elements
        for (a, b) in ((1,2), (2,3), (3,4), (4,1))
            ni = mesh.elements[element, a]
            nj = mesh.elements[element, b]
            key = (min(ni, nj), max(ni, nj))
            if haskey(edge_ownership, key)
                count, owner_i, owner_j, owner_element = edge_ownership[key]
                edge_ownership[key] = (count + 1, owner_i, owner_j, owner_element)
            else
                edge_ownership[key] = (1, ni, nj, element)
            end
        end
    end
    return edge_ownership
end

"""
    build_water_neumann_edges(mesh) -> Vector{Tuple{Int,Int}}

Boundary edges carrying a water flux (Neumann) BC.

An edge is on the geometric boundary iff exactly one element contains it. An edge carries the
flux BC iff BOTH of its nodes appear in `mesh.liquid_discharge_bc`.

# Notes
- Independent of `absolute_pressure_bc` (a gas BC) — unlike `identify_boundary_edges`.
- Returns node pairs only; edge length is computed at assembly time.
"""
function build_water_neumann_edges(mesh)::Vector{Tuple{Int,Int}}
    edges = Tuple{Int,Int}[]
    for (key, ownership) in build_boundary_edge_ownership(mesh)
        ownership[1] == 1 || continue
        ni, nj = key
        (haskey(mesh.liquid_discharge_bc, ni) &&
         haskey(mesh.liquid_discharge_bc, nj)) || continue
        push!(edges, (ni, nj))
    end
    return edges
end

"""
    assemble_water_neumann!(q_boundary_water, mesh, neumann_edges)

Integrate the prescribed flux intensity over each Neumann edge into nodal residual
contributions, using linear interpolation of q̄ along the edge:

    R_i += l_e (2 q_i + q_j) / 6
    R_j += l_e (q_i + 2 q_j) / 6

from ∫N_i² dΓ = l_e/3 and ∫N_i N_j dΓ = l_e/6. For a uniform q̄ this reduces to q̄ l_e / 2 at
each node — the standard consistent load vector.

Sign convention: q̄ > 0 is inflow and adds positively to the residual.
"""
function assemble_water_neumann!(q_boundary_water::Vector{Float64}, mesh,
                                 neumann_edges::Vector{Tuple{Int,Int}})
    fill!(q_boundary_water, 0.0)
    for (ni, nj) in neumann_edges
        qi = mesh.liquid_discharge_bc[ni]
        qj = mesh.liquid_discharge_bc[nj]
        xi = mesh.coordinates[ni, 1];  yi = mesh.coordinates[ni, 2]
        xj = mesh.coordinates[nj, 1];  yj = mesh.coordinates[nj, 2]
        l_e = sqrt((xj - xi)^2 + (yj - yi)^2)
        q_boundary_water[ni] += l_e * (2.0*qi + qj) / 6.0
        q_boundary_water[nj] += l_e * (qi + 2.0*qj) / 6.0
    end
    return nothing
end

"""
    build_free_drainage_edges(mesh) -> Vector{Tuple{Int,Int,Int}}

Find topological boundary edges whose two endpoints belong to `mesh.free_drainage_bc`.
The returned endpoint order follows the owning counter-clockwise element so
`calculate_edge_outward_normal` produces the outward normal for any edge orientation.
"""
function build_free_drainage_edges(mesh)::Vector{Tuple{Int,Int,Int}}
    edges = Tuple{Int,Int,Int}[]
    for (_, ownership) in build_boundary_edge_ownership(mesh)
        count, ni, nj, element = ownership
        count == 1 || continue
        if ni in mesh.free_drainage_bc && nj in mesh.free_drainage_bc
            push!(edges, (ni, nj, element))
        end
    end
    return edges
end

function free_drainage_edge_flux(h_curr, mesh, elem_props, e_g,
                                 edge::Tuple{Int,Int,Int})::Float64
    ni, nj, element = edge
    model = elem_props[element].model
    Kx = 0.5 * (K_h_x(model, h_curr[ni]) + K_h_x(model, h_curr[nj]))
    Ky = 0.5 * (K_h_y(model, h_curr[ni]) + K_h_y(model, h_curr[nj]))
    edge_length, outward_normal = calculate_edge_outward_normal(mesh, ni, nj)
    inward_flux = -(Kx * e_g[1] * outward_normal[1] +
                    Ky * e_g[2] * outward_normal[2])
    return inward_flux * edge_length
end

function compute_free_drainage_flux(h_curr, mesh, elem_props,
                                     e_g, free_drainage_edges)::Float64
    total_flux = 0.0
    for edge in free_drainage_edges
        total_flux += free_drainage_edge_flux(h_curr, mesh, elem_props, e_g, edge)
    end
    return total_flux
end

function compute_net_boundary_flux(reaction, q_boundary_water, free_drainage_flux)::Float64
    return -sum(reaction) + sum(q_boundary_water) + free_drainage_flux
end

"""
    assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache,
                       q_boundary_water, free_drainage_edges; reaction_out=nothing)

Assemble global sparse system matrix A and residual vector R.
Applies Dirichlet and Neumann boundary conditions.

# Algorithm
1. Assemble element contributions into A and R
2. Add prescribed Neumann and explicit free-drainage loads at all nodes
3. Capture the total residual at Dirichlet nodes
4. Apply Dirichlet masking (zero rows at prescribed nodes)

# Boundary Conditions
- Dirichlet (prescribed head): A[i,:] = 0, A[i,i] = 1, R[i] = 0
  (Prevents updates at prescribed nodes)
- Neumann (prescribed flux): Add to R before Dirichlet masking

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
- q_boundary_water: Integrated water Neumann residual contributions
- free_drainage_edges: Explicitly marked free-drainage boundary edges
- reaction_out: Optional vector receiving the total residual at Dirichlet nodes
- c_floor: Lower bound on capacity used only for matrix conditioning
- floor_bound_nodes: Optional run-level mask of nodes where the floor was active

# Notes
- Assumes sparsity pattern already built (build_richards_sparsity)
- P_boundary_water and q_boundary_water must be pre-populated
- BCs applied once per assembly (not in Picard loop)
- `reaction_out[i]` is the conventional support reaction at Dirichlet node `i`,
    including boundary loads applied there; physical inward flow is `-reaction_out[i]`.
- Global balance uses `-sum(reaction_out) + sum(q_boundary_water) + free_drainage_flux`.
    Full boundary loads are added because their Dirichlet-node portions are already
    offset inside the conventional reactions.
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
    q_boundary_water :: Vector{Float64},
    free_drainage_edges :: Vector{Tuple{Int,Int,Int}};
    reaction_out :: Union{Nothing,Vector{Float64}} = nothing,
    c_floor :: Float64 = CAPACITY_FLOOR_DEFAULT,
    floor_bound_nodes :: Union{Nothing,BitVector} = nothing
)
    global P_boundary_water

    fill!(A.nzval, 0.0)
    fill!(R, 0.0)

    Aᵉ      = zeros(Float64, 4, 4)
    Rᵉ      = zeros(Float64, 4)
    nodes    = zeros(Int, 4)
    hᵉ_curr = zeros(Float64, 4)
    hᵉ_prev = zeros(Float64, 4)
    floor_bound = fill(false, 4)

    for e in 1:mesh.num_elements
        @inbounds for a in 1:4
            nodes[a] = mesh.elements[e, a]
        end
        @inbounds for a in 1:4
            hᵉ_curr[a] = h_curr[nodes[a]]
            hᵉ_prev[a] = h_prev[nodes[a]]
        end

        element_matrices_aniso!(Aᵉ, Rᵉ, hᵉ_curr, hᵉ_prev,
                                 elem_props[e], Δt, e_g, cache, e,
                                 c_floor, floor_bound)

        @inbounds for a in 1:4
            I = nodes[a]
            if floor_bound_nodes !== nothing && floor_bound[a]
                floor_bound_nodes[I] = true
            end
            R[I] += Rᵉ[a]
            for b in 1:4
                A[I, nodes[b]] += Aᵉ[a, b]
            end
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Apply prescribed Neumann loads at all nodes before Dirichlet masking.
    # ─────────────────────────────────────────────────────────────────
    for i in 1:mesh.num_nodes
        if q_boundary_water[i] != 0.0
            R[i] += q_boundary_water[i]
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Apply solution-dependent unit-gradient free-drainage loads before masking.
    # q_in = -(Kx*e_g[1]*n_x + Ky*e_g[2]*n_y), with n outward.
    # ─────────────────────────────────────────────────────────────────
    for edge in free_drainage_edges
        ni, nj, _ = edge
        integrated_flux = free_drainage_edge_flux(h_curr, mesh, elem_props, e_g, edge)
        R[ni] += integrated_flux / 2.0
        R[nj] += integrated_flux / 2.0
    end

    # Capture conventional support reactions after all boundary loads. Since free-node
    # residuals vanish at convergence, the exact balance identity is
    # ΔW/Δt = -sum(reaction_out) + sum(q_boundary_water) + free_drainage_flux.
    if reaction_out !== nothing
        length(reaction_out) == mesh.num_nodes ||
            throw(DimensionMismatch("reaction_out must have length $(mesh.num_nodes)"))
        fill!(reaction_out, 0.0)
        @inbounds for i in 1:mesh.num_nodes
            if P_boundary_water[i] == 0
                reaction_out[i] = R[i]
            end
        end
    end

    # ─────────────────────────────────────────────────────────────────
    # Apply P_boundary_water masking: zero matrix rows and residual at BC nodes.
    # Boundary loads at prescribed nodes are discarded here and cannot move h.
    # ─────────────────────────────────────────────────────────────────
    n = mesh.num_nodes
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1] - 1)
            i = A.rowval[k]
            if P_boundary_water[i] == 0
                A.nzval[k] = (i == j) ? 1.0 : 0.0
            end
        end
    end
    for i in 1:n
        if P_boundary_water[i] == 0
            R[i] = 0.0
        end
    end

    return nothing
end




# Water Neumann BCs are integrated once at solver startup by
# assemble_water_neumann!(). The resulting q_boundary_water vector is added to
# the residual in assemble_richards!().


# Picard iteration solver (one time step)

"""
    picard_richards!(h_curr, h_prev, mesh, elem_props, Δt, e_g, A, cache,
                     q_boundary_water, free_drainage_edges; kwargs...)

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
- q_boundary_water: Integrated water Neumann residual contributions
- free_drainage_edges: Explicitly marked free-drainage boundary edges
- tol = 1e-8: Convergence tolerance
- max_iter = 100: Maximum iterations (safety limit)
- ω = 1.0: Relaxation parameter (Newton step size)
- reaction_out: Optional vector receiving final-iterate Dirichlet reactions
- c_floor: Lower bound on capacity used only for matrix conditioning
- floor_bound_nodes: Optional run-level mask of nodes where the floor was active

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
    q_boundary_water :: Vector{Float64},
    free_drainage_edges :: Vector{Tuple{Int,Int,Int}};
    tol        :: Float64 = 1e-8,
    max_iter   :: Int     = 100,
    ω          :: Float64 = 1.0,
    reaction_out :: Union{Nothing,Vector{Float64}} = nothing,
    c_floor    :: Float64 = CAPACITY_FLOOR_DEFAULT,
    floor_bound_nodes :: Union{Nothing,BitVector} = nothing
)
    N = mesh.num_nodes
    R = zeros(Float64, N)
    
    # ── Store initial Dirichlet BC values from mesh ──────────────────────────
    # Nodes where pressure_head_bc dictionary has entries are Dirichlet BCs
    dbc_nodes = collect(keys(mesh.pressure_head_bc))
    dbc_vals  = [mesh.pressure_head_bc[i] for i in dbc_nodes]
    
    res_history = Float64[]

    # Cold-start enforcement is necessary because h_prev may still contain the
    # initial value at nodes whose prescribed head differs at the new step.
    h_curr .= h_prev
    for (i, val) in zip(dbc_nodes, dbc_vals)
        h_curr[i] = val
    end

    for m in 1:max_iter
        # Step 1: Assemble nonlinear system (A and R depend on current h)
        assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache,
                   q_boundary_water, free_drainage_edges; reaction_out=reaction_out, c_floor=c_floor,
                           floor_bound_nodes=floor_bound_nodes)

        # Step 2: Solve for pressure head increments
        delta = A \ R
        delta_norm = maximum(abs.(delta))
        max_dbc_delta = isempty(dbc_nodes) ? 0.0 : maximum(abs(delta[i]) for i in dbc_nodes)
        if max_dbc_delta > 1e-10
            @warn "Non-zero Picard increment at Dirichlet nodes: max|δ| = $max_dbc_delta. " *
                  "mesh.pressure_head_bc and P_boundary_water disagree about which nodes " *
                  "are constrained."
        end
        push!(res_history, delta_norm)

        # Step 3: Check convergence
        if delta_norm < tol
            return m, res_history
        end

        # Step 4: Update pressure head
        h_curr .+= ω .* delta
        
        # Identity rows make δ_i zero at constrained nodes; this only absorbs
        # round-off from the sparse solve.
        for (i, val) in zip(dbc_nodes, dbc_vals)
            h_curr[i] = val
        end
    end

    # Non-convergence: Backward Euler is unconditionally stable so solution is accepted
    @warn "Picard did not converge in $max_iter iterations. ||δ|| = $(res_history[end])"
    if reaction_out !== nothing
        assemble_richards!(A, R, h_curr, h_prev, mesh, elem_props, Δt, e_g, cache,
                   q_boundary_water, free_drainage_edges; reaction_out=reaction_out, c_floor=c_floor,
                           floor_bound_nodes=floor_bound_nodes)
    end
    return max_iter, res_history
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
        - `picard_history`: Picard iteration count at each output
        - `mb_history`: Mass-balance ratio at each accepted time step
        - `flux_cumul`: Cumulative inflow recovered from Dirichlet reactions

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

    # ── Build sparsity pattern ────────────────────────────────────────
    # Note: Cache and elem_props are precomputed in kernel.jl and passed as arguments
    A = build_richards_sparsity(mesh)
    log_print(@sprintf("   ✓ Sparse matrix: %d nonzeros (%.2f%% density)",
                        nnz(A), 100.0 * nnz(A) / mesh.num_nodes^2))

    # ── Mark Dirichlet BCs in P_boundary_water masking matrix ─────────
    apply_water_dirichlet_bc!(mesh, materials)
    n_bc_nodes = count(P_boundary_water .== 0)
    log_print("   ✓ Dirichlet BC nodes marked: $n_bc_nodes nodes")

    # ── Integrate water Neumann boundary conditions ───────────────────
    neumann_edges = build_water_neumann_edges(mesh)
    q_boundary_water = zeros(Float64, mesh.num_nodes)
    assemble_water_neumann!(q_boundary_water, mesh, neumann_edges)
    log_print("   ✓ Water Neumann edges: $(length(neumann_edges))")
    log_print("   ✓ Water Neumann residual nodes: $(count(q_boundary_water .!= 0.0))")

    # ── Precompute explicitly marked free-drainage boundary edges ─────
    # Conductivity-dependent flux is evaluated inside each Picard iteration.
    free_drainage_edges = build_free_drainage_edges(mesh)
    log_print("   ✓ Explicit free-drainage edges: $(length(free_drainage_edges))")


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
    nonlinear_solver = get(calc_params, "nonlinear_solver", Dict{String, Any}())
    picard_tol = Float64(get(nonlinear_solver, "picard_tolerance", 1.0e-8))
    picard_max_iter = Int(get(nonlinear_solver, "picard_max_iter", 100))
    picard_relaxation = Float64(get(nonlinear_solver, "picard_relaxation", 1.0))
    capacity_floor = Float64(get(nonlinear_solver, "capacity_floor", CAPACITY_FLOOR_DEFAULT))
    capacity_floor >= 0.0 || error("nonlinear_solver.capacity_floor must be non-negative")

    log_print(@sprintf("   Picard settings: tol = %.3e, max_iter = %d, ω = %.2f",
                        picard_tol, picard_max_iter, picard_relaxation))
    log_print(@sprintf("   Capacity floor: %.3e 1/m", capacity_floor))
    step_count = 0
    picard_history = Int[]  # n_iter per output step
    reaction = zeros(Float64, mesh.num_nodes)
    W_initial = compute_total_water(h, mesh, elem_props, cache)
    W_current = W_initial
    flux_cumul = 0.0
    mb_history = Float64[]
    floor_bound_nodes = falses(mesh.num_nodes)

    while current_time < total_time - 1e-10
        step_count += 1

        # ── Clamp dt to hit the next output time or total_time ────────
        dt_step = min(dt, next_output_time - current_time, total_time - current_time)
        # Note: at_output checked AFTER the step based on time reached,
        # not whether dt was shortened (which misses exact-hit cases).

        # ── Picard iteration ──────────────────────────────────────────
        n_iter, res_history = picard_richards!(h_new, h, mesh, elem_props, dt_step, e_g,
                                   A, cache, q_boundary_water, free_drainage_edges;
                                   tol=picard_tol, max_iter=picard_max_iter,
                                   ω=picard_relaxation, reaction_out=reaction,
                                   c_floor=capacity_floor,
                                   floor_bound_nodes=floor_bound_nodes)

        # ── Accept step ───────────────────────────────────────────────
        h .= h_new
        current_time += dt_step
        enforce_water_dirichlet_bc!(mesh, materials)

        W_current = compute_total_water(h, mesh, elem_props, cache)
        free_drainage_flux = compute_free_drainage_flux(
            h, mesh, elem_props, e_g, free_drainage_edges)
        net_boundary_flux = compute_net_boundary_flux(
            reaction, q_boundary_water, free_drainage_flux)
        flux_cumul += net_boundary_flux * dt_step
        ΔW_total = W_current - W_initial
        MB = abs(flux_cumul) > eps(Float64) ? ΔW_total / flux_cumul : NaN
        push!(mb_history, MB)

        # ── Output at scheduled times and at final time ───────────────
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

    finite_mb = filter(isfinite, mb_history)
    if isempty(finite_mb)
        log_print("   Mass balance MB = ΔW/flux_cumul unavailable (zero cumulative flux)")
    else
        log_print(@sprintf("   Mass balance MB = ΔW/flux_cumul: min = %.8f, max = %.8f, mean = %.8f",
                           minimum(finite_mb), maximum(finite_mb),
                           sum(finite_mb) / length(finite_mb)))
    end
        n_floor_bound_nodes = count(floor_bound_nodes)
        log_print("   Capacity floor bound at $n_floor_bound_nodes unique nodes")

    return (current_time=current_time, output_counter=output_counter,
            next_output_time=next_output_time, picard_history=picard_history,
            mb_history=mb_history, flux_cumul=flux_cumul,
            total_water_initial=W_initial, total_water_final=W_current,
            capacity_floor=capacity_floor,
            capacity_floor_bound_nodes=n_floor_bound_nodes)
end
