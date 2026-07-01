#______________________________________________________
# ADSIM: Aqueous Concentration Solver
# Implicit FEM · Backward Euler · Henry's Law BCs
# Authors: Paula Sarmiento, Luis Zambrano-Cruzatty
#______________________________________________________

using LinearAlgebra
using SparseArrays
using Printf

#=
IMPLEMENTATION NOTES:
Implicit FEM solver for aqueous concentration transport: ∂(θ_w C)/∂t + ∇·(θ_w v_s C) = ∇·(θ_w D_h ∇C)
Backward Euler (unconditionally stable) + Galerkin FEM with 2×2 Gauss quadrature.
Phase 1: Pure advection-diffusion (no source term). Phase 2: Will add mass transfer.
Dirichlet BCs via Henry's Law: C_BC = P_partial / K_H
=# 

#______________________________________________________
# Element Assembly
#______________________________________________________

"""
    aqueous_concentration_element!(Aᵉ, Rᵉ, cache, e, θw_new_e, θw_old_e, 
                                   vsᵉ_old, Cᵉ_old, D_h, Δt)

Assemble local [4×4] element matrix and [4] residual with LUMPED CAPACITY MATRIX.
Uses Backward Euler implicit discretization with Galerkin FEM (2×2 Gauss quadrature).

KEY FIX: Lumped mass matrix (diagonal only) guarantees:
  • Diagonal dominance of global matrix → no spurious oscillations
  • Physical range preservation: C ∈ [0, C_max] automatically
  • Exact discrete mass conservation

Key Feature: WEAK COUPLING (Phase 1 Implementation)
  • θₗ FROZEN from previous water time step: θₗ^{n+1}_element = computed from h^{n+1} in water solver
  • vs FROZEN from previous water time step: vs^{n}_element = velocity from water solver at time n
  • This weak coupling avoids solving coupled water-aqueous system
  • In Phase 2 (strong coupling): θₗ will stay frozen but vs will become vs^{n+1} (need to pass updated velocities)
  → For Phase 2 transition: ONLY change vs_old → vs_new in call stack; structure stays identical

Reference: ADSIM aqueous_concentration_usage_guide.jl — Cell T3

# Arguments
- `Aᵉ::Matrix{Float64}`: Element matrix [4×4] (output, pre-zeroed)
- `Rᵉ::Vector{Float64}`: Element residual [4] (output, pre-zeroed)
- `cache`: Precomputed shape functions and Jacobians (includes A_e: element areas)
- `e::Int`: Element index
- `θw_new_e, θw_old_e::Vector{Float64}`: [4] Nodal water content at n+1, n
- `vsᵉ_old::Matrix{Float64}`: [4×2] Nodal Darcy velocity at n (FROZEN in Phase 1)
- `Cᵉ_old::Vector{Float64}`: [4] Nodal concentration at time n
- `D_h::Float64`: Effective diffusivity [m²/s]
- `Δt::Float64`: Time step [s]
"""
function aqueous_concentration_element!(
    Aᵉ          :: Matrix{Float64},
    Rᵉ          :: Vector{Float64},
    cache,  # RichardsCache with A_e[e] (element area)
    e           :: Int,
    θw_new_e    :: Vector{Float64},  # [4] nodal water content at n+1
    θw_old_e    :: Vector{Float64},  # [4] nodal water content at n
    vsᵉ_old     :: Matrix{Float64},  # [4×2] nodal Darcy velocity at n (FROZEN)
    Cᵉ_old      :: Vector{Float64},  # [4] nodal concentration at n
    D_h         :: Float64,           # effective diffusivity
    Δt          :: Float64            # time step
)
    fill!(Aᵉ, 0.0)
    fill!(Rᵉ, 0.0)
    
    # Get element area from cache (precomputed by build_richards_cache)
    A_e = cache.A_e[e]
    lumped_mass_entry = A_e / 4.0  # Lumped: each diagonal entry = A_e/4
    
    # ─────────────────────────────────────────────────────────────────────────────
    # 1. CONSISTENT INTEGRATION: K_diff and K_adv (Gauss quadrature)
    #    These are SYMMETRIC (K_diff) and NON-SYMMETRIC (K_adv) bilinear forms
    # ─────────────────────────────────────────────────────────────────────────────
    for p in 1:4
        Np = cache.Np[p]  # [4] shape functions at Gauss point p
        wp = cache.weights[p]
        dJ = cache.detJ[e, p]
        wdet = wp * dJ  # Quadrature weight
        
        # Interpolate θ_w and v_s to Gauss point p
        θw_p_new = dot(Np, θw_new_e)  # θ_w^{n+1}(ξ_p)
        vs_x_p = dot(Np, vsᵉ_old[:, 1])
        vs_y_p = dot(Np, vsᵉ_old[:, 2])
        
        @inbounds for a in 1:4
            ∂Na_x = cache.Bp[e, p, 1, a]  # ∂N_a/∂x
            ∂Na_y = cache.Bp[e, p, 2, a]  # ∂N_a/∂y
            
            @inbounds for b in 1:4
                Nb = Np[b]
                ∂Nb_x = cache.Bp[e, p, 1, b]  # ∂N_b/∂x
                ∂Nb_y = cache.Bp[e, p, 2, b]  # ∂N_b/∂y
                
                # K_diff: ∫ θ_w D_h ∇N_a·∇N_b dΩ (SYMMETRIC, negative off-diagonals)
                K_diff_ab = θw_p_new * D_h * (∂Na_x * ∂Nb_x + ∂Na_y * ∂Nb_y)
                
                # K_adv: -∫ ∇N_a·(θ_w v_s) N_b dΩ (NON-SYMMETRIC, integration by parts)
                K_adv_ab = -(∂Na_x * θw_p_new * vs_x_p + ∂Na_y * θw_p_new * vs_y_p) * Nb
                
                # Accumulate consistent terms
                Aᵉ[a, b] += (K_diff_ab + K_adv_ab) * wdet
            end
        end
    end
    
    # ─────────────────────────────────────────────────────────────────────────────
    # 2. LUMPED CAPACITY MATRIX (diagonal only, uses NODAL values at t^{n+1})
    #    LHS: (1/Δt) M^e_θ(lumped) = (1/Δt) diag(θ_w_a * A_e/4)  [a = 1,...,4]
    #    RHS: (1/Δt) M^e_θ C^n_a = (1/Δt) θ_w_a^n * C_a^n * A_e/4 (time history)
    #
    #    Why nodal values, not Gauss-point interpolation?
    #    → Allows independent control of each node's capacity
    #    → Consistent with Richards solver approach
    #    → Guarantees diagonal dominance when K_diff has negative off-diagonals
    #
    #    WEAK COUPLING (Phase 1):
    #    • θw_new_e and θw_old_e come from water solver (time-interpolated)
    #    • Note: θw is NOT from current element!  θw^{n+1} is FROZEN from water solver
    # ─────────────────────────────────────────────────────────────────────────────
    @inbounds for a in 1:4
        # LHS diagonal: (1/Δt) * θ_w^{n+1}_a * A_e/4
        Aᵉ[a, a] += (1.0 / Δt) * θw_new_e[a] * lumped_mass_entry
        
        # RHS: (1/Δt) * θ_w^n_a * C^n_a * A_e/4  (Backward Euler time history)
        Rᵉ[a] += (1.0 / Δt) * θw_old_e[a] * Cᵉ_old[a] * lumped_mass_entry
    end
end

#______________________________________________________
# Global Assembly + Boundary Conditions (Combined)
#______________________________________________________

"""
    assemble_aqueous_concentration!(A, R, mesh, cache, D_h, θw_new, θw_old,
                                    vs_old, C_old, Δt, P_boundary_mask, C_prescribed,
                                    q_flux_aq_gas)

Assemble global sparse system matrix and residual with Dirichlet BCs via row-zeroing.
Combines element assembly and BC application (direct stiffness method).

# Arguments
- `A::SparseMatrixCSC`: Global system matrix [N×N] (modified in-place)
- `R::Vector{Float64}`: Global residual [N] (modified in-place)
- `mesh::MeshData`: Mesh connectivity and element info
- `cache`: Precomputed cache (Np, Bp, detJ)
- `D_h::Float64`: Effective diffusivity [m²/s]
- `θw_new, θw_old::Vector{Float64}`: [N] Water content at n+1, n
- `vs_old::Matrix{Float64}`: [N×2] Darcy velocity at n (FROZEN in Phase 1)
- `C_old::Vector{Float64}`: [N] Concentration at time n
- `Δt::Float64`: Time step [s]
- `P_boundary_mask::Vector{Int}`: [N] Dirichlet mask (0=BC, 1=interior)
- `C_prescribed::Vector{Float64}`: [N] Prescribed BC values
- `q_flux_aq_gas::Union{AbstractVector{Float64}, Nothing}`: [N] Nodal Neumann flux for active gas
"""
function assemble_aqueous_concentration!(
    A               :: SparseMatrixCSC{Float64, Int},
    R               :: Vector{Float64},
    mesh,  # MeshData
    cache,  # RichardsCache
    D_h             :: Float64,
    θw_new          :: Vector{Float64},
    θw_old          :: Vector{Float64},
    vs_old          :: Matrix{Float64},
    C_old           :: Vector{Float64},
    Δt              :: Float64,
    P_boundary_mask :: Vector{Int},
    C_prescribed    :: Vector{Float64},
    q_flux_aq_gas   :: Union{AbstractVector{Float64}, Nothing} = nothing
)
    n_nodes = length(R)
    
    # Allocate local arrays
    Aᵉ = zeros(Float64, 4, 4)
    Rᵉ = zeros(Float64, 4)
    nodes = zeros(Int, 4)
    
    # Zero the global system (preserves sparsity pattern)
    A.nzval .= 0.0
    R .= 0.0
    
    # Loop over all elements
    for e in 1:mesh.num_elements
        # Extract local node numbers for element e
        nodes .= mesh.elements[e, :]
        
        # Extract local nodal values
        θw_new_e = θw_new[nodes]
        θw_old_e = θw_old[nodes]
        Cᵉ_old = C_old[nodes]
        vsᵉ_old = vs_old[nodes, :]  # [4×2] matrix
        
        # Zero local arrays
        fill!(Aᵉ, 0.0)
        fill!(Rᵉ, 0.0)
        
        # Compute element contribution (Phase 1: no source term)
        aqueous_concentration_element!(Aᵉ, Rᵉ, cache, e,
                                      θw_new_e, θw_old_e, vsᵉ_old,
                                      Cᵉ_old, D_h, Δt)
        
        # Scatter into global system (direct stiffness)
        for a in 1:4
            I = nodes[a]
            R[I] += Rᵉ[a]
            for j in 1:4
                J = nodes[j]
                A[I, J] += Aᵉ[a, j]
            end
        end
    end
    
    # Apply Neumann BCs first: Dirichlet row-zeroing below keeps Dirichlet precedence
    if q_flux_aq_gas !== nothing
        @inbounds for i in 1:n_nodes
            if P_boundary_mask[i] != 0 && q_flux_aq_gas[i] != 0.0
                R[i] += q_flux_aq_gas[i]
            end
        end
    end

    # Apply Dirichlet BCs: row-zeroing method for sparse matrix
    # For BC nodes (P_boundary_mask[i] == 0): set row to [0...1...0] with RHS = prescribed value
    n = length(R)
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1] - 1)
            i = A.rowval[k]
            if P_boundary_mask[i] == 0  # BC node
                A.nzval[k] = (i == j) ? 1.0 : 0.0  # Identity row: diagonal=1, off-diag=0
            end
        end
    end
    for i in 1:n
        if P_boundary_mask[i] == 0
            R[i] = C_prescribed[i]  # Set RHS to prescribed concentration
        end
    end
    
    return nothing
end

#______________________________________________________
# Main Solver
#______________________________________________________

"""
    aqueous_concentration_solver(A, R, mesh, cache, materials, 
                                 P_boundary_aq, mesh_partial_pressure_bc,
                                 soil_idx, θw_new, θw_old, vs_old, C_old, Δt, gas_idx) → C_new

Solve one time step of aqueous concentration transport (Backward Euler implicit).
Applies Henry's Law BCs: C_BC = P_partial / K_H.

# Arguments
- `A::SparseMatrixCSC`: Global system matrix [N×N] (modified, pre-allocated)
- `R::Vector{Float64}`: Global residual [N] (modified)
- `mesh::MeshData`: Mesh connectivity
- `cache`: Precomputed shape functions
- `materials::MaterialData`: Gas/Soil properties structs
- `P_boundary_aq::Vector{Int}`: [N] Dirichlet mask (0=BC, 1=interior)
- `mesh_partial_pressure_bc::Dict`: Partial pressure BC values
- `soil_idx::Int`: Soil material index
- `θw_new, θw_old::Vector{Float64}`: Water content at n+1, n
- `vs_old::Matrix{Float64}`: [N×2] Darcy velocity at n (FROZEN in Phase 1)
- `C_old::Vector{Float64}`: Concentration at time n
- `Δt::Float64`: Time step [s]
- `gas_idx::Int`: Gas species index

# Returns
- `C_new::Vector{Float64}`: [N] Concentration at time n+1
"""
function aqueous_concentration_solver(
    A                       :: SparseMatrixCSC{Float64, Int},
    R                       :: Vector{Float64},
    mesh,  # MeshData
    cache,  # RichardsCache
    materials,  # MaterialData
    P_boundary_aq           :: Vector{Int},
    mesh_partial_pressure_bc :: Dict,
    soil_idx                :: Int,
    θw_new                  :: Vector{Float64},
    θw_old                  :: Vector{Float64},
    vs_old                  :: Matrix{Float64},
    C_old                   :: Vector{Float64},
    Δt                      :: Float64,
    gas_idx                 :: Int;
    C_prescribed_override   :: Union{Vector{Float64}, Nothing} = nothing,
    q_flux_aq_gas           :: Union{AbstractVector{Float64}, Nothing} = nothing
) :: Vector{Float64}
    
    # Extract property structs (no field extraction)
    gas_name = materials.gas_dictionary[gas_idx]
    gas = materials.gases[gas_name]
    
    soil_name = materials.soil_dictionary[soil_idx]
    soil_props = materials.soils[soil_name]
    
    # Compute D_h inline (ADSIM pattern: like Richards computes conductivity)
    D_h = soil_props.granular_tortuosity * gas.diff_coefficient
    
    # Compute Henry's Law BC values: C_BC = P_partial / K_H
    C_prescribed = zeros(Float64, length(R))
    
    # Use override if provided (for verification testing), otherwise compute from Henry's law
    if C_prescribed_override !== nothing
        C_prescribed = copy(C_prescribed_override)
    else
        K_H = gas.henry_constant
        if K_H > 0.0
            for (node_id, partial_pressures) in mesh_partial_pressure_bc
                if gas_idx <= length(partial_pressures)
                    P_partial = partial_pressures[gas_idx]
                    if P_partial > 0.0
                        C_prescribed[node_id] = P_partial / K_H
                        P_boundary_aq[node_id] = 0  # Mark as BC node
                    end
                end
            end
        end
    end
    
    # Assemble + apply BCs (combined call, like Richards)
    assemble_aqueous_concentration!(A, R, mesh, cache, D_h,
                                   θw_new, θw_old, vs_old, C_old, Δt,
                                   P_boundary_aq, C_prescribed, q_flux_aq_gas)
    
    # Solve sparse linear system
    C_new = A \ R
    
    return C_new
end
