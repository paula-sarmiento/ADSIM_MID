#______________________________________________________
# ADSIM: Aqueous Concentration Solver
# Q4 Isoparametric · Backward Euler · Decoupled Phase
# Solves aqueous transport from dissolved gas species
# Authors: Paula Sarmiento, Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Implicit FEM solver for aqueous concentration transport
# Decoupled from Richards water flow (receives frozen θ_w, v_s)
# Applies Henry's Law BCs from gas partial pressure
#______________________________________________________

using LinearAlgebra
using SparseArrays

#=
IMPLEMENTATION NOTES:
=====================

This module implements an implicit finite element solver for aqueous
concentration transport of dissolved gases:

  ∂(θ_w C)/∂t + ∇·(θ_w v_s C) = ∇·(θ_w D_h ∇C) + k_g(C_eq - C)

Key features:
  • Backward Euler time discretization (unconditionally stable)
  • Advection-diffusion-reaction coupling (all implicit)
  • Phase 1: k_g = 0 (no reaction/equilibration, deferred to Phase 2)
  • Galerkin FEM with 2×2 Gauss quadrature
  • Dirichlet BCs via Henry's Law: C_BC = P_partial / K_H
  • Frozen water content θ_w^n and seepage velocity v_s at time n
  • Receives data from implicit_richards_solver output

Dependencies:
  read_materials.jl       → Gas properties with henry_constant
  shape_functions.jl      → Precomputed cache (Np, Bp, detJ, weights)
  implicit_richards_solver.jl → Water content θ_w, velocity v_s
  swrc_models.jl          → theta() function for phase 1

Architecture:
  • Element assembly: aqueous_concentration_element!()
  • Global assembly: assemble_aqueous_concentration!()
  • BC application: apply_aqueous_concentration_bc!()
  • Main solver: aqueous_concentration_solver()

References:
  • Galerkin FEM for advection-diffusion-reaction
  • Henry's Law equilibration (Phase 1)
  • Backward Euler implicit scheme for unconditional stability
=#

#______________________________________________________
# Transport Parameters and Defaults
#______________________________________________________

"""
    AqueousTransportParams

Container for aqueous concentration transport parameters.

# Fields
- `D_h::Float64`: Effective diffusivity [m²/s] (scalar, isotropic, Phase 1)
"""
mutable struct AqueousTransportParams
    D_h::Float64  # Effective diffusivity [m²/s]
    
    function AqueousTransportParams(D_h::Float64=0.0)
        new(D_h)
    end
end

#______________________________________________________
# Element Assembly
#______________________________________________________

"""
    aqueous_concentration_element!(Ae, be, cache, e, θw_new_e, θw_old_e, 
                                   vs_e, C_old_e, kg_e, Ceq_e, D_h, Δt)

Assemble local element matrix and RHS vector for aqueous concentration transport.

**Mathematical Basis (Backward Euler, Galerkin FEM)**:

Weak form of transport PDE after integration by parts:
```
∫ w ∂(θ_w C)/∂t dΩ - ∫ ∇w·(θ_w v_s C) dΩ + ∫ ∇w·θ_w D_h∇C dΩ = ∫ w k_g(C_eq - C) dΩ
```

Backward Euler time discretization:
```
[1/Δt M_θ + K_adv + K_diff + K_rxn] C^{n+1} = (1/Δt) M_θ^n C^n + f_rxn^n
```

Element matrices:
- M_θ    = ∫ θ_w N_a N_b dΩ  (capacity, uses θ_w^{n+1})
- K_adv  = -∫ ∇N_a·(θ_w v_s) N_b dΩ  (advection, minus from IBP)
- K_diff = ∫ θ_w D_h ∇N_a·∇N_b dΩ  (diffusion, symmetric)
- K_rxn  = ∫ k_g N_a N_b dΩ  (reaction, on LHS)
- RHS: (1/Δt) ∫ θ_w^n C^n N_a dΩ + ∫ k_g C_eq N_a dΩ

Quadrature: 2×2 Gauss points (4 points per element)

# Arguments
- `Ae::Matrix{Float64}`: Element matrix [4×4] (output, pre-zeroed)
- `be::Vector{Float64}`: Element RHS [4] (output, pre-zeroed)
- `cache::RichardsCache`: Precomputed shape functions, Jacobians
- `e::Int`: Element index
- `θw_new_e::Vector{Float64}`: [4] nodal water content at time n+1
- `θw_old_e::Vector{Float64}`: [4] nodal water content at time n
- `vs_e::Matrix{Float64}`: [4×2] nodal Darcy velocity (columns: x, y components)
- `C_old_e::Vector{Float64}`: [4] nodal concentration at time n
- `kg_e::Vector{Float64}`: [4] nodal mass transfer coefficient (frozen at n)
- `Ceq_e::Vector{Float64}`: [4] nodal equilibrium concentration at time n+1
- `D_h::Float64`: Effective diffusivity [m²/s]
- `Δt::Float64`: Time step [s]

# Notes
- Modified in-place for performance
- Phase 1: k_g = 0.0 typically (reaction deferred)
- C_eq = 0.0 typically (no source term in Phase 1)
"""
function aqueous_concentration_element!(
    Ae          :: Matrix{Float64},
    be          :: Vector{Float64},
    cache,  # RichardsCache
    e           :: Int,
    θw_new_e    :: Vector{Float64},
    θw_old_e    :: Vector{Float64},
    vs_e        :: Matrix{Float64},
    C_old_e     :: Vector{Float64},
    kg_e        :: Vector{Float64},
    Ceq_e       :: Vector{Float64},
    D_h         :: Float64,
    Δt          :: Float64
)
    fill!(Ae, 0.0)
    fill!(be, 0.0)
    
    # Loop over 4 Gauss points
    for p in 1:4
        # Shape functions and weights at Gauss point p
        Np = cache.Np[p]  # [4] shape functions at point p
        wp = cache.weights[p]  # Always 1.0 for 2×2 rule, but explicit for clarity
        dJ = cache.detJ[e, p]  # Jacobian determinant at point p
        wdet = wp * dJ  # Quadrature weight: w_p × det(J)
        
        # Interpolate nodal values to Gauss point p
        θw_p_new = dot(Np, θw_new_e)  # θ_w^{n+1}(ξ_p)
        θw_p_old = dot(Np, θw_old_e)  # θ_w^n(ξ_p)
        C_old_p  = dot(Np, C_old_e)   # C^n(ξ_p)
        kg_p     = dot(Np, kg_e)      # k_g^n(ξ_p) - frozen at n
        Ceq_p    = dot(Np, Ceq_e)     # C_eq^{n+1}(ξ_p)
        
        vs_x_p = dot(Np, vs_e[:, 1])  # v_s,x(ξ_p)
        vs_y_p = dot(Np, vs_e[:, 2])  # v_s,y(ξ_p)
        
        # Assemble element matrix and RHS
        for a in 1:4
            Na = Np[a]  # Shape function N_a at Gauss point p
            
            for b in 1:4
                Nb = Np[b]  # Shape function N_b at Gauss point p
                
                # Shape function gradients (physical coordinates)
                ∂Na_x = cache.Bp[e, p, 1, a]  # ∂N_a/∂x
                ∂Na_y = cache.Bp[e, p, 2, a]  # ∂N_a/∂y
                
                ∂Nb_x = cache.Bp[e, p, 1, b]  # ∂N_b/∂x
                ∂Nb_y = cache.Bp[e, p, 2, b]  # ∂N_b/∂y
                
                # ──── Capacity matrix: M_θ = ∫ θ_w N_a N_b dΩ ────
                M_ab = θw_p_new * Na * Nb
                
                # ──── Advection matrix: K_adv = -∫ ∇N_a·(θ_w v_s) N_b dΩ ────
                # After IBP: derivative acts on test function (a), minus sign from IBP
                K_adv_ab = -(∂Na_x * θw_p_new * vs_x_p + 
                            ∂Na_y * θw_p_new * vs_y_p) * Nb
                
                # ──── Diffusion matrix: K_diff = ∫ θ_w D_h ∇N_a·∇N_b dΩ ────
                # Symmetric bilinear form
                K_diff_ab = θw_p_new * D_h * (∂Na_x * ∂Nb_x + ∂Na_y * ∂Nb_y)
                
                # ──── Reaction matrix: K_rxn = ∫ k_g N_a N_b dΩ ────
                # From -k_g·C term on LHS (when moved from RHS)
                K_rxn_ab = kg_p * Na * Nb
                
                # Assemble to element matrix
                Ae[a, b] += ((1.0 / Δt) * M_ab + K_adv_ab + K_diff_ab + K_rxn_ab) * wdet
            end
            
            # ──── RHS assembly ────
            # From backward Euler: (1/Δt) ∫ θ_w^n C^n N_a dΩ + ∫ k_g C_eq N_a dΩ
            be[a] += ((1.0 / Δt) * θw_p_old * C_old_p + kg_p * Ceq_p) * Na * wdet
        end
    end
end

#______________________________________________________
# Global Assembly (Direct Stiffness Method)
#______________________________________________________

"""
    assemble_aqueous_concentration!(A, b, mesh, cache, params, θw_new, θw_old, 
                                    vs, C_old, kg, Ceq, Δt)

Assemble global sparse system matrix and RHS via direct stiffness method.

Iterates over all elements, computes local contributions, and scatters into
global sparse matrix A and vector b using element-to-global node mapping.

# Arguments
- `A::SparseMatrixCSC`: Global system matrix [N×N] (output, pre-allocated)
- `b::Vector{Float64}`: Global RHS [N] (output)
- `mesh::MeshData`: Mesh connectivity and element info
- `cache::RichardsCache`: Precomputed cache (Np, Bp, detJ)
- `params::AqueousTransportParams`: Transport parameters (D_h)
- `θw_new::Vector{Float64}`: [N] water content at time n+1
- `θw_old::Vector{Float64}`: [N] water content at time n
- `vs::Matrix{Float64}`: [N×2] Darcy velocity (columns: x, y)
- `C_old::Vector{Float64}`: [N] concentration at time n
- `kg::Vector{Float64}`: [N] mass transfer coefficient (frozen)
- `Ceq::Vector{Float64}`: [N] equilibrium concentration at time n+1
- `Δt::Float64`: Time step [s]

# Notes
- Zeroes A and b at start (non-destructive on pattern)
- Loop over mesh.num_elements
- Extract local nodes via mesh.elements[e, :]
"""
function assemble_aqueous_concentration!(
    A           :: SparseMatrixCSC{Float64, Int},
    b           :: Vector{Float64},
    mesh,  # MeshData
    cache,  # RichardsCache
    params      :: AqueousTransportParams,
    θw_new      :: Vector{Float64},
    θw_old      :: Vector{Float64},
    vs          :: Matrix{Float64},
    C_old       :: Vector{Float64},
    kg          :: Vector{Float64},
    Ceq         :: Vector{Float64},
    Δt          :: Float64
)
    n_nodes = length(b)
    
    # Allocate local arrays
    Ae = zeros(Float64, 4, 4)
    be = zeros(Float64, 4)
    nodes = zeros(Int, 4)
    
    # Zero the global system (preserves sparsity pattern)
    A.nzval .= 0.0
    b .= 0.0
    
    # Loop over all elements
    for e in 1:mesh.num_elements
        # Extract local node numbers for element e
        nodes .= mesh.elements[e, :]
        
        # Extract local nodal values
        θw_new_e = θw_new[nodes]
        θw_old_e = θw_old[nodes]
        C_old_e = C_old[nodes]
        kg_e = kg[nodes]
        Ceq_e = Ceq[nodes]
        vs_e = vs[nodes, :]  # [4×2] matrix
        
        # Zero local arrays
        fill!(Ae, 0.0)
        fill!(be, 0.0)
        
        # Compute element contribution
        aqueous_concentration_element!(Ae, be, cache, e,
                                      θw_new_e, θw_old_e, vs_e,
                                      C_old_e, kg_e, Ceq_e,
                                      params.D_h, Δt)
        
        # Scatter into global system (direct stiffness)
        for a in 1:4
            I = nodes[a]
            b[I] += be[a]
            for j in 1:4
                J = nodes[j]
                A[I, J] += Ae[a, j]
            end
        end
    end
end

#______________________________________________________
# Boundary Condition Application
#______________________________________________________

"""
    apply_aqueous_concentration_bc!(A, b, P_boundary_mask, C_prescribed)

Apply Dirichlet boundary conditions via row-zeroing method.

For each row i where P_boundary_mask[i] == 0 (BC node):
  - Zero all entries in row i except diagonal
  - Set A[i,i] = 1.0
  - Set b[i] = C_prescribed[i]

# Arguments
- `A::SparseMatrixCSC`: Global system matrix (modified in-place)
- `b::Vector{Float64}`: Global RHS (modified in-place)
- `P_boundary_mask::Vector{Int}`: [N] Dirichlet mask (0=BC, 1=interior)
- `C_prescribed::Vector{Float64}`: [N] prescribed values at BC nodes

# Notes
- Efficient for sparse matrices
- Preserves matrix sparsity pattern
- Called after global assembly
"""
function apply_aqueous_concentration_bc!(
    A               :: SparseMatrixCSC{Float64, Int},
    b               :: Vector{Float64},
    P_boundary_mask :: Vector{Int},
    C_prescribed    :: Vector{Float64}
)
    n = length(b)
    
    for i in 1:n
        if P_boundary_mask[i] == 0  # This is a BC node
            # Zero row i except diagonal, set diagonal to 1
            for k in A.colptr[i]:(A.colptr[i+1] - 1)
                if A.rowval[k] == i
                    A.nzval[k] = 1.0  # Diagonal entry
                else
                    A.nzval[k] = 0.0  # Off-diagonal entries
                end
            end
            # Set RHS to prescribed value
            b[i] = C_prescribed[i]
        end
    end
end

#______________________________________________________
# Main Solver
#______________________________________________________

"""
    aqueous_concentration_solver(A, b, mesh, cache, materials, 
                                 P_boundary_aq, mesh_partial_pressure_bc,
                                 params, θw_new, θw_old, vs, C_old, 
                                 kg, Δt, gas_idx) → C_new

Solve one time step of aqueous concentration transport for a single gas species.

**Workflow**:
1. Assemble global system (implicit FEM with Backward Euler)
2. Compute Dirichlet BC values using Henry's Law
3. Apply Dirichlet BCs via row-zeroing
4. Solve sparse linear system A·C_new = b

**Henry's Law BC**:
```
C_aq[i] = P_partial[i] / K_H
```
where:
- P_partial[i] = mesh_partial_pressure_bc[node_i, gas_idx] [Pa]
- K_H = materials.gases[gas_name].henry_constant [Pa·m³/mol]

# Arguments
- `A::SparseMatrixCSC`: Global system matrix [N×N] (modified, pre-allocated)
- `b::Vector{Float64}`: Global RHS [N] (modified)
- `mesh::MeshData`: Mesh connectivity
- `cache::RichardsCache`: Precomputed shape functions
- `materials::MaterialData`: Material properties with gas data
- `P_boundary_aq::Vector{Int}`: [N] Dirichlet mask (0=BC, 1=interior)
- `mesh_partial_pressure_bc::Dict`: Partial pressure BC data
- `params::AqueousTransportParams`: Transport parameters (D_h)
- `θw_new::Vector{Float64}`: [N] water content at time n+1
- `θw_old::Vector{Float64}`: [N] water content at time n
- `vs::Matrix{Float64}`: [N×2] Darcy velocity
- `C_old::Vector{Float64}`: [N] concentration at time n
- `kg::Vector{Float64}`: [N] mass transfer coefficient
- `Δt::Float64`: Time step [s]
- `gas_idx::Int`: Index of gas species (1, 2, 3, ...)

# Returns
- `C_new::Vector{Float64}`: [N] concentration at time n+1

# Notes
- Phase 1: k_g = 0, C_eq = 0 (no reaction/equilibration)
- Assumes materials.gas_dictionary[gas_idx] exists
- BC mask is modified in-place (reset after call if needed)
"""
function aqueous_concentration_solver(
    A                       :: SparseMatrixCSC{Float64, Int},
    b                       :: Vector{Float64},
    mesh,  # MeshData
    cache,  # RichardsCache
    materials,  # MaterialData
    P_boundary_aq           :: Vector{Int},
    mesh_partial_pressure_bc :: Dict,
    params                  :: AqueousTransportParams,
    θw_new                  :: Vector{Float64},
    θw_old                  :: Vector{Float64},
    vs                      :: Matrix{Float64},
    C_old                   :: Vector{Float64},
    kg                      :: Vector{Float64},
    Δt                      :: Float64,
    gas_idx                 :: Int
) :: Vector{Float64}
    
    # ──── Step 1: Assemble global system ────
    assemble_aqueous_concentration!(A, b, mesh, cache, params,
                                   θw_new, θw_old, vs, C_old, kg,
                                   zeros(length(C_old)), Δt)  # Ceq=0 in Phase 1
    
    # ──── Step 2: Compute Henry's Law BC values ────
    C_prescribed = zeros(Float64, length(b))
    gas_name = materials.gas_dictionary[gas_idx]
    K_H = materials.gases[gas_name].henry_constant
    
    # Apply Henry's Law at BC nodes (if K_H defined)
    if K_H > 0.0
        for (node_id, partial_pressures) in mesh_partial_pressure_bc
            if gas_idx <= length(partial_pressures)
                P_partial = partial_pressures[gas_idx]
                C_prescribed[node_id] = P_partial / K_H
                P_boundary_aq[node_id] = 0  # Mark as BC node
            end
        end
    end
    
    # ──── Step 3: Apply Dirichlet BCs ────
    apply_aqueous_concentration_bc!(A, b, P_boundary_aq, C_prescribed)
    
    # ──── Step 4: Solve sparse linear system ────
    C_new = A \ b
    
    return C_new
end

#______________________________________________________
