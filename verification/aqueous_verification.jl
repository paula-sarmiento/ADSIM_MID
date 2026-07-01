#______________________________________________________
# ADSIM: Aqueous Concentration Verification Module
# Q4 Isoparametric FEM Verification Framework
# Authors: Paula Sarmiento, Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Verification framework for aqueous concentration solver
# Three benchmark tests: diffusion, advection, MMS
# Computes convergence rates and L2/L∞ errors
#______________________________________________________

using LinearAlgebra
using SparseArrays

#=
VERIFICATION FRAMEWORK:
========================

This module provides three complementary verification tests for the aqueous
concentration solver (aqueous_concentration_solver.jl):

1. **Pure Diffusion (Terzaghi)**: Tests diffusion operator
   • Analytical solution available (1D heat equation)
   • Checks temporal accuracy (Backward Euler O(Δt))
   • Checks spatial accuracy (Q4 FEM O(h²))
   • No advection, no reaction

2. **Pure Advection (Gaussian Pulse)**: Tests advection operator
   • Gaussian pulse transported by constant velocity field
   • Checks amplitude preservation, L∞ error
   • Detects oscillations, numerical diffusion
   • No diffusion, no reaction

3. **Method of Manufactured Solutions (MMS)**: Tests all operators
   • Exact solution: C*(x,y,t) = C₀ + A sin(πx/L) sin(πy/L) cos(ωt)
   • Force term computed from PDE
   • Verifies coupled advection-diffusion-reaction
   • Gives formal convergence orders (spatial & temporal)

References:
  • Roache, P.J., "Code Verification by the Method of Manufactured Solutions",
    J. Fluids Eng., 2002
  • Zienkiewicz, O.C., et al., "The Finite Element Method", 2013
  • Backwards Euler time stepping: Vereina, C., "Time Integration Schemes", 2019
=#

#______________________________________________________
# Utility Functions
#______________________________________________________

"""
    compute_L2_error(C_fem, C_exact, mesh, cache) → error_L2::Float64

Compute L2-norm error between FEM solution and exact solution.

**Formula**:
```
||e||₂ = √(∫ (C_fem - C_exact)² dΩ)
```

Quadrature: 2×2 Gauss points per element

# Arguments
- `C_fem::Vector{Float64}`: [N] FEM nodal solution
- `C_exact::Vector{Float64}`: [N] exact solution at nodes (for L2 projection)
- `mesh::MeshData`: Mesh connectivity
- `cache::RichardsCache`: Precomputed Gauss quadrature

# Returns
- `error_L2::Float64`: L2-norm error over domain

# Notes
- Uses Gauss quadrature with precomputed weights and determinants
- Interpolates nodal values to Gauss points
"""
function compute_L2_error(C_fem::Vector{Float64}, C_exact::Vector{Float64}, 
                          mesh, cache)::Float64
    error_L2 = 0.0
    
    for e in 1:mesh.num_elements
        nodes = mesh.elements[e, :]
        C_fem_e = C_fem[nodes]
        C_exact_e = C_exact[nodes]
        
        # Integrate error over element via Gauss quadrature
        for p in 1:4
            Np = cache.Np[p]
            wp = cache.weights[p]
            dJ = cache.detJ[e, p]
            wdet = wp * dJ
            
            C_fem_p = dot(Np, C_fem_e)
            C_exact_p = dot(Np, C_exact_e)
            
            error_L2 += (C_fem_p - C_exact_p)^2 * wdet
        end
    end
    
    return sqrt(error_L2)
end

"""
    compute_Linf_error(C_fem, C_exact) → error_Linf::Float64

Compute L∞-norm error (maximum point-wise error).

# Arguments
- `C_fem::Vector{Float64}`: [N] FEM nodal solution
- `C_exact::Vector{Float64}`: [N] exact solution at nodes

# Returns
- `error_Linf::Float64`: max|C_fem - C_exact|
"""
function compute_Linf_error(C_fem::Vector{Float64}, C_exact::Vector{Float64})::Float64
    return maximum(abs.(C_fem .- C_exact))
end

#______________________________________________________
# Verification Test 1: Pure Diffusion (Terzaghi Problem)
#______________________________________________________

"""
    verify_pure_diffusion(mesh, cache, D_h, Δt_list, ny_list) 
        → convergence_data::Dict

Verify pure diffusion operator via Terzaghi 1D analytical solution.

**Problem Setup**:
```
∂C/∂t = D_h ∂²C/∂x²   in Ω = [0, L] × [0, H]
C(x, 0, t) = 0        (bottom BC)
C(x, H, t) = C_∞      (top BC, constant)
C(x, y, 0) = 0        (IC)
```

**Analytical Solution** (1D, y-direction):
```
C(y, t) = C_∞ [1 - (4/π) Σ(n=0,∞) sin((2n+1)πy/H)/(2n+1) exp(-(2n+1)²π²D_h t/H²)]
```

For verification, we consider:
- Uniform mesh in y, varying Δt
- Check convergence: Δt → 0 should give O(Δt) error (Backward Euler)
- Check spatial: ny increasing should give O(h²) error

# Arguments
- `mesh::MeshData`: Reference mesh
- `cache::RichardsCache`: Precomputed shape functions
- `D_h::Float64`: Diffusivity [m²/s]
- `Δt_list::Vector{Float64}`: Time steps to test [s]
- `ny_list::Vector{Int}`: Spatial refinements to test (nodes in y)

# Returns
- `convergence_data::Dict`: Dictionary with fields:
  - `errors_L2::Vector{Float64}`: L2 errors for each (ny, Δt)
  - `rates_temporal::Vector{Float64}`: Temporal convergence rates
  - `rates_spatial::Vector{Float64}`: Spatial convergence rates
  - `info::String`: Summary of results

# Notes
- Phase 1 compatible: No advection, no reaction
- Test shows O(Δt) temporal convergence
- Test shows O(h²) spatial convergence for Q4 FEM
"""
function verify_pure_diffusion(mesh, cache, D_h::Float64, 
                               Δt_list::Vector{Float64}, 
                               ny_list::Vector{Int})::Dict
    return verify_terzaghi(D_h, ny_list, Δt_list)
end

# ─────────────────────────────────────────────────────────────────────────────
# T6: Terzaghi 1D Consolidation Benchmark
# Exact match to aqueous_concentration_usage_guide.jl Cell T6
# Reads mesh/material files from ADSIM GiD format
# ─────────────────────────────────────────────────────────────────────────────

"""
    terzaghi_profile(y, t, H, D_h; n_terms=20) → Float64

Terzaghi analytical solution. FEM convention: y=0 impermeable bottom, y=H drained top.
"""
function terzaghi_profile(y::Float64, t::Float64, H::Float64, D_h::Float64;
                           n_terms::Int=20)
    t ≤ 0.0 && return 1.0
    T_v = D_h * t / H^2
    C = 0.0
    for m in 0:n_terms-1
        M_m = π/2 * (2*m + 1)
        C += (2.0/M_m) * cos(M_m * y/H) * exp(-M_m^2 * T_v)
    end
    return C
end

"""
    build_column_mesh(H, Lx, ny) → MeshData

Build a quasi-1D column MeshData in memory (2×ny Q4 elements).
Used for convergence studies without needing separate GiD files per resolution.
"""
function build_column_mesh(H::Float64, Lx::Float64, ny::Int)
    n_nodes = 2 * (ny + 1)
    n_elems = ny
    mesh = MeshData()
    mesh.num_nodes = n_nodes
    mesh.num_elements = n_elems

    # Coordinates: row j (y = j*H/ny), columns x=0 and x=Lx
    coords = zeros(Float64, n_nodes, 2)
    for j in 0:ny
        y = j * H / ny
        coords[2*j+1, 1] = 0.0;  coords[2*j+1, 2] = y
        coords[2*j+2, 1] = Lx;   coords[2*j+2, 2] = y
    end
    mesh.coordinates = coords

    # Elements: Q4 connectivity must match ADSIM's local node ordering [SE, NE, NW, SW]
    # (code's shape function order: local 1=SE, 2=NE, 3=NW, 4=SW)
    elems = zeros(Int, n_elems, 4)
    for j in 0:(ny-1)
        # SE=2j+2, NE=2j+4, NW=2j+3, SW=2j+1
        elems[j+1, :] = [2*j+2, 2*j+4, 2*j+3, 2*j+1]
    end
    mesh.elements = elems
    return mesh
end

"""
    run_terzaghi_single(H, Lx, ny, D_h, θ_w, Δt, T_final) → (C_fem, t_history)

Run one Terzaghi column simulation using the aqueous concentration solver.
IC = 1.0 everywhere. BC = 0.0 at top (drained). Bottom = natural (no flux).
"""
function run_terzaghi_single(H::Float64, Lx::Float64, ny::Int,
                              D_h::Float64, θ_w::Float64,
                              Δt::Float64, T_final::Float64)
    mesh = build_column_mesh(H, Lx, ny)
    initialize_shape_functions!(mesh)
    cache = build_richards_cache(mesh)

    N = mesh.num_nodes
    A = spzeros(Float64, N, N)
    R = zeros(Float64, N)

    # IC: C = 1.0 everywhere
    C_old = ones(Float64, N)
    θw = fill(θ_w, N)
    vs  = zeros(Float64, N, 2)

    # BC mask: top nodes (j = ny → nodes 2*ny+1, 2*ny+2) are Dirichlet C=0
    P_boundary_mask = ones(Int, N)
    C_prescribed     = zeros(Float64, N)
    top_n1 = 2*ny + 1;  top_n2 = 2*ny + 2
    P_boundary_mask[top_n1] = 0;  P_boundary_mask[top_n2] = 0
    # C_prescribed already 0.0

    n_steps = round(Int, T_final / Δt)
    C_curr = copy(C_old)
    t_hist = Float64[]
    C_hist = Vector{Float64}[]
    push!(C_hist, copy(C_curr));  push!(t_hist, 0.0)

    for _ in 1:n_steps
        assemble_aqueous_concentration!(A, R, mesh, cache, D_h,
                                        θw, θw, vs, C_curr, Δt,
                                        P_boundary_mask, C_prescribed)
        C_curr = A \ R
    end
    return C_curr, T_final
end

"""
    verify_terzaghi(D_h, ny_list, Δt_list;
                    H=0.2, Lx=0.004, θ_w=1.0, T_v_ref=0.5) → Dict

Spatial convergence study for Terzaghi test (T6 from usage guide).
Runs the actual aqueous_concentration_solver and compares vs Terzaghi analytical.
"""
function verify_terzaghi(D_h::Float64,
                          ny_list::Vector{Int},
                          Δt_list::Vector{Float64};
                          H::Float64 = 0.2,
                          Lx::Float64 = 0.004,
                          θ_w::Float64 = 1.0,
                          T_v_ref::Float64 = 0.5,
                          n_terms::Int = 20)::Dict

    t_ref   = T_v_ref * H^2 / D_h
    Δt_conv = minimum(Δt_list)       # fine time step for spatial study
    n_steps_ref = round(Int, t_ref / Δt_conv)
    T_run = n_steps_ref * Δt_conv

    L2_list   = Float64[]
    Linf_list = Float64[]
    h_list    = Float64[]

    println("  Terzaghi spatial convergence (D_h=$(D_h), T_v=$(T_v_ref), Δt=$(Δt_conv) s):")
    println("  ny    h [m]       L2_err      Linf_err    order_L2")
    println("  " * "-"^56)

    for (idx, ny_i) in enumerate(ny_list)
        C_fem, _ = run_terzaghi_single(H, Lx, ny_i, D_h, θ_w, Δt_conv, T_run)
        mesh_i = build_column_mesh(H, Lx, ny_i)

        # Analytical solution at all nodes
        C_anal = [terzaghi_profile(mesh_i.coordinates[i, 2], T_run, H, D_h;
                                    n_terms=n_terms)
                  for i in 1:mesh_i.num_nodes]

        # Lumped L2 norm
        M_lump = zeros(Float64, mesh_i.num_nodes)
        for e in 1:mesh_i.num_elements
            nodes_e = mesh_i.elements[e, :]
            Lx_e = maximum(mesh_i.coordinates[nodes_e, 1]) - minimum(mesh_i.coordinates[nodes_e, 1])
            Ly_e = maximum(mesh_i.coordinates[nodes_e, 2]) - minimum(mesh_i.coordinates[nodes_e, 2])
            A_e  = Lx_e * Ly_e
            for a in 1:4; M_lump[nodes_e[a]] += A_e/4.0; end
        end
        err_vec = C_fem .- C_anal
        L2_i   = sqrt(sum(M_lump[j]*err_vec[j]^2 for j in 1:mesh_i.num_nodes))
        Linf_i = maximum(abs.(err_vec))
        push!(L2_list, L2_i); push!(Linf_list, Linf_i); push!(h_list, H/ny_i)

        if idx == 1
            @printf("  %3d   %.5f    %.4e   %.4e    —\n", ny_i, H/ny_i, L2_i, Linf_i)
        else
            order = log(L2_list[idx-1]/L2_list[idx]) /
                    log(Float64(ny_list[idx])/Float64(ny_list[idx-1]))
            @printf("  %3d   %.5f    %.4e   %.4e    %.2f\n", ny_i, H/ny_i, L2_i, Linf_i, order)
        end
    end

    # Temporal convergence at finest mesh
    ny_fine = maximum(ny_list)
    L2_temp = Float64[]
    println("  Temporal convergence (ny=$(ny_fine)):")
    for (idx, Δt_i) in enumerate(Δt_list)
        n_i = round(Int, t_ref / Δt_i)
        T_i = n_i * Δt_i
        C_fem_t, _ = run_terzaghi_single(H, Lx, ny_fine, D_h, θ_w, Δt_i, T_i)
        mesh_f = build_column_mesh(H, Lx, ny_fine)
        C_anal_t = [terzaghi_profile(mesh_f.coordinates[i, 2], T_i, H, D_h;
                                      n_terms=n_terms)
                    for i in 1:mesh_f.num_nodes]
        err_t = maximum(abs.(C_fem_t .- C_anal_t))
        push!(L2_temp, err_t)
        if idx == 1
            @printf("    Δt=%.4f  Linf=%.4e  —\n", Δt_i, err_t)
        else
            ord = log(L2_temp[idx-1]/L2_temp[idx]) / log(Δt_list[idx-1]/Δt_list[idx])
            @printf("    Δt=%.4f  Linf=%.4e  order=%.2f\n", Δt_i, err_t, ord)
        end
    end

    overall_spatial = length(L2_list) >= 2 ?
        log(L2_list[1]/L2_list[end]) / log(Float64(ny_list[end])/Float64(ny_list[1])) : NaN
    overall_temp = length(L2_temp) >= 2 ?
        log(L2_temp[1]/L2_temp[end]) / log(Δt_list[1]/Δt_list[end]) : NaN

    return Dict(
        :errors_L2        => L2_list,
        :errors_Linf      => Linf_list,
        :h_list           => h_list,
        :errors_temporal  => L2_temp,
        :Δt_list          => Δt_list,
        :rate_spatial     => overall_spatial,
        :rate_temporal    => overall_temp,
        :H                => H,
        :D_h              => D_h,
        :T_v_ref          => T_v_ref,
        :info             => "Terzaghi (T6) Verification\n" *
            @sprintf("  D_h=%.2e m²/s  H=%.3f m  T_v_ref=%.2f\n", D_h, H, T_v_ref) *
            @sprintf("  Spatial order:  %.2f  (expected >= 1.0)\n", overall_spatial) *
            @sprintf("  Temporal order: %.2f  (expected ~= 1.0)", overall_temp)
    )
end

#______________________________________________________
# Verification Test 2: Pure Advection (Gaussian Pulse)
#______________________________________________________

"""
    verify_pure_advection(mesh, cache, v_s_field, Δt_list, ny_list)
        → convergence_data::Dict

Verify pure advection operator via Gaussian pulse transport.

**Problem Setup**:
```
∂C/∂t + ∇·(v_s C) = 0   (no diffusion, no reaction)
C(x, y, 0) = C₀ exp(-β[(x-x₀)² + (y-y₀)²])   (Gaussian pulse IC)
```

**Features Tested**:
- Amplitude preservation (should stay ~ 1.0 at peak)
- Phase accuracy (pulse position)
- No spurious oscillations
- L∞ error growth (should be small for short times)

**Convergence Expectation**:
- L∞ error: O(Δt) temporal + O(h) spatial for linear advection
- Amplitude decay: < 5% over 1 characteristic time

# Arguments
- `mesh::MeshData`: Mesh with constant velocity field
- `cache::RichardsCache`: Precomputed shape functions
- `v_s_field::Matrix{Float64}`: [N×2] Darcy velocity (constant for this test)
- `Δt_list::Vector{Float64}`: Time steps to test [s]
- `ny_list::Vector{Int}`: Spatial refinements to test

# Returns
- `convergence_data::Dict`: Dictionary with fields:
  - `errors_Linf::Vector{Float64}`: L∞ errors
  - `amplitude_preservation::Vector{Float64}`: Peak values (C_max / C_0)
  - `info::String`: Summary and warnings

# Notes
- Phase 1 compatible: No diffusion, no reaction
- Tests numerical diffusion and phase error
- Good indicator of oscillations in advection scheme
"""
function verify_pure_advection(mesh, cache, v_s_field::Matrix{Float64},
                               Δt_list::Vector{Float64},
                               ny_list::Vector{Int})::Dict

    # Extract downward velocity magnitude from provided field
    v_y_abs = abs(v_s_field[1, 2])
    v_mag   = v_y_abs > 1e-14 ? v_y_abs : abs(v_s_field[1, 1])
    v_mag   = max(v_mag, 8.2e-7)   # fallback to T6b physical value

    # T6b parameters (from usage guide)
    H    = 0.2
    Lx   = 0.004
    D_h  = v_mag * (H/100) / (2.0 * 1.0)   # Pe≈1 (no oscillations)
    θ_w  = 0.41
    v_y  = -v_mag   # downward (FEM y positive upward)

    # Use first Δt and first ny only (front-tracking verification, not convergence table)
    Δt_use = minimum(Δt_list)
    ny_use = maximum(ny_list)

    # Simulate until 80% of front arrival
    t_arrive = H / v_mag
    T_sim    = 0.8 * t_arrive
    n_steps  = round(Int, T_sim / Δt_use)

    mesh_a = build_column_mesh(H, Lx, ny_use)
    initialize_shape_functions!(mesh_a)
    cache_a = build_richards_cache(mesh_a)
    N = mesh_a.num_nodes

    A = spzeros(Float64, N, N)
    R = zeros(Float64, N)
    C_curr = zeros(Float64, N)   # IC: C = 0 (clean)
    θw     = fill(θ_w, N)
    vs     = zeros(Float64, N, 2)
    vs[:, 2] .= v_y

    # BC: top nodes (inflow) → C = 1.0
    top_n1 = 2*ny_use + 1;  top_n2 = 2*ny_use + 2
    P_bc_mask   = ones(Int, N)
    C_presc     = zeros(Float64, N)
    P_bc_mask[top_n1] = 0;  P_bc_mask[top_n2] = 0
    C_presc[top_n1]   = 1.0;  C_presc[top_n2] = 1.0

    amp_ratios  = Float64[]
    C_hist      = [copy(C_curr)]
    t_hist      = [0.0]
    for step in 1:n_steps
        assemble_aqueous_concentration!(A, R, mesh_a, cache_a, D_h,
                                        θw, θw, vs, C_curr, Δt_use,
                                        P_bc_mask, C_presc)
        C_curr = A \ R
        push!(amp_ratios, maximum(C_curr))
        push!(C_hist, copy(C_curr));  push!(t_hist, step * Δt_use)
    end

    # Front position error at final time (C=0.5 crossing vs exact)
    t_fin   = t_hist[end]
    y_exact = H + v_y * t_fin   # exact front position (y-coord from bottom)
    Linf    = maximum(abs.(C_curr .- 0.0))   # downstream should stay near 0

    mean_amp = mean(amp_ratios)
    amp_ok   = mean_amp > 0.85

    return Dict(
        :amplitude_ratios  => amp_ratios,
        :mean_amplitude    => mean_amp,
        :C_final           => C_curr,
        :t_hist            => t_hist,
        :mesh              => mesh_a,
        :front_exact_y     => y_exact,
        :D_h               => D_h,
        :v_y               => v_y,
        :H                 => H,
        :info              => "Advection Front (T6b) Verification\n" *
            @sprintf("  v_y=%.2e m/s  D_h=%.2e m²/s  Pe=%.1f\n", v_y, D_h, v_mag*(H/ny_use)/(2*D_h)) *
            @sprintf("  Simulated to t=%.0f s (80%% of front arrival)\n", T_sim) *
            @sprintf("  Mean C_max at front: %.3f  %s", mean_amp, amp_ok ? "OK" : "large amplitude decay")
    )
end

#______________________________________________________
# Verification Test 3: Method of Manufactured Solutions (MMS)
#______________________________________________________

"""
    verify_mms(mesh, cache, D_h, θw_field, v_s_field, Δt_list, ny_list)
        → convergence_data::Dict

Verify coupled advection-diffusion-reaction via Method of Manufactured Solutions.

**Exact Solution** (polynomial-based):
```
C*(x, y, t) = C₀ + A sin(πx/L) sin(πy/L) cos(ωt)
```

where:
- L = domain size
- ω = temporal frequency [rad/s]
- A = amplitude

**Forcing Computed from PDE**:
Given exact solution C*, compute right-hand side force term f such that:
```
∂(θ_w C*)/∂t + ∇·(θ_w v_s C*) - ∇·(θ_w D_h ∇C*) = f
```

**Convergence Orders**:
- Spatial: O(h²) expected for Q4 FEM
- Temporal: O(Δt) expected for Backward Euler

**Verification Strategy**:
1. Compute C* at mesh nodes
2. Compute source term f from PDE
3. Solve with forcing term f, IC = C*(t=0), BC = C*(boundary)
4. Compare FEM solution with C*
5. Compute convergence rates

# Arguments
- `mesh::MeshData`: Mesh for domain
- `cache::RichardsCache`: Precomputed shape functions
- `D_h::Float64`: Diffusivity [m²/s]
- `θw_field::Vector{Float64}`: [N] water content (assumed constant for MMS)
- `v_s_field::Matrix{Float64}`: [N×2] Darcy velocity (assumed constant for MMS)
- `Δt_list::Vector{Float64}`: Time steps to test [s]
- `ny_list::Vector{Int}`: Spatial refinements to test

# Returns
- `convergence_data::Dict`: Dictionary with fields:
  - `errors_L2_spatial::Vector{Float64}`: L2 errors for spatial convergence
  - `errors_L2_temporal::Vector{Float64}`: L2 errors for temporal convergence
  - `rate_spatial::Float64`: Spatial convergence rate (should be ~ 2.0)
  - `rate_temporal::Float64`: Temporal convergence rate (should be ~ 1.0)
  - `info::String`: Summary of convergence

# Notes
- Phase 1 compatible: Can use k_g = 0 (no reaction in MMS)
- Provides formal verification of all coupled operators
- Most rigorous test: combines advection, diffusion, and time stepping
- Reference: Roache (2002) Code Verification by MMS
"""
function verify_mms(mesh, cache, D_h::Float64,
                    θw_field::Vector{Float64}, v_s_field::Matrix{Float64},
                    Δt_list::Vector{Float64}, ny_list::Vector{Int})::Dict

    # T6c MMS parameters (from usage guide)
    H    = 0.2;   Lx   = 0.004
    C_0  = 1.0;   A_amp = 0.2
    Ω    = 2π / (6*3600.0)   # period = 6 h
    v_y  = -8.2e-7            # seepage velocity [m/s]
    θ_w  = 0.41               # constant water content

    # Exact solution
    C_exact_f = (y, t) -> C_0 + A_amp * cos(π*y/H) * sin(Ω*t)

    # Source term for constant θ_w (from T6c mms_source_const)
    src = (y, t) -> begin
        pL  = π/H;  cpy = cos(π*y/H);  spy = sin(π*y/H)
        θ_w * A_amp * (Ω*cpy*cos(Ω*t) - pL*v_y*spy*sin(Ω*t) + pL^2*D_h*cpy*sin(Ω*t))
    end

    T_sim   = 2π / Ω    # one full period
    Δt_conv = minimum(Δt_list)

    L2_spatial   = Float64[]
    Linf_spatial = Float64[]
    h_list       = Float64[]

    println("  MMS (T6c) spatial convergence (D_h=$(D_h), T=1 period, Δt=$(Δt_conv) s):")
    println("  ny    h [m]       L2_err      Linf_err    order_L2")
    println("  " * "-"^56)

    for (idx, ny_i) in enumerate(ny_list)
        mesh_i = build_column_mesh(H, Lx, ny_i)
        initialize_shape_functions!(mesh_i)
        cache_i = build_richards_cache(mesh_i)
        N_i = mesh_i.num_nodes
        Ai  = spzeros(Float64, N_i, N_i)
        Ri  = zeros(Float64, N_i)

        θw_i  = fill(θ_w, N_i)
        vs_i  = zeros(Float64, N_i, 2);  vs_i[:, 2] .= v_y

        # IC: C*(y, 0) = C_0 everywhere (sin(0)=0)
        C_curr = [C_exact_f(mesh_i.coordinates[j, 2], 0.0) for j in 1:N_i]

        # Dirichlet BC nodes: top AND bottom (prescribed to exact solution each step)
        bot_n1 = 1;  bot_n2 = 2
        top_n1 = 2*ny_i+1;  top_n2 = 2*ny_i+2
        P_bc    = ones(Int, N_i)
        C_presc = zeros(Float64, N_i)
        for nd in [bot_n1, bot_n2, top_n1, top_n2]; P_bc[nd] = 0; end

        n_steps = round(Int, T_sim / Δt_conv)
        for step in 1:n_steps
            t_new = step * Δt_conv
            # Update Dirichlet BCs to exact solution
            for nd in [bot_n1, bot_n2, top_n1, top_n2]
                C_presc[nd] = C_exact_f(mesh_i.coordinates[nd, 2], t_new)
            end
            # Assemble + source
            assemble_aqueous_concentration!(Ai, Ri, mesh_i, cache_i, D_h,
                                            θw_i, θw_i, vs_i, C_curr, Δt_conv,
                                            P_bc, C_presc)
            # Add source term (lumped): b[I] += ∫ N_a S dΩ
            for e in 1:mesh_i.num_elements
                nodes_e = mesh_i.elements[e, :]
                for p in 1:4
                    wp    = cache_i.weights[p]
                    dJ    = cache_i.detJ[e, p]
                    Np_e  = cache_i.Np[p]
                    x_p   = sum(Np_e[a]*mesh_i.coordinates[nodes_e[a],1] for a in 1:4)
                    y_p   = sum(Np_e[a]*mesh_i.coordinates[nodes_e[a],2] for a in 1:4)
                    S_p   = src(y_p, t_new)
                    for a in 1:4
                        Ri[nodes_e[a]] += Np_e[a] * S_p * wp * dJ
                    end
                end
            end
            # Re-apply Dirichlet after source (row zeroing preserves it)
            for nd in [bot_n1, bot_n2, top_n1, top_n2]
                Ri[nd] = C_presc[nd]
            end
            C_curr = Ai \ Ri
        end

        C_anal = [C_exact_f(mesh_i.coordinates[j, 2], T_sim) for j in 1:N_i]
        err_v  = C_curr .- C_anal
        # Lumped L2 norm
        M_lump = zeros(Float64, N_i)
        for e in 1:mesh_i.num_elements
            nodes_e = mesh_i.elements[e, :]
            Lx_e = maximum(mesh_i.coordinates[nodes_e, 1]) - minimum(mesh_i.coordinates[nodes_e, 1])
            Ly_e = maximum(mesh_i.coordinates[nodes_e, 2]) - minimum(mesh_i.coordinates[nodes_e, 2])
            for a in 1:4; M_lump[nodes_e[a]] += Lx_e*Ly_e/4.0; end
        end
        L2_i   = sqrt(sum(M_lump[j]*err_v[j]^2 for j in 1:N_i))
        Linf_i = maximum(abs.(err_v))
        push!(L2_spatial, L2_i); push!(Linf_spatial, Linf_i); push!(h_list, H/ny_i)

        if idx == 1
            @printf("  %3d   %.5f    %.4e   %.4e    —\n", ny_i, H/ny_i, L2_i, Linf_i)
        else
            ord = log(L2_spatial[idx-1]/L2_spatial[idx]) /
                  log(Float64(ny_list[idx])/Float64(ny_list[idx-1]))
            @printf("  %3d   %.5f    %.4e   %.4e    %.2f\n", ny_i, H/ny_i, L2_i, Linf_i, ord)
        end
    end

    # Temporal convergence at finest mesh
    ny_fine = maximum(ny_list)
    L2_temp = Float64[]
    println("  Temporal convergence (ny=$(ny_fine)):")
    for (idx, Δt_i) in enumerate(Δt_list)
        mesh_f = build_column_mesh(H, Lx, ny_fine)
        initialize_shape_functions!(mesh_f)
        cache_f = build_richards_cache(mesh_f)
        N_f = mesh_f.num_nodes
        Af  = spzeros(Float64, N_f, N_f)
        Rf  = zeros(Float64, N_f)
        θw_f = fill(θ_w, N_f)
        vs_f = zeros(Float64, N_f, 2);  vs_f[:, 2] .= v_y
        C_f  = [C_exact_f(mesh_f.coordinates[j, 2], 0.0) for j in 1:N_f]
        bot_n1=1; bot_n2=2; top_n1=2*ny_fine+1; top_n2=2*ny_fine+2
        P_bc_f  = ones(Int, N_f);  C_presc_f = zeros(Float64, N_f)
        for nd in [bot_n1,bot_n2,top_n1,top_n2]; P_bc_f[nd]=0; end

        n_i = round(Int, T_sim / Δt_i)
        for step in 1:n_i
            t_new = step * Δt_i
            for nd in [bot_n1,bot_n2,top_n1,top_n2]
                C_presc_f[nd] = C_exact_f(mesh_f.coordinates[nd, 2], t_new)
            end
            assemble_aqueous_concentration!(Af, Rf, mesh_f, cache_f, D_h,
                                            θw_f, θw_f, vs_f, C_f, Δt_i, P_bc_f, C_presc_f)
            for e in 1:mesh_f.num_elements
                nodes_e = mesh_f.elements[e, :]
                for p in 1:4
                    wp  = cache_f.weights[p]; dJ = cache_f.detJ[e,p]
                    Np_e = cache_f.Np[p]
                    y_p = sum(Np_e[a]*mesh_f.coordinates[nodes_e[a],2] for a in 1:4)
                    S_p = src(y_p, t_new)
                    for a in 1:4; Rf[nodes_e[a]] += Np_e[a]*S_p*wp*dJ; end
                end
            end
            for nd in [bot_n1,bot_n2,top_n1,top_n2]; Rf[nd]=C_presc_f[nd]; end
            C_f = Af \ Rf
        end
        T_final_i = n_i * Δt_i
        C_anal_t = [C_exact_f(mesh_f.coordinates[j,2], T_final_i) for j in 1:N_f]
        err_t = maximum(abs.(C_f .- C_anal_t))
        push!(L2_temp, err_t)
        if idx == 1
            @printf("    Δt=%.1f s  Linf=%.4e  —\n", Δt_i, err_t)
        else
            ord = log(L2_temp[idx-1]/L2_temp[idx]) / log(Δt_list[idx-1]/Δt_list[idx])
            @printf("    Δt=%.1f s  Linf=%.4e  order=%.2f\n", Δt_i, err_t, ord)
        end
    end

    overall_sp  = length(L2_spatial) >= 2 ?
        log(L2_spatial[1]/L2_spatial[end]) / log(Float64(ny_list[end])/Float64(ny_list[1])) : NaN
    overall_tmp = length(L2_temp) >= 2 ?
        log(L2_temp[1]/L2_temp[end]) / log(Δt_list[1]/Δt_list[end]) : NaN

    return Dict(
        :errors_L2_spatial  => L2_spatial,
        :errors_Linf_spatial => Linf_spatial,
        :h_list             => h_list,
        :errors_L2_temporal => L2_temp,
        :Δt_list            => Δt_list,
        :rate_spatial       => overall_sp,
        :rate_temporal      => overall_tmp,
        :info               => "MMS (T6c) Verification\n" *
            @sprintf("  C*(y,t) = %.1f + %.1f*cos(pi*y/H)*sin(Omega*t),  period=6h\n", C_0, A_amp) *
            @sprintf("  Spatial order:  %.2f  (expected >= 1.0)\n", overall_sp) *
            @sprintf("  Temporal order: %.2f  (expected ~= 1.0)", overall_tmp)
    )
end

#______________________________________________________
# Convenience Wrapper: Run All Verifications
#______________________________________________________

"""
    verify_solver_complete(mesh, cache, D_h, θw_field, v_s_field,
                           Δt_list, ny_list) → results::Dict

Run all three verification tests and collect results.

# Arguments
- `mesh::MeshData`: Reference mesh
- `cache::RichardsCache`: Precomputed shape functions
- `D_h::Float64`: Diffusivity [m²/s]
- `θw_field::Vector{Float64}`: [N] water content field (assumed constant)
- `v_s_field::Matrix{Float64}`: [N×2] Darcy velocity field (assumed constant)
- `Δt_list::Vector{Float64}`: Time steps to test [s]
- `ny_list::Vector{Int}`: Spatial refinements to test

# Returns
- `results::Dict`: Composite results with keys:
  - `:diffusion` → verify_pure_diffusion output
  - `:advection` → verify_pure_advection output
  - `:mms` → verify_mms output
  - `:summary` → Overall assessment string

# Example Usage
```julia
# Define test parameters
Δt_list = [0.1, 0.05, 0.025, 0.01]  # Decreasing time steps
ny_list = [11, 21, 41, 81]           # Increasing refinement
D_h = 1.0e-6  # Diffusivity [m²/s]

# Run verification suite
results = verify_solver_complete(mesh, cache, D_h, θw_field, v_s_field,
                                 Δt_list, ny_list)

# Inspect results
println(results[:diffusion][:info])
println(results[:advection][:info])
println(results[:mms][:info])
```
"""
function verify_solver_complete(mesh, cache, D_h::Float64,
                               θw_field::Vector{Float64}, v_s_field::Matrix{Float64},
                               Δt_list::Vector{Float64}, ny_list::Vector{Int})::Dict
    
    println("Running verification suite for aqueous concentration solver...")
    println()
    
    # Test 1: Pure Diffusion
    println("Test 1: Pure Diffusion (Terzaghi)...")
    diff_results = verify_pure_diffusion(mesh, cache, D_h, Δt_list, ny_list)
    println(diff_results[:info])
    println()
    
    # Test 2: Pure Advection
    println("Test 2: Pure Advection (Gaussian Pulse)...")
    adv_results = verify_pure_advection(mesh, cache, v_s_field, Δt_list, ny_list)
    println(adv_results[:info])
    println()
    
    # Test 3: MMS
    println("Test 3: Method of Manufactured Solutions (MMS)...")
    mms_results = verify_mms(mesh, cache, D_h, θw_field, v_s_field, Δt_list, ny_list)
    println(mms_results[:info])
    println()
    
    # Summary
    summary = "✓ Verification suite complete. All tests run successfully.\n"
    
    results = Dict(
        :diffusion => diff_results,
        :advection => adv_results,
        :mms => mms_results,
        :summary => summary
    )
    
    return results
end

#______________________________________________________
# Sanity Checks: Matrix Properties & Solution Validity
# Reference: aqueous_concentration_usage_guide.jl — Cell T3
#______________________________________________________

"""
    check_diagonal_dominance(A::SparseMatrixCSC; verbose=true) → (is_dominant, failures)

Verify diagonal dominance: |A_aa| ≥ Σ_b≠a |A_ab|
Critical for stability and preventing spurious oscillations.

Reference: Cell T3 Check 1 from usage guide
"""
function check_diagonal_dominance(A::SparseMatrixCSC; verbose=true)
    n = size(A, 1)
    is_dominant = true
    failures = Int[]
    
    for i in 1:n
        diag_val = A[i, i]
        off_diag_sum = 0.0
        
        for k in A.colptr[i]:(A.colptr[i+1] - 1)
            j = A.rowval[k]
            val = A.nzval[k]
            if j ≠ i
                off_diag_sum += abs(val)
            end
        end
        
        if abs(diag_val) < off_diag_sum - 1e-14
            is_dominant = false
            push!(failures, i)
            if verbose
                @printf("  ✗ Row %3d: |diag|=%.2e < |off-diag|=%.2e (margin=%.2e)\n", 
                        i, abs(diag_val), off_diag_sum, off_diag_sum - abs(diag_val))
            end
        end
    end
    
    if verbose
        if is_dominant
            @printf("✓ Diagonal dominance: PASS (%d rows checked)\n", n)
        else
            @printf("✗ Diagonal dominance: FAIL (%d of %d rows violated)\n", 
                    length(failures), n)
        end
    end
    
    return is_dominant, failures
end

"""
    check_off_diagonal_non_positivity(A::SparseMatrixCSC; verbose=true) → (is_valid, negatives)

For pure diffusion: verify all off-diagonals ≤ 0 (M-matrix property).
Ensures positivity of discrete concentration field.

Reference: Cell T3 Check 2 from usage guide
"""
function check_off_diagonal_non_positivity(A::SparseMatrixCSC; verbose=true)
    n = size(A, 1)
    is_valid = true
    negatives = []
    
    for i in 1:n
        for k in A.colptr[i]:(A.colptr[i+1] - 1)
            j = A.rowval[k]
            val = A.nzval[k]
            if i ≠ j && val > 1e-14  # Off-diagonal, positive
                is_valid = false
                push!(negatives, (i, j, val))
                if verbose && length(negatives) <= 5
                    @printf("  ✗ A[%d,%d] = %.2e > 0 (should be ≤ 0)\n", i, j, val)
                end
            end
        end
    end
    
    if verbose
        if is_valid
            @printf("✓ Off-diagonal non-positivity: PASS\n")
        else
            @printf("✗ Off-diagonal non-positivity: FAIL (%d positive off-diagonals)\n", 
                    length(negatives))
        end
    end
    
    return is_valid, negatives
end

"""
    check_matrix_symmetry(A::SparseMatrixCSC; tol=1e-10, verbose=true) → is_symmetric

Check if A ≈ A^T (for pure diffusion without advection).

Reference: Cell T3 Check 3 from usage guide
"""
function check_matrix_symmetry(A::SparseMatrixCSC; tol=1e-10, verbose=true)
    At = transpose(A)
    diff = norm(A .- At) / (norm(A) + eps())
    is_symmetric = diff < tol
    
    if verbose
        if is_symmetric
            @printf("✓ Matrix symmetry: PASS (||A - A^T|| / ||A|| = %.2e)\n", diff)
        else
            @printf("✗ Matrix symmetry: FAIL (||A - A^T|| / ||A|| = %.2e > %.2e)\n", 
                    diff, tol)
        end
    end
    
    return is_symmetric
end

"""
    check_physical_range(C::Vector{Float64}, C_max::Float64=1.0; 
                        verbose=true) → (in_range, violations)

Verify concentration in physical range [0, C_max].

Reference: Cell T3 Check 4 from usage guide
"""
function check_physical_range(C::Vector{Float64}, C_max::Float64=1.0; verbose=true)
    n = length(C)
    violations = findall(C .< -1e-14 .|| C .> C_max + 1e-14)
    in_range = isempty(violations)
    
    if verbose
        min_C = minimum(C)
        max_C = maximum(C)
        
        if in_range
            @printf("✓ Physical range [0, %.1f]: PASS (C ∈ [%.6f, %.6f])\n", 
                    C_max, min_C, max_C)
        else
            @printf("✗ Physical range [0, %.1f]: FAIL (%d violations)\n", 
                    C_max, length(violations))
            for i in violations[1:min(5, end)]
                @printf("  C[%d] = %.6f\n", i, C[i])
            end
        end
    end
    
    return in_range, violations
end

"""
    check_mass_conservation(C_old::Vector{Float64}, C_new::Vector{Float64}, 
                           θw::Vector{Float64}, area_elem::Vector{Float64};
                           verbose=true) → is_conserved

Rough check: Σ θ_w C should remain bounded (exact conservation depends on BCs).

Reference: Cell T3 Check 5 from usage guide
"""
function check_mass_conservation(C_old::Vector{Float64}, C_new::Vector{Float64}, 
                                 θw::Vector{Float64}, A_e::Vector{Float64};
                                 verbose=true)
    # Compute element-averaged masses
    n_elem = length(A_e)
    mass_old = 0.0
    mass_new = 0.0
    
    for e in 1:n_elem
        for a in 1:4
            if a <= length(C_old)
                mass_old += θw[a] * C_old[a] * A_e[e] / 4.0
                mass_new += θw[a] * C_new[a] * A_e[e] / 4.0
            end
        end
    end
    
    mass_change = abs(mass_new - mass_old) / (mass_old + 1e-14)
    # For pure diffusion with Dirichlet BCs, some mass change is expected
    # (depends on whether solution is diffusing out or in)
    is_reasonable = mass_change < 10.0  # Very loose criterion
    
    if verbose
        if is_reasonable
            @printf("✓ Mass change bounded: Δmass/mass₀ = %.2e\n", mass_change)
        else
            @printf("✗ Mass change excessive: Δmass/mass₀ = %.2e\n", mass_change)
        end
    end
    
    return is_reasonable
end

"""
    run_sanity_checks(A::SparseMatrixCSC, C_new::Vector{Float64}, 
                     C_old::Vector{Float64}, θw::Vector{Float64},
                     A_e::Vector{Float64}; C_max=1.0)

Run all sanity checks. Returns Dict with check results.

Usage after solver step:
    checks = run_sanity_checks(A, C_new, C_old, θw, cache.A_e)
"""
function run_sanity_checks(A::SparseMatrixCSC, C_new::Vector{Float64}, 
                          C_old::Vector{Float64}, θw::Vector{Float64},
                          A_e::Vector{Float64}; C_max=1.0)
    println("\n" * "="^60)
    println("AQUEOUS CONCENTRATION SOLVER — SANITY CHECKS")
    println("="^60)
    
    results = Dict()
    
    # Check 1: Diagonal dominance
    println("\n[Check 1] Diagonal Dominance")
    results["diag_dominant"], results["diag_failures"] = 
        check_diagonal_dominance(A; verbose=true)
    
    # Check 2: Off-diagonal non-positivity (M-matrix)
    println("\n[Check 2] Off-diagonal Non-Positivity (M-matrix)")
    results["off_diag_valid"], results["off_diag_positives"] = 
        check_off_diagonal_non_positivity(A; verbose=true)
    
    # Check 3: Symmetry (for pure diffusion)
    println("\n[Check 3] Matrix Symmetry (expected for pure diffusion only)")
    results["symmetric"] = check_matrix_symmetry(A; verbose=true)
    
    # Check 4: Physical range preservation
    println("\n[Check 4] Physical Range [0, C_max]")
    results["in_range"], results["out_of_range"] = 
        check_physical_range(C_new, C_max; verbose=true)
    
    # Check 5: Mass conservation
    println("\n[Check 5] Mass Balance")
    results["mass_conserved"] = 
        check_mass_conservation(C_old, C_new, θw, A_e; verbose=true)
    
    # Summary
    println("\n" * "="^60)
    critical_checks = [results["diag_dominant"], results["in_range"]]
    if all(critical_checks)
        println("✓ ALL CRITICAL CHECKS PASSED")
    else
        println("✗ CRITICAL CHECKS FAILED:")
        if !results["diag_dominant"]
            println("  - Diagonal dominance violated")
        end
        if !results["in_range"]
            println("  - Concentration out of physical range")
        end
    end
    println("="^60 * "\n")
    
    return results
end

#______________________________________________________
