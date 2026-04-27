#______________________________________________________
# ADSIM: Celia Verification
# Implicit FEM for the Richards equation — three SWRC models
# Following Celia et al. (1990), WRR 26(7):1483–1496
# Author: Paula Sarmiento — April 2026
#______________________________________________________

#______________________________________________________
# Mixed-form Richards equation (gravity ON):
#   ∂θ/∂t = ∇·(K(h)∇h) − ∇·(K(h)ê_g)
#
# Three SWRC models compared on identical domain:
#   VanGenuchten   — α=2.0, n=2.0, K_s=0.1
#   LinearSoil     — C=const, K∝Se, K_s=0.1
#   CavalcanteZornberg — δ=2.0, K_s=0.1
#
# Domain: 1 m × 1 m
# IC: h = −0.5 m (uniform)
# BC: h(y=0) = −1.0 m (dry bottom)
#     h(y=1) =  0.0 m (wet top, ponding)
#     sides: no-flow (natural Neumann = 0)
#______________________________________________________

using LinearAlgebra, SparseArrays, Printf, Plots, Statistics, TOML, Dates

const PROJECT_ROOT = @__DIR__
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")

include(joinpath(SRC_DIR, "swrc_models.jl"))
include(joinpath(SRC_DIR, "shape_functions.jl"))
using .ShapeFunctions: initialize_shape_functions!, build_richards_cache, RichardsCache

#______________________________________________________
# Minimal mesh struct
#______________________________________________________
struct MinMesh
    num_nodes    :: Int
    num_elements :: Int
    coordinates  :: Matrix{Float64}
    elements     :: Matrix{Int}
end

#______________________________________________________
# Mesh reader
#______________________________________________________
function read_adsim_mesh(filename)
    lines = readlines(filename)
    idx = 1
    while idx <= length(lines) && !startswith(strip(lines[idx]), "MESH"); idx += 1; end
    parts   = split(strip(lines[idx])); idx += 1
    N_nodes = parse(Int, parts[2])
    N_elem  = parse(Int, parts[3])

    coords   = zeros(Float64, N_nodes, 2)
    elements = zeros(Int,     N_elem,  4)
    pressure_head_bc      = Dict{Int,Float64}()
    initial_pressure_head = Dict{Int,Float64}()

    while idx <= length(lines)
        line = strip(lines[idx]); idx += 1
        if line == "coordinates"
            for i in 1:N_nodes
                p = split(strip(lines[idx])); idx += 1
                coords[i,1] = parse(Float64, p[1])
                coords[i,2] = parse(Float64, p[2])
            end
        elseif line == "elements"
            for e in 1:N_elem
                p = split(strip(lines[idx])); idx += 1
                for a in 1:4; elements[e,a] = parse(Int, p[a]); end
            end
        elseif line == "pressure_head_bc"
            n_bc = parse(Int, strip(lines[idx])); idx += 1
            for _ in 1:n_bc
                p = split(strip(lines[idx])); idx += 1
                pressure_head_bc[parse(Int,p[1])] = parse(Float64,p[2])
            end
        elseif line == "initial_pressure_head"
            n_ic = parse(Int, strip(lines[idx])); idx += 1
            for _ in 1:n_ic
                p = split(strip(lines[idx])); idx += 1
                initial_pressure_head[parse(Int,p[1])] = parse(Float64,p[2])
            end
        end
    end
    return (N_nodes=N_nodes, N_elem=N_elem, coordinates=coords, elements=elements,
            pressure_head_bc=pressure_head_bc, initial_pressure_head=initial_pressure_head)
end

#______________________________________________________
# SWRC model factory (reads from TOML soil section)
#______________________________________________________
function create_model(soil_toml, K_s::Float64)
    model_name = soil_toml["swrc_model"]
    theta_s    = Float64(soil_toml["porosity"])
    theta_r    = Float64(soil_toml["residual_water_content"])

    if model_name == "Van_Genuchten"
        return VanGenuchten(
            theta_s = theta_s,
            theta_r = theta_r,
            alpha   = Float64(soil_toml["swrc_vg_alpha"]),
            n_vg    = Float64(soil_toml["swrc_vg_n"]),
            K_s     = K_s
        )
    elseif model_name == "LinearSoil"
        return LinearSoil(
            theta_r = theta_r,
            theta_s = theta_s,
            h_min   = -Float64(soil_toml["swrc_vg_alpha"]),
            K_s     = K_s
        )
    elseif model_name == "CavalcanteZornberg"
        return CavalcanteZornberg(
            theta_s = theta_s,
            theta_r = theta_r,
            delta   = Float64(soil_toml["swrc_cav_delta"]),
            K_s     = K_s
        )
    else
        error("Unknown SWRC model in TOML: $model_name")
    end
end

#______________________________________________________
# FEM element matrices (gravity in residual only)
#______________________________________________________
function element_matrices_local!(
    Aᵉ      :: Matrix{Float64},
    Rᵉ      :: Vector{Float64},
    hᵉ_curr :: Vector{Float64},
    hᵉ_prev :: Vector{Float64},
    model   :: SWRCModel,
    Δt      :: Float64,
    e_g     :: Vector{Float64},
    cache   :: RichardsCache,
    e       :: Int
)
    fill!(Aᵉ, 0.0); fill!(Rᵉ, 0.0)

    for p in 1:4
        Np = cache.Np[p]
        w  = cache.weights[p] * cache.detJ[e, p]

        hp = Np[1]*hᵉ_curr[1] + Np[2]*hᵉ_curr[2] +
             Np[3]*hᵉ_curr[3] + Np[4]*hᵉ_curr[4]
        Kx = K_h_x(model, hp)
        Ky = K_h_y(model, hp)

        @inbounds for a in 1:4, b in 1:4
            Aᵉ[a,b] += (Kx*cache.Bp[e,p,1,a]*cache.Bp[e,p,1,b] +
                         Ky*cache.Bp[e,p,2,a]*cache.Bp[e,p,2,b]) * w
        end

        gx = sum(cache.Bp[e,p,1,a]*hᵉ_curr[a] for a in 1:4)
        gy = sum(cache.Bp[e,p,2,a]*hᵉ_curr[a] for a in 1:4)

        @inbounds for a in 1:4
            Rᵉ[a] += (-(cache.Bp[e,p,1,a]*Kx*gx + cache.Bp[e,p,2,a]*Ky*gy)
                      + (cache.Bp[e,p,1,a]*Kx*e_g[1] + cache.Bp[e,p,2,a]*Ky*e_g[2])) * w
        end
    end

    ML = cache.A_e[e] / 4.0
    @inbounds for a in 1:4
        Ca = C_moist(model, hᵉ_curr[a])
        Aᵉ[a,a] += Ca / Δt * ML
        Rᵉ[a]   -= ML / Δt * (theta(model, hᵉ_curr[a]) - theta(model, hᵉ_prev[a]))
    end
end

function build_sparsity(mesh::MinMesh, N::Int)
    I_vec, J_vec = Int[], Int[]
    for e in 1:mesh.num_elements
        nd = mesh.elements[e,:]
        for a in 1:4, b in 1:4; push!(I_vec, nd[a]); push!(J_vec, nd[b]); end
    end
    return sparse(I_vec, J_vec, zeros(length(I_vec)), N, N)
end

function assemble!(A, R, h_curr, h_prev, mesh::MinMesh, model, Δt, e_g, cache, P_bc)
    fill!(A.nzval, 0.0); fill!(R, 0.0)
    Aᵉ = zeros(4,4); Rᵉ = zeros(4)
    hc  = zeros(4);  hp  = zeros(4);  nd = zeros(Int,4)

    for e in 1:mesh.num_elements
        @inbounds for a in 1:4
            nd[a] = mesh.elements[e,a]; hc[a] = h_curr[nd[a]]; hp[a] = h_prev[nd[a]]
        end
        element_matrices_local!(Aᵉ, Rᵉ, hc, hp, model, Δt, e_g, cache, e)
        @inbounds for a in 1:4
            R[nd[a]] += Rᵉ[a]
            for b in 1:4; A[nd[a], nd[b]] += Aᵉ[a,b]; end
        end
    end

    for j in 1:mesh.num_nodes
        for k in A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[k]
            P_bc[i] == 0 && (A.nzval[k] = (i == j) ? 1.0 : 0.0)
        end
    end
    for i in 1:mesh.num_nodes; P_bc[i] == 0 && (R[i] = 0.0); end
end

function picard!(h_curr, h_prev, mesh::MinMesh, model, Δt, e_g, A, cache, P_bc;
                 tol=1e-8, max_iter=100, ω=1.0)
    R = zeros(mesh.num_nodes)
    h_curr .= h_prev
    for m in 1:max_iter
        assemble!(A, R, h_curr, h_prev, mesh, model, Δt, e_g, cache, P_bc)
        δ = A \ R
        h_curr .+= ω .* δ
        maximum(abs.(δ)) < tol && return m
    end
    @warn "Picard did not converge in $max_iter iterations"
    return max_iter
end

#______________________________________________________
# Mass balance utility
#______________________________________________________
function total_water(h::Vector{Float64}, mesh::MinMesh,
                     model::SWRCModel, cache::RichardsCache)
    W = 0.0
    for e in 1:mesh.num_elements
        ML = cache.A_e[e] / 4.0
        for a in 1:4
            W += ML * theta(model, h[mesh.elements[e,a]])
        end
    end
    return W
end

#______________________________________________________
# Run one full simulation (returns h_hist, t_hist, W_hist)
#______________________________________________________
function run_simulation(mesh, model, h_ic, P_bc, Δt, T, e_g, A, cache, log_print)
    N       = mesh.num_nodes
    n_steps = round(Int, T / Δt)
    h       = copy(h_ic)
    h_new   = copy(h)
    h_hist  = [copy(h)]
    t_hist  = [0.0]
    W_hist  = [total_water(h, mesh, model, cache)]

    log_print("-"^64)
    let t = 0.0
        for step in 1:n_steps
            n_iter = picard!(h_new, h, mesh, model, Δt, e_g, A, cache, P_bc)
            h .= h_new
            t += Δt
            log_print(@sprintf("   Step %3d | t = %.4f s | Picard = %d iters", step, t, n_iter))
            push!(h_hist, copy(h))
            push!(t_hist, t)
            push!(W_hist, total_water(h, mesh, model, cache))
        end
    end
    log_print("-"^64)
    return h_hist, t_hist, W_hist
end

#______________________________________________________
# Main execution
#______________________________________________________
function main()
    if length(ARGS) < 1
        println("Error: No project name provided")
        println("Usage: julia run_celia_verification.jl <project_name>")
        println("Example: julia run_celia_verification.jl Celia")
        exit(1)
    end

    project_name = ARGS[1]
    data_dir     = joinpath(SRC_DIR, "data")
    output_dir   = joinpath(PROJECT_ROOT, "output")
    mesh_file    = joinpath(data_dir, "$(project_name).mesh")
    mat_file     = joinpath(data_dir, "$(project_name)_mat.toml")
    calc_file    = joinpath(data_dir, "$(project_name)_calc.toml")

    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file"); exit(1)
    end

    isdir(output_dir) || mkdir(output_dir)

    log_file_path = joinpath(output_dir, "$(project_name)_verification.log")
    isfile(log_file_path) && rm(log_file_path)
    log_file = open(log_file_path, "w")

    function log_print(msg::String)
        println(msg); println(log_file, msg); flush(log_file)
    end

    try
        start_time = now()

        log_print("="^64)
        log_print("ADSIM: Celia Verification — Richards Equation (3 SWRC Models)")
        log_print("Following Celia et al. (1990), WRR 26(7):1483–1496")
        log_print("Project: $project_name")
        log_print("="^64)

        # ── [1/7] Read mesh ──────────────────────────────────────────
        log_print("\n[1/7] Reading mesh file: $mesh_file")
        msh  = read_adsim_mesh(mesh_file)
        mesh = MinMesh(msh.N_nodes, msh.N_elem, msh.coordinates, msh.elements)
        log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")

        # ── [2/7] Read material properties ──────────────────────────
        log_print("\n[2/7] Reading material properties: $mat_file")
        mat_data   = TOML.parsefile(mat_file)
        soil_names = mat_data["soil_dictionary_"]
        rho_w      = Float64(mat_data["liquid"]["density"])
        mu_w       = Float64(mat_data["liquid"]["dynamic_viscosity"])
        log_print("   ✓ Found $(length(soil_names)) soil models: $(join(soil_names, ", "))")

        # ── [3/7] Read calculation parameters ───────────────────────
        log_print("\n[3/7] Reading calculation parameters: $calc_file")
        calc_data = TOML.parsefile(calc_file)
        g_mag     = Float64(calc_data["gravity"]["gravity_magnitude"])
        e_g       = [Float64(calc_data["gravity"]["gravity_x_component"]),
                     Float64(calc_data["gravity"]["gravity_y_component"])]
        T_total   = Float64(calc_data["time_stepping"]["total_simulation_time"])
        dt        = Float64(calc_data["time_stepping"]["time_per_step"])
        n_steps   = round(Int, T_total / dt)

        log_print("   ✓ Gravity: e_g = [$(e_g[1]), $(e_g[2])]  " *
                  (all(e_g .== 0) ? "(OFF)" : "(ON — downward ✓)"))
        log_print(@sprintf("   ✓ T = %.2f s,  Δt = %.3f s,  %d steps", T_total, dt, n_steps))

        # ── [4/7] Initialize IC, BCs, shape functions ────────────────
        log_print("\n[4/7] Applying initial and boundary conditions")
        N    = mesh.num_nodes
        P_bc = ones(Int, N)
        for node_id in keys(msh.pressure_head_bc); P_bc[node_id] = 0; end

        h_ic_val = isempty(msh.initial_pressure_head) ? -0.5 :
                   first(values(msh.initial_pressure_head))
        h_ic = fill(h_ic_val, N)
        for (nid, hval) in msh.pressure_head_bc; h_ic[nid] = hval; end

        log_print(@sprintf("   ✓ IC: h = %.3f m (uniform)", h_ic_val))
        log_print("   ✓ Dirichlet BC: $(count(P_bc .== 0)) nodes")

        # ── [5/7] Shape functions & Richards cache ───────────────────
        log_print("\n[5/7] Initializing shape functions and Richards cache")
        initialize_shape_functions!(mesh)
        cache = build_richards_cache(mesh)
        A     = build_sparsity(mesh, N)
        log_print("   ✓ RichardsCache built: $(mesh.num_elements) elements × 4 Gauss points")
        log_print("   ✓ Sparsity pattern: $(nnz(A)) nonzeros")

        # ── [6/7] Run three simulations ──────────────────────────────
        log_print("\n[6/7] Running three Picard simulations (backward Euler)")

        h_histories = Vector{Vector{Vector{Float64}}}()
        t_histories = Vector{Vector{Float64}}()
        W_histories = Vector{Vector{Float64}}()
        models      = SWRCModel[]
        labels      = String[]

        for soil_name in soil_names
            soil_toml  = mat_data["soil"][soil_name]
            k_int      = Float64(soil_toml["intrinsic_permeability"])
            K_s        = k_int * rho_w * g_mag / mu_w
            model      = create_model(soil_toml, K_s)

            log_print("\n   ── $(soil_name)  (K_s = $(round(K_s, digits=4)) m/s) ──")

            h_hist, t_hist, W_hist = run_simulation(
                mesh, model, h_ic, copy(P_bc), dt, T_total, e_g, A, cache, log_print)

            push!(h_histories, h_hist)
            push!(t_histories, t_hist)
            push!(W_histories, W_hist)
            push!(models, model)
            push!(labels, soil_name)
        end

        log_print("\n   ✓ All three simulations complete")

        # ── [7/7] Post-processing + plots ────────────────────────────
        log_print("\n[7/7] Post-processing: profiles and mass balance")

        # Extract centre-column nodes (x ≈ Lx/2)
        Lx    = maximum(mesh.coordinates[:,1])
        x_mid = Lx / 2.0
        tol_x = Lx / (2 * 5)
        col_nodes = sort(findall(i -> abs(mesh.coordinates[i,1] - x_mid) < tol_x, 1:N),
                         by = i -> mesh.coordinates[i,2])
        y_col = mesh.coordinates[col_nodes, 2]
        log_print("   ✓ Extraction column: x ≈ $(round(x_mid,digits=2)) m  ($(length(col_nodes)) nodes)")

        colors_m = [:steelblue, :firebrick, :seagreen]
        colors_t = [:navy, :dodgerblue, :deepskyblue, :skyblue, :lightblue]
        times_snap = [0.0, 0.5, 1.0, 1.5, 2.0]

        # ── Per-model profile plots (h and θ) ─────────────────────────
        profile_panels = Plots.Plot[]
        for (mi, (h_hist, t_hist, model, lbl)) in enumerate(zip(h_histories, t_histories, models, labels))
            p_h = plot(xlabel="h [m]", ylabel="y [m]", title="$lbl — h(y)",
                       legend=:topright, grid=true, gridalpha=0.3)
            p_θ = plot(xlabel="θ [-]", ylabel="y [m]", title="$lbl — θ(y)",
                       legend=:topright, grid=true, gridalpha=0.3)

            for (ki, t_plot) in enumerate(times_snap)
                snap = argmin(abs.(t_hist .- t_plot))
                t_act = t_hist[snap]
                col   = colors_t[ki]
                lbl_t = @sprintf("t=%.1f s", t_act)
                h_col = h_hist[snap][col_nodes]
                θ_col = [theta(model, h_col[j]) for j in eachindex(h_col)]
                plot!(p_h, h_col, y_col, label=lbl_t, color=col, lw=2, marker=:circle, ms=3)
                plot!(p_θ, θ_col, y_col, label=lbl_t, color=col, lw=2, marker=:circle, ms=3)
            end
            push!(profile_panels, p_h, p_θ)
        end

        p_profiles = plot(profile_panels..., layout=(3,2), size=(900, 1100),
                          plot_title="Celia Verification — h and θ profiles (gravity ON)")
        out_prof = joinpath(output_dir, "$(project_name)_profiles.png")
        savefig(p_profiles, out_prof)
        log_print("   ✓ Profile plot saved: $out_prof")

        # ── Mass balance W(t) ─────────────────────────────────────────
        p_mass = plot(xlabel="t [s]", ylabel="W = ∫θ dΩ [m²]",
                      title="Total stored water W(t)",
                      legend=:topright, grid=true, gridalpha=0.3)
        for (mi, (W_hist, t_hist, lbl)) in enumerate(zip(W_histories, t_histories, labels))
            plot!(p_mass, t_hist, W_hist, label=lbl, color=colors_m[mi], lw=2, marker=:circle, ms=2)
        end

        p_dW = plot(xlabel="t [s]", ylabel="|ΔW| [m²]",
                    title="Change in water content per step",
                    legend=:topright, grid=true, gridalpha=0.3, yscale=:log10)
        for (mi, (W_hist, t_hist, lbl)) in enumerate(zip(W_histories, t_histories, labels))
            dW = [abs(W_hist[i+1] - W_hist[i]) for i in 1:length(W_hist)-1]
            plot!(p_dW, t_hist[2:end], dW, label=lbl, color=colors_m[mi], lw=2, marker=:circle, ms=2)
        end

        p_mb = plot(p_mass, p_dW, layout=(1,2), size=(1000, 400),
                    plot_title="Mass Balance Check")
        out_mb = joinpath(output_dir, "$(project_name)_mass_balance.png")
        savefig(p_mb, out_mb)
        log_print("   ✓ Mass balance plot saved: $out_mb")

        # ── Final comparison (all models at t=T) ─────────────────────
        p_hf = plot(xlabel="h [m]", ylabel="y [m]", title="Final h (t=$(T_total) s)",
                    legend=:topright, grid=true, gridalpha=0.3)
        p_θf = plot(xlabel="θ [-]", ylabel="y [m]", title="Final θ (t=$(T_total) s)",
                    legend=:topright, grid=true, gridalpha=0.3)

        for (mi, (h_hist, model, lbl)) in enumerate(zip(h_histories, models, labels))
            h_fin = h_hist[end][col_nodes]
            θ_fin = [theta(model, h_fin[j]) for j in eachindex(h_fin)]
            plot!(p_hf, h_fin, y_col, label=lbl, color=colors_m[mi], lw=2, marker=:circle, ms=4)
            plot!(p_θf, θ_fin, y_col, label=lbl, color=colors_m[mi], lw=2, marker=:circle, ms=4)
        end

        p_cmp = plot(p_hf, p_θf, layout=(1,2), size=(900, 450),
                     plot_title="Model comparison at t = $(T_total) s  (gravity ON, e_g=[$(e_g[1]),$(e_g[2])])")
        out_cmp = joinpath(output_dir, "$(project_name)_comparison.png")
        savefig(p_cmp, out_cmp)
        log_print("   ✓ Comparison plot saved: $out_cmp")

        # ── SWRC curves ───────────────────────────────────────────────
        h_range = range(-1.5, 0.0, length=300)
        p_θh = plot(xlabel="h [m]", ylabel="θ [-]", title="θ(h)",
                    legend=:bottomright, grid=true, gridalpha=0.3)
        p_Ch = plot(xlabel="h [m]", ylabel="C(h) [1/m]", title="C(h)",
                    legend=:topright, grid=true, gridalpha=0.3)
        p_Kh = plot(xlabel="h [m]", ylabel="K(h) [m/s]", title="K(h)",
                    legend=:bottomright, grid=true, gridalpha=0.3, yscale=:log10)

        for (mi, (model, lbl)) in enumerate(zip(models, labels))
            θv = [theta(model, Float64(h)) for h in h_range]
            Cv = [C_moist(model, Float64(h)) for h in h_range]
            Kv = [max(K_h_x(model, Float64(h)), 1e-12) for h in h_range]
            plot!(p_θh, h_range, θv, label=lbl, color=colors_m[mi], lw=2)
            plot!(p_Ch, h_range, Cv, label=lbl, color=colors_m[mi], lw=2)
            plot!(p_Kh, h_range, Kv, label=lbl, color=colors_m[mi], lw=2)
        end

        p_swrc = plot(p_θh, p_Ch, p_Kh, layout=(1,3), size=(1200, 380),
                      plot_title="SWRC Comparison")
        out_swrc = joinpath(output_dir, "$(project_name)_swrc.png")
        savefig(p_swrc, out_swrc)
        log_print("   ✓ SWRC plot saved: $out_swrc")

        # ── Summary ───────────────────────────────────────────────────
        log_print("\n" * "-"^64)
        log_print("   Final water content W(t=T) per model:")
        for (mi, (W_hist, lbl)) in enumerate(zip(W_histories, labels))
            ΔW = W_hist[end] - W_hist[1]
            log_print(@sprintf("   %-30s  W₀=%.4f  Wₜ=%.4f  ΔW=%.4f m²",
                               lbl, W_hist[1], W_hist[end], ΔW))
        end

        end_time   = now()
        total_time = (end_time - start_time).value / 1000.0
        log_print("\nTotal run time: $(total_time) s")
        log_print("\n" * "="^64)
        log_print("Celia verification completed successfully")
        log_print("="^64)

    catch e
        log_print("\n" * "="^64)
        log_print("FATAL ERROR OCCURRED")
        log_print("="^64)
        log_print("\nError Type: $(typeof(e))")
        log_print("\nError Message:")
        log_print(sprint(showerror, e))
        log_print("\nStack Trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(log_file, exc, bt); println(log_file); flush(log_file)
        end
        log_print("\n" * "="^64)
        log_print("Program terminated due to error")
        log_print("="^64)
        println("\nFATAL ERROR — check log: $log_file_path")
        close(log_file); exit(1)

    finally
        isopen(log_file) && close(log_file)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
