#!/usr/bin/env julia
"""
extract_benchmark_comparison.jl

Extract vertical water content profiles from all benchmark runs
and generate side-by-side comparison plots.

Compares:
  - Final profiles at t=24h
  - Time snapshots at t=0h, 6h, 12h, 18h, 24h (finest run)
  - Quantitative metrics (front penetration, saturation levels)

Run from: verification/scripts with julia extract_benchmark_comparison.jl
"""

using Printf
using Plots
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "swrc_models.jl"))

println("="^70)
println("  BENCHMARK PROFILE EXTRACTION & COMPARISON")
println("="^70)

# VanGenuchten model for Celia 1990 benchmark (Haverkamp sand, SI units)
# K_sat computed from intrinsic_permeability=9.4e-12 m², μ=1e-3 Pa·s, ρg=9810 N/m³
# K_sat = k·ρg/μ = 9.4e-12 * 9810 / 1e-3 ≈ 9.22e-5 m/s
vg_celia = VanGenuchten(
    theta_s=0.368, theta_r=0.102,
    alpha=3.35, n_vg=2.0,
    K_s=9.22e-5,
    K_s_x=9.22e-5, K_s_y=9.22e-5,
    L=0.5
)

# Verification benchmark configurations (verification-only cases)
benchmark_cases = [
    (prefix="CeliaCol_dt144",  label="Δt=2.4 min (144s)",  color=:black),
    (prefix="CeliaCol_dt720",  label="Δt=12 min (720s)",   color=:steelblue),
    (prefix="CeliaCol_dt3600", label="Δt=60 min (3600s)",  color=:firebrick)
]

output_dir = joinpath(@__DIR__, "..", "..", "output")
plots_dir = joinpath(@__DIR__, "..", "plots")
mkpath(plots_dir)

# ──────────────────────────────────────────────────────────────────────
# VTK parsing function
# ──────────────────────────────────────────────────────────────────────

function read_vtk_points_field(filename; field_name::String="Matric_Head")
    open(filename, "r") do f
        lines = readlines(f)

        # Find POINTS and parse all 3*n_pts coordinate values
        pt_idx = findfirst(l -> startswith(l, "POINTS"), lines)
        if pt_idx === nothing
            error("No POINTS found in $filename")
        end
        n_pts = parse(Int, split(lines[pt_idx])[2])

        coords_flat = Float64[]
        for i in (pt_idx + 1):length(lines)
            line = strip(lines[i])
            if isempty(line)
                continue
            end
            startswith(line, "CELLS") && break
            startswith(line, "POINT_DATA") && break

            for tok in split(line)
                try
                    push!(coords_flat, parse(Float64, tok))
                catch
                end
                length(coords_flat) >= 3 * n_pts && break
            end
            length(coords_flat) >= 3 * n_pts && break
        end

        if length(coords_flat) < 3 * n_pts
            error("Incomplete POINTS data in $filename: expected $(3*n_pts), got $(length(coords_flat))")
        end

        points = reshape(coords_flat[1:3*n_pts], 3, n_pts)'

        # Find POINT_DATA and requested scalar field
        pd_idx = findfirst(l -> startswith(l, "POINT_DATA"), lines)
        pd_idx === nothing && error("No POINT_DATA found in $filename")

        sc_idx = nothing
        for i in (pd_idx + 1):length(lines)
            line = strip(lines[i])
            if startswith(line, "SCALARS")
                parts = split(line)
                if length(parts) >= 2 && parts[2] == field_name
                    sc_idx = i
                    break
                end
            end
        end

        # Fallback: first scalar field if requested one is absent
        if sc_idx === nothing
            sc_idx = findnext(l -> startswith(l, "SCALARS"), lines, pd_idx)
            sc_idx === nothing && error("No SCALARS field found in $filename")
        end

        lk_idx = findnext(l -> startswith(l, "LOOKUP_TABLE"), lines, sc_idx)
        lk_idx === nothing && error("No LOOKUP_TABLE after SCALARS in $filename")

        values = Float64[]
        for i in (lk_idx + 1):length(lines)
            line = strip(lines[i])
            isempty(line) && continue
            startswith(line, "SCALARS") && break

            for tok in split(line)
                try
                    push!(values, parse(Float64, tok))
                catch
                end
                length(values) >= n_pts && break
            end
            length(values) >= n_pts && break
        end

        if length(values) < n_pts
            error("Incomplete SCALARS data in $filename: expected $n_pts, got $(length(values))")
        end

        return points[:, 1:2], values[1:n_pts]
    end
end

# ──────────────────────────────────────────────────────────────────────
# Extract column profile
# ──────────────────────────────────────────────────────────────────────

function extract_profile(points, values, x_tol=0.01)
    col_idx = findall(i -> abs(points[i,1]) < x_tol, 1:size(points,1))
    if isempty(col_idx) 
        return Float64[], Float64[] 
    end
    
    y = points[col_idx, 2]
    v = values[col_idx]
    perm = sortperm(y)
    return y[perm], v[perm]
end

# ──────────────────────────────────────────────────────────────────────
# Extract all profiles for a case
# ──────────────────────────────────────────────────────────────────────

function extract_case_profiles(prefix, output_dir, vg_model)
    profiles = Dict{Float64, Tuple}()  # time_hours => (y, h, θ)
    
    for step in 0:24
        t_hours = step * 1.0  # Each step is 1 hour
        file = joinpath(output_dir, @sprintf("%s_water_%06d.vtk", prefix, step))
        
        if isfile(file)
            try
                points, h_vals = read_vtk_points_field(file)
                y, h = extract_profile(points, h_vals)
                
                if !isempty(y)
                    θ = [theta(vg_model, hi) for hi in h]
                    profiles[t_hours] = (y, h, θ)
                end
            catch e
                @warn "Could not read step $step: $e"
            end
        end
    end
    
    return profiles
end

# ══════════════════════════════════════════════════════════════════════
# MAIN EXTRACTION
# ══════════════════════════════════════════════════════════════════════

println("\nExtracting profiles from all benchmark runs...")

all_profiles = Dict()
for case in benchmark_cases
    println("  • $(case.label)...")
    profiles = extract_case_profiles(case.prefix, output_dir, vg_celia)
    if !isempty(profiles)
        all_profiles[case.label] = (profiles=profiles, color=case.color)
        println("    ✓ Extracted $(length(profiles)) profiles")
    else
        @warn "No profiles extracted for $(case.label)"
    end
end

if isempty(all_profiles)
    error("No profiles could be extracted. Check that benchmark runs completed.")
end

println("✓ Profile extraction complete\n")

# ══════════════════════════════════════════════════════════════════════
# PLOT 1: Final profiles at t=24h
# ══════════════════════════════════════════════════════════════════════

println("Generating comparison plots...")

p1 = plot(xlabel="h [m]", ylabel="Depth from surface [m]", 
          title="Final Pressure Head Profile (t=24h)", 
          legend=:bottomleft, grid=true, gridalpha=0.3,
          size=(700, 500), xlims=(-10, 0))

p2 = plot(xlabel="θ [-]", ylabel="Depth from surface [m]",
          title="Final Water Content Profile (t=24h)",
          legend=:bottomright, grid=true, gridalpha=0.3,
          size=(700, 500), xlims=(0.1, 0.4))

for case in benchmark_cases
    if haskey(all_profiles, case.label)
        profs = all_profiles[case.label].profiles
        if haskey(profs, 24.0)  # t=24h
            y, h, θ = profs[24.0]
            depth = 1.0 .- y ./ 1.0  # Convert to depth from surface
            
            plot!(p1, h, depth, label=case.label, 
                  color=case.color, lw=2.5, marker=:circle, ms=3)
            plot!(p2, θ, depth, label=case.label,
                  color=case.color, lw=2.5, marker=:circle, ms=3)
        end
    end
end

p_final = plot(p1, p2, layout=(1,2), size=(1400, 550),
               plot_title="Celia et al. (1990) Fig 6B: Convergence Study")
display(p_final)
savefig(p_final, joinpath(plots_dir, "extract_benchmark_final_profiles.png"))
println("✓ Saved: $(joinpath(plots_dir, "extract_benchmark_final_profiles.png"))")

# ══════════════════════════════════════════════════════════════════════
# PLOT 2: Time evolution for finest timestep (dt=144s)
# ══════════════════════════════════════════════════════════════════════

println("Generating time evolution plots...")

finest_label = benchmark_cases[1].label
if haskey(all_profiles, finest_label)
    profs_finest = all_profiles[finest_label].profiles
    
    p_t_h = plot(xlabel="h [m]", ylabel="Depth from surface [m]",
               title="Time Evolution (Δt=2.4 min) — Matric Head",
               legend=:bottomleft, grid=true, gridalpha=0.3,
               size=(700, 500), xlims=(-10, 0))
    
    p_t_theta = plot(xlabel="θ [-]", ylabel="Depth from surface [m]",
                   title="Time Evolution (Δt=2.4 min) — Water Content",
                   legend=:bottomright, grid=true, gridalpha=0.3,
                   size=(700, 500), xlims=(0.1, 0.4))
    
    colors_t = [:gray70, :steelblue, :seagreen, :darkorange, :firebrick]
    for (t_idx, t_hour) in enumerate([0.0, 6.0, 12.0, 18.0, 24.0])
        if haskey(profs_finest, t_hour)
            y, h, θ = profs_finest[t_hour]
            depth = 1.0 .- y ./ 1.0
            
            plot!(p_t_h, h, depth, label="t=$(Int(t_hour))h",
                  color=colors_t[t_idx], lw=2, marker=:circle, ms=2)
            plot!(p_t_theta, θ, depth, label="t=$(Int(t_hour))h",
                  color=colors_t[t_idx], lw=2, marker=:circle, ms=2)
        end
    end
    
    p_evol = plot(p_t_h, p_t_theta, layout=(1,2), size=(1400, 550),
                  plot_title="Time Evolution (Finest: Δt=2.4 min)")
    display(p_evol)
    savefig(p_evol, joinpath(plots_dir, "extract_benchmark_time_evolution.png"))
    println("✓ Saved: $(joinpath(plots_dir, "extract_benchmark_time_evolution.png"))")
end

# ══════════════════════════════════════════════════════════════════════
# QUANTITATIVE COMPARISON
# ══════════════════════════════════════════════════════════════════════

println("\n" * "="^70)
println("QUANTITATIVE METRICS AT t=24h")
println("="^70)

theta_threshold = vg_celia.theta_r + 0.05 * (vg_celia.theta_s - vg_celia.theta_r)

println(@sprintf("%-35s %15s %15s", "Configuration", "Front Depth", "S_e (top 20cm)"))
println("-"^65)

for case in benchmark_cases
    if haskey(all_profiles, case.label)
        profs = all_profiles[case.label].profiles
        if haskey(profs, 24.0)
            y, h, θ = profs[24.0]
            
            # Front penetration (where θ exceeds threshold)
            wet_idx = findall(theta_val -> theta_val > theta_threshold, θ)
            if !isempty(wet_idx)
                front_depth = maximum(1.0 .- y[wet_idx] ./ 1.0) * 100  # cm
            else
                front_depth = 0.0
            end
            
            # Average saturation in top 20cm (y > 0.8m)
            top_20_idx = findall(i -> y[i] > 0.8, 1:length(y))
            if !isempty(top_20_idx)
                S_e_top = mean([(θ[i] - vg_celia.theta_r) / (vg_celia.theta_s - vg_celia.theta_r) 
                               for i in top_20_idx])
            else
                S_e_top = 0.0
            end
            
            @printf("%-35s %12.2f cm %12.4f\n",
                    case.label, front_depth, S_e_top)
        end
    end
end

println("-"^65)
println("\nExpected from Celia et al. (1990) Figure 6B:")
println("  • Front penetration at t=24h: ~60–70 cm")
println("  • Average S_e in top 20cm: 0.60–0.80")
println("  • Profile shape: Sharp wetting front with saturated zone above")

println("\n" * "="^70)
println("✓ Comparison complete")
println("="^70)
