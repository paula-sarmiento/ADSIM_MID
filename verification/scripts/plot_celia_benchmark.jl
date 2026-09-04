#!/usr/bin/env julia
"""
Plot Celia et al. (1990) Figure 6B benchmark verification from VTK output

Reads VTK files from implicit Richards solver and creates:
  1. Vertical pressure head profiles vs. time (at domain center x=0.005m)
  2. Vertical water content profiles vs. time
  3. Wetting front evolution over 24 hours
  4. Comparison between different time step factors

Usage:
  julia plot_celia_benchmark.jl

Output: verification/plots/ directory with PNG files

Dependencies: Plots, StatsPlots, Printf, glob
"""

using Plots
using Printf
using Statistics
using LinearAlgebra


# ============================================================================
# UTILITIES
# ============================================================================

function glob_pattern(dir::String, pattern::String)::Vector{String}
    """Find files matching pattern in directory"""
    if !isdir(dir)
        return String[]
    end
    try
        files = readdir(dir, join=true)
        prefix = split(pattern, "*", limit=2)[1]
        return filter(f -> isfile(f) && startswith(basename(f), prefix), files)
    catch
        return String[]
    end
end


# ============================================================================
# VTK READER
# ============================================================================

function read_vtk_unstructured(filename::String)
    """
    Read unstructured VTK file (ADSIM output format).
    
    Returns: (coords, data_dict, time_val)
      coords: (n_nodes, 3) matrix of node coordinates
      data_dict: Dict of arrays {variable_name => vector}
      time_val: Scalar time value from VTK file
    """
    
    content = read(filename, String)
    
    lines = split(content, '\n')
    
    # Parse time value from either legacy TIME_VALUE metadata or the file title line
    time_val = nothing
    for line in lines
        if contains(line, "TIME_VALUE")
            time_val = parse(Float64, split(line)[end])
            break
        elseif contains(line, "Time =")
            m = match(r"Time\s*=\s*([0-9eE+\-.]+)", line)
            if m !== nothing
                time_val = parse(Float64, m.captures[1])
                break
            end
        end
    end
    
    # Find POINTS section
    points_idx = findfirst(line -> startswith(line, "POINTS"), lines)
    if points_idx === nothing
        error("No POINTS section found in $filename")
    end
    
    # Parse number of points
    n_points = parse(Int, split(lines[points_idx])[2])
    
    # Extract coordinates
    coords = zeros(n_points, 3)
    idx = 1
    for i in (points_idx + 1):length(lines)
        line = strip(lines[i])
        if line == "" || !isdigit(line[1]) && line[1] != '-'
            break
        end
        try
            vals = parse.(Float64, split(line))
            if idx <= n_points && length(vals) >= 3
                coords[idx, :] .= vals[1:3]
                idx += 1
            end
        catch
            continue
        end
    end
    
    # Find POINT_DATA section and extract variables
    data_dict = Dict{String, Vector{Float64}}()
    pointdata_idx = findfirst(line -> startswith(line, "POINT_DATA"), lines)
    
    if pointdata_idx !== nothing
        i = pointdata_idx + 1
        while i <= length(lines)
            line = strip(lines[i])
            
            if startswith(line, "SCALARS")
                parts = split(line)
                var_name = parts[2]
                n_components = length(parts) > 3 ? parse(Int, parts[4]) : 1
                
                # Skip LOOKUP_TABLE line if present
                if i + 1 <= length(lines) && contains(lines[i+1], "LOOKUP_TABLE")
                    i += 2
                else
                    i += 1
                end
                
                # Read data
                data = Float64[]
                for _ in 1:(n_components * n_points)
                    if i > length(lines)
                        break
                    end
                    try
                        val = parse(Float64, split(strip(lines[i]))[1])
                        push!(data, val)
                        i += 1
                    catch
                        i += 1
                    end
                end
                
                if length(data) == n_components * n_points
                    data_dict[var_name] = data
                end
            else
                i += 1
            end
        end
    end
    
    return coords, data_dict, time_val
end


# ============================================================================
# PROFILE EXTRACTION
# ============================================================================

function extract_vertical_profile(coords::Matrix, h_data::Vector, theta_data::Vector;
                                  x_target::Float64=0.005, tol::Float64=0.001)
    """
    Extract vertical profile at x ≈ x_target.
    
    Returns: (y_vals, h_vals, theta_vals) sorted by y (ascending)
    """
    
    # Find nodes near x_target
    x_diffs = abs.(coords[:, 1] .- x_target)
    x_indices = findall(x_diffs .< tol)
    
    if isempty(x_indices)
        # Use closest x
        x_closest = argmin(x_diffs)
        x_actual = coords[x_closest, 1]
        x_indices = findall(abs.(coords[:, 1] .- x_actual) .< tol)
    end
    
    if isempty(x_indices)
        error("No nodes found near x=$x_target")
    end
    
    # Extract y and data
    y_vals = coords[x_indices, 2]
    h_vals = h_data[x_indices]
    theta_vals = theta_data[x_indices]
    
    # Sort by y (ascending)
    sort_idx = sortperm(y_vals)
    
    return y_vals[sort_idx], h_vals[sort_idx], theta_vals[sort_idx]
end


function convert_height_to_depth(y_vals::Vector; y_max::Float64=1.0)
    """Convert y-coordinate to depth from surface (z = y_max - y)"""
    return y_max .- y_vals
end


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

function plot_pressure_head_profiles(vtk_files::Vector{String}; output_dir::String="verification/plots")
    """
    Plot pressure head profiles at selected times.
    Mimics Celia et al. (1990) Figure 6B: h(z) at t = 0, 6, 12, 24 hours
    """
    
    mkpath(output_dir)
    
    # Extract time and data for each file
    times = Float64[]
    profiles_h = Tuple{Vector{Float64}, Vector{Float64}}[]
    profiles_theta = Tuple{Vector{Float64}, Vector{Float64}}[]
    
    for filename in sort(vtk_files)
        try
            coords, data_dict, time_val = read_vtk_unstructured(filename)
            
            if time_val === nothing
                continue
            end
            
            h_key = haskey(data_dict, "Matric_Head") ? "Matric_Head" : (haskey(data_dict, "h") ? "h" : nothing)
            theta_key = haskey(data_dict, "Water_Content") ? "Water_Content" : (haskey(data_dict, "theta_w") ? "theta_w" : nothing)

            if h_key === nothing || theta_key === nothing
                println("Skipping $(basename(filename)): Missing required fields")
                continue
            end
            
            y, h, theta = extract_vertical_profile(coords, data_dict[h_key], data_dict[theta_key])
            
            push!(times, time_val)
            push!(profiles_h, (y, h))
            push!(profiles_theta, (y, theta))
            
            @printf("✓ Loaded: %-40s t=%8.1f s (%.2f h)\n", 
                    basename(filename), time_val, time_val/3600)
            
        catch e
            println("✗ Failed to load $(basename(filename)): $e")
        end
    end
    
    if isempty(times)
        println("No valid VTK files found!")
        return
    end
    
    # ── Figure 1: Pressure Head Profiles ────────────────────────────────
    p1 = plot(xlabel="Pressure Head h (cm)", ylabel="Depth from surface (m)",
              title="Celia et al. (1990) Figure 6B: Pressure Head Evolution",
              legend=:bottomright, size=(900, 750), dpi=150, grid=true, gridalpha=0.3)
    
    for (i, (t, (y, h))) in enumerate(zip(times, profiles_h))
        depth = convert_height_to_depth(y)
        plot!(p1, h, depth, marker=:circle, markersize=4, linewidth=2.5,
              label="t = $(round(t/3600, digits=2)) h", alpha=0.8)
    end
    
    yaxis!(p1, :flip)  # Depth increases downward
    
    savefig(p1, joinpath(output_dir, "celia_pressure_head_profiles.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_pressure_head_profiles.png"))")
    
    # ── Figure 2: Water Content Profiles ───────────────────────────────
    p2 = plot(xlabel="Volumetric Water Content θ (-)", ylabel="Depth from surface (m)",
              title="Celia et al. (1990): Water Content Evolution",
              legend=:bottomright, size=(900, 750), dpi=150, grid=true, gridalpha=0.3)
    
    for (i, (t, (y, theta))) in enumerate(zip(times, profiles_theta))
        depth = convert_height_to_depth(y)
        plot!(p2, theta, depth, marker=:square, markersize=4, linewidth=2.5,
              label="t = $(round(t/3600, digits=2)) h", alpha=0.8)
    end
    
    yaxis!(p2, :flip)
    
    savefig(p2, joinpath(output_dir, "celia_water_content_profiles.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_water_content_profiles.png"))")
    
    # ── Figure 3: Wetting Front Evolution ───────────────────────────────
    theta_r = 0.102
    theta_s = 0.368
    theta_front = (theta_r + theta_s) / 2.0
    
    wetting_front_depths = Float64[]
    wetting_front_times = Float64[]
    
    for (t, (y, theta)) in zip(times, profiles_theta)
        depth = convert_height_to_depth(y)
        
        # Find deepest point where theta > theta_front
        indices = findall(theta .> theta_front)
        if !isempty(indices)
            wetting_front_depth = depth[last(indices)]
            push!(wetting_front_depths, wetting_front_depth)
            push!(wetting_front_times, t / 3600.0)
        end
    end
    
    p3 = plot(xlabel="Time (hours)", ylabel="Depth from surface (m)",
              title="Celia Benchmark: Wetting Front Propagation",
              legend=:topleft, size=(900, 750), dpi=150, grid=true, gridalpha=0.3)
    
    if !isempty(wetting_front_times)
        plot!(p3, wetting_front_times, wetting_front_depths, marker=:circle, markersize=7,
              linewidth=2.5, label="Wetting front", color=:darkblue, alpha=0.8)
        
        # Fill between
        plot!(p3, wetting_front_times, wetting_front_depths, fillrange=0,
              alpha=0.2, label="Wet zone", color=:blue, linewidth=0)
    end
    
    xlims!(p3, 0, 24)
    yaxis!(p3, :flip)
    
    savefig(p3, joinpath(output_dir, "celia_wetting_front.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_wetting_front.png"))")
    
    println("\n" * "="^70)
    println("SUMMARY: Generated 3 publication-quality figures")
    println("  • Pressure head profiles ($(length(times)) time points)")
    println("  • Water content profiles ($(length(times)) time points)")
    println("  • Wetting front evolution over 24 hours")
    println("Output directory: $(abspath(output_dir))")
    println("="^70)
end


function compare_timesteps(test_dirs::Dict{String, String}; output_dir::String="verification/plots")
    """
    Compare solutions with different time step factors.
    
    Example:
      compare_timesteps(Dict(
        "Δt=0.54s (F=100k)" => "output",
        "Δt=2.71s (F=500k)" => "output"
      ))
    """
    
    mkpath(output_dir)
    
    p_h = plot(xlabel="Pressure Head h (cm)", ylabel="Depth from surface (m)",
               title="Celia Benchmark: Time Step Convergence (Pressure Head, t=24h)",
               legend=:bottomright, size=(900, 750), dpi=150, grid=true, gridalpha=0.3)
    
    p_theta = plot(xlabel="Volumetric Water Content θ (-)", ylabel="Depth from surface (m)",
                   title="Celia Benchmark: Time Step Convergence (Water Content, t=24h)",
                   legend=:bottomright, size=(900, 750), dpi=150, grid=true, gridalpha=0.3)
    
    for (label, test_prefix) in test_dirs
        # Find final VTK file
        dir = dirname(test_prefix)
        base = basename(test_prefix)
        files = glob_pattern(dir, base * "*_water_*.vtk")
        if isempty(files)
            println("✗ No VTK files found matching $(test_prefix)")
            continue
        end
        
        # Extract step numbers and keep only files that match the suffix pattern
        step_pairs = Int[]
        valid_files = String[]
        for f in files
            m = match(r"_(\d+)\.vtk$", basename(f))
            if m !== nothing
                push!(step_pairs, parse(Int, m.captures[1]))
                push!(valid_files, f)
            end
        end

        if isempty(valid_files)
            println("✗ No step-numbered VTK files found in $(test_prefix)")
            continue
        end

        final_file = valid_files[argmax(step_pairs)]
        
        try
            coords, data_dict, time_val = read_vtk_unstructured(final_file)
            
            h_key = haskey(data_dict, "Matric_Head") ? "Matric_Head" : (haskey(data_dict, "h") ? "h" : nothing)
            theta_key = haskey(data_dict, "Water_Content") ? "Water_Content" : (haskey(data_dict, "theta_w") ? "theta_w" : nothing)

            if h_key === nothing || theta_key === nothing
                continue
            end
            
            y, h, theta = extract_vertical_profile(coords, data_dict[h_key], data_dict[theta_key])
            depth = convert_height_to_depth(y)
            
            plot!(p_h, h, depth, marker=:circle, markersize=5, linewidth=2.5,
                  label=label, alpha=0.8)
            plot!(p_theta, theta, depth, marker=:square, markersize=5, linewidth=2.5,
                  label=label, alpha=0.8)
            
            println("✓ Loaded $label from $(basename(final_file))")
            
        catch e
            println("✗ Failed to load $label: $e")
        end
    end
    
    yaxis!(p_h, :flip)
    yaxis!(p_theta, :flip)
    
    p_combined = plot(p_h, p_theta, layout=(1,2), size=(1600, 700))
    savefig(p_combined, joinpath(output_dir, "celia_timestep_comparison.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_timestep_comparison.png"))")
end


# ============================================================================
# CELIA FIG 6B STYLE: head [cm] vs depth [cm], selected snapshots
# ============================================================================

"""
    plot_celia_fig6b_cm(vtk_files; target_hours, output_dir)

Reproduce Celia et al. (1990) Figure 6B style:
  - x-axis: matric head h [cm]  (-1000 → 0)
  - y-axis: depth from surface  [cm]  (0 at top → 100 at bottom)
  - Selected snapshots at `target_hours` (default 0, 6, 12, 18, 24 h)

Also generates matching θ(z) [cm] figure.
All input VTK coordinates and h values are assumed to be in metres (SI).
"""
function plot_celia_fig6b_cm(vtk_files::Vector{String};
                              target_hours::Vector{<:Real}=[0, 6, 12, 18, 24],
                              output_dir::String="verification/plots",
                              domain_height_m::Float64=1.0)

    mkpath(output_dir)

    # Collect only the snapshots closest to target_hours
    tol_s = 60.0  # within 1 min is a match
    selected = Dict{Float64, NamedTuple}()  # target_hour => (y, h, θ)

    for filename in sort(vtk_files)
        try
            coords, data_dict, time_val = read_vtk_unstructured(filename)
            time_val === nothing && continue

            h_key     = haskey(data_dict, "Matric_Head")   ? "Matric_Head"   :
                        haskey(data_dict, "h")              ? "h"             : nothing
            theta_key = haskey(data_dict, "Water_Content") ? "Water_Content" :
                        haskey(data_dict, "theta_w")        ? "theta_w"       : nothing
            (h_key === nothing || theta_key === nothing) && continue

            y, h, θ = extract_vertical_profile(coords, data_dict[h_key], data_dict[theta_key])

            t_h = time_val / 3600.0
            for th in target_hours
                if abs(t_h - th) * 3600.0 < tol_s
                    if !haskey(selected, Float64(th)) || abs(t_h - th) < abs(selected[Float64(th)].t_h - th)
                        selected[Float64(th)] = (y=y, h=h, θ=θ, t_h=t_h)
                    end
                end
            end
        catch e
            println("✗ $(basename(filename)): $e")
        end
    end

    if isempty(selected)
        println("No snapshots matched target hours $target_hours")
        return
    end

    colors_t = [:gray60, :steelblue, :seagreen, :darkorange, :firebrick]
    sorted_hours = sort(collect(keys(selected)))

    # ── Figure A: Matric head [cm] vs Depth [cm] ───────────────────────
    pA = plot(xlabel="Matric Head h (cm)", ylabel="Depth from surface (cm)",
              title="Celia et al. (1990) Fig. 6B — Matric Head",
              legend=:bottomleft, size=(700, 700), dpi=150,
              grid=true, gridalpha=0.3,
              xlims=(-1050, 50), ylims=(0, domain_height_m * 100))
    yaxis!(pA, :flip)

    for (ci, th) in enumerate(sorted_hours)
        d = selected[th]
        h_cm    = d.h .* 100.0                          # m → cm
        depth_cm = (domain_height_m .- d.y) .* 100.0   # m → cm, 0=surface
        lbl = th == 0 ? "t = 0 h (IC)" : "t = $(Int(th)) h"
        plot!(pA, h_cm, depth_cm; label=lbl,
              color=colors_t[mod1(ci, length(colors_t))], lw=2.5,
              marker=:circle, ms=3, alpha=0.9)
    end

    # Annotate BCs
    annotate!(pA, [(-950, 2,  text("h_top = −75 cm",  8, :left, :black)),
                   (-950, 98, text("h_bot = −1000 cm", 8, :left, :black))])

    savefig(pA, joinpath(output_dir, "celia_fig6b_head_cm.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_fig6b_head_cm.png"))")

    # ── Figure B: Water content θ vs Depth [cm] ────────────────────────
    pB = plot(xlabel="Volumetric Water Content θ (−)", ylabel="Depth from surface (cm)",
              title="Celia et al. (1990) Fig. 6B — Water Content",
              legend=:bottomright, size=(700, 700), dpi=150,
              grid=true, gridalpha=0.3,
              xlims=(0.08, 0.40), ylims=(0, domain_height_m * 100))
    yaxis!(pB, :flip)

    for (ci, th) in enumerate(sorted_hours)
        d = selected[th]
        depth_cm = (domain_height_m .- d.y) .* 100.0
        lbl = th == 0 ? "t = 0 h (IC)" : "t = $(Int(th)) h"
        plot!(pB, d.θ, depth_cm; label=lbl,
              color=colors_t[mod1(ci, length(colors_t))], lw=2.5,
              marker=:circle, ms=3, alpha=0.9)
    end

    savefig(pB, joinpath(output_dir, "celia_fig6b_theta_cm.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_fig6b_theta_cm.png"))")

    # ── Figure C: Side-by-side ──────────────────────────────────────────
    pC = plot(pA, pB, layout=(1,2), size=(1400, 700), dpi=150,
              plot_title="Celia et al. (1990) Figure 6B Verification")
    savefig(pC, joinpath(output_dir, "celia_fig6b_combined_cm.png"))
    println("✓ Saved: $(joinpath(output_dir, "celia_fig6b_combined_cm.png"))")

    println("\n  Snapshots included:")
    for th in sorted_hours
        @printf("    t = %5.1f h → file t_h = %.4f h\n", th, selected[th].t_h)
    end
end


# ============================================================================
# MAIN
# ============================================================================

println("\n" * "="^70)
println("CELIA BENCHMARK VISUALIZATION")
println("="^70)

# Get project root (script now lives in verification/scripts)
project_root = normpath(joinpath(@__DIR__, "..", ".."))
output_base = joinpath(project_root, "verification", "output")
figures_dir = joinpath(project_root, "verification", "plots")
mkpath(figures_dir)

# Benchmark sweep cases (Celia et al. 1990 Fig 6B, aligned with Colab)
sweep_cases = [
    (prefix="CeliaCol_dt144",  label="Δt=2.4 min (144 s)",  color=:black),
    (prefix="CeliaCol_dt720",  label="Δt=12 min (720 s)",   color=:steelblue),
    (prefix="CeliaCol_dt3600", label="Δt=60 min (3600 s)",  color=:firebrick),
]

# [1] Profiles for finest timestep (dt=144s) at t = 0, 6, 12, 18, 24 h
finest = sweep_cases[1]
finest_files = sort(glob_pattern(output_base, "$(finest.prefix)_water_"))
if !isempty(finest_files)
    println("\n[1] Plotting time-evolution profiles for finest run ($(finest.label))...")
    plot_pressure_head_profiles(finest_files, output_dir=figures_dir)
else
    println("✗ No VTK files found for $(finest.prefix)")
end

# [2] Timestep convergence comparison at t=24h
available = [(c.label => joinpath(output_base, c.prefix))
             for c in sweep_cases
             if !isempty(glob_pattern(output_base, "$(c.prefix)_water_"))]

if length(available) >= 2
    println("\n[2] Comparing $(length(available)) timestep cases at t=24h...")
    compare_timesteps(Dict(available...), output_dir=figures_dir)
elseif length(available) == 1
    println("⚠ Only one case found; skipping timestep comparison.")
else
    println("✗ No benchmark VTK files found in $output_base")
end

# [3] Celia Fig 6B style: h [cm] vs depth [cm], snapshots at 0/6/12/18/24 h
if !isempty(finest_files)
    println("\n[3] Celia Fig 6B style figure (h in cm, depth in cm)...")
    plot_celia_fig6b_cm(finest_files, output_dir=figures_dir)
end

println("\n✓ Visualization complete!\n")
