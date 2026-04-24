"""
    plot_celia_comparison.jl

Extract water content (θ) and matric head (h) profiles from ADSIM VTK output
and compare quantitatively against Celia et al. (1990) Figure 6B benchmark.

Reads VTK files, extracts left column (x ≈ 0) profiles at key times,
converts h to θ using VanGenuchten SWRC, and plots for visual verification.

Usage:
    include("plot_celia_comparison.jl")
    plot_celia_comparison(project_name="CeliaCol", output_dir="output")
"""

using Printf
using Plots
using Statistics

# Include SWRC models for h→θ conversion
include("swrc_models.jl")


"""
    read_vtk_points_and_field(filename) → (points, field_values, field_name)

Parse ASCII VTK file and extract POINTS and first POINT_DATA scalar field.

Returns:
    - points :: Matrix{Float64}  (N×3, x y z coordinates)
    - field_values :: Vector{Float64}  (N scalars at each point)
    - field_name :: String  (name of scalar field, e.g., "Matric_Head")
"""
function read_vtk_points_and_field(filename::String)
    open(filename, "r") do f
        lines = readlines(f)

        # Find POINTS line
        points_idx = findfirst(l -> startswith(l, "POINTS"), lines)
        if points_idx === nothing
            error("No POINTS section found in $filename")
        end

        # Parse POINTS header: "POINTS N float"
        header_parts = split(lines[points_idx])
        N_points = parse(Int, header_parts[2])

        # Read N_points coordinate triplets
        points = zeros(Float64, N_points, 3)
        line_idx = points_idx + 1
        pt = 0
        while pt < N_points && line_idx <= length(lines)
            coords = split(lines[line_idx])
            for i in 1:min(3, length(coords))
                pt += 1
                if pt > N_points
                    break
                end
                points[pt, mod1(i, 3)] = parse(Float64, coords[i])
            end
            line_idx += 1
        end

        # Find POINT_DATA and first scalar field
        pdata_idx = findfirst(l -> startswith(l, "POINT_DATA"), lines)
        if pdata_idx === nothing
            error("No POINT_DATA section found in $filename")
        end

        # Look for SCALARS line
        scalars_idx = findnext(l -> startswith(l, "SCALARS"), lines, pdata_idx)
        if scalars_idx === nothing
            error("No SCALARS field found in POINT_DATA")
        end

        # Parse SCALARS header: "SCALARS field_name float"
        scalars_parts = split(lines[scalars_idx])
        field_name = scalars_parts[2]

        # Skip LOOKUP_TABLE line
        lookup_idx = findnext(l -> startswith(l, "LOOKUP_TABLE"), lines, scalars_idx)
        data_start = (lookup_idx !== nothing) ? lookup_idx + 1 : scalars_idx + 1

        # Read N_points scalar values
        field_values = Float64[]
        line_idx = data_start
        while length(field_values) < N_points && line_idx <= length(lines)
            values_str = split(lines[line_idx])
            for val_str in values_str
                if length(field_values) < N_points
                    push!(field_values, parse(Float64, val_str))
                end
            end
            line_idx += 1
        end

        if length(field_values) != N_points
            error("Expected $N_points scalar values, got $(length(field_values))")
        end

        return points, field_values, field_name
    end
end


"""
    extract_column_profile(points, field_values, x_target, tol_x) → (y_col, field_col)

Extract profile along vertical line at x ≈ x_target ± tol_x.
Returns: y coordinates (sorted ascending), field values at those points.
"""
function extract_column_profile(points::Matrix{Float64}, field_values::Vector{Float64},
                                 x_target::Float64, tol_x::Float64)
    # Find all nodes with x ≈ x_target
    col_indices = findall(i -> abs(points[i, 1] - x_target) < tol_x, 1:size(points, 1))

    if isempty(col_indices)
        @warn "No nodes found with x ≈ $(x_target)"
        return Float64[], Float64[]
    end

    # Extract y and field values
    y_col = points[col_indices, 2]
    field_col = field_values[col_indices]

    # Sort by y (ascending = bottom to top)
    perm = sortperm(y_col)
    y_col = y_col[perm]
    field_col = field_col[perm]

    return y_col, field_col
end


"""
    read_all_vtk_files(output_dir, project_name) → Dict{Float64, String}

Find all VTK files matching pattern and return time → filename mapping.
"""
function read_all_vtk_files(output_dir::String, project_name::String)
    time_file_map = Dict{Float64, String}()

    # Pattern: project_name_water_NNNNNN.vtk
    pattern = project_name * "_water_"

    for file in readdir(output_dir)
        if endswith(file, ".vtk") && contains(file, pattern)
            # Extract step number from filename
            # e.g., "CeliaCol_water_000001.vtk" → step=1
            parts = split(file, "_")
            if length(parts) >= 3
                step_str = parts[end][1:6]  # first 6 digits before .vtk
                try
                    step = parse(Int, step_str)
                    # For now, store by step number (will map to time later)
                    time_file_map[Float64(step)] = joinpath(output_dir, file)
                catch
                    # Skip if step number doesn't parse
                end
            end
        end
    end

    return time_file_map
end


"""
    plot_celia_comparison(project_name, output_dir; times_to_extract=[1,6,12,24])

Main function: extract profiles at key times, convert to θ, and plot comparison.

Arguments:
    - project_name: Prefix of VTK files (e.g., "CeliaCol")
    - output_dir: Directory containing VTK files (default: "output")
    - times_to_extract: Hours to extract (default: [1,6,12,24])
    - vg_model: VanGenuchten SWRC model (optional, default: Celia 1990 parameters)
"""
function plot_celia_comparison(project_name::String="CeliaCol";
                                output_dir::String="output",
                                times_to_extract::Vector{Int}=[1, 6, 12, 24],
                                vg_model=nothing)

    # Use default Celia 1990 parameters if not provided
    if vg_model === nothing
        vg_model = VanGenuchten(
            θ_r=0.102, θ_s=0.368,
            α=0.335, n_vg=2.0,
            K_sx=9.4e-12, K_sy=9.4e-12,
            L=0.5
        )
    end

    # Find all VTK files
    println("Searching for VTK files in: $output_dir")
    file_map = read_all_vtk_files(output_dir, project_name)

    if isempty(file_map)
        @warn "No VTK files found matching pattern: $(project_name)_water_*.vtk"
        return nothing
    end

    println("Found $(length(file_map)) VTK files")

    # Convert step numbers to times (hours)
    # Assuming: step 0 = t=0, step 1 = t=1h, step N = t=Nh
    times_hours = sort(collect(keys(file_map)))
    filenames = [file_map[t] for t in times_hours]

    println("\nAvailable times (hours):")
    for t in times_hours
        println("  Step $(lpad(Int(t), 3)): t=$(t) hours")
    end

    # Create profiles dictionary
    profiles = Dict{Int, NamedTuple}()  # time_hours → (y, h, θ)

    for tp in times_to_extract
        if tp ∉ times_hours
            @warn "Time $tp hours not found in output"
            continue
        end

        filename = file_map[Float64(tp)]
        println("\nExtracting profile from: $(basename(filename))")

        try
            points, h_values, field_name = read_vtk_points_and_field(filename)
            println("  ✓ Loaded $(size(points, 1)) nodes, field: $field_name")

            # Extract left column (x ≈ 0)
            x_min = minimum(points[:, 1])
            tol_x = (maximum(points[:, 1]) - x_min) / 20.0  # 5% tolerance
            y_col, h_col = extract_column_profile(points, h_values, x_min, tol_x)

            # Convert h to θ
            θ_col = [theta(vg_model, h) for h in h_col]

            profiles[tp] = (y=y_col, h=h_col, θ=θ_col)

            # Compute statistics
            θ_min = minimum(θ_col)
            θ_max = maximum(θ_col)
            h_min = minimum(h_col)
            h_max = maximum(h_col)
            @printf("  Column profile: %d nodes\n", length(y_col))
            @printf("    h range: [%.3f, %.3f] m\n", h_min, h_max)
            @printf("    θ range: [%.4f, %.4f] [-]\n", θ_min, θ_max)

        catch e
            @warn "Failed to extract profile from $filename: $e"
            continue
        end
    end

    # Create comparison plots
    if !isempty(profiles)
        p_h = plot_profiles_h(profiles, vg_model)
        p_θ = plot_profiles_theta(profiles, vg_model)
        p_combined = plot(p_h, p_θ, layout=(1,2), size=(1000, 500),
                         plot_title="Celia 1990 Fig 6B Benchmark: FEM Results")
        display(p_combined)

        # Save plots
        try
            savefig(p_h, "output/celia_profiles_h.png")
            savefig(p_θ, "output/celia_profiles_theta.png")
            savefig(p_combined, "output/celia_comparison.png")
            println("\n✓ Saved comparison plots to output/")
        catch e
            @warn "Could not save plots: $e"
        end
    end

    return profiles
end


"""
    plot_profiles_h(profiles, vg_model) → plot object

Plot matric head h(y) at different times.
"""
function plot_profiles_h(profiles::Dict{Int, NamedTuple}, vg_model)
    p = plot(xlabel="h [m]", ylabel="Depth from surface [m]",
             title="Matric Head Profiles",
             legend=:bottomright, grid=true, gridalpha=0.3,
             size=(500, 400))

    times = sort(collect(keys(profiles)))
    colors = [:gray70, :steelblue, :seagreen, :firebrick, :darkorange]

    for (i, tp) in enumerate(times)
        prof = profiles[tp]
        y_depth = maximum(prof.y) .- prof.y  # Depth from surface (positive downward)
        color = colors[mod1(i, length(colors))]
        plot!(p, prof.h, y_depth, label="t=$(tp)h",
              color=color, lw=2, marker=:circle, ms=3)
    end

    # Flip y-axis so depth increases downward
    plot!(p, yflip=true)

    return p
end


"""
    plot_profiles_theta(profiles, vg_model) → plot object

Plot water content θ(y) at different times.
"""
function plot_profiles_theta(profiles::Dict{Int, NamedTuple}, vg_model)
    p = plot(xlabel="θ [-]", ylabel="Depth from surface [m]",
             title="Water Content Profiles",
             legend=:bottomright, grid=true, gridalpha=0.3,
             size=(500, 400))

    times = sort(collect(keys(profiles)))
    colors = [:gray70, :steelblue, :seagreen, :firebrick, :darkorange]

    for (i, tp) in enumerate(times)
        prof = profiles[tp]
        y_depth = maximum(prof.y) .- prof.y  # Depth from surface
        color = colors[mod1(i, length(colors))]
        plot!(p, prof.θ, y_depth, label="t=$(tp)h",
              color=color, lw=2, marker=:circle, ms=3)
    end

    # Flip y-axis so depth increases downward
    plot!(p, yflip=true)

    return p
end


"""
    quantitative_comparison(profiles, vg_model) → summary statistics

Compute front penetration, saturation level, and other metrics for comparison.
"""
function quantitative_comparison(profiles::Dict{Int, NamedTuple}, vg_model)
    println("\n" * "="^65)
    println("QUANTITATIVE COMPARISON")
    println("="*65)

    times = sort(collect(keys(profiles)))

    for tp in times
        prof = profiles[tp]
        θ_init = vg_model.θ_r + 0.05 * (vg_model.θ_s - vg_model.θ_r)  # "wet" threshold
        y_depth = maximum(prof.y) .- prof.y

        # Front penetration: depth where θ first exceeds threshold
        wet_indices = findall(θ -> θ > θ_init, prof.θ)
        if !isempty(wet_indices)
            front_depth = y_depth[wet_indices[1]]
        else
            front_depth = 0.0
        end

        # Average saturation in top 20 cm
        near_surface = findall(d -> d <= 0.20, y_depth)
        if !isempty(near_surface)
            Se_avg = mean((prof.θ[near_surface] .- vg_model.θ_r) ./ (vg_model.θ_s - vg_model.θ_r))
        else
            Se_avg = 0.0
        end

        @printf("\n  t = %2d hours:\n", tp)
        @printf("    Front penetration:  %.3f m (%.1f cm)\n", front_depth, front_depth*100)
        @printf("    Avg saturation (top 20cm): %.4f (%.1f%%)\n", Se_avg, 100*Se_avg)
    end

    println("="*65)
end


# ══════════════════════════════════════════════════════════════════════════════
# USAGE
# ══════════════════════════════════════════════════════════════════════════════

if abspath(PROGRAM_FILE) == @__FILE__
    # Run as script
    println("Celia 1990 Fig 6B Comparison Script")
    println("="*65)

    # Extract profiles
    profiles = plot_celia_comparison("CeliaCol"; output_dir="output")

    # Print statistics
    if profiles !== nothing && !isempty(profiles)
        quantitative_comparison(profiles, VanGenuchten(
            θ_r=0.102, θ_s=0.368,
            α=0.335, n_vg=2.0,
            K_sx=9.4e-12, K_sy=9.4e-12,
            L=0.5
        ))
    end
end
