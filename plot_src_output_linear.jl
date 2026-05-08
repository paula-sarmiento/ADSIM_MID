# Plot src/output Linear Diffusion VTK Data
# Extracts and visualizes Matric_Head profiles from
# src/output/LinDiffusion_water_*.vtk files

using Printf, Plots, Statistics

const PROJECT_ROOT = pwd()
const SRC_DIR = joinpath(PROJECT_ROOT, "src")

# Try to include version.jl if available, otherwise skip
try
    include(joinpath(SRC_DIR, "version.jl"))
    using .ADSIMVersion: get_version
catch
    # Version not available, continue anyway
end

# VTK Parser: Extract Matric_Head and coordinates
function parse_vtk_linear_data(filepath::String)
    """
    Parse VTK file and extract spatial points and Matric_Head field.
    Returns: (y_coords, matric_head, time)
    """
    lines = readlines(filepath)
    
    # Extract time from line 2
    time = 0.0
    if length(lines) >= 2
        line2 = lines[2]
        if occursin("Time", line2)
            m = match(r"Time\s*=\s*([0-9eE+\.-]+)", line2)
            if m !== nothing
                time = parse(Float64, m.captures[1])
            end
        end
    end
    
    # Find POINT_DATA and Matric_Head
    pd_idx = findfirst(l -> startswith(l, "POINT_DATA"), lines)
    if pd_idx === nothing
        error("POINT_DATA not found in $filepath")
    end
    
    mh_search = findfirst(l -> occursin("SCALARS Matric_Head", l), lines[pd_idx:end])
    if mh_search === nothing
        error("Matric_Head not found in $filepath")
    end
    
    mh_idx = pd_idx + mh_search - 1
    data_start = mh_idx + 2  # Skip SCALARS and LOOKUP_TABLE lines
    
    # Find POINTS section for y coordinates
    points_idx = findfirst(l -> startswith(l, "POINTS"), lines)
    if points_idx === nothing
        error("POINTS not found in $filepath")
    end
    
    m = match(r"POINTS\s+(\d+)", lines[points_idx])
    n_points = parse(Int, m.captures[1])
    
    y_coords = Float64[]
    for i in 1:n_points
        parts = split(lines[points_idx + i])
        push!(y_coords, parse(Float64, parts[2]))  # second coordinate is y
    end
    
    # Read Matric_Head values (stop at next SCALARS or end of data)
    matric_head = Float64[]
    for i in data_start:length(lines)
        line = strip(lines[i])
        if startswith(line, "SCALARS") || startswith(line, "VECTORS")
            break
        end
        if line != ""
            try
                push!(matric_head, parse(Float64, line))
            catch
                break
            end
        end
    end
    
    if length(matric_head) != length(y_coords)
        error("Mismatch in data: $(length(y_coords)) points but $(length(matric_head)) Matric_Head values")
    end
    
    # Sort by y coordinate
    perm = sortperm(y_coords)
    return (y_coords[perm], matric_head[perm], time)
end

# Main plotting routine
function main()
    output_dir = "src/output"
    
    if !isdir(output_dir)
        error("Output directory not found: $output_dir")
    end
    
    # Collect all Linear VTK files
    vtk_files = filter(f -> occursin(r"^LinDiffusion_water_\d{6}\.vtk$", f), readdir(output_dir))
    sort!(vtk_files)
    
    if isempty(vtk_files)
        error("No LinDiffusion VTK files found in $output_dir")
    end
    
    println("Found $(length(vtk_files)) VTK files")
    
    # Select representative times for plotting
    indices = [1, div(length(vtk_files), 4), div(length(vtk_files), 2), div(3*length(vtk_files), 4), length(vtk_files)]
    unique!(indices)
    
    # Parse VTK data at selected times
    plot_data = []
    times = []
    
    for idx in indices
        filepath = joinpath(output_dir, vtk_files[idx])
        println("Processing: $(vtk_files[idx])")
        
        try
            y_coords, matric_head, time = parse_vtk_linear_data(filepath)
            push!(plot_data, (y_coords, matric_head))
            push!(times, time)
        catch e
            println("Error processing $(vtk_files[idx]): $e")
        end
    end
    
    if isempty(plot_data)
        error("No valid data extracted from VTK files")
    end
    
    # Create plot
    p = plot(xlabel="Matric Head h (m)", ylabel="Depth y (m)", title="src/output Linear Diffusion - Matric Head Profiles",
             legend=:bottomright, size=(900, 600), margin=5Plots.mm)
    
    for i in eachindex(plot_data)
        y, h = plot_data[i]
        if length(y) == length(h)
            label_str = @sprintf("t = %.3f s", times[i])
            plot!(p, h, y, label=label_str, linewidth=2, marker=:o, markersize=3, alpha=0.7)
        else
            println("Warning: length mismatch for index $i (y=$(length(y)), h=$(length(h)))")
        end
    end
    
    # Save plot
    output_plot = joinpath(output_dir, "src_output_linear_profiles.png")
    savefig(p, output_plot)
    println("Plot saved: $output_plot")
    
    # Summary statistics
    println("\n" * "-"^60)
    println("src/output Linear Diffusion VTK Data Summary")
    println("-"^60)
    println("Total VTK files: $(length(vtk_files))")
    println("Time range: $(times[begin]) to $(times[end]) s")
    println("Number of spatial profiles plotted: $(length(plot_data))")
    for i in eachindex(plot_data)
        y, h = plot_data[i]
        println(@sprintf("  t=%.3f s: y ∈ [%.3f, %.3f] m, h ∈ [%.3f, %.3f] m", 
                        times[i], minimum(y), maximum(y), minimum(h), maximum(h)))
    end
    println("-"^60)
end

if !isinteractive()
    main()
end
