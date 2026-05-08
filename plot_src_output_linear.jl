# Plot Linear Diffusion VTK Data - Overlay src/output vs verification/output
# Extracts and visualizes Matric_Head profiles from both output directories

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
    src_dir  = "src/output"
    ver_dir  = "verification/output"

    for d in (src_dir, ver_dir)
        isdir(d) || error("Output directory not found: $d")
    end

    vtk_pattern = r"^LinDiffusion_water_\d{6}\.vtk$"
    src_files = sort(filter(f -> occursin(vtk_pattern, f), readdir(src_dir)))
    ver_files = sort(filter(f -> occursin(vtk_pattern, f), readdir(ver_dir)))

    isempty(src_files) && error("No VTK files found in $src_dir")
    isempty(ver_files) && error("No VTK files found in $ver_dir")

    println("src/output:          $(length(src_files)) VTK files")
    println("verification/output: $(length(ver_files)) VTK files")

    # Select representative snapshot indices (t=0, 25%, 50%, 75%, 100%)
    n = length(src_files)
    indices = unique([1, div(n,4), div(n,2), div(3n,4), n])

    # Colours — one per time snapshot; src = solid, ver = dashed
    colours = [:royalblue, :darkorange, :green3, :crimson, :purple]

    p = plot(
        xlabel = "Matric Head h (m)",
        ylabel = "Elevation y (m)",
        title  = "Linear Diffusion — Kernel vs Verification",
        legend = :outertopright,
        size   = (1000, 650),
        margin = 5Plots.mm,
    )

    for (ci, idx) in enumerate(indices)
        src_path = joinpath(src_dir, src_files[idx])
        ver_path = joinpath(ver_dir, ver_files[idx])

        y_s, h_s, t_s = parse_vtk_linear_data(src_path)
        y_v, h_v, t_v = parse_vtk_linear_data(ver_path)

        tstr = @sprintf("t = %.3f s", t_s)
        c    = colours[ci]

        plot!(p, h_s, y_s; label="Kernel — $tstr",       color=c, linewidth=2,
              linestyle=:solid,  marker=:circle, markersize=3, alpha=0.85)
        plot!(p, h_v, y_v; label="Verification — $tstr", color=c, linewidth=2,
              linestyle=:dash,   marker=:diamond, markersize=4, alpha=0.65)
    end

    output_plot = joinpath(src_dir, "src_vs_verification_profiles.png")
    savefig(p, output_plot)
    println("\nPlot saved: $output_plot")
    println("-"^60)
end

if !isinteractive()
    main()
end
