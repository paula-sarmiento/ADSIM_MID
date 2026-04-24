#!/usr/bin/env julia
"""
Inspect benchmark VTK outputs to diagnose profile extraction issues
"""

using Printf

include("swrc_models.jl")

function read_vtk_points_field(filename)
    open(filename, "r") do f
        lines = readlines(f)
        
        # Find POINTS
        pt_idx = findfirst(l -> startswith(l, "POINTS"), lines)
        if pt_idx === nothing
            error("No POINTS found in $filename")
        end
        n_pts = parse(Int, split(lines[pt_idx])[2])
        
        points = zeros(Float64, n_pts, 3)
        idx = 0
        for i in (pt_idx+1):length(lines)
            coords = split(lines[i])
            for c in coords
                idx += 1
                if idx > n_pts break end
                points[idx, mod1(idx, 3)] = parse(Float64, c)
            end
            idx > n_pts && break
        end
        
        # Find POINT_DATA scalar
        pd_idx = findfirst(l -> startswith(l, "POINT_DATA"), lines)
        sc_idx = findnext(l -> startswith(l, "SCALARS"), lines, pd_idx)
        lk_idx = findnext(l -> startswith(l, "LOOKUP_TABLE"), lines, sc_idx)
        
        values = Float64[]
        for i in (lk_idx+1):length(lines)
            for v in split(lines[i])
                try
                    push!(values, parse(Float64, v))
                    length(values) == n_pts && break
                catch end
            end
            length(values) == n_pts && break
        end
        
        return points[:, 1:2], values
    end
end

vg_celia = VanGenuchten(
    theta_s=0.368, theta_r=0.102,
    alpha=0.335, n_vg=2.0,
    K_s=9.4e-12
)

println("="^70)
println("BENCHMARK VTK INSPECTION")
println("="^70)

file = "../output/CeliaCol_dt144_water_000024.vtk"
println("\nReading: $file")

points, h_vals = read_vtk_points_field(file)

println("\nPoints shape: $(size(points))")
println("h_vals length: $(length(h_vals))")

# Extract column at x≈0
x_tol = 0.01
col_idx = findall(i -> abs(points[i,1]) < x_tol, 1:size(points,1))
y = points[col_idx, 2]
h = h_vals[col_idx]

perm = sortperm(y)
y_sorted = y[perm]
h_sorted = h[perm]

println("\nVertical profile (x ≈ 0):")
println("  y[m]    |    h[m]    |    θ[-]    |  Depth[cm]")
println("  " * "-"^50)

theta_vals = [theta(vg_celia, hi) for hi in h_sorted]

for i in 1:length(y_sorted)
    depth_cm = (1.0 - y_sorted[i]) * 100
    @printf("  %.4f   |   %.6f   |   %.6f   |  %.2f\n", 
            y_sorted[i], h_sorted[i], theta_vals[i], depth_cm)
end

# Check metrics
theta_threshold = vg_celia.theta_r + 0.05 * (vg_celia.theta_s - vg_celia.theta_r)
wet_idx = findall(ϑ -> ϑ > theta_threshold, theta_vals)

println("\nMetrics:")
println("  theta_threshold = $theta_threshold")
println("  Wet nodes (θ > threshold): $(length(wet_idx))")

if !isempty(wet_idx)
    front_depth = maximum(1.0 .- y_sorted[wet_idx] ./ 1.0) * 100
    println("  Front penetration: $front_depth cm")
end

top_20_idx = findall(i -> y_sorted[i] > 0.8, 1:length(y_sorted))
if !isempty(top_20_idx)
    S_e_top = mean([(theta_vals[i] - vg_celia.theta_r) / (vg_celia.theta_s - vg_celia.theta_r) 
                    for i in top_20_idx])
    println("  S_e (top 20cm): $S_e_top")
end

println("\n" * "="^70)
