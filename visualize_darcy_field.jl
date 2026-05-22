#!/usr/bin/env julia
#=
Compute and visualize Darcy velocity field using compute_darcy_nodes_extrapolation
for the constant soil verification test case.

This script:
1. Loads the checkpoint data from LinVerif_diffusion simulation
2. Extracts mesh, pressure head field, and element properties
3. Computes nodal Darcy velocities using Gauss-point extrapolation
4. Visualizes the velocity field with a quiver plot
=#

const PROJECT_ROOT = @__DIR__
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")

include("src/version.jl")
include("src/read_mesh.jl")
include("src/read_materials.jl")
include("src/read_calc_params.jl")
include("src/swrc_models.jl")
include("src/shape_functions.jl")
include("src/initialize_variables.jl")
include("src/write_vtk.jl")
include("src/read_checkpoint.jl")

using JLD2
using LinearAlgebra
using Statistics

# Try to use Plots.jl with PyPlot backend for visualization
has_plots = false
try
    using Plots
    pyplot()
    global has_plots = true
catch
    println("⚠ Warning: Plots.jl not available, will save data to CSV instead")
end

println("\n" * "="^80)
println("DARCY VELOCITY VISUALIZATION: LinVerif_diffusion (Constant Soil)")
println("="^80)

# ──────────────────────────────────────────────────────────────────────────────
# Load checkpoint data and reconstruct mesh/materials
# ──────────────────────────────────────────────────────────────────────────────

println("\n[SETUP] Loading checkpoint data...")

checkpoint_file = "output/LinVerif_diffusion_stage1.jld2"
if !isfile(checkpoint_file)
    error("Checkpoint file not found: $checkpoint_file")
end

# Load checkpoint state
checkpoint_data = jldopen(checkpoint_file, "r") do f
    Dict(key => f[key] for key in keys(f))
end
h_global = checkpoint_data["P"]  # Pressure head stored as "P"
Nnodes = checkpoint_data["Nnodes"]
Nelements = checkpoint_data["Nelements"]

println("✓ Checkpoint loaded:")
println("  - Nodes: $Nnodes, Elements: $Nelements")
println("  - Pressure head range: [$(@sprintf "%.3f" minimum(h_global)), $(@sprintf "%.3f" maximum(h_global))] m")

# Load mesh from mesh file (must exist for the project)
println("\n[SETUP] Loading mesh from data file...")

mesh_file = joinpath(SRC_DIR, "data", "LinVerif_diffusion.mesh")
if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

mesh = read_mesh_file(mesh_file)
println("✓ Mesh loaded from: $mesh_file")
println("  - Nodes: $(mesh.num_nodes), Elements: $(mesh.num_elements)")

# Verify mesh matches checkpoint
if mesh.num_nodes != Nnodes || mesh.num_elements != Nelements
    error("Mesh dimensions don't match checkpoint: " *
          "mesh($Nnodes, $Nelements) vs checkpoint($Nnodes, $Nelements)")
end

# Extract pressure head vector (h = P)
h = h_global[1:Nnodes]

# Load materials
println("\n[SETUP] Loading materials...")
mat_file = joinpath(SRC_DIR, "data", "LinVerif_diffusion_mat.toml")
if !isfile(mat_file)
    error("Materials file not found: $mat_file")
end

materials = read_materials_file(mat_file)
println("✓ Materials loaded")

# Load calculation parameters
println("\n[SETUP] Loading calculation parameters...")
calc_file = joinpath(SRC_DIR, "data", "LinVerif_diffusion_calc.toml")
if !isfile(calc_file)
    error("Calculation parameters file not found: $calc_file")
end

calc_params = get_all_calc_params(calc_file)
println("✓ Calculation parameters loaded")

# Initialize SWRC models
println("\n[SETUP] Initializing SWRC models...")
compute_K_sat_runtime!(materials, calc_params)
println("✓ SWRC models initialized")

# DEBUG: Print actual K_val being used
println("\n[DEBUG] Material properties:")
for (soil_name, soil) in materials.soils
    swrc_model = soil.water.swrc_model
    K_sat = soil.water.K_sat
    println("  Soil: $soil_name | SWRC: $swrc_model | K_sat: $K_sat m/s")
    
    if swrc_model == "ConstantSoil"
        model = soil.water.swrc_model_instance
        println("    ConstantSoil K_val: $(model.K_val) m/s")
        println("    ConstantSoil h_min: $(model.h_min) m")
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Build element properties from materials
# ──────────────────────────────────────────────────────────────────────────────

println("\n[SETUP] Building element properties...")

elem_props = Vector{ElementWaterProps}(undef, mesh.num_elements)
for e in 1:mesh.num_elements
    elem_props[e] = get_element_water_props(mesh, materials, e)
end

# Get gravity from calculation parameters
gravity_x = calc_params["gravity"]["x_component"]
gravity_y = calc_params["gravity"]["y_component"]
e_g = [gravity_x, gravity_y]  # Gravity vector from parameters

println("✓ Element properties reconstructed:")
println("  - Number of elements: $(mesh.num_elements)")
println("  - Gravity vector: $e_g (from calculation parameters)")

# DEBUG: Check a sample element's K value
if mesh.num_elements > 0
    sample_model = elem_props[1].model
    sample_h = 0.0
    K_sample = K_h_y(sample_model, sample_h)
    println("  - Sample element K_h_y(h=0): $K_sample m/s")
end

# ──────────────────────────────────────────────────────────────────────────────
# Initialize shape functions and cache
# ──────────────────────────────────────────────────────────────────────────────

println("\n[SETUP] Initializing shape functions and cache...")

ShapeFunctions.initialize_shape_functions!(mesh)
cache = ShapeFunctions.build_richards_cache(mesh)

println("✓ Shape functions and cache initialized")

# ──────────────────────────────────────────────────────────────────────────────
# Compute Darcy velocities using extrapolation method
# ──────────────────────────────────────────────────────────────────────────────

println("\n[COMPUTING] Darcy velocity field via Gauss-point extrapolation...")

q_x, q_y = compute_darcy_nodes_extrapolation(h, mesh, elem_props, cache, e_g)

vel_mag = sqrt.(q_x.^2 .+ q_y.^2)

println("✓ Velocities computed:")
println("  - Velocity magnitude range: [$(@sprintf "%.6f" minimum(vel_mag)), $(@sprintf "%.6f" maximum(vel_mag))] m/s")
println("  - Mean velocity: $(@sprintf "%.6f" mean(vel_mag)) m/s")
println("  - Std deviation: $(@sprintf "%.6f" std(vel_mag)) m/s")

# DEBUG: Analyze gradient and theoretical velocity
println("\n[DEBUG] Velocity analysis:")
model = elem_props[1].model
if isa(model, ConstantSoil)
    K_val = model.K_val
    println("  K_val from ConstantSoil: $K_val m/s")
    
    # Estimate theoretical max velocity
    local grad_h_max = 0.0
    for e in 1:mesh.num_elements
        h_e = [h[mesh.elements[e, a]] for a in 1:4]
        for p in 1:4
            grad_h = [sum(cache.Bp[e, p, d, a] * h_e[a] for a in 1:4) for d in 1:2]
            grad_mag = sqrt(grad_h[1]^2 + grad_h[2]^2)
            grad_h_max = max(grad_h_max, grad_mag)
        end
    end
    
    v_theory = K_val * grad_h_max
    v_obs = maximum(vel_mag)
    
    println("  Max |∇h|: $(@sprintf "%.6f" grad_h_max) m/m")
    println("  Theoretical max q: $(@sprintf "%.6f" v_theory) m/s (K × max|∇h|)")
    println("  Observed max q: $(@sprintf "%.6f" v_obs) m/s")
    if v_theory > 0
        println("  Ratio (obs/theory): $(@sprintf "%.2f" v_obs/v_theory)×")
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Analyze velocity field structure
# ──────────────────────────────────────────────────────────────────────────────

println("\n[ANALYSIS] Velocity field structure:")

# Count nodes with negligible vs significant flow
tol_negligible = 1e-6
n_negligible = count(vel_mag .< tol_negligible)
n_significant = mesh.num_nodes - n_negligible

println("  - Nodes with negligible flow (<1e-6): $n_negligible")
println("  - Nodes with significant flow: $n_significant")

# Find max/min velocity nodes
idx_max = argmax(vel_mag)
idx_min = argmin(vel_mag)

println("  - Max velocity: $(@sprintf "%.6f" vel_mag[idx_max]) m/s at node $idx_max")
println("    Position: $(@sprintf "[%.4f, %.4f]" mesh.coordinates[idx_max, 1] mesh.coordinates[idx_max, 2])")
println("    Components: (q_x=$(@sprintf "%.6f" q_x[idx_max]), q_y=$(@sprintf "%.6f" q_y[idx_max]))")

println("  - Min velocity: $(@sprintf "%.6f" vel_mag[idx_min]) m/s at node $idx_min")

# ──────────────────────────────────────────────────────────────────────────────
# Save velocity field data to CSV for visualization
# ──────────────────────────────────────────────────────────────────────────────

println("\n[OUTPUT] Saving velocity field to CSV...")

using DelimitedFiles

output_file = "output/darcy_velocity_field.csv"

open(output_file, "w") do io
    # Write header
    write(io, "node_id,x,y,h,q_x,q_y,velocity_magnitude\n")
    # Write data rows
    for i in 1:mesh.num_nodes
        write(io, "$i,$(@sprintf "%.6f" mesh.coordinates[i,1]),$(@sprintf "%.6f" mesh.coordinates[i,2]),$(@sprintf "%.6f" h[i]),$(@sprintf "%.6f" q_x[i]),$(@sprintf "%.6f" q_y[i]),$(@sprintf "%.6f" vel_mag[i])\n")
    end
end

println("✓ Velocity field saved to: $output_file")

# ──────────────────────────────────────────────────────────────────────────────
# Create visualization if Plots.jl is available
# ──────────────────────────────────────────────────────────────────────────────

if has_plots
    println("\n[VISUALIZATION] Creating quiver plot...")
    
    # Create figure with multiple subplots
    fig = plot(layout=(2,2), size=(1200, 1000))
    
    # Plot 1: Pressure head contour
    scatter!(fig[1], mesh.coordinates[:, 1], mesh.coordinates[:, 2], 
             markersize=5, markerstrokewidth=0, c=h, 
             title="Pressure Head (h) [m]", xlabel="x [m]", ylabel="y [m]",
             legend=false, colorbar=true, palette=:viridis)
    
    # Plot 2: Velocity magnitude
    scatter!(fig[2], mesh.coordinates[:, 1], mesh.coordinates[:, 2],
             markersize=5, markerstrokewidth=0, c=vel_mag,
             title="Velocity Magnitude [m/s]", xlabel="x [m]", ylabel="y [m]",
             legend=false, colorbar=true, clim=(0, maximum(vel_mag)*0.9), palette=:hot)
    
    # Plot 3: Velocity field (quiver plot)
    # Subsample if too many nodes for clarity
    n_skip = max(1, div(mesh.num_nodes, 30))
    idx_subsample = 1:n_skip:mesh.num_nodes
    
    quiver!(fig[3], mesh.coordinates[idx_subsample, 1], mesh.coordinates[idx_subsample, 2],
            quiver=(q_x[idx_subsample], q_y[idx_subsample]),
            title="Velocity Field (Gauss-Point Extrapolation)", xlabel="x [m]", ylabel="y [m]",
            arrow=arrow(:closed, 1), linewidth=1.5, color=:darkblue, legend=false,
            aspect_ratio=:equal)
    
    # Plot 4: Velocity histogram
    histogram!(fig[4], vel_mag, bins=min(30, mesh.num_nodes÷2), 
               title="Velocity Magnitude Distribution",
               xlabel="Velocity [m/s]", ylabel="Count", legend=false, color=:steelblue)
    
    # Save figure
    plot_file = "output/darcy_velocity_field.png"
    savefig(fig, plot_file)
    println("✓ Visualization saved to: $plot_file")
else
    println("\n[INFO] Plots.jl not available. Velocity data saved to CSV:")
    println("       File: $output_file")
    println("       Use external tool (Python, Excel, etc.) to visualize")
end

# ──────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ──────────────────────────────────────────────────────────────────────────────

println("\n" * "="^80)
println("SUMMARY STATISTICS")
println("="^80)

println("\nPressure Head Field (h):")
println("  Range:     [$(@sprintf "%.4f" minimum(h)), $(@sprintf "%.4f" maximum(h))] m")
println("  Mean:      $(@sprintf "%.4f" mean(h)) m")
println("  Std Dev:   $(@sprintf "%.4f" std(h)) m")

println("\nVelocity Field (Darcy - Extrapolation Method):")
println("  |q| range: [$(@sprintf "%.6f" minimum(vel_mag)), $(@sprintf "%.6f" maximum(vel_mag))] m/s")
println("  |q| mean:  $(@sprintf "%.6f" mean(vel_mag)) m/s")
println("  |q| std:   $(@sprintf "%.6f" std(vel_mag)) m/s")

println("\nVelocity Components:")
println("  q_x range: [$(@sprintf "%.6f" minimum(q_x)), $(@sprintf "%.6f" maximum(q_x))] m/s")
println("  q_x mean:  $(@sprintf "%.6f" mean(q_x)) m/s")
println("  q_y range: [$(@sprintf "%.6f" minimum(q_y)), $(@sprintf "%.6f" maximum(q_y))] m/s")
println("  q_y mean:  $(@sprintf "%.6f" mean(q_y)) m/s")

println("\n" * "="^80)
println("✓ COMPLETE - Darcy velocity field computed and visualized!")
println("="^80 * "\n")
