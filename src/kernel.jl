#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.x.x
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Kernel functions for ADSIM
# Calls procedures in order
#______________________________________________________

# Include data reading modules
include("read_mesh.jl")

#______________________________________________________
# Main execution script
#______________________________________________________

#starts timer
using Dates
start_time = now()

# Check if mesh file is provided as command-line argument
if length(ARGS) < 1
    println("Error: No mesh file provided")
    println("Usage: julia kernel.jl <mesh_file_path>")
    println("Example: julia kernel.jl ../data/problem.mesh")
    exit(1)
end

mesh_file = ARGS[1]

# Verify file exists
if !isfile(mesh_file)
    println("Error: Mesh file not found: ", mesh_file)
    exit(1)
end

println("="^64)
println("ADSIM: Advection-Diffusion for Soil Improvement and Modification")
println("="^64)

# Step 1: Read mesh data
println("\n[1/N] Reading mesh file: ", mesh_file)
mesh = read_mesh_file(mesh_file)
println("   ✓ Loaded ", mesh.num_nodes, " nodes and ", mesh.num_elements, " elements")
println("   ✓ Loaded initial and boundary conditions")

# Step 2: Additional processing steps (to be implemented)
# - Material properties
# - Assembly of system matrices
# - Time stepping
# - Solution
# - Post-processing


#Print total run Time
end_time = now()
total_time = end_time - start_time
println("Total run time: ", total_time)

println("\n", "="^64)
println("Mesh reading completed successfully")
println("="^64)

