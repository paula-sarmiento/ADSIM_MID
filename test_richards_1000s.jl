"""
Richards Solver Test: t = 1000 s, Δt = 10 s
Tests water infiltration over 1000 seconds with 10-second time steps
"""

println("="^80)
println("RICHARDS SOLVER TEST: t = 1000 s, Δt = 10 s")
println("="^80)

# Use the existing ADSIM kernel with modified parameters
# We'll create a test data file and run through the main kernel

# First, create a modified calc file for this test
calc_file = "src/data/Test1000s_calc.toml"

# Read original calc params
orig_file = "src/data/Test04102026_calc.toml"
calc_content = read(orig_file, String)

# Modify parameters: t=1000s, dt=10s (already correct in file)
modified_content = replace(calc_content, "total_simulation_time = 1000" => "total_simulation_time = 1000")
modified_content = replace(modified_content, "time_per_step = 10" => "time_per_step = 10")

# Write modified calc file
write(calc_file, modified_content)
println("\n✓ Created test configuration: $calc_file")
println("  - Total simulation time: 1000 s")
println("  - Time per step: 10 s")
println("  - Expected steps: 100")

# Copy mesh and materials files for test project
cp("src/data/Test04102026.mesh", "src/data/Test1000s.mesh", force=true)
cp("src/data/Test04102026_mat.toml", "src/data/Test1000s_mat.toml", force=true)
println("✓ Copied mesh and material files")

# Set ARGS to run with Test1000s
empty!(ARGS)
push!(ARGS, "Test1000s")

# Now run the main kernel
println("\n" * "="^80)
println("Running Richards solver with Test1000s configuration...")
println("="^80 * "\n")

# Include and run the kernel
include("src/kernel.jl")
main()
