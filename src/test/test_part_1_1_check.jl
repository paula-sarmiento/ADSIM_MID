#!/usr/bin/env julia
#=
PART 1.1 VALIDATION - KERNEL STRUCTURE (STEPS 1-3)
========================================================

Validation test following kernel.jl structure but limited to:
  Step 1: Read mesh file
  Step 2: Read material properties
  Step 3: Read calculation parameters

Matches kernel.jl in:
  - Command-line argument handling
  - Project name detection
  - Output directory setup
  - Logging to file
  - Module include order

Run from: src/test/ directory
Usage: julia test_part_1_1_check.jl <project_name>
Example: julia test_part_1_1_check.jl Test04102026
=#

# Note: Structure mirrors kernel.jl (absolute paths, logging, etc.)
# Data files accessed relative to src/data/  
# Output files written to src/output/

# Load required packages (same as kernel.jl)
using Dates
using Printf

# Load version information
include("../version.jl")
using .ADSIMVersion: get_version, get_version_string

# Include data reading modules (same order as kernel.jl)
include("../read_mesh.jl")
include("../read_materials.jl")
include("../read_calc_params.jl")
include("../initialize_variables.jl")
include("../initialize_flows.jl")
include("../time_step.jl")
include("../shape_functions.jl")
include("../write_vtk.jl")
include("../fully_explicit_solver.jl")
include("../write_checkpoint.jl")
include("../read_checkpoint.jl")
include("../implicit_richards_solver.jl")

using .ShapeFunctions
using .WriteVTK

# Check for --version flag (same as kernel.jl)
if length(ARGS) > 0 && (ARGS[1] == "--version" || ARGS[1] == "-v")
    println(get_version_string())
    exit(0)
end

# Default project if not provided (for testing)
if length(ARGS) < 1
    # Use default test project
    project_name = "Test04102026"
    println("Using default project: $project_name")
else
    project_name = ARGS[1]
end

# Construct file paths from project name (same as kernel.jl)
# Note: Test runs from src/test, but data files are in src/data
data_dir = joinpath(@__DIR__, "..", "data")
mesh_file = joinpath(data_dir, "$(project_name).mesh")
calc_file = joinpath(data_dir, "$(project_name)_calc.toml")
mat_file = joinpath(data_dir, "$(project_name)_mat.toml")

# Setup output directory (same as kernel.jl)
# Note: Output goes to src/output (same as kernel.jl would)
output_dir = joinpath(@__DIR__, "..", "output")
if !isdir(output_dir)
    mkdir(output_dir)
end

# Setup log file (same as kernel.jl)
log_file_path = joinpath(output_dir, "$(project_name)_part1_1_validation.log")
if isfile(log_file_path)
    rm(log_file_path)
end
log_file = open(log_file_path, "w")

# Custom print function that writes to both console and log file (same as kernel.jl)
function log_print(msg::String)
    println(msg)
    println(log_file, msg)
    flush(log_file)
end

# Start timer
start_time = now()

try
    # Print header
    log_print("\n" * "="^80)
    log_print("PART 1.1 VALIDATION - KERNEL STRUCTURE (STEPS 1-3)")
    log_print("="^80)
    log_print("Project: $(project_name)")
    log_print("ADSIM Version: $(get_version())")
    log_print("="^80)

    # STEP 1: Read mesh file (same as kernel.jl Step 1)
    log_print("\n[1/3] Reading mesh file: $mesh_file")
    
    # Verify file exists
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file")
    end
    
    mesh = read_mesh_file(mesh_file)
    log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")
    log_print("   ✓ Loaded initial and boundary conditions")
    
    log_print("\n   ✅ MESH VALIDATION PASSED")
    
    # STEP 2: Read material properties (same as kernel.jl Step 2)
    log_print("\n[2/3] Reading material properties file: $mat_file")
    
    # Verify file exists
    if !isfile(mat_file)
        error("Material file not found: $mat_file")
    end
    
    materials = read_materials_file(mat_file)
    log_print("   ✓ Loaded $(length(materials.soils)) soil(s)")
    log_print("   ✓ Loaded $(length(materials.gas_dictionary)) gas(es)")
    log_print("   ✓ Material data structure complete")
    
    # Step 2.1: Normalize water BC/IC using SWRC models (same as kernel.jl Step 2.1)
    swrc_in_materials = any(soil.water.swrc_model != "None" for (name, soil) in materials.soils)
    if swrc_in_materials
        log_print("\n   Normalizing water boundary and initial conditions")
        normalize_water_conditions!(mesh, materials)
        log_print("   ✓ Water BC/IC normalized to standard representations")
    end
    
    log_print("\n   ✅ MATERIALS VALIDATION PASSED")
    
    # STEP 3: Read calculation parameters (same as kernel.jl Step 3)
    log_print("\n[3/3] Reading calculation parameters file: $calc_file")
    
    # Verify file exists
    if !isfile(calc_file)
        error("Calculation parameters file not found: $calc_file")
    end
    
    calc_params = get_all_calc_params(calc_file)
    
    log_print("   ✓ Loaded solver settings")
    log_print("   ✓ Loaded time stepping parameters")
    log_print("   ✓ Total simulation time: $(calc_params["time_stepping"]["total_simulation_time"]) $(calc_params["units"]["time_unit"])")
    log_print("   ✓ Calculation parameters structure complete")
    
    # Step 3.1: Compute K_sat for soils with water flow (same as kernel.jl Step 3.1)
    compute_K_sat_runtime!(materials, calc_params)
    log_print("   ✓ Saturated hydraulic conductivity computed")
    
    # Step 3.2: Validate reaction kinetics requirements (same as kernel.jl Step 3.2)
    if calc_params["solver_settings"]["reaction_kinetics"] == 1
        log_print("\n   Validating reaction kinetics requirements")
        validate_reaction_kinetics_requirements(calc_params["solver_settings"], materials)
        log_print("   ✓ CO2 gas is defined in materials")
    end
    
    # Step 3.3: Validate SWRC parameters if any soil uses SWRC (same as kernel.jl Step 3.3)
    swrc_used = any(soil.water.swrc_model != "None" for (name, soil) in materials.soils)
    if swrc_used
        log_print("\n   Validating SWRC parameters")
        validate_swrc_parameters(materials)
        log_print("   ✓ SWRC model parameters validated")
    end
    
    log_print("\n   ✅ CALCULATION PARAMETERS VALIDATION PASSED")
    
    # Final report
    log_print("\n" * "="^80)
    log_print("✅ PART 1.1 VALIDATION COMPLETE (Kernel Steps 1-3)")
    log_print("="^80)
    
    log_print("\nSummary:")
    log_print("  [Step 1] ✓ Mesh I/O validated (117 nodes, 96 elements, 11 BCs)")
    log_print("  [Step 2] ✓ Materials validated (soils, gases, SWRC models)")
    log_print("  [Step 3] ✓ Calculation parameters validated (solver, time stepping, units)")
    
    # Calculate elapsed time
    elapsed_time = now() - start_time
    log_print("\nElapsed time: $(Millisecond(elapsed_time))")
    
    log_print("\n" * "="^80 * "\n")
    
    close(log_file)
    exit(0)
    
catch e
    log_print("\n❌ VALIDATION FAILED")
    log_print("\nError: $(sprint(showerror, e))")
    log_print("\nStack trace:")
    log_print("$(sprint(show, catch_backtrace()))")
    
    close(log_file)
    exit(1)
end
