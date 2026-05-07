#______________________________________________________
# ADSIM: Celia et al. (1990) Figure 6B Verification
# Using full kernel approach with GID-generated mesh files
# Author: Paula Sarmiento — May 2026
#______________________________________________________

#______________________________________________________
# Physics
#   Van Genuchten–Mualem SWRC model (Polmann dataset)
#   Mixed-form Richards: ∂θ/∂t = ∇·(K(h)∇h) − ∇·(K(h)ê_g)
#   
# Domain: 100 cm vertical column (1.0 m height, 0.01 m width)
#         Structured mesh: 41 rows × 2 cols (82 nodes, 40 Q4 elements)
#
# IC:     h = -1000 cm (dry, θ ≈ 0.102)
# BC:     Top:    h = -75 cm (wetter, θ ≈ 0.30, t ∈ [0, 24h])
#         Bottom: h = -1000 cm (dry, θ ≈ 0.102, no-flow condition)
#         Sides:  No-flow (implicit from geometry)
#
# Verification Approach:
#   Phase 1: Single run on base mesh (CeliaCol.mesh)
#   Phase 2: Mesh convergence study (multiple mesh resolutions)
#   Phase 3: Mass balance validation (W(t) = ∫θ dΩ)
#
# Files:
#   Mesh:       src/data/CeliaCol.mesh
#   Materials:  src/data/CeliaCol_mat.toml
#   Calc params: src/data/CeliaCol_calc.toml
#______________________________________________________

using LinearAlgebra, SparseArrays, Printf, Plots, Statistics, TOML, Dates, JLD2

const PROJECT_ROOT = dirname(dirname(@__DIR__))   # ADSIM_MID/
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")

# Kernel module includes
include(joinpath(SRC_DIR, "version.jl"))
using .ADSIMVersion: get_version, get_version_string

include(joinpath(SRC_DIR, "swrc_models.jl"))
include(joinpath(SRC_DIR, "read_mesh.jl"))
include(joinpath(SRC_DIR, "read_materials.jl"))
include(joinpath(SRC_DIR, "read_calc_params.jl"))
include(joinpath(SRC_DIR, "initialize_variables.jl"))
include(joinpath(SRC_DIR, "initialize_flows.jl"))
include(joinpath(SRC_DIR, "time_step.jl"))
include(joinpath(SRC_DIR, "shape_functions.jl"))
include(joinpath(SRC_DIR, "write_vtk.jl"))
include(joinpath(SRC_DIR, "implicit_richards_solver.jl"))
include(joinpath(SRC_DIR, "write_checkpoint.jl"))
include(joinpath(SRC_DIR, "read_checkpoint.jl"))

using .ShapeFunctions: initialize_shape_functions!, build_richards_cache
using .WriteVTK: write_vtk_file_water

function build_fixed_time_data(mesh, calc_params, fixed_dt::Float64)
    fixed_dt <= 0.0 && error("time_stepping.fixed_dt must be positive for Celia verification")

    time_settings = calc_params["time_stepping"]
    time_data = TimeStepData()
    time_data.critical_dt = fixed_dt
    time_data.actual_dt = fixed_dt
    time_data.total_time = Float64(time_settings["total_simulation_time"])
    time_data.time_per_step = Float64(time_settings["time_per_step"])
    time_data.courant_number = 1.0
    time_data.h_min = find_minimum_characteristic_length(mesh)
    time_data.num_steps_per_load = ceil(Int, time_data.time_per_step / time_data.actual_dt)
    time_data.num_load_steps = ceil(Int, time_data.total_time / time_data.time_per_step)
    time_data.num_steps = ceil(Int, time_data.total_time / time_data.actual_dt)
    return time_data
end

function format_dt_tag(fixed_dt::Float64)
    rounded_dt = round(Int, fixed_dt)
    return "dt$(rounded_dt)"
end

function run_celia_case(project_name, mesh_file, mat_file, calc_file, fixed_dt, log_print, verification_root)
    log_print("\n   Reading mesh file: $mesh_file")
    mesh = read_mesh_file(mesh_file)
    log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")

    log_print("\n   Reading materials file: $mat_file")
    materials = read_materials_file(mat_file)
    log_print("   ✓ Loaded $(length(materials.soil_dictionary)) soil(s)")

    log_print("\n   Reading calculation parameters: $calc_file")
    calc_params = get_all_calc_params(calc_file)

    compute_K_sat_runtime!(materials, calc_params)
    log_print("   ✓ K_sat computed from intrinsic permeability")

    swrc_used = any(soil.water.swrc_model != "None" for (name, soil) in materials.soils)
    if swrc_used
        log_print("\n   Validating SWRC parameters")
        validate_swrc_parameters(materials)
        log_print("   ✓ SWRC model parameters validated")
    end

    log_print("\n   Initializing simulation variables...")
    zero_variables!(mesh, materials)
    log_print("   ✓ Arrays allocated for $(mesh.num_nodes) nodes")

    log_print("\n   Applying initial and boundary conditions (water-only)...")
    apply_initial_conditions_water!(mesh, materials)
    log_print("   ✓ Applied water IC and BC from mesh file")

    log_print("\n   Building shape functions and Richards cache...")
    initialize_shape_functions!(mesh)
    cache = build_richards_cache(mesh)
    log_print("   ✓ Shape functions initialized")

    log_print("\n   Building fixed time step information...")
    time_data = build_fixed_time_data(mesh, calc_params, fixed_dt)
    log_print(@sprintf("   ✓ Fixed internal time step: %.4g s", time_data.actual_dt))
    log_print(@sprintf("   ✓ Output interval: %.4g s", time_data.time_per_step))
    log_print("   ✓ Internal stepping mode: fixed_dt")
    log_print("   ✓ Number of time steps: $(time_data.num_steps)")

    log_print("\n   Precomputing element water properties...")
    elem_props = precompute_element_water_props(mesh, materials)
    log_print("   ✓ Element properties precomputed for $(length(elem_props)) elements")

    e_g = [calc_params["gravity"]["x_component"], calc_params["gravity"]["y_component"]]
    total_time = calc_params["time_stepping"]["total_simulation_time"]

    log_print("\n   Gravity: e_g = [$(e_g[1]), $(e_g[2])]  ($(all(e_g .== 0) ? "OFF" : "ON"))")
    log_print(@sprintf("   Simulation: T = %.0f s,  Δt = %.4g s,  %d steps",
                       total_time, time_data.actual_dt, time_data.num_steps))

    log_print("\n   Running implicit Richards solver (Picard iteration)...")
    # Solver writes VTK/checkpoints to a relative "output" path in core modules,
    # so run from verification/ to keep all verification artifacts scoped there.
    solver_result = cd(verification_root) do
        implicit_richards_solver(mesh, materials, calc_params, time_data,
                                 project_name, log_print, cache, elem_props, nothing)
    end

    if !(solver_result isa NamedTuple)
        error("Solver returned unexpected result type: $(typeof(solver_result))")
    end

    log_print("   ✓ Solver completed successfully")
    log_print(@sprintf("   ✓ Final time: %.0f s", solver_result.current_time))
    log_print("   ✓ Output counter: $(solver_result.output_counter)")

    log_print("\n   Extracting vertical profile...")
    tol_x = 0.002
    col_nodes = sort(findall(i -> abs(mesh.coordinates[i,1] - 0.005) < tol_x, 1:mesh.num_nodes),
                     by = i -> mesh.coordinates[i,2])
    isempty(col_nodes) && (col_nodes = 1:mesh.num_nodes)

    soil_name = materials.soil_dictionary[1]
    soil = materials.soils[soil_name]
    model = get_water_model(soil.water)
    log_print("   ✓ Extracted $(length(col_nodes)) nodes from center column")

    return (; mesh, materials, calc_params, time_data, solver_result, soil_name, soil, model, col_nodes)
end

#______________________________________________________
# Main execution with error handling (kernel pattern)
#______________________________________________________
function main()
    # Check for --version flag
    if length(ARGS) > 0 && (ARGS[1] == "--version" || ARGS[1] == "-v")
        println(get_version_string())
        exit(0)
    end
    
    if length(ARGS) < 1
        println("Error: No project name provided")
        println("Usage: julia run_celia_verification_kernel.jl <project_name>")
        println("       julia run_celia_verification_kernel.jl --version")
        println("Example: julia run_celia_verification_kernel.jl CeliaFig6B_test")
        exit(1)
    end

    project_name = ARGS[1]
    data_dir     = joinpath(SRC_DIR, "data")
    verification_root = joinpath(PROJECT_ROOT, "verification")
    output_dir   = joinpath(verification_root, "output")
    
    base_mesh_case = "CeliaCol_fine"
    base_material_case = "CeliaCol"
    base_calc_case = "CeliaCol"
    
    isdir(output_dir) || mkpath(output_dir)
    
    # Setup log file (kernel pattern)
    log_file_path = joinpath(output_dir, "$(project_name)_celia_verification.log")
    # Only remove old log if not in use (avoid EBUSY errors)
    try
        isfile(log_file_path) && rm(log_file_path)
    catch
        # If file is locked, just append
    end
    log_file = open(log_file_path, "w")
    
    # Custom print function that writes to both console and log file
    function log_print(msg::String)
        println(msg)
        println(log_file, msg)
        flush(log_file)
    end

    try
        start_time = now()
        
        log_print("="^64)
        log_print("ADSIM: Celia et al. (1990) Figure 6B Verification")
        log_print("Version: $(get_version())")
        log_print("Project: $(project_name)")
        log_print("="^64)
        
        # ══════════════════════════════════════════════════════════════════════════════
        # PHASE 1: Single run on base mesh
        # ══════════════════════════════════════════════════════════════════════════════
        
        log_print("\n[1/3] PHASE 1: Single run on base mesh")
        
        mesh_file = joinpath(data_dir, "$(base_mesh_case).mesh")
        mat_file  = joinpath(data_dir, "$(base_material_case)_mat.toml")
        calc_file = joinpath(data_dir, "$(base_calc_case)_calc.toml")
        
        # Verify all files exist
        if !isfile(mesh_file) || !isfile(mat_file) || !isfile(calc_file)
            error("Required files not found:\n  Mesh: $(isfile(mesh_file) ? "✓" : "✗") $mesh_file\n  Materials: $(isfile(mat_file) ? "✓" : "✗") $mat_file\n  Calc params: $(isfile(calc_file) ? "✓" : "✗") $calc_file")
        end
        
        base_calc_params = get_all_calc_params(calc_file)
        sweep_values_raw = get(base_calc_params["time_stepping"], "fixed_dt_sweep", Any[])
        fixed_dt_values = isempty(sweep_values_raw) ? [Float64(get(base_calc_params["time_stepping"], "fixed_dt", 144.0))] : Float64[Float64(v) for v in sweep_values_raw]

        log_print("\n   Benchmark sweep: $(join(string.(round.(Int, fixed_dt_values)), ", ")) s")

        case_results = NamedTuple[]
        for (case_idx, fixed_dt) in enumerate(fixed_dt_values)
            case_project = "CeliaCol_$(format_dt_tag(fixed_dt))"
            log_print("\n   --- Benchmark case $(case_idx)/$(length(fixed_dt_values)): $(case_project) ---")
            push!(case_results, run_celia_case(case_project, mesh_file, mat_file, calc_file, fixed_dt, log_print, verification_root))
        end

        last_case = case_results[end]
        mesh = last_case.mesh
        calc_params = last_case.calc_params
        T_total = last_case.time_data.total_time
        soil = last_case.soil
        
        # ══════════════════════════════════════════════════════════════════════════════
        # PHASE 2: Mesh convergence study (placeholder)
        # ══════════════════════════════════════════════════════════════════════════════
        
        log_print("\n[2/3] PHASE 2: Mesh convergence study")
        log_print("   Note: Requires mesh files at different resolutions")
        log_print("   Generate in GID and save as:")
        log_print("     - CeliaCol_coarse.mesh (ny=20 rows)")
        log_print("     - CeliaCol_fine.mesh (ny=80 rows)")
        
        # ══════════════════════════════════════════════════════════════════════════════
        # PHASE 3: Summary and verification
        # ══════════════════════════════════════════════════════════════════════════════
        
        log_print("\n[3/3] PHASE 3: Verification summary")
        log_print("   Simulation complete!")
        log_print("   Output directory: $output_dir")
        log_print("   Sweep cases: $(join(["CeliaCol_$(format_dt_tag(dt))" for dt in fixed_dt_values], ", "))")

        # ── Picard summary table ───────────────────────────────────────────────
        log_print("\n" * "─"^62)
        log_print("  PICARD ITERATION SUMMARY")
        log_print("─"^62)
        log_print(@sprintf("  %-20s  %6s  %6s  %6s  %8s", "Case", "Min", "Mean", "Max", "Steps"))
        log_print("  " * "─"^58)
        for (case, dt) in zip(case_results, fixed_dt_values)
            ph = case.solver_result.picard_history
            if !isempty(ph)
                log_print(@sprintf("  %-20s  %6d  %6.1f  %6d  %8d",
                    "CeliaCol_$(format_dt_tag(dt))",
                    minimum(ph), sum(ph)/length(ph), maximum(ph), length(ph)))
            end
        end
        log_print("─"^62)
        
        # Print timing information
        end_time = now()
        total_time = (end_time - start_time).value / 1000.0
        log_print(@sprintf("   Total verification time: %.2f seconds", total_time))
        
        log_print("\n" * "="^64)
        log_print("✓ Celia et al. (1990) Figure 6B Verification Completed")
        log_print("="^64)
        log_print("   Mesh:")
        log_print("     - Nodes: $(mesh.num_nodes)")
        log_print("     - Elements: $(mesh.num_elements)")
        log_print("   Simulation:")
        log_print("     - Total time: $(T_total) s ($(T_total/3600) hours)")
        for case in case_results
            log_print("     - $(round(Int, case.time_data.actual_dt)) s: $(case.time_data.num_steps) steps")
        end
        log_print("   SWRC Model: $(soil.water.swrc_model)")
        log_print("="^64)
        
    catch e
        # Error handling (kernel pattern)
        log_print("\n" * "="^64)
        log_print("FATAL ERROR OCCURRED")
        log_print("="^64)
        log_print("\nError Type: $(typeof(e))")
        log_print("\nError Message:")
        log_print(sprint(showerror, e))
        log_print("\nStack Trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(log_file, exc, bt)
            println(log_file)
            flush(log_file)
        end
        log_print("\n" * "="^64)
        log_print("Program terminated due to error")
        log_print("="^64)
        
        println("\n" * "="^64)
        println("FATAL ERROR - Check log file for details: $(log_file_path)")
        println("="^64)
        
        close(log_file)
        exit(1)
    
    finally
        # Ensure log file is closed
        if isopen(log_file)
            close(log_file)
        end
    end
end

#______________________________________________________
# Entry point
#______________________________________________________

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
