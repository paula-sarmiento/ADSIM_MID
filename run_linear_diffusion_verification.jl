#______________________________________________________
# ADSIM: Linear Diffusion Verification
# Pure-diffusion test — Richards equation with ConstantSoil
# reduces to the heat equation (D = K/C = const)
# Author: Paula Sarmiento — April 2026
#______________________________________________________

#______________________________________________________
# Physics
#   ConstantSoil: C = (θ_s − θ_r)/|h_min|,  K = K_val
#   gravity OFF  →  ∂h/∂t = D ∇²h,  D = K/C
#
# Domain: 1 m × 1 m  (1-D in y by symmetry)
# IC:     h = −0.5 m (uniform)
# BC:     h(y=0) = −1.0 m,  h(y=1) = 0.0 m,  sides no-flow
#
# Analytical solution (Fourier sine series, N_modes terms):
#   h_ss(y) = h_bot + (h_top − h_bot)·y/L
#   Bₙ = (1 + (−1)ⁿ) / (nπ)
#   h(y,t) = h_ss(y) + Σ Bₙ sin(nπy/L) exp(−D(nπ/L)²t)
#______________________________________________________

using LinearAlgebra, SparseArrays, Printf, Plots, Statistics, TOML, Dates

const PROJECT_ROOT = @__DIR__
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")

include(joinpath(SRC_DIR, "version.jl"))
using .ADSIMVersion: get_version

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

#______________________________________________________
# Analytical solution (Fourier sine series)
#______________________________________________________
function h_analytical(y, t, D; L=1.0, h_bot=-1.0, h_top=0.0, N_modes=500)
    h_ss = h_bot + (h_top - h_bot) * y / L
    s = sum((1.0 + (-1.0)^n) / (n * π) * sin(n * π * y / L) *
            exp(-D * (n * π / L)^2 * t) for n in 1:N_modes)
    return h_ss + s
end

#______________________________________________________
# Main execution
#______________________________________________________
function main()
    if length(ARGS) < 1
        println("Error: No project name provided")
        println("Usage: julia run_linear_diffusion_verification.jl <project_name>")
        println("Example: julia run_linear_diffusion_verification.jl LinVerif_diffusion")
        exit(1)
    end

    project_name = ARGS[1]
    data_dir     = joinpath(SRC_DIR, "data")
    output_dir   = joinpath(PROJECT_ROOT, "output")
    mesh_file    = joinpath(data_dir, "$(project_name).mesh")
    mat_file     = joinpath(data_dir, "$(project_name)_mat.toml")
    calc_file    = joinpath(data_dir, "$(project_name)_calc.toml")

    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        exit(1)
    end

    isdir(output_dir) || mkdir(output_dir)

    # ── Determine current stage by checking for existing checkpoint files (kernel pattern) ──
    if isdir(output_dir)
        all_files = readdir(output_dir)
        checkpoint_files = filter(f -> occursin(r"_stage\d+\.jld2$", f) && startswith(f, project_name), all_files)
        
        if isempty(checkpoint_files)
            current_stage = 1
        else
            stage_numbers = Int[]
            for filename in checkpoint_files
                m = match(r"_stage(\d+)\.jld2$", filename)
                if m !== nothing
                    push!(stage_numbers, parse(Int, m.captures[1]))
                end
            end
            current_stage = isempty(stage_numbers) ? 1 : maximum(stage_numbers) + 1
        end
    else
        current_stage = 1
    end

    log_file_path = joinpath(output_dir, "$(project_name)_stage$(current_stage).log")
    isfile(log_file_path) && rm(log_file_path)
    log_file = open(log_file_path, "w")

    function log_print(msg::String)
        println(msg)
        println(log_file, msg)
        flush(log_file)
    end

    try
        start_time = now()

        log_print("="^64)
        log_print("ADSIM: Linear Diffusion Verification")
        log_print("Project: $project_name")
        log_print("Stage: $(current_stage)")
        log_print("="^64)

        # ── [1/8] Read mesh ──────────────────────────────────────────
        log_print("\n[1/8] Reading mesh file: $mesh_file")
        mesh = read_mesh_file(mesh_file)
        log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")
        log_print("   ✓ Loaded initial and boundary conditions")

        # ── [2/8] Read material properties ────────────────────────
        log_print("\n[2/8] Reading material properties: $mat_file")
        materials = read_materials_file(mat_file)
        log_print("   ✓ Loaded $(length(materials.soil_dictionary)) soils")

        # ── [3/8] Read calculation parameters ─────────────────────
        log_print("\n[3/8] Reading calculation parameters: $calc_file")
        calc_params = get_all_calc_params(calc_file)
        log_print(log_analysis_type(calc_params["solver_settings"]))
        log_print("   ✓ Total simulation time: $(calc_params["time_stepping"]["total_simulation_time"]) s")
        
        # Step 3.1: Compute K_sat (using kernel pattern - enables any SWRC model)
        compute_K_sat_runtime!(materials, calc_params)
        log_print("   ✓ K_sat computed for all soils from intrinsic permeability")

        # ── Step 3.7: Check for existing checkpoint from previous stage (kernel pattern) ──
        checkpoint_file, prev_stage = find_latest_checkpoint(project_name, output_dir)
        checkpoint_loaded = false
        initial_state = nothing
        
        if checkpoint_file !== nothing && current_stage > 1
            log_print("\n    Loading checkpoint from previous stage")
            log_print("   Found checkpoint: $(basename(checkpoint_file)) (Stage $(prev_stage))")
            
            # Initialize arrays first (dimensions only)
            zero_variables!(mesh, materials)
            
            # Load checkpoint data
            checkpoint_result = load_checkpoint(checkpoint_file, mesh, materials)
            
            if checkpoint_result.success
                checkpoint_loaded = true
                initial_state = (current_time=checkpoint_result.current_time, 
                               output_counter=checkpoint_result.output_counter,
                               next_output_time=checkpoint_result.next_output_time)
                
                checkpoint_size = get_checkpoint_file_size(checkpoint_file)
                log_print("   ✓ Checkpoint loaded successfully ($(checkpoint_size))")
                log_print("   ✓ Restored state at time: $(checkpoint_result.current_time) s")
                log_print("   ✓ Continuing from output counter: $(checkpoint_result.output_counter)")
            else
                log_print("   ⚠ Warning: $(checkpoint_result.message)")
                log_print("   ⚠ Proceeding with normal initialization instead")
            end
        end

        # ── [4/8] Initialize simulation variables (conditional on checkpoint) ─────────────────
        if !checkpoint_loaded
            log_print("\n[4/8] Initializing simulation variables")
            zero_variables!(mesh, materials)
            log_print("   ✓ Allocated arrays for $(mesh.num_nodes) nodes")
        else
            log_print("\n[4/8] Simulation variables initialized from checkpoint")
            log_print("   ✓ Using $(mesh.num_nodes) nodes from checkpoint")
        end

        # ── [5/8] Apply initial conditions and initialize flows (conditional on checkpoint) ────
        if !checkpoint_loaded
            log_print("\n[5/8] Applying initial conditions and initializing flows")
            
            # Get soil and extract parameters for verification reporting (kernel pattern)
            soil_name = materials.soil_dictionary[1]
            soil = materials.soils[soil_name]
            model = get_water_model(soil.water)
            
            if model === nothing
                error("No SWRC model found for soil '$soil_name'. Ensure swrc_model != 'None' in materials file.")
            end
            
            # Extract parameters via polymorphic dispatch (works for all SWRC models)
            liquid_props = get_liquid_properties(materials)
            g_mag = calc_params["gravity"]["magnitude"]
            
            # Polymorphic: K_h at saturation (h=0.0) gives K_sat equivalent for all models
            K_val = K_h(model, 0.0)
            C_val = C_moist(model, 0.0)
            D_val = K_val / C_val
            
            # Apply water-specific initial and boundary conditions using kernel pattern
            apply_initial_conditions_water!(mesh, materials)
            
            # For verification test: override with uniform IC (h = -0.5 m everywhere)
            # This ensures all nodes are initialized regardless of mesh IC format
            for node_id in 1:mesh.num_nodes
                if !haskey(mesh.pressure_head_bc, node_id)
                    h[node_id] = -0.5
                    theta_w[node_id] = theta(model, -0.5)
                else
                    h[node_id] = mesh.pressure_head_bc[node_id]
                    theta_w[node_id] = theta(model, mesh.pressure_head_bc[node_id])
                end
            end
            
            # Initialize flows using kernel pattern (handles boundary edges, pressure BC flows)
            initialize_all_flows!(mesh, materials, mesh.num_nodes, 0)
            log_print("   ✓ Initial and boundary conditions applied (water-only)")
            log_print("   ✓ Flow vectors initialized")
            log_print("   ✓ SWRC Model: $(typeof(model).name)")
        else
            log_print("\n[5/8] Re-initializing flows from checkpoint state")
            # Get soil model (needed for D_val calculation in plotting)
            soil_name = materials.soil_dictionary[1]
            soil = materials.soils[soil_name]
            model = get_water_model(soil.water)
            K_val = K_h(model, 0.0)
            C_val = C_moist(model, 0.0)
            D_val = K_val / C_val
            
            # Reapply boundary conditions from mesh file (may have changed between stages)
            apply_pressure_bc!(mesh)
            initialize_all_flows!(mesh, materials, mesh.num_nodes, 0)
            log_print("   ✓ Boundary conditions reapplied")
            log_print("   ✓ Flow arrays initialized")
        end

        # ── [6/8] Initialize shape functions ──────────────────────
        log_print("\n[6/8] Initializing shape functions")
        initialize_shape_functions!(mesh)
        log_print("   ✓ Shape functions and Jacobians precomputed")
        
        # ── [7/8] Calculate time step information (kernel pattern) ──────────────────
        log_print("\n[7/8] Calculating time step information")
        time_data, limiting_scale = calculate_time_step_info(mesh, materials, calc_params)
        dt_out = calc_params["data_saving_interval"]
        
        log_print(@sprintf("   ✓ Minimum characteristic length: %.3g %s", time_data.h_min, calc_params["units"]["geometry_unit"]))
        log_print(@sprintf("   ✓ Critical time step: %.4g %s", time_data.critical_dt, calc_params["units"]["time_unit"]))
        log_print("   ✓ Limiting time scale: $(limiting_scale)")
        log_print("   ✓ Courant number: $(time_data.courant_number)")
        log_print(@sprintf("   ✓ Actual time step: %.4g %s", time_data.actual_dt, calc_params["units"]["time_unit"]))
        log_print("   ✓ Number of time steps: $(time_data.num_steps)")
        
        # Build Richards cache and element properties
        cache = build_richards_cache(mesh)
        elem_props = [ElementWaterProps(model) for _ in 1:mesh.num_elements]
        log_print("   ✓ RichardsCache built")
        log_print("   ✓ Element properties created for $(mesh.num_elements) elements")
        
        # ── [8/8] Call production implicit solver ──────────────────
        log_print("\n[8/8] Running implicit Richards solver (production wrapper)")
        log_print("   ✓ Solver: implicit_richards_solver()")
        log_print("   ✓ Pattern: identical to ADSIM kernel")
        log_print(@sprintf("   ✓ Configuration: D = %.6f m²/s, K = %.4f m/s, C = %.4f m⁻¹", D_val, K_val, C_val))
        log_print("-"^64)
        
        final_state = implicit_richards_solver(
            mesh, materials, calc_params, time_data,
            project_name, log_print, cache, elem_props, initial_state
        )
        
        log_print("-"^64)
        log_print(@sprintf("   ✓ Solver completed at t = %.4f s", final_state.current_time))
        log_print(@sprintf("   ✓ Output steps: %d", final_state.output_counter))

        # ── Post-processing: Extract solution for verification ────────
        log_print("\nPost-processing: Extracting solution for verification")
        N = mesh.num_nodes
        
        # Collect all VTK files generated (more robust matching)
        vtk_pattern = "$(project_name)_water_*.vtk"
        vtk_files = sort(filter(f -> occursin(r"_water_\d{6}\.vtk$", f),
                                 readdir(output_dir)))
        
        # Parse VTK files to recover solution history
        function parse_vtk_head(filepath)
            try
                lines = readlines(filepath)
                h_vec = zeros(N)
                idx = findfirst(l -> startswith(l, "SCALARS Matric_Head"), lines)
                if idx !== nothing
                    data_start = idx + 2
                    for i in 1:N
                        if data_start + i - 1 <= length(lines)
                            h_vec[i] = parse(Float64, strip(lines[data_start + i - 1]))
                        end
                    end
                    return h_vec
                end
            catch e
                return zeros(N)
            end
            return zeros(N)
        end
        
        # Read all VTK files for complete history
        t_hist = Float64[]
        h_hist = Vector{Float64}[]
        
        for (i, vtkfile) in enumerate(vtk_files)
            filepath = joinpath(output_dir, vtkfile)
            h_recovered = parse_vtk_head(filepath)
            # All files should have valid data (solver wrote them)
            push!(h_hist, h_recovered)
            push!(t_hist, (i - 1) * dt_out)  # Time based on output interval
        end
        
        log_print(@sprintf("   ✓ Recovered %d solution snapshots from VTK files", length(h_hist)))
        log_print(@sprintf("   ✓ Time history range: t = [%.1f, %.1f] s", t_hist[1], t_hist[end]))
        log_print(@sprintf("   ✓ Time step between snapshots: %.2f s", length(t_hist) > 1 ? t_hist[2] - t_hist[1] : 0.0))

        # ── Post-processing: L2 errors + plots ────────────────────────
        log_print("\nPost-processing: L2 error vs. analytical solution")

        # Extract LEFT column only (x = 0) — effective 1D verification
        tol_x = 1.0e-6
        col_nodes = sort(findall(i -> mesh.coordinates[i,1] < tol_x, 1:N),
                         by = i -> mesh.coordinates[i,2])
        y_col = mesh.coordinates[col_nodes, 2]

        log_print(@sprintf("   ✓ Extraction column: x = 0.0 m  (%d nodes)", length(col_nodes)))
        log_print("-"^64)

        L2_errors = Float64[]
        for (snap_idx, t_val) in enumerate(t_hist)
            h_fem = h_hist[snap_idx][col_nodes]
            h_ana = [h_analytical(y, t_val, D_val) for y in y_col]
            err   = sqrt(mean((h_fem .- h_ana).^2))
            push!(L2_errors, err)
            log_print(@sprintf("   t = %.3f s  |  L2 error = %.2e m", t_val, err))
        end
        log_print("-"^64)
        log_print(@sprintf("   Max L2 error : %.2e m", maximum(L2_errors)))
        log_print(@sprintf("   Final L2 error (t = %.2f s): %.2e m", t_hist[end], L2_errors[end]))

        # Plots with debug output — Colab style (3 columns)
        colors      = [:steelblue, :firebrick, :seagreen, :darkorange, :purple]
        y_fine      = range(0.0, 1.0, length=300)
        
        # Select specific time points to plot
        times_plot = [0.0, 0.1, 0.25, 0.5, 1.0]
        
        log_print("\n   Plotting time profiles (3 columns):")

        # COLUMN 1: h(y) profiles FEM vs Analytical
        p_h = plot(xlabel="h [m]", ylabel="y [m]",
                   title="Pressure Head h(y)",
                   legend=:bottomright, grid=true, gridalpha=0.3)
        
        # COLUMN 2: θ(y) profiles
        p_theta = plot(xlabel="θ [-]", ylabel="y [m]",
                       title="Water Content θ(y)",
                       legend=:bottomright, grid=true, gridalpha=0.3)

        # COLUMN 3: L2 error
        p_error = plot(t_hist, L2_errors,
                       xlabel="t [s]", ylabel="L2 error [m]",
                       title="L2 Error Convergence",
                       yscale=:log10, lw=2.5, color=:firebrick,
                       marker=:circle, ms=5, legend=false, grid=true, gridalpha=0.3)

        for (k, t_plot) in enumerate(times_plot)
            snap_idx = argmin(abs.(t_hist .- t_plot))
            t_act    = t_hist[snap_idx]
            col      = colors[mod1(k, length(colors))]
            lbl      = @sprintf("t=%.2f s", t_act)
            
            log_print(@sprintf("     Requested t=%.2f s → Found snap_idx=%d (t=%.3f s)", t_plot, snap_idx, t_act))
            
            # Analytical h profile
            h_ana_ln = [h_analytical(y, t_act, D_val) for y in y_fine]
            plot!(p_h, h_ana_ln, collect(y_fine), color=col, lw=3, ls=:dash, label="Ana $lbl", alpha=0.7)
            scatter!(p_h, h_hist[snap_idx][col_nodes], y_col,
                     color=col, ms=6, marker=:circle, label="FEM $lbl", alpha=0.85)
            
            # θ profile (only FEM)
            h_fem_snapshot = h_hist[snap_idx][col_nodes]
            theta_fem = [theta(model, h_val) for h_val in h_fem_snapshot]
            plot!(p_theta, theta_fem, y_col,
                  color=col, lw=3, marker=:circle, ms=6, label=lbl, alpha=0.85)
        end

        # Combine: 3 columns layout
        p_all = plot(p_h, p_theta, p_error, layout=(1,3), size=(1350,450),
                     plot_title="Linear Diffusion Verification — ConstantSoil (D = $(round(D_val, digits=4)) m²/s)")
        display(p_all)

        out_png = joinpath(output_dir, "$(project_name)_verification.png")
        savefig(p_all, out_png)
        log_print("\n   ✓ Plot saved: $out_png")

        # Write checkpoint file for multi-stage calculations (kernel pattern)
        log_print("\nWriting checkpoint file for stage $(current_stage)...")
        checkpoint_file = write_checkpoint(project_name, current_stage,
                                          final_state.current_time,
                                          final_state.output_counter,
                                          final_state.next_output_time)
        checkpoint_size = get_checkpoint_file_size(checkpoint_file)
        log_print("   ✓ Checkpoint saved: $(basename(checkpoint_file)) ($(checkpoint_size))")

        # Summary
        end_time   = now()
        total_time = (end_time - start_time).value / 1000.0
        log_print("\nTotal run time: $(total_time) s")
        log_print("\n" * "="^64)
        log_print("Verification completed successfully")
        log_print("="^64)

    catch e
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
        println("\nFATAL ERROR — check log: $log_file_path")
        close(log_file)
        exit(1)

    finally
        isopen(log_file) && close(log_file)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
