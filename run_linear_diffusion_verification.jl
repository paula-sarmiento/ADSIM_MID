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

    log_file_path = joinpath(output_dir, "$(project_name)_verification.log")
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
        log_print("   ✓ Total simulation time: $(calc_params["time_stepping"]["total_simulation_time"]) s")
        
        # Step 3.1: Compute K_sat (kernel pattern, skip for ConstantSoil)
        # compute_K_sat_runtime!(materials, calc_params)
        log_print("   ✓ (K_sat computation skipped for ConstantSoil verification)")

        # ── [4/8] Initialize simulation variables ─────────────────
        log_print("\n[4/8] Initializing simulation variables")
        zero_variables!(mesh, materials)
        log_print("   ✓ Allocated arrays for $(mesh.num_nodes) nodes")

        # ── [5/8] Apply initial conditions and initialize flows ────
        log_print("\n[5/8] Applying initial conditions and initializing flows")
        
        # For ConstantSoil verification: extract parameters and create model
        soil_props = get_soil_properties(materials, materials.soil_dictionary[1])
        liquid_props = get_liquid_properties(materials)
        
        theta_s = soil_props.porosity
        theta_r = soil_props.water.theta_r
        h_min   = -soil_props.water.swrc_vg_alpha  # Repurposed for ConstantSoil
        k_int   = soil_props.intrinsic_permeability
        rho_w   = liquid_props.density
        mu_w    = liquid_props.dynamic_viscosity
        
        g_mag = calc_params["gravity"]["magnitude"]
        K_val = k_int * rho_w * g_mag / mu_w
        C_val = (theta_s - theta_r) / abs(h_min)
        D_val = K_val / C_val
        
        model = ConstantSoil(theta_r=theta_r, theta_s=theta_s, h_min=h_min, K_val=K_val)
        
        # Apply initial and boundary conditions manually for ConstantSoil
        for node_id in 1:mesh.num_nodes
            if !haskey(mesh.pressure_head_bc, node_id)
                h[node_id] = -0.5
                theta_w[node_id] = theta(model, -0.5)
            else
                h[node_id] = mesh.pressure_head_bc[node_id]
                theta_w[node_id] = theta(model, mesh.pressure_head_bc[node_id])
            end
        end
        
        # Initialize flow vectors
        zero_flow_vectors_water!(mesh.num_nodes)
        log_print("   ✓ Initial conditions applied (h = -0.5 m, BC from mesh)")
        log_print("   ✓ Flow vectors initialized")

        # ── [6/8] Initialize shape functions ──────────────────────
        log_print("\n[6/8] Initializing shape functions")
        initialize_shape_functions!(mesh)
        log_print("   ✓ Shape functions and Jacobians precomputed")
        
        # ── [7/8] Calculate time step information ──────────────────
        log_print("\n[7/8] Setting up solver infrastructure")
        T_total = calc_params["time_stepping"]["total_simulation_time"]
        dt      = calc_params["time_stepping"]["time_per_step"]
        dt_out  = calc_params["data_saving_interval"]
        
        time_data = TimeStepData()
        time_data.total_time = T_total
        time_data.time_per_step = dt
        time_data.actual_dt = dt
        time_data.num_steps = round(Int, T_total / dt)
        time_data.h_min = dt
        log_print(@sprintf("   ✓ Time stepping: %d steps × %.4f s", time_data.num_steps, dt))
        log_print(@sprintf("   ✓ Output interval: %.3f s", dt_out))
        
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
            project_name, log_print, cache, elem_props, nothing
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

        # ── Post-processing: L2 errors + plots ────────────────────────
        log_print("\nPost-processing: L2 error vs. analytical solution")

        Lx    = maximum(mesh.coordinates[:,1])
        x_mid = Lx / 2.0
        tol_x = Lx / (2 * 5)
        col_nodes = sort(findall(i -> abs(mesh.coordinates[i,1] - x_mid) < tol_x, 1:N),
                         by = i -> mesh.coordinates[i,2])
        y_col = mesh.coordinates[col_nodes, 2]

        log_print(@sprintf("   ✓ Extraction column: x ≈ %.2f m  (%d nodes)", x_mid, length(col_nodes)))
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

        # Plots
        colors      = [:steelblue, :firebrick, :seagreen, :darkorange, :purple, :black]
        y_fine      = range(0.0, 1.0, length=300)
        times_plot  = [0.0, 0.1, 0.2, 0.3, 0.5]

        p1 = plot(xlabel="h [m]", ylabel="y [m]",
                  title="h(y) profiles — FEM vs. Analytical",
                  legend=:topright, grid=true, gridalpha=0.3)
        for (k, t_plot) in enumerate(times_plot)
            snap_idx = argmin(abs.(t_hist .- t_plot))
            t_act    = t_hist[snap_idx]
            col      = colors[mod1(k, length(colors))]
            lbl      = @sprintf("t = %.2f s", t_act)
            h_ana_ln = [h_analytical(y, t_act, D_val) for y in y_fine]
            plot!(p1,  h_ana_ln, collect(y_fine), color=col, lw=2, ls=:dash, label="Ana $lbl")
            scatter!(p1, h_hist[snap_idx][col_nodes], y_col,
                     color=col, ms=4, marker=:circle, label="FEM $lbl")
        end

        p2 = plot(t_hist, L2_errors,
                  xlabel="t [s]", ylabel="L2 error [m]",
                  title="L2 error vs. analytical solution",
                  yscale=:log10, lw=2, color=:steelblue,
                  marker=:circle, ms=4, legend=false, grid=true, gridalpha=0.3)

        p3 = plot(xlabel="h [m]", ylabel="y [m]",
                  title="Steady state (t → ∞): FEM vs. exact",
                  legend=:topright, grid=true, gridalpha=0.3)
        h_ss_line = collect(-1.0 .+ y_fine)
        plot!(p3, h_ss_line, collect(y_fine), color=:black, lw=2, ls=:dash, label="h_ss exact")
        scatter!(p3, h_hist[end][col_nodes], y_col,
                 color=:steelblue, ms=5, marker=:circle,
                 label=@sprintf("FEM  t = %.2f s", t_hist[end]))

        p_all = plot(p1, p2, p3, layout=(1,3), size=(1400,500),
                     plot_title="Linear Diffusion Verification — ConstantSoil  " *
                                "(D = K/C = $(round(D_val, digits=4)) m²/s)")
        display(p_all)

        out_png = joinpath(output_dir, "$(project_name)_verification.png")
        savefig(p_all, out_png)
        log_print("\n   ✓ Plot saved: $out_png")

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
