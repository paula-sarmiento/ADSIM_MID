#______________________________________________________
# ADSIM: Aqueous Concentration Solver Verification
# Kernel-pattern implementation
# Pure Diffusion Benchmark — Analytical Fourier Solution
# Input: GID mesh file
#______________________________________________________

#______________________________________________________
# Physics
#   Pure diffusion: ∂C/∂t = ∇·(D_h ∇C)
#   No advection, no source terms (Phase 1 only)
#   Domain: User-supplied via GID mesh file
#   IC: C = 0, BC: C(y=0)=1, C(y=1)=0
#   D_h = 1e-5 m²/s, θ_w = 0.3, Δt = 0.1 s
#
#   Boundary Conditions (Option B):
#   - Specify mole_fraction_bc in GID (yi = composition)
#   - Compute: P_partial = yi × P_total
#   - Apply Henry's Law: C = P_partial / K_H
#______________________________________________________

using LinearAlgebra, SparseArrays, Printf, Statistics, Dates

const PROJECT_ROOT = dirname(dirname(@__DIR__))
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")

include(joinpath(SRC_DIR, "version.jl"))
using .ADSIMVersion: get_version, get_version_string

include(joinpath(SRC_DIR, "read_mesh.jl"))
include(joinpath(SRC_DIR, "read_materials.jl"))
include(joinpath(SRC_DIR, "shape_functions.jl"))
include(joinpath(SRC_DIR, "aqueous_concentration_solver.jl"))

using .ShapeFunctions: initialize_shape_functions!, build_richards_cache

#______________________________________________________
# Analytical Fourier Solution
#______________________________________________________

function C_analytical(y, t, D_h; L=1.0, C_bot=1.0, C_top=0.0, N_modes=200)
    C_ss = C_bot + (C_top - C_bot) * y / L
    transient = sum(((-2.0) / (n * π)) * sin(n * π * y / L) *
                    exp(-D_h * (n * π / L)^2 * t) for n in 1:N_modes)
    return C_ss + transient
end

#______________________________________________________
# Error Metrics
#______________________________________________________

function compute_errors(C_num, C_analytical_vec)
    error = C_num .- C_analytical_vec
    L2 = sqrt(mean(error.^2))
    Linf = maximum(abs.(error))
    return L2, Linf
end

#______________________________________________________
# Main Verification Run
#______________________________________________________

function run_aqueous_verification(project_name, mesh_file, D_h, θw, Δt, n_steps, P_total, log_print)
    log_print("\n[SETUP]")
    log_print("  Reading mesh file: $(mesh_file)")
    mesh = read_mesh_file(mesh_file)
    log_print(@sprintf("  ✓ Mesh: %d nodes, %d elements", mesh.num_nodes, mesh.num_elements))
    
    initialize_shape_functions!(mesh)
    cache = build_richards_cache(mesh)
    log_print("  ✓ Shape functions initialized (2×2 Gauss quadrature)")
    
    log_print("  Reading materials...")
    materials = read_materials_file(joinpath(SRC_DIR, "data", "standard_materials.toml"))
    log_print("  ✓ Materials loaded (CO2, water, soil)")
    
    log_print("\n[INITIALIZATION]")
    C_old = zeros(Float64, mesh.num_nodes)
    θw_arr = fill(θw, mesh.num_nodes)
    vs = zeros(Float64, mesh.num_nodes, 2)  # No advection
    
    # Identify BC nodes (y ≈ 0 and y ≈ 1)
    y_coords = mesh.coordinates[:, 2]
    y_min = minimum(y_coords)
    y_max = maximum(y_coords)
    y_range = y_max - y_min
    
    P_boundary_mask = ones(Int, mesh.num_nodes)
    for i in 1:mesh.num_nodes
        y_normalized = (mesh.coordinates[i, 2] - y_min) / y_range
        if y_normalized < 0.01 || y_normalized > 0.99
            P_boundary_mask[i] = 0
        end
    end
    
    A = spzeros(Float64, mesh.num_nodes, mesh.num_nodes)
    R = zeros(Float64, mesh.num_nodes)
    
    times = Float64[]; L2_errs = Float64[]; Linf_errs = Float64[]
    
    log_print(@sprintf("  ✓ Initialized: D_h=%.2e, θ_w=%.2f, Δt=%.2e, P_total=%.0f Pa", D_h, θw, Δt, P_total))
    
    log_print("\n[TIME-STEPPING]")
    log_print(@sprintf("  n_steps=%d, T_total=%.2e s", n_steps, Δt * n_steps))
    
    t = 0.0
    for step in 1:n_steps
        C_analytical_current = [C_analytical(
            (mesh.coordinates[i, 2] - y_min) / y_range, t, D_h
        ) for i in 1:mesh.num_nodes]
        
        # For verification: directly prescribe concentrations matching analytical solution
        # Skip Henry's law conversion to isolate FEM solver accuracy
        mesh_partial_pressure_bc = Dict{Int, Vector{Float64}}()
        
        C_prescribed = zeros(Float64, mesh.num_nodes)
        for i in 1:mesh.num_nodes
            y_normalized = (mesh.coordinates[i, 2] - y_min) / y_range
            if y_normalized < 0.01
                C_prescribed[i] = 1.0
            elseif y_normalized > 0.99
                C_prescribed[i] = 0.0
            end
        end
        
        C_new = aqueous_concentration_solver(A, R, mesh, cache, materials,
                                             P_boundary_mask, mesh_partial_pressure_bc, 1,
                                             θw_arr, θw_arr, vs, C_old, Δt, 1;
                                             C_prescribed_override=C_prescribed)
        C_old = C_new
        t += Δt
        
        if step % max(1, div(n_steps, 10)) == 0
            L2, Linf = compute_errors(C_new, C_analytical_current)
            push!(times, t)
            push!(L2_errs, L2)
            push!(Linf_errs, Linf)
            log_print(@sprintf("  Step %4d: t=%.2e, L₂=%.2e, L∞=%.2e", step, t, L2, Linf))
        end
    end
    
    log_print("\n[RESULTS]")
    final_L2 = L2_errs[end]
    final_Linf = Linf_errs[end]
    pass_L2 = final_L2 < 1e-3
    pass_Linf = final_Linf < 5e-3
    
    log_print("─"^60)
    log_print(@sprintf("  L₂ error:     %.2e mol/L  (threshold: 1e-3)  %s", 
                      final_L2, pass_L2 ? "✓" : "✗"))
    log_print(@sprintf("  L∞ error:     %.2e mol/L  (threshold: 5e-3)  %s", 
                      final_Linf, pass_Linf ? "✓" : "✗"))
    log_print("─"^60)
    
    return (pass=pass_L2 && pass_Linf, L2=final_L2, Linf=final_Linf)
end

#______________________________________________________
# Main Entry Point
#______________________________________________________

function main()
    if length(ARGS) > 0 && (ARGS[1] == "--version" || ARGS[1] == "-v")
        println(get_version_string())
        exit(0)
    end
    
    if length(ARGS) < 2
        println("Error: Missing arguments")
        println("Usage: julia run_aqueous_verification_kernel.jl <project_name> <mesh_file>")
        println("Example: julia run_aqueous_verification_kernel.jl AqueousVerif_test data/AqueousMesh.mesh")
        println("         julia run_aqueous_verification_kernel.jl --version")
        exit(1)
    end
    
    project_name = ARGS[1]
    mesh_file = ARGS[2]
    data_dir = joinpath(SRC_DIR, "data")
    
    # Resolve mesh file path
    mesh_path = if isfile(mesh_file)
        mesh_file
    elseif isfile(joinpath(data_dir, mesh_file))
        joinpath(data_dir, mesh_file)
    else
        nothing
    end
    
    if mesh_path === nothing || !isfile(mesh_path)
        println("Error: Mesh file not found: $mesh_file")
        println("Searched in:")
        println("  - $(mesh_file)")
        println("  - $(joinpath(data_dir, mesh_file))")
        exit(1)
    end
    
    output_dir = joinpath(dirname(@__DIR__), "output")
    isdir(output_dir) || mkpath(output_dir)
    
    log_file_path = joinpath(output_dir, "$(project_name)_aqueous.log")
    try
        isfile(log_file_path) && rm(log_file_path)
    catch
    end
    log_file = open(log_file_path, "w")
    
    function log_print(msg::String)
        println(msg)
        println(log_file, msg)
        flush(log_file)
    end
    
    try
        start_time = now()
        
        log_print("="^64)
        log_print("ADSIM: Aqueous Concentration Solver Verification")
        log_print("Version: $(get_version())")
        log_print("Pure Diffusion — Analytical Fourier Solution (Kernel)")
        log_print("Option B: mole_fraction_bc (yi composition-based)")
        log_print("="^64)
        log_print("\nProject: $(project_name)")
        log_print("Mesh: $(mesh_path)")
        
        # Verification parameters
        D_h = 1e-5              # Diffusivity [m²/s]
        θw = 0.3                # Water content (constant)
        Δt = 0.1                # Time step [s]
        n_steps = 1000          # Number of steps
        P_total = 101325.0      # Total gas pressure [Pa] (1 atm)
        
        result = run_aqueous_verification(project_name, mesh_path, D_h, θw, Δt, n_steps, P_total, log_print)
        
        end_time = now()
        total_time = (end_time - start_time).value / 1000.0
        
        log_print("\n" * "="^64)
        log_print(result.pass ? "✓ VERIFICATION PASSED" : "✗ VERIFICATION FAILED")
        log_print("="^64)
        log_print(@sprintf("Verification time: %.2f seconds", total_time))
        log_print("Output: $(log_file_path)")
        log_print("="^64)
        
        close(log_file)
        exit(result.pass ? 0 : 1)
        
    catch e
        log_print("\n" * "="^64)
        log_print("FATAL ERROR")
        log_print("="^64)
        log_print("Error: $(e)")
        for (exc, bt) in Base.catch_stack()
            showerror(log_file, exc, bt)
            println(log_file)
        end
        close(log_file)
        exit(1)
    finally
        isopen(log_file) && close(log_file)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
