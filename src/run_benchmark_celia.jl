#!/usr/bin/env julia
"""
run_benchmark_celia.jl

Execute Celia et al. (1990) benchmark with multiple timesteps via ADSIM kernel.
Each run goes through full kernel workflow (mesh, BCs, element properties, etc.)

Configurations:
  - Δt = 144s (2.4 min)  — finest
  - Δt = 720s (12 min)   — medium
  - Δt = 3600s (60 min)  — coarse (current default)

Stores outputs with dt-specific prefixes for comparison tracking.
"""

using Printf

println("="^70)
println("  CELIA BENCHMARK — MULTIPLE TIMESTEPS (via ADSIM kernel)")
println("="^70)

dt_configs = [
    (dt=144.0,  label="Δt=2.4 min   (144s)",  idx=1),
    (dt=720.0,  label="Δt=12 min    (720s)",  idx=2),
    (dt=3600.0, label="Δt=60 min   (3600s)",  idx=3)
]

project_dir = joinpath(@__DIR__, "data")
project_root = dirname(@__DIR__)
calc_file = joinpath(project_dir, "CeliaCol_calc.toml")
calc_backup = joinpath(project_dir, "CeliaCol_calc.orig")

# Save original
cp(calc_file, calc_backup, force=true)
original_content = read(calc_file, String)

results = Dict{String, Dict}()

for config in dt_configs
    dt = config.dt
    label = config.label
    cfg_idx = config.idx
    
    println("\n" * "-"^70)
    println("Configuration $cfg_idx / $(length(dt_configs)): $label")
    println("-"^70)
    
    try
        # Delete previous checkpoint to force fresh start
        output_dir = joinpath(project_root, "output")
        checkpoint_files = filter(f -> startswith(f, "CeliaCol_stage") && endswith(f, ".jld2"), readdir(output_dir))
        println("🔍 Found $(length(checkpoint_files)) checkpoint file(s) to delete...")
        for cf in checkpoint_files
            rm(joinpath(output_dir, cf), force=true)
            println("  ✗ Deleted: $cf")
        end
        
        # Modify calc.toml
        modified_content = replace(original_content, 
            r"time_per_step\s*=\s*\d+" => "time_per_step = $(Int(dt))")
        write(calc_file, modified_content)
        println("✓ Modified calc.toml: time_per_step = $(Int(dt))s")
        
        # Run via ADSIM kernel (same as run_celia_test.jl)
        println("Running ADSIM kernel (full workflow)...")
        
        # Push project name to ARGS
        empty!(ARGS)
        push!(ARGS, "CeliaCol")
        
        # Load and run ADSIM
        cd(project_root)
        include(joinpath(@__DIR__, "ADSIM.jl"))
        using .ADSIM
        
        ADSIM.main()
        
        println("✓ Kernel execution complete")
        
        # Copy outputs with dt-specific prefix for tracking
        output_base = joinpath(project_root, "output")
        new_prefix = "CeliaCol_dt$(Int(dt))"
        
        for step in 0:24
            old_file = joinpath(output_base, @sprintf("CeliaCol_water_%06d.vtk", step))
            new_file = joinpath(output_base, @sprintf("%s_water_%06d.vtk", new_prefix, step))
            if isfile(old_file)
                cp(old_file, new_file, force=true)
            end
        end
        println("✓ Outputs saved with prefix: $new_prefix")
        
        results[label] = Dict(
            "dt" => dt,
            "output_prefix" => new_prefix,
            "status" => "success"
        )
        
    catch e
        @warn "Error in configuration $cfg_idx: $(sprint(showerror, e))"
        results[label] = Dict("status" => "failed", "error" => "$(e)")
    finally
        # Always restore original before next iteration
        write(calc_file, original_content)
    end
end

# Final restore
write(calc_file, original_content)

println("\n" * "="^70)
println("✓ BENCHMARK EXECUTION COMPLETE (all via kernel)")
println("="^70)
println("\nResults summary:")
for (label, result) in results
    status = result["status"]
    println("  • $label: $status")
end
println("\nNext: Run extract_benchmark_comparison.jl to visualize results")
println("  cd src")
println("  julia --project=.. extract_benchmark_comparison.jl")
println("="^70)
