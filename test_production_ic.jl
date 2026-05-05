#!/usr/bin/env julia
#=
Test script to verify fresh initialization produces correct L2 error
This script deletes the checkpoint and runs a fresh initialization test
=#

using Printf

const PROJECT_ROOT = dirname(@__FILE__)
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "output")
const CHECKPOINT_FILE = joinpath(OUTPUT_DIR, "LinVerif_diffusion_stage1.jld2")
const LOG_FILE = joinpath(OUTPUT_DIR, "LinVerif_diffusion_verification.log")

println("\n" * "="^70)
println("Testing Fresh Initialization (ConstantSoil with Manual IC/BC)")
println("="^70)

# Step 1: Delete checkpoint to force fresh initialization
println("\n[Step 1] Deleting checkpoint file to force fresh initialization...")
if isfile(CHECKPOINT_FILE)
    rm(CHECKPOINT_FILE)
    println("   ✓ Deleted: $CHECKPOINT_FILE")
else
    println("   ℹ Checkpoint not found (already deleted)")
end

# Step 2: Delete old log
println("\n[Step 2] Clearing previous log file...")
if isfile(LOG_FILE)
    rm(LOG_FILE)
    println("   ✓ Deleted old log")
else
    println("   ℹ No previous log")
end

# Step 3: Run the verification script
println("\n[Step 3] Running verification script with fresh initialization...")
println("   Command: julia run_linear_diffusion_verification.jl LinVerif_diffusion")
println()

# Change to project directory
cd(PROJECT_ROOT)

# Set ARGS BEFORE including the script (so main() sees the arguments)
empty!(ARGS)
push!(ARGS, "LinVerif_diffusion")

# Include the script module (it defines main() function)
include("run_linear_diffusion_verification.jl")

# Now call main() with ARGS set correctly
try
    main()
catch e
    println("\n⚠ Script execution completed or encountered error")
end

# Step 4: Check logs to see if initialization worked correctly
println("\n[Step 4] Analyzing results...")
println()

if isfile(LOG_FILE)
    log_content = read(LOG_FILE, String)
    
    # Extract and show final L2 error
    l2_error = nothing
    final_line = ""
    for line in split(log_content, "\n")
        if contains(line, "Final L2 error")
            final_line = line
            # Try to extract the numeric value
            if contains(line, "e-04")
                l2_error = "e-04 (✓ CORRECT)"
            elseif contains(line, "e-01") || contains(line, "e-00")
                l2_error = "e-01/e-00 (❌ WRONG - indicates bad IC/BC)"
            end
            break
        end
    end
    
    if !isnothing(l2_error)
        println("═══════════════════════════════════════════════════════════════════")
        println("PERFORMANCE METRIC:")
        println("═══════════════════════════════════════════════════════════════════")
        println(final_line)
        println("   → L2 error magnitude: $l2_error")
        println("═══════════════════════════════════════════════════════════════════")
    end
    
    # Show Step 5 from log
    println("\n[Step 5] from log file:")
    println("─" * "─"^69)
    lines = split(log_content, "\n")
    in_step5 = false
    for line in lines
        if contains(line, "[5/8]")
            in_step5 = true
        end
        if in_step5
            if contains(line, "[6/8]")
                break
            end
            println(line)
        end
    end
    println("─" * "─"^69)
    
else
    println("❌ ERROR: Log file not created")
end

println("\n" * "="^70)
println("Test Complete")
println("="^70 * "\n")
