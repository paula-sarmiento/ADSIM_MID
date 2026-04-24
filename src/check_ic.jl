#!/usr/bin/env julia
"""
Check initial and final conditions
"""

using Printf

include("swrc_models.jl")

function read_vtk_field(filename)
    open(filename, "r") do f
        lines = readlines(f)
        
        # Find POINT_DATA scalar
        pd_idx = findfirst(l -> startswith(l, "POINT_DATA"), lines)
        sc_idx = findnext(l -> startswith(l, "SCALARS"), lines, pd_idx)
        lk_idx = findnext(l -> startswith(l, "LOOKUP_TABLE"), lines, sc_idx)
        
        values = Float64[]
        for i in (lk_idx+1):length(lines)
            for v in split(lines[i])
                try
                    push!(values, parse(Float64, v))
                catch end
            end
        end
        return values
    end
end

vg_celia = VanGenuchten(
    theta_s=0.368, theta_r=0.102,
    alpha=0.335, n_vg=2.0,
    K_s=9.4e-12
)

println("="^70)
println("INITIAL vs FINAL CONDITION CHECK")
println("="^70)

# Check step 0
println("\nStep 0 (Initial Condition):")
h_init = read_vtk_field("../output/CeliaCol_dt144_water_000000.vtk")
theta_init = [theta(vg_celia, h) for h in h_init]
println("  h: min=$(minimum(h_init))  max=$(maximum(h_init))")
println("  θ: min=$(minimum(theta_init))  max=$(maximum(theta_init))")
println("  Expected: θ ≈ 0.178 (dry)")

# Check step 24
println("\nStep 24 (Final at t=24h):")
h_final = read_vtk_field("../output/CeliaCol_dt144_water_000024.vtk")
theta_final = [theta(vg_celia, h) for h in h_final]
println("  h: min=$(minimum(h_final))  max=$(maximum(h_final))")
println("  θ: min=$(minimum(theta_final))  max=$(maximum(theta_final))")

println("\nChange:")
dtheta = mean(theta_final) - mean(theta_init)
println("  Δθ = $(mean(theta_final)) - $(mean(theta_init)) = $dtheta")

if dtheta < 0.001
    println("  WARNING: Almost no infiltration - simulation may not have run!")
end

println("\n" * "="^70)
