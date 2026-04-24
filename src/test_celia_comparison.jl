#!/usr/bin/env julia
# Test script for profile comparison

include("swrc_models.jl")
include("plot_celia_comparison.jl")

# Extract profiles
profiles = plot_celia_comparison("CeliaCol"; output_dir="../output")

# Quantitative analysis
if profiles !== nothing && !isempty(profiles)
    quantitative_comparison(profiles, VanGenuchten(
        θ_r=0.102, θ_s=0.368,
        α=0.335, n_vg=2.0,
        K_sx=9.4e-12, K_sy=9.4e-12,
        L=0.5
    ))
end

println("\n✓ Profile extraction complete")
