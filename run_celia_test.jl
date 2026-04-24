#!/usr/bin/env julia
# Quick test script to run CeliaCol with cold start

# Push project name to ARGS if not already there
if isempty(ARGS)
    push!(ARGS, "CeliaCol")
end

include("src/ADSIM.jl")
using .ADSIM

# Run CeliaCol project from root directory
# (kernel.jl handles cd into src internally)
ADSIM.main()
