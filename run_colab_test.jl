#!/usr/bin/env julia
# Run ColabTest: Colab implicit Richards benchmark
# VG: alpha=2.0, Ks=0.1 m/s, IC h=-0.5, Bot BC h=-1 (Dirichlet), Top BC h=0 (Dirichlet)
# T=2.0 s, dt=0.05 s, 40 steps, output every 0.5 s

if isempty(ARGS)
    push!(ARGS, "ColabTest")
end

include("src/ADSIM.jl")
using .ADSIM

ADSIM.main()
