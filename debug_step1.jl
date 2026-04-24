using Printf, TOML

# Set: 3 steps at dt=10s, save every step
data = TOML.parsefile("src/data/CeliaCol_calc.toml")
data["time_stepping"]["total_simulation_time"] = 30.0
data["time_stepping"]["time_per_step"] = 10.0
data["data_saving"]["data_saving_interval"] = 10.0
open("src/data/CeliaCol_calc.toml", "w") do f
    TOML.print(f, data)
end

# Run the solver
empty!(ARGS)
push!(ARGS, "CeliaCol")
include("src/ADSIM.jl")
using .ADSIM
ADSIM.main()

# Check profiles at each step
function get_h(filepath)
    lines = readlines(filepath)
    idx = findfirst(l -> startswith(l, "SCALARS Matric_Head"), lines)
    [parse(Float64, strip(lines[idx+2+i])) for i in 0:81]
end

println("\n=== Profile at top 5 nodes (z=0.875 to 1.0) ===")
left = collect(1:2:81)
y    = [0.025*(i-1) for i in 1:41]
for step in 0:3
    fname = @sprintf("output/CeliaCol_water_%06d.vtk", step)
    isfile(fname) || continue
    h = get_h(fname)
    top5_idx = left[end-4:end]
    top5_y   = y[end-4:end]
    top5_h   = h[top5_idx]
    @printf "step %d  " step
    for (yi, hi) in zip(top5_y, top5_h)
        @printf "z=%.3f:%.4f  " yi hi
    end
    println()
end
