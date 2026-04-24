using Printf, TOML

include("src/ADSIM.jl")
using .ADSIM

function run_dt(dt_s, label)
    # Update calc.toml
    data = TOML.parsefile("src/data/CeliaCol_calc.toml")
    data["time_stepping"]["total_simulation_time"] = 86400.0
    data["time_stepping"]["time_per_step"] = Float64(dt_s)
    data["data_saving"]["data_saving_interval"] = 21600.0   # save at 6h, 12h, 18h, 24h
    open("src/data/CeliaCol_calc.toml", "w") do f
        TOML.print(f, data)
    end

    # Remove old checkpoint
    for f in readdir("output")
        if startswith(f, "CeliaCol_stage") && endswith(f, ".jld2")
            rm(joinpath("output", f))
        end
    end

    # Rename output files by copying final vtk
    empty!(ARGS); push!(ARGS, "CeliaCol")
    ADSIM.main()

    # Copy the 24h VTK to a labeled file
    src = "output/CeliaCol_water_000004.vtk"
    dst = "output/CeliaCol_dt$(dt_s)s_24h.vtk"
    isfile(src) && cp(src, dst; force=true)
    println(">>> Saved 24h profile for dt=$(dt_s)s → $dst")
end

run_dt(144,  "2.4min")
run_dt(720,  "12min")
run_dt(3600, "60min")
println("\nAll three runs complete.")
