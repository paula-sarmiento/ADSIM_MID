# Plot CeliaCol kernel output — reproduces Celia et al. (1990) Figure 6B
# Reads src/output/CeliaCol_water_*.vtk and plots h(z) and θ(z) profiles.

using Printf, Plots

const TARGET_HOURS = [0.0, 6.0, 12.0, 18.0, 24.0]

function parse_vtk_celia(filepath::String)
    lines = readlines(filepath)

    # Time from header line 2
    time_s = 0.0
    if length(lines) >= 2
        m = match(r"Time\s*=\s*([0-9eE+.\-]+)", lines[2])
        m !== nothing && (time_s = parse(Float64, m.captures[1]))
    end

    # POINTS → y-coordinates
    pi = findfirst(l -> startswith(l, "POINTS"), lines)
    pi === nothing && error("POINTS not found in $filepath")
    n = parse(Int, match(r"POINTS\s+(\d+)", lines[pi]).captures[1])
    y = [parse(Float64, split(lines[pi + i])[2]) for i in 1:n]

    function read_scalar(name)
        idx = findfirst(l -> occursin("SCALARS $name", l), lines)
        idx === nothing && return nothing
        vals = Float64[]
        for i in (idx + 2):length(lines)
            l = strip(lines[i])
            (startswith(l, "SCALARS") || startswith(l, "VECTORS")) && break
            l == "" && continue
            try push!(vals, parse(Float64, l)) catch; break end
        end
        length(vals) == n ? vals : nothing
    end

    h     = read_scalar("Matric_Head")
    theta = read_scalar("Water_Content")
    h === nothing && error("Matric_Head not found in $filepath")
    theta === nothing && error("Water_Content not found in $filepath")

    perm  = sortperm(y)
    return (y[perm], h[perm], theta[perm], time_s)
end

function main()
    src_dir = "src/output"
    vtk_files = sort(filter(f -> occursin(r"^CeliaCol_water_\d{6}\.vtk$", f), readdir(src_dir)))
    isempty(vtk_files) && error("No CeliaCol VTK files found in $src_dir")
    println("Found $(length(vtk_files)) CeliaCol VTK files")

    # Load all snapshots
    snaps = []
    for f in vtk_files
        y, h, θ, t = parse_vtk_celia(joinpath(src_dir, f))
        push!(snaps, (y=y, h=h, θ=θ, t_h=t/3600.0))
    end

    # Select snapshot closest to each target time
    selected = Dict{Float64, NamedTuple}()
    for th in TARGET_HOURS
        _, best = findmin(s -> abs(s.t_h - th), snaps)
        selected[th] = snaps[best]
        @printf("  Target t=%4.0f h → found t=%.4f h\n", th, snaps[best].t_h)
    end

    colors  = [:black, :royalblue, :darkorange, :green3, :crimson]
    domain_m = maximum(snaps[1].y)   # column height in metres

    # ── Plot A: Matric Head h [cm] vs Depth [cm] ─────────────────────────────
    pA = plot(
        xlabel = "Pressure Head h (cm)",
        ylabel = "Depth from surface (cm)",
        title  = "Celia et al. (1990) Fig. 6B — Pressure Head",
        legend = :bottomright,
        size   = (700, 700), dpi = 150,
        grid = true, gridalpha = 0.3,
        yflip  = true,
        ylims  = (0, domain_m * 100),
    )
    for (ci, th) in enumerate(TARGET_HOURS)
        d = selected[th]
        depth_cm = (domain_m .- d.y) .* 100.0
        h_cm     = d.h .* 100.0
        lbl = th == 0 ? "t = 0 h (IC)" : "t = $(Int(th)) h"
        plot!(pA, h_cm, depth_cm; label=lbl, color=colors[ci], lw=2.5,
              marker=:circle, ms=3, alpha=0.9)
    end

    # ── Plot B: Water Content θ vs Depth [cm] ────────────────────────────────
    pB = plot(
        xlabel = "Volumetric Water Content θ (−)",
        ylabel = "Depth from surface (cm)",
        title  = "Celia et al. (1990) Fig. 6B — Water Content",
        legend = :bottomright,
        size   = (700, 700), dpi = 150,
        grid = true, gridalpha = 0.3,
        yflip  = true,
        ylims  = (0, domain_m * 100),
    )
    for (ci, th) in enumerate(TARGET_HOURS)
        d = selected[th]
        depth_cm = (domain_m .- d.y) .* 100.0
        lbl = th == 0 ? "t = 0 h (IC)" : "t = $(Int(th)) h"
        plot!(pB, d.θ, depth_cm; label=lbl, color=colors[ci], lw=2.5,
              marker=:circle, ms=3, alpha=0.9)
    end

    # ── Combined figure ───────────────────────────────────────────────────────
    pC = plot(pA, pB, layout=(1,2), size=(1400, 700), dpi=150,
              plot_title="Celia et al. (1990) Figure 6B — Kernel Output")

    out = joinpath(src_dir, "celia_fig6b.png")
    savefig(pC, out)
    println("\nPlot saved: $out")
end

main()
