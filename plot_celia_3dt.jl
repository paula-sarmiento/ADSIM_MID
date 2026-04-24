using Plots, Printf

function parse_vtk_water(filepath)
    lines = readlines(filepath)
    n = 82
    pt_start = findfirst(l -> startswith(l, "POINTS"), lines) + 1
    coords = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        parts = split(strip(lines[pt_start+i-1]))
        coords[i, 1] = parse(Float64, parts[1])
        coords[i, 2] = parse(Float64, parts[2])
    end
    h_idx    = findfirst(l -> startswith(l, "SCALARS Matric_Head"), lines)
    theta_idx = findfirst(l -> startswith(l, "SCALARS Water_Content"), lines)
    h     = [parse(Float64, strip(lines[h_idx+2+i]))     for i in 0:n-1]
    theta = [parse(Float64, strip(lines[theta_idx+2+i])) for i in 0:n-1]
    return coords, h, theta
end

left = collect(1:2:81)   # left-column node indices (x=0), 41 nodes

dt_vals  = [144, 720, 3600]
dt_labels = ["Δt = 144 s (2.4 min)", "Δt = 720 s (12 min)", "Δt = 3600 s (60 min)"]
colors    = [:blue, :red, :green]

# Parse initial condition from dt=144 run (same IC for all)
coords0, _, _ = parse_vtk_water("output/CeliaCol_dt144s_24h.vtk")
z = coords0[left, 2]

# ── θ at t = 24 h, three Δt ────────────────────────────────────────────────
p1 = plot(title="Water Content at t=24 h — Celia 1990 Fig 6B\n(three time steps)",
          xlabel="θ [-]", ylabel="z [m]",
          legend=:bottomright, size=(600, 750), dpi=150,
          framestyle=:box, grid=true, xlims=(0.15, 0.42), ylims=(0.0, 1.0))

# IC profile from initial VTK (step 0 of dt=144 run)
coords_ic, _, theta_ic = parse_vtk_water("output/CeliaCol_water_000000.vtk")
plot!(p1, theta_ic[left], coords_ic[left, 2];
      label="t = 0 (IC)", color=:black, lw=2, ls=:dash)

for (dt, lbl, col) in zip(dt_vals, dt_labels, colors)
    fname = "output/CeliaCol_dt$(dt)s_24h.vtk"
    _, _, theta = parse_vtk_water(fname)
    plot!(p1, theta[left], z; label=lbl, color=col, lw=2)
end
vline!(p1, [0.178];  label="IC: θ=0.178",    ls=:dot,  color=:gray,    lw=1)
vline!(p1, [0.3601]; label="Top BC: θ=0.360", ls=:dot, color=:darkred, lw=1)
savefig(p1, "output/Celia3dt_theta_24h.png")
println("Saved: output/Celia3dt_theta_24h.png")

# ── h at t = 24 h, three Δt ────────────────────────────────────────────────
p2 = plot(title="Matric Head at t=24 h — Celia 1990 Fig 6B\n(three time steps)",
          xlabel="h [m]", ylabel="z [m]",
          legend=:bottomright, size=(600, 750), dpi=150,
          framestyle=:box, grid=true, xlims=(-12.0, 2.0), ylims=(0.0, 1.0))

coords_ic2, h_ic, _ = parse_vtk_water("output/CeliaCol_water_000000.vtk")
plot!(p2, h_ic[left], coords_ic2[left, 2];
      label="t = 0 (IC)", color=:black, lw=2, ls=:dash)

for (dt, lbl, col) in zip(dt_vals, dt_labels, colors)
    fname = "output/CeliaCol_dt$(dt)s_24h.vtk"
    _, h, _ = parse_vtk_water(fname)
    plot!(p2, h[left], z; label=lbl, color=col, lw=2)
end
hline!(p2, [0.0];  label="",                ls=:dot, color=:gray, lw=1)
scatter!(p2, [-10.0], [0.0];  label="Bot BC: h=-10 m", ms=7, marker=:diamond, color=:purple)
scatter!(p2, [-0.75], [1.0];  label="Top BC: h=-0.75 m", ms=7, marker=:star5, color=:darkred)
savefig(p2, "output/Celia3dt_head_24h.png")
println("Saved: output/Celia3dt_head_24h.png")
