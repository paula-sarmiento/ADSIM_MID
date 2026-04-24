using Plots, Printf

function parse_vtk_h(filepath)
    lines = readlines(filepath)
    n = 42
    pt_start = findfirst(l -> startswith(l, "POINTS"), lines) + 1
    coords = [parse(Float64, split(strip(lines[pt_start+i-1]))[2]) for i in 1:n]  # y only
    idx = findfirst(l -> startswith(l, "SCALARS Matric_Head"), lines)
    h = [parse(Float64, strip(lines[idx+2+i])) for i in 0:n-1]
    theta_idx = findfirst(l -> startswith(l, "SCALARS Water_Content"), lines)
    theta = [parse(Float64, strip(lines[theta_idx+2+i])) for i in 0:n-1]
    return coords, h, theta
end

left = 1:2:41   # left-column nodes (x=0)

files = ["output/CeliaCol_water_$(lpad(i,6,'0')).vtk" for i in [0,1,4,8,12,24]]
labels = ["t = 0 h", "t = 1 h", "t = 4 h", "t = 8 h", "t = 12 h", "t = 24 h"]
colors = [:black, :blue, :cyan, :green, :orange, :red]

y_ref, _, _ = parse_vtk_h(files[1])
y = y_ref[left]

# ── Matric head ────────────────────────────────────────────────
p1 = plot(title="Matric Head — Celia 1990 Setup\n(VG model, Δt=360 s)",
          xlabel="h [m]", ylabel="z [m]",
          legend=:right, size=(550,700), dpi=150, framestyle=:box, grid=true)

for (i, (f, lbl, col)) in enumerate(zip(files, labels, colors))
    _, h, _ = parse_vtk_h(f)
    plot!(p1, h[left], y; label=lbl, color=col, lw=2)
end
# BC markers
scatter!(p1, [-10.0], [0.0];  label="BC bottom (h=-10 m)", ms=8, marker=:diamond, color=:purple)
scatter!(p1, [-0.75], [0.60]; label="BC top (h=-0.75 m)",  ms=8, marker=:star5,   color=:darkred)
savefig(p1, "output/CeliaFixed_head.png")
println("Saved: output/CeliaFixed_head.png")

# ── Water content ──────────────────────────────────────────────
p2 = plot(title="Water Content — Celia 1990 Setup\n(VG model, Δt=360 s)",
          xlabel="θ [-]", ylabel="z [m]",
          legend=:right, size=(550,700), dpi=150, framestyle=:box, grid=true)

for (i, (f, lbl, col)) in enumerate(zip(files, labels, colors))
    _, _, theta = parse_vtk_h(f)
    plot!(p2, theta[left], y; label=lbl, color=col, lw=2)
end
vline!(p2, [0.178];  label="IC: θ=0.178", ls=:dash, color=:gray, lw=1)
vline!(p2, [0.3601]; label="Top BC: θ=0.360", ls=:dot, color=:darkred, lw=1)
savefig(p2, "output/CeliaFixed_theta.png")
println("Saved: output/CeliaFixed_theta.png")
