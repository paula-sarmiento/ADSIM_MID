using Plots, Printf

function parse_vtk_water(filepath)
    lines = readlines(filepath)
    n = 42

    # Extract coordinates (lines after "POINTS 42 float")
    coords = zeros(n, 2)
    pt_start = findfirst(l -> startswith(l, "POINTS"), lines) + 1
    for i in 1:n
        parts = split(strip(lines[pt_start + i - 1]))
        coords[i, 1] = parse(Float64, parts[1])  # x
        coords[i, 2] = parse(Float64, parts[2])  # y
    end

    function read_scalar(label)
        idx = findfirst(l -> startswith(l, "SCALARS $label"), lines)
        idx === nothing && return fill(NaN, n)
        data_start = idx + 2  # skip LOOKUP_TABLE line
        return [parse(Float64, strip(lines[data_start + i - 1])) for i in 1:n]
    end

    h     = read_scalar("Matric_Head")
    theta = read_scalar("Water_Content")
    sat   = read_scalar("Saturation")

    return coords, h, theta, sat
end

output_dir = "output"
files = [joinpath(output_dir, @sprintf("ColabTest_water_%06d.vtk", i)) for i in 0:4]
times = [0.0, 0.5, 1.0, 1.5, 2.0]

# Use left-column nodes (x=0): nodes 1,3,5,...,41 → indices 1,3,5,...,41
left_idx = 1:2:41

# Gather profiles
coords0, _, _, _ = parse_vtk_water(files[1])
all_y   = coords0[left_idx, 2]
all_h   = Vector{Vector{Float64}}()
all_th  = Vector{Vector{Float64}}()
all_sat = Vector{Vector{Float64}}()

for f in files
    coords, h, theta, sat = parse_vtk_water(f)
    push!(all_h,   h[left_idx])
    push!(all_th,  theta[left_idx])
    push!(all_sat, sat[left_idx])
end

colors = [:black, :blue, :green, :orange, :red]
labels = ["t = 0.0 s", "t = 0.5 s", "t = 1.0 s", "t = 1.5 s", "t = 2.0 s"]

# ── Figure 1: Matric Head ──────────────────────────────────────────────────
p1 = plot(title="Matric Head Profile — ColabTest",
          xlabel="Matric Head h [m]", ylabel="Depth y [m]",
          legend=:topright, size=(600,700), dpi=150,
          framestyle=:box, grid=true)
for i in 1:5
    plot!(p1, all_h[i], all_y; label=labels[i], color=colors[i], lw=2)
end
# Mark BCs
scatter!(p1, [-1.0], [0.0]; label="BC: h=-1 (bottom)", marker=:diamond, ms=8, color=:purple)
scatter!(p1, [ 0.0], [1.0]; label="BC: h=0 (top)",    marker=:star5,   ms=8, color=:red)
savefig(p1, "output/ColabTest_matric_head.png")
println("Saved: output/ColabTest_matric_head.png")

# ── Figure 2: Water Content ────────────────────────────────────────────────
p2 = plot(title="Water Content Profile — ColabTest",
          xlabel="Volumetric Water Content θ [-]", ylabel="Depth y [m]",
          legend=:topright, size=(600,700), dpi=150,
          framestyle=:box, grid=true)
for i in 1:5
    plot!(p2, all_th[i], all_y; label=labels[i], color=colors[i], lw=2)
end
vline!(p2, [0.2975]; label="IC: θ=0.2975", ls=:dash, color=:gray, lw=1)
savefig(p2, "output/ColabTest_water_content.png")
println("Saved: output/ColabTest_water_content.png")

# ── Figure 3: Both side by side ───────────────────────────────────────────
p3 = plot(p1, p2; layout=(1,2), size=(1200,700), dpi=150)
savefig(p3, "output/ColabTest_profiles.png")
println("Saved: output/ColabTest_profiles.png")

println("\nAll figures saved to output/")
