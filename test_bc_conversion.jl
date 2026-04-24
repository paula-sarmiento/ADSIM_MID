#!/usr/bin/env julia

include("src/read_mesh.jl")
include("src/read_materials.jl")
include("src/swrc_models.jl")

# Read mesh
mesh = read_mesh("src/data/Test1000s.mesh")
println("Before normalization:")
println("  volumetric_content_bc entries: $(length(mesh.volumetric_content_bc))")
println("  pressure_head_bc entries: $(length(mesh.pressure_head_bc))")

# Read materials
materials = read_materials("src/data/Test1000s_mat.toml")

# Create SWRC models after reading materials
for (soil_name, soil) in materials.soils
    if soil.water.swrc_model == "Van_Genuchten"
        soil.water.swrc_model_instance = VanGenuchten(
            theta_s = soil.water.theta_s,
            theta_r = soil.water.theta_r,
            alpha = soil.water.swrc_vg_alpha,
            n_vg = soil.water.swrc_vg_n,
            K_sat = 1e-5
        )
    end
end

# Normalize (this should convert volumetric_content_bc to pressure_head_bc)
normalize_water_conditions!(mesh, materials)

println("\nAfter normalization:")
println("  volumetric_content_bc entries: $(length(mesh.volumetric_content_bc))")
println("  pressure_head_bc entries: $(length(mesh.pressure_head_bc))")
if length(mesh.pressure_head_bc) > 0
    println("  pressure_head_bc values:")
    for (node_id, h_val) in sort(collect(mesh.pressure_head_bc))
        println("    Node $node_id: h = $h_val m")
    end
end
