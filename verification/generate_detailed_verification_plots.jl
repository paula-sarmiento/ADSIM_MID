"""
Generate detailed verification plots with mesh, parameters, and 6-timestep comparison
Format: Main plot + separate legend boxes + equation legend
Saved in SVG format for publication
"""

using Plots, Plots.PlotMeasures, Printf, LaTeXStrings
include("../src/read_mesh.jl")
include("../src/read_materials.jl")
include("../src/shape_functions.jl")
include("../src/aqueous_concentration_solver.jl")
include("aqueous_verification.jl")

using .ShapeFunctions

# Helper: Create a legend text box (separate from plot)
function create_legend_textbox(title::String, info_lines::Vector{String})
    p = plot(legend=false, axis=false, grid=false, framestyle=:box, 
             margin=8mm, xlims=(0, 100), ylims=(0, 100))
    
    y_pos = 90
    annotate!(p, 5, y_pos, text(title, 12, :left, :top, weight=:bold))
    y_pos -= 12
    
    for line in info_lines
        annotate!(p, 5, y_pos, text(line, 9, :left, :top))
        y_pos -= 10
    end
    
    return p
end

# Helper: Create equation legend with term explanations
function create_equation_legend(eq_str::String, terms_dict::Dict{String, String})
    p = plot(legend=false, axis=false, grid=false, framestyle=:box,
             margin=8mm, xlims=(0, 100), ylims=(0, 100))
    
    y_pos = 95
    annotate!(p, 5, y_pos, text("EQUATION:", 11, :left, :top, weight=:bold))
    y_pos -= 12
    
    annotate!(p, 5, y_pos, text(eq_str, 10, :left, :top))
    y_pos -= 15
    
    annotate!(p, 5, y_pos, text("TERMS:", 10, :left, :top, weight=:bold))
    y_pos -= 12
    
    for (term, meaning) in sort(collect(terms_dict))
        line_str = term * " = " * meaning
        annotate!(p, 5, y_pos, text(line_str, 8, :left, :top))
        y_pos -= 9
    end
    
    return p
end

# Skip helper function - will use inline code instead

# Start main script
println("Loading materials...")
materials = read_materials_file("../src/data/AqueousVerif_T6_mat.toml")
soil = materials.soils["VerifSoil"]
gas = materials.gases["VerifGas"]
D_h = soil.granular_tortuosity * gas.diff_coefficient
H = 0.2; Lx = 0.004

# ─────────────────────────────────────────────────────────────────────────────
# T6: TERZAGHI
# ─────────────────────────────────────────────────────────────────────────────

println("Generating T6 Terzaghi detailed plot...")

T_final_t6 = 0.5 * H^2 / D_h
ny_t6 = 50
# Timesteps as multiples of 5 hours
t_steps_t6 = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0] .* 3600.0  # Convert to seconds  # 6 timesteps, exclude t=0

mesh_t6 = build_column_mesh(H, Lx, ny_t6)
y_t6 = mesh_t6.coordinates[:, 2]
initialize_shape_functions!(mesh_t6)
cache_t6 = build_richards_cache(mesh_t6)

# Run simulation step by step, capturing solutions at target times
N_t6 = mesh_t6.num_nodes
A_t6 = spzeros(Float64, N_t6, N_t6)
R_t6 = zeros(Float64, N_t6)

C_t6 = ones(Float64, N_t6)  # IC: C = 1.0
θw_t6 = fill(1.0, N_t6)
vs_t6 = zeros(Float64, N_t6, 2)

P_bc_t6 = ones(Int, N_t6)
C_presc_t6 = zeros(Float64, N_t6)
top_n1_t6 = 2*ny_t6 + 1;  top_n2_t6 = 2*ny_t6 + 2
P_bc_t6[top_n1_t6] = 0;  P_bc_t6[top_n2_t6] = 0

Δt_t6 = 25.0
C_fem_at_times_t6 = Dict()  # Store solutions at each target time

t_current_t6 = 0.0
next_target_idx_t6 = 1

while t_current_t6 < T_final_t6 && next_target_idx_t6 <= length(t_steps_t6)
    global C_t6, t_current_t6, next_target_idx_t6
    assemble_aqueous_concentration!(A_t6, R_t6, mesh_t6, cache_t6, D_h,
                                    θw_t6, θw_t6, vs_t6, C_t6, Δt_t6,
                                    P_bc_t6, C_presc_t6)
    C_t6 = A_t6 \ R_t6
    t_current_t6 = t_current_t6 + Δt_t6
    
    # Store if close to target time
    target_t = t_steps_t6[next_target_idx_t6]
    if t_current_t6 >= target_t - 0.5*Δt_t6
        C_fem_at_times_t6[target_t] = copy(C_t6)
        next_target_idx_t6 = next_target_idx_t6 + 1
    end
end

# Get analytical at each t_step
C_anal_t6 = [terzaghi_profile.(y_t6, t, H, D_h) for t in t_steps_t6]

# Create LaTeX info boxes
nx_t6 = 2  # Fixed for quasi-1D column
info_mesh_t6 = [
    @sprintf("Mesh Grid:"),
    @sprintf("  nx (x-direction) = %d", nx_t6),
    @sprintf("  ny (y-direction) = %d", ny_t6),
    @sprintf("  Total Elements = %d", ny_t6 * nx_t6),
    @sprintf("  Total Nodes = %d", N_t6),
    @sprintf("  Element Type: Q4 (Bilinear Quadrilateral)"),
    @sprintf("  Nodes per Element: 4"),
    @sprintf("  Gauss Quadrature: 2×2 (4 integration points)"),
    @sprintf("Domain Dimensions:"),
    @sprintf("  H (vertical) = 0.20 m"),
    @sprintf("  L_x (horizontal) = 0.004 m"),
    @sprintf("Materials:"),
    @sprintf("  Soil: VerifSoil (granular_tortuosity = %.4f)", soil.granular_tortuosity),
    @sprintf("  Gas: VerifGas (diffusion = %.2e m²/s)", gas.diff_coefficient)
]

info_params_t6 = [
    @sprintf("Physical Parameters:"),
    @sprintf("  D_h (hydrodynamic diffusivity) = %.2e m²/s", D_h),
    @sprintf("  θ_w (water content) = 1.0 [-]"),
    @sprintf("  v_y (advective velocity) = 0.0 m/s (pure diffusion)"),
    @sprintf("Time Integration:"),
    @sprintf("  Δt (timestep) = 25 s"),
    @sprintf("  T_final = %.1f s (≈ %.2f h)", T_final_t6, T_final_t6/3600),
    @sprintf("Boundary Conditions:"),
    @sprintf("  Top (y=H): C = 0 (Dirichlet)"),
    @sprintf("  Bottom (y=0): Natural BC (no flux)")
]

# Generate SINGLE plot with all timesteps overlaid (different colors)
colors_t6 = [:blue, :red, :green, :purple, :orange, :brown]

p_prof_t6 = plot(legend=:bottomright, grid=true, gridalpha=0.3, 
                 xlabel="C", ylabel="y (m)", title="T6 Terzaghi - All Timesteps",
                 size=(700, 500), margin=4mm)

for i in 1:6
    t = t_steps_t6[i]
    C_anal = C_anal_t6[i]
    C_fem = C_fem_at_times_t6[t]
    color = colors_t6[i]
    t_hours = t / 3600.0
    
    # Analytical: solid line
    plot!(p_prof_t6, C_anal, y_t6, color=color, linewidth=2.5, linestyle=:solid, 
          label="Analytical t=$(round(t_hours, digits=1))h", alpha=0.9)
    
    # FEM: scatter points
    scatter!(p_prof_t6, C_fem, y_t6, color=color, markersize=2, markerstrokewidth=0, 
             label="FEM t=$(round(t_hours, digits=1))h", alpha=0.7)
end

p_prof_t6 = plot(legend=:bottomright, grid=true, gridalpha=0.3, 
                 xlabel="C", ylabel="y (m)", title="T6 Terzaghi - All Timesteps",
                 size=(900, 500), margin=4mm)

for i in 1:6
    t = t_steps_t6[i]
    C_anal = C_anal_t6[i]
    C_fem = C_fem_at_times_t6[t]
    color = colors_t6[i]
    t_hours = t / 3600.0
    
    # Analytical: solid line
    plot!(p_prof_t6, C_anal, y_t6, color=color, linewidth=2.5, linestyle=:solid, 
          label="Analytical t=$(round(t_hours, digits=1))h", alpha=0.9)
    
    # FEM: scatter points
    scatter!(p_prof_t6, C_fem, y_t6, color=color, markersize=2, markerstrokewidth=0, 
             label="FEM t=$(round(t_hours, digits=1))h", alpha=0.7)
end

xlims!(p_prof_t6, -0.05, 1.05)
ylims!(p_prof_t6, 0, H)

# Create LEGEND BOXES (separate from plot)
legend_info_t6 = [
    "Domain: Elements: $ny_t6  |  Nodes: $N_t6  |  H = 0.20 m  |  L_x = 0.004 m",
    "Parameters: D_h = $(Printf.@sprintf("%.2e", D_h)) m²/s  |  θ_w = 1.0  |  Δt = 25 s  |  BC: C=0 at y=H",
    "Numerics: Timesteps: 6  |  Time range: 5-30 h  |  Scheme: Backward Euler"
]

p_legend_t6 = create_legend_textbox("LEGEND", legend_info_t6)

# Create EQUATION LEGEND
eq_t6_str = "∂(θ_w C)/∂t = D_h ∂²C/∂y²"
terms_t6 = Dict(
    "C" => "Aqueous solute concentration [mol/m³]",
    "θ_w" => "Water content [-]",
    "t" => "Time [s]",
    "D_h" => "Hydrodynamic diffusivity [m²/s]",
    "y" => "Vertical coordinate [m]"
)

# Save just the plot to SVG
savefig(p_prof_t6, "../output/T6_Terzaghi_Detailed_Verification.svg")
println("✓ Saved: ../output/T6_Terzaghi_Detailed_Verification.svg")

# Save legend and equation info to TXT file
open("../output/T6_Terzaghi_LEGEND.txt", "w") do f
    write(f, "="^80 * "\n")
    write(f, "T6 TERZAGHI VERIFICATION BENCHMARK\n")
    write(f, "="^80 * "\n\n")
    
    write(f, "MESH STRUCTURE\n")
    write(f, "-"^80 * "\n")
    write(f, "Grid Dimensions:\n")
    write(f, "  nx (x-direction) = 2 elements\n")
    write(f, "  ny (y-direction) = 50 elements\n")
    write(f, "  Total Elements = 100 (quasi-1D column: 2 × 50)\n")
    write(f, "  Total Nodes = 102 (Q4 bilinear elements)\n")
    write(f, "  Element Type: Q4 (Bilinear Quadrilateral)\n")
    write(f, "  Nodes per Element: 4 corner nodes\n")
    write(f, "  Integration Rule: 2×2 Gauss Quadrature (4 integration points per element)\n\n")
    write(f, "Physical Domain:\n")
    write(f, "  Height H = 0.20 m\n")
    write(f, "  Width Lx = 0.004 m\n")
    write(f, "  Aspect Ratio H/Lx = 50:1 (quasi-1D column)\n\n")
    
    write(f, "MATERIALS\n")
    write(f, "-"^80 * "\n")
    write(f, "Soil: VerifSoil\n")
    write(f, "  Granular Tortuosity τ = $(Printf.@sprintf("%.4f", soil.granular_tortuosity))\n")
    write(f, "Gas: VerifGas\n")
    write(f, "  Molecular Diffusion D_mol = $(Printf.@sprintf("%.2e", gas.diff_coefficient)) m²/s\n")
    write(f, "  Hydrodynamic Diffusivity:\n")
    write(f, "    D_h = τ × D_mol = $(Printf.@sprintf("%.2e", D_h)) m²/s\n\n")
    
    write(f, "PHYSICAL PARAMETERS\n")
    write(f, "-"^80 * "\n")
    write(f, "Transport Properties:\n")
    write(f, "  D_h (hydrodynamic diffusivity) = $(Printf.@sprintf("%.2e", D_h)) m²/s\n")
    write(f, "  θ_w (water content) = 1.0 (fully saturated)\n")
    write(f, "  v_y (advective velocity) = 0.0 m/s (PURE DIFFUSION)\n")
    write(f, "  Initial Condition: C₀ = 1.0 mol/m³ (uniform)\n\n")
    
    write(f, "NUMERICAL SCHEME\n")
    write(f, "-"^80 * "\n")
    write(f, "Time Integration:\n")
    write(f, "  Scheme: Backward Euler (implicit)\n")
    write(f, "  Δt (timestep) = 25 s\n")
    write(f, "  T_final = $(Printf.@sprintf("%.1f", T_final_t6)) s (≈ $(Printf.@sprintf("%.1f", T_final_t6/3600)) hours)\n")
    write(f, "  Stability: Unconditionally stable (any Δt allowed)\n")
    write(f, "  Accuracy: O(Δt) first-order in time\n\n")
    write(f, "Spatial Discretization:\n")
    write(f, "  Method: Q4 Finite Elements (bilinear shape functions)\n")
    write(f, "  Assembly: Lumped mass matrix (diagonal)\n")
    write(f, "  Linear Solver: Direct LU decomposition\n\n")
    
    write(f, "BOUNDARY CONDITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "Top (y = H = 0.20 m):\n")
    write(f, "  Type: Dirichlet (essential boundary condition)\n")
    write(f, "  Value: C = 0.0 mol/m³ (drained condition)\n")
    write(f, "Bottom (y = 0):\n")
    write(f, "  Type: Natural (Neumann, no flux)\n")
    write(f, "  Value: ∂C/∂y = 0 (impermeable base)\n\n")
    
    write(f, "VERIFICATION TIMESTEPS\n")
    write(f, "-"^80 * "\n")
    write(f, "Capture times (normalized by T_final = 0.5H²/D_h):\n")
    for (i, t) in enumerate(t_steps_t6)
        write(f, "  Step $(i): t = $(Printf.@sprintf("%.1f", t)) s = $(Printf.@sprintf("%.1f", t/3600)) h " *
                  "(T/T_final = $(Printf.@sprintf("%.2f", t/T_final_t6)))\n")
    end
    write(f, "\n")
    
    write(f, "GOVERNING EQUATION\n")
    write(f, "-"^80 * "\n")
    write(f, "Transport Equation (unsteady pure diffusion):\n\n")
    write(f, "  ∂(θ_w C)/∂t = D_h ∂²C/∂y²\n\n")
    write(f, "In this case (θ_w=1.0, no advection):\n")
    write(f, "  ∂C/∂t = D_h ∂²C/∂y²\n\n")
    write(f, "Standard heat/diffusion equation with:\n")
    write(f, "  - Initial condition: C(y,0) = 1.0 ∀y ∈ [0,H]\n")
    write(f, "  - Dirichlet BC: C(H,t) = 0 (top boundary)\n")
    write(f, "  - Natural BC: ∂C/∂y|_{y=0} = 0 (bottom boundary)\n")
    write(f, "  - Analytical Solution: Terzaghi-type exponential decay with sine basis\n\n")
    
    write(f, "TERM DEFINITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "  C            = Aqueous solute concentration [mol/m³]\n")
    write(f, "  θ_w          = Water content (saturation) [-]\n")
    write(f, "  D_h          = Hydrodynamic (effective) diffusivity [m²/s]\n")
    write(f, "  t            = Time [s]\n")
    write(f, "  y            = Vertical spatial coordinate [m]\n")
    write(f, "  τ            = Granular tortuosity [-]\n")
    write(f, "  D_mol        = Molecular diffusion coefficient [m²/s]\n")
    write(f, "  T_final      = Final simulation time [s]\n")
    write(f, "  H            = Domain height [m]\n")
    write(f, "  Lx           = Domain width [m]\n")
    write(f, "="^80 * "\n")
end
println("✓ Saved: ../output/T6_Terzaghi_LEGEND.txt")

# ─────────────────────────────────────────────────────────────────────────────
# T6b: ADVECTION
# ─────────────────────────────────────────────────────────────────────────────

println("Generating T6b Advection detailed plot...")

v_y_adv = -8.2e-7
θ_w_adv = 0.41
D_h_adv = 8.2e-10
ny_adv = 50
t_arrive = H / abs(v_y_adv)
T_final_adv = 0.8 * t_arrive
# Timesteps as multiples of 5 hours
t_steps_adv = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0] .* 3600.0  # Convert to seconds

mesh_adv = build_column_mesh(H, Lx, ny_adv)
y_adv = mesh_adv.coordinates[:, 2]

# Run advection step by step, capturing solutions at target times
initialize_shape_functions!(mesh_adv)
cache_adv = build_richards_cache(mesh_adv)
N_adv = mesh_adv.num_nodes
A_adv = spzeros(Float64, N_adv, N_adv)
R_adv = zeros(Float64, N_adv)

C_adv = zeros(Float64, N_adv)  # IC: C = 0
θw_adv = fill(θ_w_adv, N_adv)
vs_adv = zeros(Float64, N_adv, 2); vs_adv[:, 2] .= v_y_adv

P_bc_adv = ones(Int, N_adv)
top_n1_adv = 2*ny_adv + 1; top_n2_adv = 2*ny_adv + 2
P_bc_adv[top_n1_adv] = 0; P_bc_adv[top_n2_adv] = 0
C_presc_adv = zeros(Float64, N_adv)
C_presc_adv[top_n1_adv] = 1.0; C_presc_adv[top_n2_adv] = 1.0

Δt_adv = 100.0
C_fem_at_times_adv = Dict()  # Store solutions at each target time

t_current_adv = 0.0
next_target_idx_adv = 1

while t_current_adv < T_final_adv && next_target_idx_adv <= length(t_steps_adv)
    global C_adv, t_current_adv, next_target_idx_adv
    assemble_aqueous_concentration!(A_adv, R_adv, mesh_adv, cache_adv, D_h_adv,
                                    θw_adv, θw_adv, vs_adv, C_adv, Δt_adv,
                                    P_bc_adv, C_presc_adv)
    C_adv = A_adv \ R_adv
    t_current_adv = t_current_adv + Δt_adv
    
    # Store if close to target time
    target_t_adv = t_steps_adv[next_target_idx_adv]
    if t_current_adv >= target_t_adv - 0.5*Δt_adv
        C_fem_at_times_adv[target_t_adv] = copy(C_adv)
        next_target_idx_adv = next_target_idx_adv + 1
    end
end

# Create LaTeX info boxes
nx_adv = 2  # Fixed for quasi-1D column
info_mesh_adv = [
    @sprintf("Mesh Grid:"),
    @sprintf("  nx (x-direction) = %d", nx_adv),
    @sprintf("  ny (y-direction) = %d", ny_adv),
    @sprintf("  Total Elements = %d", ny_adv * nx_adv),
    @sprintf("  Total Nodes = %d", N_adv),
    @sprintf("  Element Type: Q4 (Bilinear Quadrilateral)"),
    @sprintf("  Nodes per Element: 4"),
    @sprintf("  Gauss Quadrature: 2×2 (4 integration points)"),
    @sprintf("Domain Dimensions:"),
    @sprintf("  H (vertical) = 0.20 m"),
    @sprintf("  L_x (horizontal) = 0.004 m"),
    @sprintf("Materials:"),
    @sprintf("  Soil: VerifSoil (granular_tortuosity = %.4f)", soil.granular_tortuosity),
    @sprintf("  Gas: VerifGas (diffusion = %.2e m²/s)", gas.diff_coefficient)
]

info_params_adv = [
    @sprintf("Transport Parameters:"),
    @sprintf("  v_y (advective velocity) = %.2e m/s (downward)", v_y_adv),
    @sprintf("  θ_w (water content) = %.2f [-]", θ_w_adv),
    @sprintf("  D_h (hydrodynamic diffusivity) = %.2e m²/s", D_h_adv),
    @sprintf("Time Integration:"),
    @sprintf("  Δt (timestep) = 100 s"),
    @sprintf("  T_final = %.1f s (≈ %.2f h)", T_final_adv, T_final_adv/3600),
    @sprintf("Boundary Conditions:"),
    @sprintf("  Top (y=H): C = 1.0 (Dirichlet - inlet)"),
    @sprintf("  Bottom (y=0): Natural BC (no flux)")
]

# Generate SINGLE plot with all timesteps overlaid (different colors)
colors_adv = [:blue, :red, :green, :purple, :orange, :brown]

p_prof_adv = plot(legend=:bottomright, grid=true, gridalpha=0.3,
                  xlabel="C", ylabel="y (m)", title="T6b Advection - All Timesteps",
                  size=(900, 500), margin=4mm)

for i in 1:6
    t_adv = t_steps_adv[i]
    C_fem_adv = C_fem_at_times_adv[t_adv]
    
    # Analytical: Moving front (solid step)
    y_front = H - abs(v_y_adv) * t_adv
    C_anal_adv = [y > y_front ? 1.0 : 0.0 for y in y_adv]
    
    color = colors_adv[i]
    t_hours = t_adv / 3600.0
    
    # Analytical: solid line (step function)
    plot!(p_prof_adv, C_anal_adv, y_adv, color=color, linewidth=2.5, linestyle=:solid,
          label="Analytical t=$(round(t_hours, digits=1))h", alpha=0.9)
    
    # FEM: scatter points
    scatter!(p_prof_adv, C_fem_adv, y_adv, color=color, markersize=2, markerstrokewidth=0,
             label="FEM t=$(round(t_hours, digits=1))h", alpha=0.7)
end

xlims!(p_prof_adv, -0.05, 1.1)
ylims!(p_prof_adv, 0, H)

# Create LEGEND BOXES (separate from plot)
legend_info_adv = [
    "Domain: Elements: $ny_adv  |  Nodes: $N_adv  |  H = 0.20 m  |  L_x = 0.004 m",
    "Parameters: v_y = $(Printf.@sprintf("%.2e", v_y_adv)) m/s  |  D_h = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s  |  θ_w = 0.41  |  Δt = 100 s",
    "Numerics: Timesteps: 6  |  Time range: 5-30 h  |  Scheme: Backward Euler"
]

p_legend_adv = create_legend_textbox("LEGEND", legend_info_adv)

# Create EQUATION LEGEND
eq_adv_str = "∂(θ_w C)/∂t + θ_w v_y ∂C/∂y = D_h ∂²C/∂y²"
terms_adv = Dict(
    "C" => "Aqueous solute concentration [mol/m³]",
    "θ_w" => "Water content [-]",
    "t" => "Time [s]",
    "v_y" => "Advective velocity in y-direction [m/s]",
    "D_h" => "Hydrodynamic diffusivity [m²/s]",
    "y" => "Vertical coordinate [m]"
)

# Save just the plot to SVG
savefig(p_prof_adv, "../output/T6b_Advection_Detailed_Verification.svg")
println("✓ Saved: ../output/T6b_Advection_Detailed_Verification.svg")

# Save legend and equation info to TXT file
open("../output/T6b_Advection_LEGEND.txt", "w") do f
    write(f, "="^80 * "\n")
    write(f, "T6b ADVECTION VERIFICATION BENCHMARK\n")
    write(f, "="^80 * "\n\n")
    
    write(f, "MESH STRUCTURE\n")
    write(f, "-"^80 * "\n")
    write(f, "Grid Dimensions:\n")
    write(f, "  nx (x-direction) = 2 elements\n")
    write(f, "  ny (y-direction) = 50 elements\n")
    write(f, "  Total Elements = 100 (quasi-1D column: 2 × 50)\n")
    write(f, "  Total Nodes = 102 (Q4 bilinear elements)\n")
    write(f, "  Element Type: Q4 (Bilinear Quadrilateral)\n")
    write(f, "  Nodes per Element: 4 corner nodes\n")
    write(f, "  Integration Rule: 2×2 Gauss Quadrature (4 integration points per element)\n\n")
    write(f, "Physical Domain:\n")
    write(f, "  Height H = 0.20 m\n")
    write(f, "  Width Lx = 0.004 m\n")
    write(f, "  Aspect Ratio H/Lx = 50:1 (quasi-1D column)\n\n")
    
    write(f, "MATERIALS\n")
    write(f, "-"^80 * "\n")
    write(f, "Soil: VerifSoil\n")
    write(f, "  Granular Tortuosity τ = $(Printf.@sprintf("%.4f", soil.granular_tortuosity))\n")
    write(f, "Gas: VerifGas\n")
    write(f, "  Molecular Diffusion D_mol = $(Printf.@sprintf("%.2e", gas.diff_coefficient)) m²/s\n")
    write(f, "  Hydrodynamic Diffusivity:\n")
    write(f, "    D_h = τ × D_mol = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s\n\n")
    
    write(f, "TRANSPORT PARAMETERS\n")
    write(f, "-"^80 * "\n")
    write(f, "Velocity:\n")
    write(f, "  v_y (advective velocity) = $(Printf.@sprintf("%.2e", v_y_adv)) m/s (downward)\n")
    write(f, "  Advective timescale: t_adv = H/|v_y| = $(Printf.@sprintf("%.1f", H/abs(v_y_adv))) s\n\n")
    write(f, "Diffusivity:\n")
    write(f, "  D_h (hydrodynamic diffusivity) = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s\n")
    write(f, "  Péclet number Pe = |v_y|H/D_h = $(Printf.@sprintf("%.1f", abs(v_y_adv)*H/D_h_adv))\n")
    write(f, "    → Pe >> 1: ADVECTION-DOMINATED (sharp front)\n\n")
    write(f, "Saturation:\n")
    write(f, "  θ_w (water content) = $(Printf.@sprintf("%.2f", θ_w_adv)) (unsaturated)\n")
    write(f, "  Initial Condition: C₀ = 0.0 (clean)\n\n")
    
    write(f, "NUMERICAL SCHEME\n")
    write(f, "-"^80 * "\n")
    write(f, "Time Integration:\n")
    write(f, "  Scheme: Backward Euler (implicit)\n")
    write(f, "  Δt (timestep) = 100 s\n")
    write(f, "  Courant Number: Cr = |v_y| Δt / H = $(Printf.@sprintf("%.4f", abs(v_y_adv)*100.0/H))\n")
    write(f, "  T_final = $(Printf.@sprintf("%.1f", T_final_adv)) s (≈ $(Printf.@sprintf("%.1f", T_final_adv/3600)) hours)\n")
    write(f, "  Stability: Unconditionally stable (any Δt allowed)\n\n")
    write(f, "Spatial Discretization:\n")
    write(f, "  Method: Q4 Finite Elements (bilinear shape functions)\n")
    write(f, "  Assembly: Lumped mass matrix (diagonal)\n")
    write(f, "  Stabilization: None (implicit scheme provides stability)\n\n")
    
    write(f, "BOUNDARY CONDITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "Top (y = H = 0.20 m):\n")
    write(f, "  Type: Dirichlet (essential boundary condition)\n")
    write(f, "  Value: C = 1.0 mol/m³ (contaminated inlet)\n")
    write(f, "Bottom (y = 0):\n")
    write(f, "  Type: Natural (Neumann, no flux)\n")
    write(f, "  Value: ∂C/∂y = 0 (impermeable base)\n\n")
    
    write(f, "VERIFICATION TIMESTEPS\n")
    write(f, "-"^80 * "\n")
    write(f, "Capture times (propagation of front):\n")
    for (i, t) in enumerate(t_steps_adv)
        y_front_calc = H - abs(v_y_adv) * t
        write(f, "  Step $(i): t = $(Printf.@sprintf("%.1f", t)) s = $(Printf.@sprintf("%.1f", t/3600)) h " *
                  "(Front y = $(Printf.@sprintf("%.4f", y_front_calc)) m)\n")
    end
    write(f, "\n")
    
    write(f, "GOVERNING EQUATION\n")
    write(f, "-"^80 * "\n")
    write(f, "Full Advection-Diffusion Equation:\n\n")
    write(f, "  ∂(θ_w C)/∂t + θ_w v_y ∂C/∂y = D_h ∂²C/∂y²\n\n")
    write(f, "In this benchmark (constant θ_w = 0.41):\n")
    write(f, "  ∂C/∂t + v_y ∂C/∂y = D_h/θ_w ∂²C/∂y²\n\n")
    write(f, "Physical Interpretation:\n")
    write(f, "  - Transport dominated by advective velocity (Péclet >> 1)\n")
    write(f, "  - Front moves downward with velocity ≈ v_y (sharp step)\n")
    write(f, "  - Diffusion causes minor smoothing at front edges\n")
    write(f, "  - Analytical solution: step function at y_front(t) = H - |v_y|·t\n\n")
    write(f, "Boundary Conditions:\n")
    write(f, "  - Dirichlet: C(H,t) = 1.0 (top boundary, all t > 0)\n")
    write(f, "  - Natural: ∂C/∂y|_{y=0} = 0 (bottom boundary)\n")
    write(f, "  - Initial: C(y,0) = 0.0 (clean domain)\n\n")
    
    write(f, "TERM DEFINITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "  C            = Aqueous solute concentration [mol/m³]\n")
    write(f, "  θ_w          = Water content (saturation) [-]\n")
    write(f, "  v_y          = Advective velocity in y-direction [m/s]\n")
    write(f, "  D_h          = Hydrodynamic (effective) diffusivity [m²/s]\n")
    write(f, "  t            = Time [s]\n")
    write(f, "  y            = Vertical spatial coordinate [m]\n")
    write(f, "  y_front      = Position of concentration front [m]\n")
    write(f, "  Pe           = Péclet number (advection/diffusion ratio) [-]\n")
    write(f, "  Cr           = Courant number [−]\n")
    write(f, "  H            = Domain height [m]\n")
    write(f, "  Lx           = Domain width [m]\n")
    write(f, "="^80 * "\n")
end
println("✓ Saved: ../output/T6b_Advection_LEGEND.txt")

# ─────────────────────────────────────────────────────────────────────────────
# T6c: MMS
# ─────────────────────────────────────────────────────────────────────────────

println("Generating T6c MMS detailed plot...")

C_0_mms = 1.0
A_amp = 0.2
Ω_mms = 2π / (6.0*3600.0)
v_y_mms = -8.2e-7
θ_w_mms = 0.41

C_exact_f = (y, t) -> C_0_mms + A_amp * cos(π*y/H) * sin(Ω_mms*t)

ny_mms = 80
T_period_mms = 6.0 * 3600.0
# Timesteps as multiples of 0.5 hours (30 minutes)
t_steps_mms = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0] .* 3600.0  # Convert to seconds

mesh_mms = build_column_mesh(H, Lx, ny_mms)
y_mms = mesh_mms.coordinates[:, 2]

initialize_shape_functions!(mesh_mms)
cache_mms = build_richards_cache(mesh_mms)
N_mms = mesh_mms.num_nodes
A_mms = spzeros(Float64, N_mms, N_mms)
R_mms = zeros(Float64, N_mms)

C_mms = [C_exact_f(y_mms[j], 0.0) for j in 1:N_mms]
θw_mms = fill(θ_w_mms, N_mms)
vs_mms = zeros(Float64, N_mms, 2); vs_mms[:, 2] .= v_y_mms

P_bc_mms = ones(Int, N_mms)
for nd in [1, 2, 2*ny_mms+1, 2*ny_mms+2]; P_bc_mms[nd] = 0; end

Δt_mms = 5.0
C_fem_at_times_mms = Dict()  # Store solutions at each target time

# Define source term for MMS
pL = π/H  # for spatial derivatives
src_mms = (y, t) -> begin
    cpy = cos(π*y/H)
    spy = sin(π*y/H)
    θ_w_mms * A_amp * (Ω_mms*cpy*cos(Ω_mms*t) - pL*v_y_mms*spy*sin(Ω_mms*t) + pL^2*D_h_adv*cpy*sin(Ω_mms*t))
end

t_current_mms = 0.0
next_target_idx_mms = 1
n_steps_mms = round(Int, T_period_mms / Δt_mms)

for step in 1:n_steps_mms
    global C_mms, t_current_mms, next_target_idx_mms
    t_new = step * Δt_mms
    C_presc_mms = [C_exact_f(y_mms[j], t_new) for j in 1:N_mms]
    assemble_aqueous_concentration!(A_mms, R_mms, mesh_mms, cache_mms, D_h_adv,
                                    θw_mms, θw_mms, vs_mms, C_mms, Δt_mms,
                                    P_bc_mms, C_presc_mms)
    # Add source term (lumped)
    for e in 1:mesh_mms.num_elements
        nodes_e = mesh_mms.elements[e, :]
        for p in 1:4
            wp = cache_mms.weights[p]
            dJ = cache_mms.detJ[e, p]
            Np_e = cache_mms.Np[p]
            y_p = sum(Np_e[a]*mesh_mms.coordinates[nodes_e[a],2] for a in 1:4)
            S_p = src_mms(y_p, t_new)
            for a in 1:4
                R_mms[nodes_e[a]] += Np_e[a] * S_p * wp * dJ
            end
        end
    end
    # Re-apply Dirichlet BCs
    for nd in [1, 2, 2*ny_mms+1, 2*ny_mms+2]
        R_mms[nd] = C_presc_mms[nd]
    end
    C_mms = A_mms \ R_mms
    t_current_mms = t_new
    
    # Store if close to target time
    if next_target_idx_mms <= length(t_steps_mms)
        target_t_mms = t_steps_mms[next_target_idx_mms]
        if t_current_mms >= target_t_mms - 0.5*Δt_mms
            C_fem_at_times_mms[target_t_mms] = copy(C_mms)
            next_target_idx_mms = next_target_idx_mms + 1
        end
    end
end

# Create LaTeX info boxes
nx_mms = 2  # Fixed for quasi-1D column
info_mesh_mms = [
    @sprintf("Mesh Grid:"),
    @sprintf("  nx (x-direction) = %d", nx_mms),
    @sprintf("  ny (y-direction) = %d", ny_mms),
    @sprintf("  Total Elements = %d", ny_mms * nx_mms),
    @sprintf("  Total Nodes = %d", N_mms),
    @sprintf("  Element Type: Q4 (Bilinear Quadrilateral)"),
    @sprintf("  Nodes per Element: 4"),
    @sprintf("  Gauss Quadrature: 2×2 (4 integration points)"),
    @sprintf("Domain Dimensions:"),
    @sprintf("  H (vertical) = 0.20 m"),
    @sprintf("  L_x (horizontal) = 0.004 m"),
    @sprintf("Materials:"),
    @sprintf("  Soil: VerifSoil (granular_tortuosity = %.4f)", soil.granular_tortuosity),
    @sprintf("  Gas: VerifGas (diffusion = %.2e m²/s)", gas.diff_coefficient)
]

info_params_mms = [
    @sprintf("Transport Parameters:"),
    @sprintf("  v_y (advective velocity) = %.2e m/s", v_y_mms),
    @sprintf("  θ_w (water content) = %.2f [-]", θ_w_mms),
    @sprintf("  D_h (hydrodynamic diffusivity) = %.2e m²/s", D_h_adv),
    @sprintf("Temporal Behavior:"),
    @sprintf("  T_period (oscillation period) = 6.0 h"),
    @sprintf("  Ω (angular frequency) = 2π/T_period = %.4e rad/s", Ω_mms),
    @sprintf("  C₀ (mean concentration) = %.2f mol/m³", C_0_mms),
    @sprintf("  A (amplitude) = %.2f mol/m³", A_amp),
    @sprintf("Time Integration:"),
    @sprintf("  Δt (timestep) = 5 s"),
    @sprintf("Boundary Conditions:"),
    @sprintf("  Dirichlet BCs at corners: C = C_exact(y,t)")
]


# Generate SINGLE plot with all timesteps overlaid (different colors)
colors_mms = [:blue, :red, :green, :purple, :orange, :brown]

p_prof_mms = plot(legend=:bottomright, grid=true, gridalpha=0.3,
                  xlabel="C", ylabel="y (m)", title="T6c MMS - All Timesteps",
                  size=(900, 500), margin=4mm)

for i in 1:6
    t_mms = t_steps_mms[i]
    C_exact = C_exact_f.(y_mms, t_mms)
    C_fem_mms = C_fem_at_times_mms[t_mms]
    
    color = colors_mms[i]
    t_hours = t_mms / 3600.0
    
    # Exact solution: solid line
    plot!(p_prof_mms, C_exact, y_mms, color=color, linewidth=2.5, linestyle=:solid,
          label="Exact t=$(round(t_hours, digits=2))h", alpha=0.9)
    
    # FEM: scatter points
    scatter!(p_prof_mms, C_fem_mms, y_mms, color=color, markersize=2, markerstrokewidth=0,
             label="FEM t=$(round(t_hours, digits=2))h", alpha=0.7)
end

xlims!(p_prof_mms, 0.7, 1.3)
ylims!(p_prof_mms, 0, H)

# Create LEGEND BOXES (separate from plot)
legend_info_mms = [
    "Domain: Elements: $ny_mms  |  Nodes: $N_mms  |  H = 0.20 m  |  L_x = 0.004 m",
    "Parameters: D_h = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s  |  θ_w = 0.41  |  Period = 6 h  |  Δt = 5 s",
    "Numerics: Timesteps: 6  |  Time range: 0.5-3 h  |  Scheme: Backward Euler  |  Source: S(y,t) included"
]

p_legend_mms = create_legend_textbox("LEGEND", legend_info_mms)

# Create EQUATION LEGEND
eq_mms_str = "∂(θ_w C)/∂t + θ_w v_y ∂C/∂y = D_h ∂²C/∂y² + S(y,t)"
terms_mms = Dict(
    "C" => "Aqueous solute concentration [mol/m³]",
    "θ_w" => "Water content [-]",
    "t" => "Time [s]",
    "v_y" => "Advective velocity [m/s]",
    "D_h" => "Hydrodynamic diffusivity [m²/s]",
    "y" => "Vertical coordinate [m]",
    "S(y,t)" => "Manufactured source term [mol/(m³·s)]"
)

# Save just the plot to SVG
savefig(p_prof_mms, "../output/T6c_MMS_Detailed_Verification.svg")
println("✓ Saved: ../output/T6c_MMS_Detailed_Verification.svg")

# Save legend and equation info to TXT file
open("../output/T6c_MMS_LEGEND.txt", "w") do f
    write(f, "="^80 * "\n")
    write(f, "T6c METHOD OF MANUFACTURED SOLUTIONS (MMS) VERIFICATION BENCHMARK\n")
    write(f, "="^80 * "\n\n")
    
    write(f, "MESH STRUCTURE\n")
    write(f, "-"^80 * "\n")
    write(f, "Grid Dimensions:\n")
    write(f, "  nx (x-direction) = 2 elements\n")
    write(f, "  ny (y-direction) = 80 elements  [REFINED for accuracy]\n")
    write(f, "  Total Elements = 160 (quasi-1D column: 2 × 80)\n")
    write(f, "  Total Nodes = 162 (Q4 bilinear elements)\n")
    write(f, "  Element Type: Q4 (Bilinear Quadrilateral)\n")
    write(f, "  Nodes per Element: 4 corner nodes\n")
    write(f, "  Integration Rule: 2×2 Gauss Quadrature (4 integration points per element)\n\n")
    write(f, "Physical Domain:\n")
    write(f, "  Height H = 0.20 m\n")
    write(f, "  Width Lx = 0.004 m\n")
    write(f, "  Aspect Ratio H/Lx = 50:1 (quasi-1D column)\n")
    write(f, "  Element size: Δy ≈ H/ny = $(Printf.@sprintf("%.2e", H/ny_mms)) m\n\n")
    
    write(f, "MATERIALS\n")
    write(f, "-"^80 * "\n")
    write(f, "Soil: VerifSoil\n")
    write(f, "  Granular Tortuosity τ = $(Printf.@sprintf("%.4f", soil.granular_tortuosity))\n")
    write(f, "Gas: VerifGas\n")
    write(f, "  Molecular Diffusion D_mol = $(Printf.@sprintf("%.2e", gas.diff_coefficient)) m²/s\n")
    write(f, "  Hydrodynamic Diffusivity:\n")
    write(f, "    D_h = τ × D_mol = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s\n\n")
    
    write(f, "TRANSPORT PARAMETERS\n")
    write(f, "-"^80 * "\n")
    write(f, "Advection:\n")
    write(f, "  v_y (advective velocity) = $(Printf.@sprintf("%.2e", v_y_mms)) m/s (downward)\n\n")
    write(f, "Diffusion:\n")
    write(f, "  D_h (hydrodynamic diffusivity) = $(Printf.@sprintf("%.2e", D_h_adv)) m²/s\n")
    write(f, "  Péclet number Pe = |v_y|H/D_h = $(Printf.@sprintf("%.1f", abs(v_y_mms)*H/D_h_adv))\n")
    write(f, "    → Pe >> 1: ADVECTION-DOMINATED\n\n")
    write(f, "Saturation:\n")
    write(f, "  θ_w (water content) = $(Printf.@sprintf("%.2f", θ_w_mms)) (unsaturated)\n\n")
    
    write(f, "MANUFACTURED SOLUTION\n")
    write(f, "-"^80 * "\n")
    write(f, "Exact Solution: C*(y,t) = C₀ + A·cos(πy/H)·sin(Ω·t)\n\n")
    write(f, "Parameters:\n")
    write(f, "  C₀ (mean value) = $(Printf.@sprintf("%.2f", C_0_mms)) mol/m³\n")
    write(f, "  A (amplitude) = $(Printf.@sprintf("%.2f", A_amp)) mol/m³\n")
    write(f, "  T_period (period) = 6.0 hours\n")
    write(f, "  Ω (angular frequency) = 2π/T_period = $(Printf.@sprintf("%.4e", Ω_mms)) rad/s\n\n")
    write(f, "Spatial Component: cos(πy/H)\n")
    write(f, "  - Nodal boundary conditions: C(0,t) = C_exact(0,t), C(H,t) = C_exact(H,t)\n")
    write(f, "  - Interior oscillation in y-direction\n\n")
    write(f, "Temporal Component: sin(Ω·t)\n")
    write(f, "  - Harmonic oscillation at period T = 6 hours\n")
    write(f, "  - Phase lag in FEM solution relative to exact (test convergence)\n\n")
    
    write(f, "SOURCE TERM DERIVATION\n")
    write(f, "-"^80 * "\n")
    write(f, "To satisfy the PDE exactly, the source term is:\n")
    write(f, "  S(y,t) = ∂(θ_w C*)/∂t + θ_w v_y ∂C*/∂y - D_h ∂²C*/∂y²\n\n")
    write(f, "Components:\n")
    write(f, "  ∂C*/∂t = A·cos(πy/H)·Ω·cos(Ω·t)\n")
    write(f, "  ∂C*/∂y = -A·π/H·sin(πy/H)·sin(Ω·t)\n")
    write(f, "  ∂²C*/∂y² = -A·(π/H)²·cos(πy/H)·sin(Ω·t)\n\n")
    write(f, "Result:\n")
    write(f, "  S(y,t) = θ_w [A·Ω·cos(πy/H)·cos(Ω·t) - v_y·A·π/H·sin(πy/H)·sin(Ω·t)\n")
    write(f, "           + D_h·A·(π/H)²·cos(πy/H)·sin(Ω·t)]\n\n")
    
    write(f, "NUMERICAL SCHEME\n")
    write(f, "-"^80 * "\n")
    write(f, "Time Integration:\n")
    write(f, "  Scheme: Backward Euler (implicit)\n")
    write(f, "  Δt (timestep) = 5 s\n")
    write(f, "  Courant Number: Cr = |v_y| Δt / H = $(Printf.@sprintf("%.4f", abs(v_y_mms)*5.0/H))\n")
    write(f, "  Total Steps per Period: N_steps = T_period/Δt = $(Printf.@sprintf("%.0f", T_period_mms/5.0))\n\n")
    write(f, "Spatial Discretization:\n")
    write(f, "  Method: Q4 Finite Elements (bilinear shape functions)\n")
    write(f, "  Assembly: Lumped mass matrix (diagonal)\n")
    write(f, "  Source: Lumped assembly at Gauss points\n")
    write(f, "  Linear Solver: Direct LU decomposition\n\n")
    
    write(f, "BOUNDARY CONDITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "Corners and edges (essential Dirichlet):\n")
    write(f, "  C(0, t) = C_exact(0, t) = C₀ + A·cos(0)·sin(Ω·t) = C₀ + A·sin(Ω·t)\n")
    write(f, "  C(H, t) = C_exact(H, t) = C₀ + A·cos(π)·sin(Ω·t) = C₀ - A·sin(Ω·t)\n\n")
    write(f, "Interior points:\n")
    write(f, "  Governed by advection-diffusion equation with source term\n\n")
    
    write(f, "VERIFICATION TIMESTEPS\n")
    write(f, "-"^80 * "\n")
    write(f, "Capture times (phases of oscillation):\n")
    for (i, t) in enumerate(t_steps_mms)
        phase = (Ω_mms * t) * 180 / π  # Phase in degrees
        write(f, "  Step $(i): t = $(Printf.@sprintf("%.1f", t)) s = $(Printf.@sprintf("%.2f", t/3600)) h " *
                  "(Phase = $(Printf.@sprintf("%.0f", phase))°)\n")
    end
    write(f, "\n")
    
    write(f, "GOVERNING EQUATION\n")
    write(f, "-"^80 * "\n")
    write(f, "Advection-Diffusion with Source:\n\n")
    write(f, "  ∂(θ_w C)/∂t + θ_w v_y ∂C/∂y = D_h ∂²C/∂y² + S(y,t)\n\n")
    write(f, "In this benchmark (constant θ_w = 0.41):\n")
    write(f, "  ∂C/∂t + v_y ∂C/∂y = D_h/θ_w ∂²C/∂y² + S(y,t)/θ_w\n\n")
    write(f, "Purpose:\n")
    write(f, "  - Verify FEM solver on problem with KNOWN EXACT SOLUTION\n")
    write(f, "  - Test combined advection-diffusion-source coupling\n")
    write(f, "  - Measure spatial and temporal convergence rates\n")
    write(f, "  - Validate source term assembly in solver\n\n")
    
    write(f, "TERM DEFINITIONS\n")
    write(f, "-"^80 * "\n")
    write(f, "  C            = Aqueous solute concentration [mol/m³]\n")
    write(f, "  C*           = Manufactured (exact) solution [mol/m³]\n")
    write(f, "  θ_w          = Water content (saturation) [-]\n")
    write(f, "  v_y          = Advective velocity in y-direction [m/s]\n")
    write(f, "  D_h          = Hydrodynamic (effective) diffusivity [m²/s]\n")
    write(f, "  S(y,t)       = Manufactured source/sink term [mol/(m³·s)]\n")
    write(f, "  t            = Time [s]\n")
    write(f, "  y            = Vertical spatial coordinate [m]\n")
    write(f, "  C₀           = Mean concentration (base value) [mol/m³]\n")
    write(f, "  A            = Amplitude of spatial oscillation [mol/m³]\n")
    write(f, "  T_period     = Period of temporal oscillation [s]\n")
    write(f, "  Ω            = Angular frequency (2π/T_period) [rad/s]\n")
    write(f, "  Pe           = Péclet number (advection/diffusion ratio) [-]\n")
    write(f, "  H            = Domain height [m]\n")
    write(f, "  Lx           = Domain width [m]\n")
    write(f, "="^80 * "\n")
end
println("✓ Saved: ../output/T6c_MMS_LEGEND.txt")

println("\n" * "="^70)
println("✓ ALL DETAILED VERIFICATION PLOTS GENERATED (SVG format)")
println("="^70)
println("  Layout: [Main Plot] | [Legend] | [Equation Legend]")
println("  1. ../output/T6_Terzaghi_Detailed_Verification.svg")
println("  2. ../output/T6b_Advection_Detailed_Verification.svg")
println("  3. ../output/T6c_MMS_Detailed_Verification.svg")
println("="^70)
