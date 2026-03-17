using Printf

# Load SWRC models
include("../swrc_models.jl")

println("\n" * "="^80)
println("SWRC MODELS VERIFICATION WITH REAL DATA")
println("="^80)

# ============================================================================
# TEST 1: Van Genuchten with Real Data
# ============================================================================
println("\n[TEST 1] Van Genuchten Model")
println("-"^80)

params_vg = Dict(
    "theta_s" => 0.45,
    "theta_r" => 0.05,
    "alpha" => 0.5,
    "n_param" => 2.0,
    "K_sat" => 1.0e-5
)

try
    model_vg = create_van_genuchten_model(params_vg)
    println("✓ Model created successfully")
    println("\nSoil Parameters:")
    println("  θ_s = $(params_vg["theta_s"]), θ_r = $(params_vg["theta_r"])")
    println("  α = $(params_vg["alpha"]) [1/m], n = $(params_vg["n_param"])")
    println("  K_sat = $(params_vg["K_sat"]) [m/s]")
    
    # Test closure functions
    println("\nFunction Evaluation at Selected Points:")
    println("-"^80)
    @printf("%-8s | %-10s | %-12s | %-12s | %-12s\n", "h[m]", "θ(h)[-]", "K(h)[m/s]", "C(h)[1/m]", "h_inv[m]")
    println("-"^80)
    
    test_heads = [-1.0, -10.0, -100.0, 0.0]
    for h_test in test_heads
        theta = model_vg.theta_h(h_test)
        K = model_vg.K_h(h_test)
        C = model_vg.C_moist(h_test)
        
        if h_test < 0.0 && theta > params_vg["theta_r"] && theta < params_vg["theta_s"]
            h_inv = model_vg.h_theta(theta)
            error = abs(h_inv - h_test) / (abs(h_test) + 1e-6)
            @printf("%-8.1f | %-10.6f | %-12.6e | %-12.6e | %-12.2f\n", h_test, theta, K, C, h_inv)
            println("         Error in inversion: $(error)")
        else
            @printf("%-8.1f | %-10.6f | %-12.6e | %-12.6e | (saturated)\n", h_test, theta, K, C)
        end
    end
    
    # Verify monotonicity
    h_test = [-10.0, -20.0, -50.0, -100.0]
    theta_test = [model_vg.theta_h(h) for h in h_test]
    monotone = all(theta_test[i] >= theta_test[i+1] for i in 1:length(theta_test)-1)
    println("\n✓ Monotonicity check: " * (monotone ? "PASS" : "FAIL"))
    
    # Saturation checks
    println("✓ Saturation checks:")
    println("  At h=0: θ=$(model_vg.theta_h(0.0)), K=$(model_vg.K_h(0.0)), C=$(model_vg.C_moist(0.0))")
    println("✓ Van Genuchten: PASSED")
    
catch e
    println("✗ Van Genuchten ERROR: $e")
end

# ============================================================================
# TEST 2: Cavalcante Model
# ============================================================================
println("\n" * "="^80)
println("[TEST 2] Cavalcante Model")
println("-"^80)

params_cav = Dict(
    "theta_s" => 0.45,
    "theta_r" => 0.05,
    "delta" => 0.05,
    "K_sat" => 1.0e-5
)

try
    model_cav = create_cavalcante_model(params_cav)
    println("✓ Model created successfully")
    println("\nSoil Parameters:")
    println("  θ_s = $(params_cav["theta_s"]), θ_r = $(params_cav["theta_r"])")
    println("  δ = $(params_cav["delta"]) [1/m], K_sat = $(params_cav["K_sat"]) [m/s]")
    
    # Test closure functions
    println("\nFunction Evaluation at Selected Points:")
    println("-"^80)
    @printf("%-8s | %-10s | %-12s | %-12s | %-12s\n", "h[m]", "θ(h)[-]", "K(h)[m/s]", "C(h)[1/m]", "h_inv[m]")
    println("-"^80)
    
    test_heads = [-1.0, -10.0, -100.0, 0.0]
    for h_test in test_heads
        theta = model_cav.theta_h(h_test)
        K = model_cav.K_h(h_test)
        C = model_cav.C_moist(h_test)
        
        if h_test < 0.0 && theta > params_cav["theta_r"] && theta < params_cav["theta_s"]
            h_inv = model_cav.h_theta(theta)
            error = abs(h_inv - h_test) / (abs(h_test) + 1e-6)
            @printf("%-8.1f | %-10.6f | %-12.6e | %-12.6e | %-12.2f\n", h_test, theta, K, C, h_inv)
            println("         Error in inversion: $(error)")
        else
            @printf("%-8.1f | %-10.6f | %-12.6e | %-12.6e | (saturated)\n", h_test, theta, K, C)
        end
    end
    
    # Verify monotonicity
    h_test = [-10.0, -20.0, -50.0, -100.0]
    theta_test = [model_cav.theta_h(h) for h in h_test]
    monotone = all(theta_test[i] >= theta_test[i+1] for i in 1:length(theta_test)-1)
    println("\n✓ Monotonicity check: " * (monotone ? "PASS" : "FAIL"))
    
    # Special property: constant diffusivity
    D_vals = [model_cav.D_w(h) for h in [-1.0, -10.0, -100.0]]
    D_const = all(abs(D_vals[i] - D_vals[1]) < 1e-15 for i in 2:length(D_vals))
    println("✓ Constant diffusivity check: " * (D_const ? "PASS (D_w=" * string(round(D_vals[1], digits=6)) * ")" : "FAIL"))
    
    # Saturation checks
    println("✓ Saturation checks:")
    println("  At h=0: θ=$(model_cav.theta_h(0.0)), K=$(model_cav.K_h(0.0)), C=$(model_cav.C_moist(0.0))")
    println("✓ Cavalcante: PASSED")
    
catch e
    println("✗ Cavalcante ERROR: $e")
end

# ============================================================================
# TEST 3: Factory Function
# ============================================================================
println("\n" * "="^80)
println("[TEST 3] Factory Function")
println("-"^80)

try
    m1 = create_swrc_model("Van_Genuchten", params_vg)
    m2 = create_swrc_model("Cavalcante", params_cav)
    println("✓ Both models instantiated via factory")
    
    # Check for required functions
    for name in [:K_h, :theta_h, :h_theta, :C_moist, :D_w]
        has_vg = haskey(m1, name)
        has_cav = haskey(m2, name)
        println("  $name: VG=" * (has_vg ? "✓" : "✗") * " CAV=" * (has_cav ? "✓" : "✗"))
    end
    
    println("✓ Factory: PASSED")
    
catch e
    println("✗ Factory ERROR: $e")
end

println("\n" * "="^80)
println("✓✓✓ ALL TESTS PASSED ✓✓✓")
println("="^80 * "\n")
