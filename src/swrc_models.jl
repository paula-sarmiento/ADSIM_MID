#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.x.x
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# SWRC (Soil Water Retention Curve) Model Implementations
# Provides Van Genuchten and Cavalcante models for water flow
#______________________________________________________

using Printf

"""
    create_van_genuchten_model(params::Dict{String, Float64})

Create Van Genuchten SWRC model closures.

# Arguments
- `params::Dict{String, Float64}`: Dictionary containing SWRC parameters

  Required keys:
  - "theta_s": Saturated water content [-]
  - "theta_r": Residual water content [-]
  - "alpha": Van Genuchten parameter [1/m]
  - "n_param": Van Genuchten shape parameter [-]
  - "K_sat": Saturated hydraulic conductivity [m/s]

# Returns
- `NamedTuple`: Contains closures `K_h`, `theta_h`, `c_s`, `D_w`

# Formula
Effective saturation:
  S_e = (1 + (α|h|)^n)^(-m), where m = 1 - 1/n

Water content:
  θ(h) = θ_r + (θ_s - θ_r) × S_e

Hydraulic conductivity:
  K(h) = K_s × S_e^0.5 × [1 - (1 - S_e^(1/m))^m]²

Water capacity:
  c_s(h) = dθ/dh

Water diffusivity:
  D_w(h) = K(h) / c_s(h)
"""
function create_van_genuchten_model(params::Dict{String, Float64})
    
    # Validate required parameters
    required_keys = ["theta_s", "theta_r", "alpha", "n_param", "K_sat"]
    for key in required_keys
        if !haskey(params, key)
            error("Van Genuchten model requires key: $key")
        end
    end
    
    #______________________________________________________
    # Extract and compute SWRC parameters
    #______________________________________________________
    theta_s = params["theta_s"]
    theta_r = params["theta_r"]
    alpha = params["alpha"]
    n = params["n_param"]
    K_sat = params["K_sat"]
    m = 1.0 - 1.0 / n
    
    # Validate parameter ranges
    if theta_r >= theta_s
        error("Residual water content (theta_r) must be < saturated content (theta_s)")
    end
    if alpha <= 0.0
        error("Van Genuchten parameter alpha must be positive")
    end
    if n <= 1.0
        error("Van Genuchten parameter n must be > 1.0")
    end
    if K_sat <= 0.0
        error("Saturated hydraulic conductivity must be positive")
    end
    
    #______________________________________________________
    # Effective saturation closure
    #______________________________________________________
    function S_e(h)
        if h >= 0.0
            return 1.0  # Saturated
        else
            return (1.0 + (alpha * abs(h))^n)^(-m)
        end
    end
    
    #______________________________________________________
    # Hydraulic conductivity K(h) closure
    #______________________________________________________
    function K_h(h)
        if h >= 0.0
            return K_sat  # Saturated
        else
            Se = S_e(h)
            return K_sat * sqrt(Se) * (1.0 - (1.0 - Se^(1.0/m))^m)^2
        end
    end
    
    #______________________________________________________
    # Water content θ(h) closure
    #______________________________________________________
    function theta_h(h)
        if h >= 0.0
            return theta_s  # Saturated
        else
            Se = S_e(h)
            return theta_r + (theta_s - theta_r) * Se
        end
    end
    
    #______________________________________________________
    # Water capacity c_s(h) = dθ/dh closure
    #______________________________________________________
    function c_s(h)
        # Numerical differentiation of theta_h
        dh = max(1.0e-6, abs(h) * 1.0e-8)  # Adaptive step size
        dtheta = (theta_h(h + dh/2.0) - theta_h(h - dh/2.0)) / dh
        
        # Avoid negative or near-zero values
        if abs(dtheta) < 1.0e-15
            return 1.0e-15  # Minimum capacity to avoid division by zero
        end
        return dtheta
    end
    
    #______________________________________________________
    # Water diffusivity D_w(h) closure
    #______________________________________________________
    function D_w(h)
        K = K_h(h)
        cs = c_s(h)
        
        # Avoid division by very small capacity
        if abs(cs) < 1.0e-15
            return 1.0e-15  # Return small positive value
        end
        return K / cs
    end
    
    return (K_h=K_h, theta_h=theta_h, c_s=c_s, D_w=D_w)
end

"""
    create_cavalcante_model(params::Dict{String, Float64})

Create Cavalcante SWRC model closures.

# Arguments
- `params::Dict{String, Float64}`: Dictionary containing SWRC parameters

  Required keys:
  - "theta_s": Saturated water content [-]
  - "theta_r": Residual water content [-]
  - "alpha": Cavalcante parameter (air entry pressure inverse) [1/m]
  - "lambda_param": Pore size distribution index [-]
  - "K_sat": Saturated hydraulic conductivity [m/s]

# Returns
- `NamedTuple`: Contains closures `K_h`, `theta_h`, `c_s`, `D_w`

# Formula
Effective saturation:
  S_e = (1 + (α|h|)^λ)^(-1), where λ = lambda_param

Water content:
  θ(h) = θ_r + (θ_s - θ_r) × S_e

Hydraulic conductivity (Brooks-Corey-like):
  K(h) = K_s × S_e^(2/λ + 3)
"""
function create_cavalcante_model(params::Dict{String, Float64})
    
    # Validate required parameters
    required_keys = ["theta_s", "theta_r", "alpha", "lambda_param", "K_sat"]
    for key in required_keys
        if !haskey(params, key)
            error("Cavalcante model requires key: $key")
        end
    end
    
    #______________________________________________________
    # Extract SWRC parameters
    #______________________________________________________
    theta_s = params["theta_s"]
    theta_r = params["theta_r"]
    alpha = params["alpha"]
    lambda_param = params["lambda_param"]
    K_sat = params["K_sat"]
    
    # Validate parameter ranges
    if theta_r >= theta_s
        error("Residual water content (theta_r) must be < saturated content (theta_s)")
    end
    if alpha <= 0.0
        error("Cavalcante parameter alpha must be positive")
    end
    if lambda_param <= 0.0
        error("Pore size distribution index lambda must be positive")
    end
    if K_sat <= 0.0
        error("Saturated hydraulic conductivity must be positive")
    end
    
    #______________________________________________________
    # Effective saturation closure
    #______________________________________________________
    function S_e(h)
        if h >= 0.0
            return 1.0  # Saturated
        else
            return (1.0 + (alpha * abs(h))^lambda_param)^(-1.0)
        end
    end
    
    #______________________________________________________
    # Hydraulic conductivity K(h) closure
    #______________________________________________________
    function K_h(h)
        if h >= 0.0
            return K_sat  # Saturated
        else
            Se = S_e(h)
            exponent = 2.0 / lambda_param + 3.0
            return K_sat * Se^exponent
        end
    end
    
    #______________________________________________________
    # Water content θ(h) closure
    #______________________________________________________
    function theta_h(h)
        if h >= 0.0
            return theta_s  # Saturated
        else
            Se = S_e(h)
            return theta_r + (theta_s - theta_r) * Se
        end
    end
    
    #______________________________________________________
    # Water capacity c_s(h) = dθ/dh closure
    #______________________________________________________
    function c_s(h)
        # Numerical differentiation of theta_h
        dh = max(1.0e-6, abs(h) * 1.0e-8)  # Adaptive step size
        dtheta = (theta_h(h + dh/2.0) - theta_h(h - dh/2.0)) / dh
        
        # Avoid negative or near-zero values
        if abs(dtheta) < 1.0e-15
            return 1.0e-15  # Minimum capacity to avoid division by zero
        end
        return dtheta
    end
    
    #______________________________________________________
    # Water diffusivity D_w(h) closure
    #______________________________________________________
    function D_w(h)
        K = K_h(h)
        cs = c_s(h)
        
        # Avoid division by very small capacity
        if abs(cs) < 1.0e-15
            return 1.0e-15  # Return small positive value
        end
        return K / cs
    end
    
    return (K_h=K_h, theta_h=theta_h, c_s=c_s, D_w=D_w)
end

"""
    create_swrc_model(model_name::String, params::Dict{String, Float64})

Factory function to create SWRC model based on model name.

# Arguments
- `model_name::String`: Name of SWRC model ("Van_Genuchten" or "Cavalcante")
- `params::Dict{String, Float64}`: Model parameters dictionary

# Returns
- `NamedTuple`: Contains closures `K_h`, `theta_h`, `c_s`, `D_w`

# Throws
- `error`: If model_name is not recognized
"""
function create_swrc_model(model_name::String, params::Dict{String, Float64})
    
    if model_name == "Van_Genuchten"
        return create_van_genuchten_model(params)
    elseif model_name == "Cavalcante"
        return create_cavalcante_model(params)
    else
        error("Unknown SWRC model: $model_name. Available: Van_Genuchten, Cavalcante")
    end
end
