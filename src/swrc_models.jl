#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and
# Modification
# v0.x.x
# Author: Paula Sarmiento
#______________________________________________________

#______________________________________________________
# SWRC (Soil Water Retention Curve) Model Implementations
# Provides Van Genuchten and Cavalcante models for water flow
#
# All models return a NamedTuple with the following closures:
#   K_h(h)      -> hydraulic conductivity           [L/T]
#   theta_h(h)  -> volumetric water content         [-]
#   h_theta(θ)  -> pressure head from water content [L]
#   C_moist(h)  -> specific moisture capacity dθ/dh [1/L]
#   D_w(h)      -> hydraulic diffusivity K/C        [L²/T]
#
# Convention: h < 0 in the unsaturated zone, h >= 0 at saturation.
#______________________________________________________

using Printf


"""
    create_van_genuchten_model(params::Dict{String, Float64})

Van Genuchten (1980) SWRC with Mualem (1976) conductivity model.

# Arguments
- `params::Dict{String, Float64}`: Dictionary containing model parameters.

  Required keys:
  - `"theta_s"` : Saturated volumetric water content       [-]
  - `"theta_r"` : Residual volumetric water content        [-]
  - `"alpha"`   : Inverse air-entry pressure               [1/L]
  - `"n_param"` : Pore-size distribution index (> 1)       [-]
  - `"K_sat"`   : Saturated hydraulic conductivity         [L/T]

  Optional keys:
  - `"L_mualem"` : Mualem pore-connectivity exponent (default = 0.5) [-]

# Returns
- `NamedTuple` with closures: `K_h`, `theta_h`, `h_theta`, `C_moist`, `D_w`

# Governing equations

Effective saturation (h < 0):
```
S_e(h) = [1 + (α|h|)^n]^(−m),   m = 1 − 1/n
```

Volumetric water content:
```
θ(h) = θ_r + (θ_s − θ_r) S_e(h)
```

Specific moisture capacity — ANALYTICAL derivative dθ/dh:
```
C(h) = (θ_s − θ_r) m n α^n |h|^(n−1) / [1 + (α|h|)^n]^(m+1)   for h < 0
C(h) = 0                                                           for h ≥ 0
```

Hydraulic conductivity — Mualem (1976):
```
K(h) = K_s S_e^L [1 − (1 − S_e^(1/m))^m]²
```
where L = `L_mualem` (default 0.5).

Inverse retention (analytical, exact):
```
h(θ) = −(1/α) [(S_e^(−1/m) − 1)^(1/n)],   S_e = (θ − θ_r)/(θ_s − θ_r)
```

Hydraulic diffusivity:
```
D_w(h) = K(h) / C(h)
```
"""
function create_van_genuchten_model(params::Dict{String, Float64})

    #______________________________________________________
    # Validate required parameters
    #______________________________________________________
    required_keys = ["theta_s", "theta_r", "alpha", "n_param", "K_sat"]
    for key in required_keys
        if !haskey(params, key)
            error("Van Genuchten model requires key: $key")
        end
    end

    #______________________________________________________
    # Extract parameters
    #______________________________________________________
    theta_s = params["theta_s"]
    theta_r = params["theta_r"]
    alpha   = params["alpha"]
    n       = params["n_param"]
    K_sat   = params["K_sat"]
    L       = get(params, "L_mualem", 0.5)   # Mualem exponent, default 0.5
    m       = 1.0 - 1.0 / n

    #______________________________________________________
    # Validate parameter ranges
    #______________________________________________________
    theta_r >= theta_s  && error("theta_r must be < theta_s")
    alpha   <= 0.0      && error("alpha must be positive")
    n       <= 1.0      && error("n_param must be > 1.0")
    K_sat   <= 0.0      && error("K_sat must be positive")
    L       <  0.0      && error("L_mualem must be non-negative")

    #______________________________________________________
    # Effective saturation  S_e(h)  [internal helper]
    #______________________________________________________
    function S_e(h)
        h >= 0.0 && return 1.0
        return (1.0 + (alpha * abs(h))^n)^(-m)
    end

    #______________________________________________________
    # Volumetric water content  θ(h)
    #______________________________________________________
    function theta_h(h)
        h >= 0.0 && return theta_s
        return theta_r + (theta_s - theta_r) * S_e(h)
    end

    #______________________________________________________
    # Hydraulic conductivity  K(h)  — Mualem (1976)
    # K(h) = K_s S_e^L [1 − (1 − S_e^(1/m))^m]²
    #______________________________________________________
    function K_h(h)
        h >= 0.0 && return K_sat
        Se    = S_e(h)
        inner = clamp(1.0 - Se^(1.0 / m), 0.0, 1.0)
        return K_sat * Se^L * (1.0 - inner^m)^2
    end

    #______________________________________________________
    # Specific moisture capacity  C(h) = dθ/dh  — ANALYTICAL
    #
    # Starting from:  θ(h) = θ_r + (θ_s − θ_r)[1 + (α|h|)^n]^(−m)
    #
    # dθ/dh = (θ_s − θ_r) × d/dh {[1 + (α|h|)^n]^(−m)}
    #
    # Let u = (α|h|)^n.  For h < 0, |h| = −h, so d|h|/dh = −1.
    # du/dh = n α^n |h|^(n−1) × (−1)
    #
    # d/dh {(1+u)^(−m)} = −m(1+u)^(−m−1) du/dh
    #                    = −m(1+u)^(−m−1) × n α^n |h|^(n−1) × (−1)
    #                    =  m n α^n |h|^(n−1) / (1 + (α|h|)^n)^(m+1)
    #
    # Therefore:
    # C(h) = (θ_s − θ_r) × m n α^n |h|^(n−1) / [1 + (α|h|)^n]^(m+1)
    #______________________________________________________
    function C_moist(h)
        h >= 0.0 && return 0.0
        h_abs = abs(h)
        # Guard against h = 0 approached from below
        h_abs < 1.0e-14 && return 0.0
        numerator   = m * n * alpha^n * h_abs^(n - 1.0)
        denominator = (1.0 + (alpha * h_abs)^n)^(m + 1.0)
        return (theta_s - theta_r) * numerator / denominator
    end

    #______________________________________________________
    # Inverse retention  h(θ)  — exact analytical inversion
    #
    # From S_e = (θ − θ_r)/(θ_s − θ_r) and
    #      S_e = [1 + (α|h|)^n]^(−m):
    #
    # S_e^(−1/m) = 1 + (α|h|)^n
    # (α|h|)^n   = S_e^(−1/m) − 1
    # |h|         = (1/α)[S_e^(−1/m) − 1]^(1/n)
    # h           = −(1/α)[S_e^(−1/m) − 1]^(1/n)   (negative in unsat. zone)
    #______________________________________________________
    function h_theta(theta)
        theta >= theta_s - 1.0e-12 && return 0.0
        theta <= theta_r + 1.0e-12 && return -1.0e8   # very dry limit

        Se_val = (theta - theta_r) / (theta_s - theta_r)
        base   = Se_val^(-1.0 / m) - 1.0
        base < 1.0e-15 && return 0.0
        return -(1.0 / alpha) * base^(1.0 / n)
    end

    #______________________________________________________
    # Hydraulic diffusivity  D_w(h) = K(h) / C(h)
    #______________________________________________________
    function D_w(h)
        h >= 0.0 && return 0.0   # saturated: capillary diffusivity undefined
        C = C_moist(h)
        C < 1.0e-15 && return 0.0
        return K_h(h) / C
    end

    return (K_h=K_h, theta_h=theta_h, h_theta=h_theta, c_s=C_moist, D_w=D_w)
end


"""
    create_cavalcante_model(params::Dict{String, Float64})

Cavalcante & Zornberg (2017) exponential SWRC model.

# Arguments
- `params::Dict{String, Float64}`: Dictionary containing model parameters.

  Required keys:
  - `"theta_s"` : Saturated volumetric water content  [-]
  - `"theta_r"` : Residual volumetric water content   [-]
  - `"delta"`   : Exponential desaturation parameter  [1/L]
  - `"K_sat"`   : Saturated hydraulic conductivity    [L/T]

# Returns
- `NamedTuple` with closures: `K_h`, `theta_h`, `h_theta`, `C_moist`, `D_w`

# Governing equations

Volumetric water content (exponential retention):
```
θ(h) = θ_r + (θ_s − θ_r) exp(−δ|h|)
```

Hydraulic conductivity — Cavalcante & Zornberg (2017), Eq. (14):
```
K(h) = K_s exp(−δ|h|)
```
Note: the relative permeability is K_r = exp(−δ|h|), i.e. exponent −δ|h|,
NOT −2δ|h|.  The quadratic form would apply to a Burdine-type model, which
is not the formulation in the original paper.

Specific moisture capacity — ANALYTICAL derivative dθ/dh:
```
C(h) = (θ_s − θ_r) δ exp(−δ|h|)   for h < 0
C(h) = 0                            for h ≥ 0
```

Derivation:  θ(h) = θ_r + (θ_s − θ_r) exp(−δ|h|).
For h < 0, |h| = −h, so d|h|/dh = −1.
dθ/dh = (θ_s − θ_r) exp(−δ|h|) × (−δ) × (−1) = (θ_s − θ_r) δ exp(−δ|h|)

Inverse retention (analytical, exact):
```
h(θ) = −(1/δ) ln[(θ − θ_r)/(θ_s − θ_r)]
```

Hydraulic diffusivity:
```
D_w(h) = K(h) / C(h) = K_s / [(θ_s − θ_r) δ]   (constant in h!)
```
"""
function create_cavalcante_model(params::Dict{String, Float64})

    #______________________________________________________
    # Validate required parameters
    #______________________________________________________
    required_keys = ["theta_s", "theta_r", "delta", "K_sat"]
    for key in required_keys
        if !haskey(params, key)
            error("Cavalcante model requires key: $key")
        end
    end

    #______________________________________________________
    # Extract parameters
    #______________________________________________________
    theta_s = params["theta_s"]
    theta_r = params["theta_r"]
    delta   = params["delta"]
    K_sat   = params["K_sat"]

    #______________________________________________________
    # Validate parameter ranges
    #______________________________________________________
    theta_r >= theta_s && error("theta_r must be < theta_s")
    delta   <= 0.0     && error("delta must be positive")
    K_sat   <= 0.0     && error("K_sat must be positive")

    # Pre-compute the constant diffusivity for efficiency
    # D_w = K_s / [(θ_s − θ_r) δ]  (follows from K/C with exponential forms)
    D_w_const = K_sat / ((theta_s - theta_r) * delta)

    #______________________________________________________
    # Volumetric water content  θ(h)
    # θ(h) = θ_r + (θ_s − θ_r) exp(−δ|h|)
    #______________________________________________________
    function theta_h(h)
        h >= 0.0 && return theta_s
        return theta_r + (theta_s - theta_r) * exp(-delta * abs(h))
    end

    #______________________________________________________
    # Hydraulic conductivity  K(h)
    # K(h) = K_s exp(−δ|h|)   [Cavalcante & Zornberg 2017, Eq. 14]
    #______________________________________________________
    function K_h(h)
        h >= 0.0 && return K_sat
        return K_sat * exp(-delta * abs(h))
    end

    #______________________________________________________
    # Specific moisture capacity  C(h) = dθ/dh  — ANALYTICAL
    # C(h) = (θ_s − θ_r) δ exp(−δ|h|)
    #______________________________________________________
    function C_moist(h)
        h >= 0.0 && return 0.0
        return (theta_s - theta_r) * delta * exp(-delta * abs(h))
    end

    #______________________________________________________
    # Inverse retention  h(θ)  — exact analytical inversion
    #
    # From θ = θ_r + (θ_s − θ_r) exp(−δ|h|):
    # exp(−δ|h|) = (θ − θ_r)/(θ_s − θ_r)
    # −δ|h|      = ln[(θ − θ_r)/(θ_s − θ_r)]
    # |h|        = −(1/δ) ln[(θ − θ_r)/(θ_s − θ_r)]
    #
    # Since h < 0 in unsaturated zone:
    # h = −|h| = (1/δ) ln[(θ − θ_r)/(θ_s − θ_r)]
    #______________________________________________________
    function h_theta(theta)
        theta >= theta_s - 1.0e-12 && return 0.0
        theta <= theta_r + 1.0e-12 && return -1.0e8   # very dry limit

        theta_norm = clamp((theta - theta_r) / (theta_s - theta_r),
                           1.0e-15, 1.0 - 1.0e-15)
        return (1.0 / delta) * log(theta_norm)
    end

    #______________________________________________________
    # Hydraulic diffusivity  D_w(h) = K(h) / C(h)
    #
    # Because both K and C share the same exponential exp(−δ|h|),
    # the ratio is a constant:
    #   D_w = K_s exp(−δ|h|) / [(θ_s − θ_r) δ exp(−δ|h|)]
    #       = K_s / [(θ_s − θ_r) δ]
    #______________________________________________________
    function D_w(h)
        h >= 0.0 && return 0.0   # saturated: capillary diffusivity undefined
        return D_w_const
    end

    return (K_h=K_h, theta_h=theta_h, h_theta=h_theta, c_s=C_moist, D_w=D_w)
end


"""
    create_swrc_model(model_name::String, params::Dict{String, Float64})

Factory function to instantiate an SWRC model by name.

# Arguments
- `model_name::String` : `"Van_Genuchten"` or `"Cavalcante"`
- `params`             : Parameter dictionary (see individual model docs)

# Returns
- `NamedTuple` with closures: `K_h`, `theta_h`, `h_theta`, `C_moist`, `D_w`

# Throws
- `error` if `model_name` is not recognised
"""
function create_swrc_model(model_name::String, params::Dict{String, Float64})
    if model_name == "Van_Genuchten"
        return create_van_genuchten_model(params)
    elseif model_name == "Cavalcante"
        return create_cavalcante_model(params)
    else
        error("Unknown SWRC model: \"$model_name\". " *
              "Available options: \"Van_Genuchten\", \"Cavalcante\"")
    end
end