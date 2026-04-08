#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and
# Modification
# v0.x.x
# Author: Paula Sarmiento
#______________________________________________________

#______________________________________________________
# SWRC (Soil Water Retention Curve) Model Implementations
# Provides Van Genuchten, Cavalcante, LinearSoil, and ConstantSoil models (Latter two for verification)
#
# All models are subtypes of abstract SWRCModel and implement:
#   Se(model, h)        -> effective saturation         [-]
#   theta(model, h)     -> volumetric water content     [-]
#   h_inv(model, theta) -> pressure head from theta     [L]
#   C_moist(model, h)   -> specific moisture capacity   [1/L]
#   K_h(model, h)       -> hydraulic conductivity       [L/T]
#   D_w(model, h)       -> hydraulic diffusivity K/C    [L²/T]
#
# Convention: h < 0 in the unsaturated zone, h >= 0 at saturation.
#______________________________________________________

using Printf

# Abstract type for dispatch
abstract type SWRCModel end


# ══════════════════════════════════════════════════════════════════════════════
# Van Genuchten–Mualem Model
# ══════════════════════════════════════════════════════════════════════════════

"""
    struct VanGenuchten <: SWRCModel

Van Genuchten (1980) SWRC with Mualem (1976) conductivity model.

# Fields
- `theta_s::Float64` : Saturated volumetric water content       [-]
- `theta_r::Float64` : Residual volumetric water content        [-]
- `alpha::Float64`   : Inverse air-entry pressure               [1/L]
- `n_vg::Float64`    : Pore-size distribution index (> 1)       [-]
- `m::Float64`       : Complementary parameter m = 1 - 1/n      [-]
- `K_s::Float64`     : Saturated hydraulic conductivity         [L/T]
- `K_s_x::Float64`   : Directional x-component of K_s           [L/T]
- `K_s_y::Float64`   : Directional y-component of K_s           [L/T]
- `L::Float64`       : Mualem pore-connectivity exponent        [-]

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
struct VanGenuchten <: SWRCModel
    theta_s::Float64
    theta_r::Float64
    alpha::Float64
    n_vg::Float64
    m::Float64
    K_s::Float64
    K_s_x::Float64
    K_s_y::Float64
    L::Float64
end

function VanGenuchten(;
    theta_s::Float64,
    theta_r::Float64,
    alpha::Float64,
    n_vg::Float64,
    K_s::Float64,
    K_s_x::Float64 = K_s,
    K_s_y::Float64 = K_s,
    L::Float64 = 0.5
)
    # Validate parameter ranges
    theta_r >= theta_s  && error("theta_r must be < theta_s")
    alpha   <= 0.0      && error("alpha must be positive")
    n_vg    <= 1.0      && error("n_vg must be > 1.0")
    K_s     <= 0.0      && error("K_s must be positive")
    K_s_x   <= 0.0      && error("K_s_x must be positive (or omit for isotropic)")
    K_s_y   <= 0.0      && error("K_s_y must be positive (or omit for isotropic)")
    L       <  0.0      && error("L must be non-negative")
    
    m = 1.0 - 1.0 / n_vg
    
    return VanGenuchten(theta_s, theta_r, alpha, n_vg, m, K_s, K_s_x, K_s_y, L)
end

"""
    Se(model::VanGenuchten, h::Float64) -> Float64

Compute effective saturation for Van Genuchten model.
"""
function Se(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return 1.0
    return (1.0 + (model.alpha * abs(h))^model.n_vg)^(-model.m)
end

"""
    theta(model::VanGenuchten, h::Float64) -> Float64

Compute volumetric water content for Van Genuchten model.
"""
function theta(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return model.theta_s
    return model.theta_r + (model.theta_s - model.theta_r) * Se(model, h)
end

"""
    C_moist(model::VanGenuchten, h::Float64) -> Float64

Compute specific moisture capacity (analytical dθ/dh) for Van Genuchten model.
"""
function C_moist(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return 0.0
    h_abs = abs(h)
    # Guard against h = 0 approached from below
    h_abs < 1.0e-14 && return 0.0
    numerator   = model.m * model.n_vg * model.alpha^model.n_vg * h_abs^(model.n_vg - 1.0)
    denominator = (1.0 + (model.alpha * h_abs)^model.n_vg)^(model.m + 1.0)
    return (model.theta_s - model.theta_r) * numerator / denominator
end

"""
    K_h(model::VanGenuchten, h::Float64) -> Float64

Compute hydraulic conductivity (isotropic) for Van Genuchten model using Mualem (1976).
"""
function K_h(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return model.K_s
    Se_val  = Se(model, h)
    inner   = clamp(1.0 - Se_val^(1.0 / model.m), 0.0, 1.0)
    return model.K_s * Se_val^model.L * (1.0 - inner^model.m)^2
end

"""
    K_h_x(model::VanGenuchten, h::Float64) -> Float64

Compute directional hydraulic conductivity (x-component) for Van Genuchten model.
"""
function K_h_x(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return model.K_s_x
    Se_val  = Se(model, h)
    inner   = clamp(1.0 - Se_val^(1.0 / model.m), 0.0, 1.0)
    return model.K_s_x * Se_val^model.L * (1.0 - inner^model.m)^2
end

"""
    K_h_y(model::VanGenuchten, h::Float64) -> Float64

Compute directional hydraulic conductivity (y-component) for Van Genuchten model.
"""
function K_h_y(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return model.K_s_y
    Se_val  = Se(model, h)
    inner   = clamp(1.0 - Se_val^(1.0 / model.m), 0.0, 1.0)
    return model.K_s_y * Se_val^model.L * (1.0 - inner^model.m)^2
end

"""
    h_inv(model::VanGenuchten, theta::Float64) -> Float64

Compute pressure head from volumetric water content (inverse retention curve).
Uses exact analytical inversion for Van Genuchten model.
"""
function h_inv(model::VanGenuchten, theta::Float64)::Float64
    theta >= model.theta_s - 1.0e-12 && return 0.0
    theta <= model.theta_r + 1.0e-12 && return -1.0e8   # very dry limit

    Se_val = (theta - model.theta_r) / (model.theta_s - model.theta_r)
    base   = Se_val^(-1.0 / model.m) - 1.0
    base < 1.0e-15 && return 0.0
    return -(1.0 / model.alpha) * base^(1.0 / model.n_vg)
end

"""
    D_w(model::VanGenuchten, h::Float64) -> Float64

Compute hydraulic diffusivity D_w(h) = K(h) / C(h) for Van Genuchten model.
"""
function D_w(model::VanGenuchten, h::Float64)::Float64
    h >= 0.0 && return 0.0   # saturated: capillary diffusivity undefined
    C = C_moist(model, h)
    C < 1.0e-15 && return 0.0
    return K_h(model, h) / C
end

# ══════════════════════════════════════════════════════════════════════════════
# Cavalcante-Zornberg Model
# ══════════════════════════════════════════════════════════════════════════════
"""
    struct CavalcanteZornberg <: SWRCModel

Cavalcante & Zornberg (2017) exponential SWRC model.

# Fields
- `theta_s::Float64` : Saturated volumetric water content  [-]
- `theta_r::Float64` : Residual volumetric water content   [-]
- `delta::Float64`   : Exponential desaturation parameter  [1/L]
- `K_s::Float64`     : Saturated hydraulic conductivity    [L/T]
- `K_s_x::Float64`   : Directional x-component of K_s      [L/T]
- `K_s_y::Float64`   : Directional y-component of K_s      [L/T]

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

Inverse retention (analytical, exact):
```
h(θ) = −(1/δ) ln[(θ − θ_r)/(θ_s − θ_r)]
```

Hydraulic diffusivity:
```
D_w(h) = K(h) / C(h) = K_s / [(θ_s − θ_r) δ]   (constant in h!)
```
"""
struct CavalcanteZornberg <: SWRCModel
    theta_s::Float64
    theta_r::Float64
    delta::Float64
    K_s::Float64
    K_s_x::Float64
    K_s_y::Float64
end

function CavalcanteZornberg(;
    theta_s::Float64,
    theta_r::Float64,
    delta::Float64,
    K_s::Float64,
    K_s_x::Float64 = K_s,
    K_s_y::Float64 = K_s
)
    # Validate parameter ranges
    theta_r >= theta_s && error("theta_r must be < theta_s")
    delta   <= 0.0     && error("delta must be positive")
    K_s     <= 0.0     && error("K_s must be positive")
    K_s_x   <= 0.0     && error("K_s_x must be positive (or omit for isotropic)")
    K_s_y   <= 0.0     && error("K_s_y must be positive (or omit for isotropic)")
    
    return CavalcanteZornberg(theta_s, theta_r, delta, K_s, K_s_x, K_s_y)
end

"""
    Se(model::CavalcanteZornberg, h::Float64) -> Float64

Compute effective saturation for Cavalcante model.
"""
function Se(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return 1.0
    return exp(-model.delta * abs(h))
end

"""
    theta(model::CavalcanteZornberg, h::Float64) -> Float64

Compute volumetric water content for Cavalcante model.
"""
function theta(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return model.theta_s
    return model.theta_r + (model.theta_s - model.theta_r) * Se(model, h)
end

"""
    C_moist(model::CavalcanteZornberg, h::Float64) -> Float64

Compute specific moisture capacity (analytical dθ/dh) for Cavalcante model.
"""
function C_moist(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return 0.0
    return (model.theta_s - model.theta_r) * model.delta * exp(-model.delta * abs(h))
end

"""
    K_h(model::CavalcanteZornberg, h::Float64) -> Float64

Compute hydraulic conductivity (isotropic) for Cavalcante model.
"""
function K_h(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return model.K_s
    return model.K_s * exp(-model.delta * abs(h))
end

"""
    K_h_x(model::CavalcanteZornberg, h::Float64) -> Float64

Compute directional hydraulic conductivity (x-component) for Cavalcante model.
"""
function K_h_x(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return model.K_s_x
    return model.K_s_x * exp(-model.delta * abs(h))
end

"""
    K_h_y(model::CavalcanteZornberg, h::Float64) -> Float64

Compute directional hydraulic conductivity (y-component) for Cavalcante model.
"""
function K_h_y(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return model.K_s_y
    return model.K_s_y * exp(-model.delta * abs(h))
end

"""
    h_inv(model::CavalcanteZornberg, theta::Float64) -> Float64

Compute pressure head from volumetric water content (inverse retention curve).
Uses exact analytical inversion for Cavalcante model.
"""
function h_inv(model::CavalcanteZornberg, theta::Float64)::Float64
    theta >= model.theta_s - 1.0e-12 && return 0.0
    theta <= model.theta_r + 1.0e-12 && return -1.0e8   # very dry limit

    theta_norm = clamp((theta - model.theta_r) / (model.theta_s - model.theta_r),
                       1.0e-15, 1.0 - 1.0e-15)
    return (1.0 / model.delta) * log(theta_norm)
end

"""
    D_w(model::CavalcanteZornberg, h::Float64) -> Float64

Compute hydraulic diffusivity D_w(h) = K(h) / C(h) for Cavalcante model.
Because both K and C share the same exponential exp(−δ|h|), the ratio is constant.
"""
function D_w(model::CavalcanteZornberg, h::Float64)::Float64
    h >= 0.0 && return 0.0   # saturated: capillary diffusivity undefined
    D_w_const = model.K_s / ((model.theta_s - model.theta_r) * model.delta)
    return D_w_const
end

# ══════════════════════════════════════════════════════════════════════════════
# Linear Soil (verification baseline)
# ══════════════════════════════════════════════════════════════════════════════

"""
    struct LinearSoil <: SWRCModel

Piecewise-linear SWRC model for verification and testing.

# Fields
- `theta_r::Float64` : Residual water content [-]
- `theta_s::Float64` : Saturated water content [-]
- `h_min::Float64`   : Lower pressure head bound [L] (typically negative)
- `K_s::Float64`     : Saturated hydraulic conductivity [L/T]

# Notes
- Effective saturation: ``S_e(h) = \\text{clamp}\\left(\\frac{h - h_{\\min}}{-h_{\\min}}, 0, 1\\right)``
- Water content: ``\\theta(h) = \\theta_r + (\\theta_s - \\theta_r) S_e(h)``
- Capacity (constant): ``C_{\\text{moist}} = \\frac{\\theta_s - \\theta_r}{|h_{\\min}|}``
- Hydraulic conductivity: ``K_h(h) = K_{\\text{sat}} \\cdot S_e(h)``
"""
struct LinearSoil <: SWRCModel
    theta_r::Float64
    theta_s::Float64
    h_min::Float64
    K_s::Float64
end

function LinearSoil(;
    theta_r::Float64,
    theta_s::Float64,
    h_min::Float64,
    K_s::Float64
)
    # Validate parameter ranges
    theta_r >= theta_s && error("theta_r must be < theta_s")
    h_min >= 0.0       && error("h_min must be negative")
    K_s <= 0.0         && error("K_s must be positive")
    
    return LinearSoil(theta_r, theta_s, h_min, K_s)
end

"""
    Se(model::LinearSoil, h::Float64) -> Float64

Compute normalized effective saturation for LinearSoil model.
"""
function Se(model::LinearSoil, h::Float64)::Float64
    return clamp((h - model.h_min) / (-model.h_min), 0.0, 1.0)
end

"""
    theta(model::LinearSoil, h::Float64) -> Float64

Compute volumetric water content for LinearSoil model.
"""
function theta(model::LinearSoil, h::Float64)::Float64
    return model.theta_r + (model.theta_s - model.theta_r) * Se(model, h)
end

"""
    C_moist(model::LinearSoil, h::Float64) -> Float64

Compute specific moisture capacity (constant) for LinearSoil model.
"""
function C_moist(model::LinearSoil, h::Float64)::Float64
    return (model.theta_s - model.theta_r) / abs(model.h_min)
end

"""
    K_h(model::LinearSoil, h::Float64) -> Float64

Compute hydraulic conductivity for LinearSoil model.
"""
function K_h(model::LinearSoil, h::Float64)::Float64
    Se_val = clamp(Se(model, h), 1.0e-6, 1.0)
    return model.K_s * Se_val
end

"""
    h_inv(model::LinearSoil, theta::Float64) -> Float64

Compute pressure head from volumetric water content for LinearSoil model.
"""
function h_inv(model::LinearSoil, theta::Float64)::Float64
    theta_clamped = clamp(theta, model.theta_r, model.theta_s)
    Se_val = (theta_clamped - model.theta_r) / (model.theta_s - model.theta_r)
    return model.h_min * (1.0 - Se_val)
end

"""
    D_w(model::LinearSoil, h::Float64) -> Float64

Compute hydraulic diffusivity for LinearSoil model.
"""
function D_w(model::LinearSoil, h::Float64)::Float64
    return K_h(model, h) / C_moist(model, h)
end

# ══════════════════════════════════════════════════════════════════════════════
# Constant Soil (pure diffusion verification)
# ══════════════════════════════════════════════════════════════════════════════

"""
    struct ConstantSoil <: SWRCModel

Constant-coefficient model for verification and testing.

Richards equation reduces to the heat equation: ``\\frac{∂\\theta}{∂t} = D \\nabla^2 h``,
where ``D = K / C`` is a true constant diffusivity.

# Fields
- `theta_r::Float64` : Residual water content [-]
- `theta_s::Float64` : Saturated water content [-]
- `h_min::Float64`   : Lower pressure head bound [L] (typically negative)
- `K_val::Float64`   : Constant hydraulic conductivity [L/T]

# Notes
- Effective saturation: ``S_e(h) = \\text{clamp}\\left(\\frac{h - h_{\\min}}{-h_{\\min}}, 0, 1\\right)``
- Water content: ``\\theta(h) = \\theta_r + (\\theta_s - \\theta_r) S_e(h)``
- Capacity (constant): ``C_{\\text{moist}} = \\frac{\\theta_s - \\theta_r}{|h_{\\min}|}``
- Hydraulic conductivity (constant): ``K_h(h) = K``
- Diffusivity (constant): ``D_w = K / C``
"""
struct ConstantSoil <: SWRCModel
    theta_r::Float64
    theta_s::Float64
    h_min::Float64
    K_val::Float64
end

function ConstantSoil(;
    theta_r::Float64,
    theta_s::Float64,
    h_min::Float64,
    K_val::Float64
)
    # Validate parameter ranges
    theta_r >= theta_s && error("theta_r must be < theta_s")
    h_min >= 0.0       && error("h_min must be negative")
    K_val <= 0.0       && error("K_val must be positive")
    
    return ConstantSoil(theta_r, theta_s, h_min, K_val)
end

"""
    Se(model::ConstantSoil, h::Float64) -> Float64

Compute normalized effective saturation for ConstantSoil model.
"""
function Se(model::ConstantSoil, h::Float64)::Float64
    return clamp((h - model.h_min) / (-model.h_min), 0.0, 1.0)
end

"""
    theta(model::ConstantSoil, h::Float64) -> Float64

Compute volumetric water content for ConstantSoil model.
"""
function theta(model::ConstantSoil, h::Float64)::Float64
    return model.theta_r + (model.theta_s - model.theta_r) * Se(model, h)
end

"""
    C_moist(model::ConstantSoil, h::Float64) -> Float64

Compute specific moisture capacity (constant) for ConstantSoil model.
"""
function C_moist(model::ConstantSoil, h::Float64)::Float64
    return (model.theta_s - model.theta_r) / abs(model.h_min)
end

"""
    K_h(model::ConstantSoil, h::Float64) -> Float64

Compute hydraulic conductivity (constant) for ConstantSoil model.
"""
function K_h(model::ConstantSoil, h::Float64)::Float64
    return model.K_val
end

"""
    h_inv(model::ConstantSoil, theta::Float64) -> Float64

Compute pressure head from volumetric water content for ConstantSoil model.
"""
function h_inv(model::ConstantSoil, theta::Float64)::Float64
    theta_clamped = clamp(theta, model.theta_r, model.theta_s)
    Se_val = (theta_clamped - model.theta_r) / (model.theta_s - model.theta_r)
    return model.h_min * (1.0 - Se_val)
end

"""
    D_w(model::ConstantSoil, h::Float64) -> Float64

Compute hydraulic diffusivity (constant) for ConstantSoil model.
"""
function D_w(model::ConstantSoil, h::Float64)::Float64
    C = C_moist(model, h)
    return model.K_val / C
end


# ══════════════════════════════════════════════════════════════════════════════
# Factory Function for Creating SWRC Models
# ══════════════════════════════════════════════════════════════════════════════

"""
    create_swrc_model(name::String, params::Dict{String, Float64}) -> SWRCModel

Factory function to instantiate an SWRC model by name and parameter dictionary.

# Arguments
- `name::String` : Model identifier. Valid options:
  - `"Van_Genuchten"` : Van Genuchten–Mualem model
  - `"Cavalcante"` : Cavalcante & Zornberg exponential model
  - `"LinearSoil"` : Piecewise-linear model (verification)
  - `"ConstantSoil"` : Constant-coefficient model (verification)
  
- `params::Dict{String, Float64}` : Parameter dictionary (see individual model constructors)

# Returns
- `SWRCModel` : A concrete model instance (VanGenuchten, CavalcanteZornberg, LinearSoil, or ConstantSoil)

# Throws
- `error` if `name` is not recognized or required parameters are missing

# Example
```julia
params = Dict(
    "theta_s" => 0.45,
    "theta_r" => 0.05,
    "alpha" => 0.02,
    "n_param" => 2.0,
    "K_sat" => 0.1
)
model = create_swrc_model("Van_Genuchten", params)
theta_val = theta(model, -100.0)  # Get water content at h = -100 cm
```
"""
function create_swrc_model(name::String, params::Dict{String, Float64})::SWRCModel
    if name == "Van_Genuchten"
        # Van Genuchten model requires: theta_s, theta_r, alpha, n_param (or n_vg), K_sat
        n_key = haskey(params, "n_param") ? "n_param" : "n_vg"
        return VanGenuchten(
            theta_s = params["theta_s"],
            theta_r = params["theta_r"],
            alpha = params["alpha"],
            n_vg = params[n_key],
            K_s = params["K_sat"],
            K_s_x = get(params, "K_sat_x", params["K_sat"]),
            K_s_y = get(params, "K_sat_y", params["K_sat"]),
            L = get(params, "L_mualem", 0.5)
        )
    
    elseif name == "Cavalcante"
        # Cavalcante model requires: theta_s, theta_r, delta, K_sat
        return CavalcanteZornberg(
            theta_s = params["theta_s"],
            theta_r = params["theta_r"],
            delta = params["delta"],
            K_s = params["K_sat"],
            K_s_x = get(params, "K_sat_x", params["K_sat"]),
            K_s_y = get(params, "K_sat_y", params["K_sat"])
        )
    
    elseif name == "LinearSoil"
        # LinearSoil model requires: theta_s, theta_r, h_min, K_sat
        return LinearSoil(
            theta_r = params["theta_r"],
            theta_s = params["theta_s"],
            h_min = params["h_min"],
            K_s = params["K_sat"]
        )
    
    elseif name == "ConstantSoil"
        # ConstantSoil model requires: theta_s, theta_r, h_min, K_val
        return ConstantSoil(
            theta_r = params["theta_r"],
            theta_s = params["theta_s"],
            h_min = params["h_min"],
            K_val = params["K_val"]
        )
    
    else
        error("Unknown SWRC model: \"$name\". " *
              "Available options: \"Van_Genuchten\", \"Cavalcante\", \"LinearSoil\", \"ConstantSoil\"")
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Public API Exports
# ══════════════════════════════════════════════════════════════════════════════

export SWRCModel
export VanGenuchten, CavalcanteZornberg, LinearSoil, ConstantSoil
export Se, theta, C_moist, K_h, K_h_x, K_h_y, h_inv, D_w
export create_swrc_model