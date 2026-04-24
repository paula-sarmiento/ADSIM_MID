include("src/version.jl")
using .ADSIMVersion: get_version

include("src/read_mesh.jl")
include("src/read_materials.jl")
include("src/read_calc_params.jl")
include("src/initialize_variables.jl")
include("src/initialize_flows.jl")
include("src/time_step.jl")
include("src/shape_functions.jl")
include("src/swrc_models.jl")
include("src/write_vtk.jl")
include("src/fully_explicit_solver.jl")
include("src/write_checkpoint.jl")
include("src/read_checkpoint.jl")
include("src/implicit_richards_solver.jl")

using .ShapeFunctions
println("✓ Richards solver compiled successfully")
