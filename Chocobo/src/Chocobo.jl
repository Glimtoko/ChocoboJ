module Chocobo

using Distributed

include("data.jl")
include("user_input.jl")
include("mesh_reader.jl")
include("set_initial_conditions.jl")
include("old_hydro.jl")
include("new_hydro.jl")
include("mainloop.jl")


end # module
