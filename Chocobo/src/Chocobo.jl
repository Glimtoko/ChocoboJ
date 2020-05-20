module Chocobo

using Distributed
using ArgParse

include("data.jl")
include("user_input.jl")
include("mesh_reader.jl")
include("set_initial_conditions.jl")
include("old_geometry.jl")
include("new_geometry.jl")
include("eos.jl")
include("new_hydro.jl")
include("output_text.jl")
include("output_silo.jl")
include("mainloop.jl")

function run_chocobo()
settings = ArgParseSettings()
@add_arg_table! settings begin
    "--control", "-c"
        help = "Name of control (input) file"
        default = "input.in"

    "--eos", "-e"
        help = "Name of EoS data file"
        default = "eos_data.jl"

    "--mesh", "-m"
        help = "Name of mesh file"
        default = "sod.key"

    "--problem", "-p"
        help = "Problem name - used as directory name where run files will be found"
        default = "."

    "--nprocs", "-n"
        help = "Number of processors to run on"
        arg_type = Int
        default = 1
end

parsed_args = parse_args(ARGS, settings)
problem = parsed_args["problem"]
controlf = "$problem/$(parsed_args["control"])"
meshf = "$problem/$(parsed_args["mesh"])"
eosf = "$problem/$(parsed_args["eos"])"
nprocs = parsed_args["nprocs"]

# Check these files exist
isfile(controlf) || error("Cannot open control file: $controlf")
isfile(meshf) || error("Cannot open mesh file: $meshf")
isfile(eosf) || error("Cannot open eos file: $eosf")

# gascodeloc="/home/nick/Programming/git/gascode2/Chocobo/src/"
# push!(LOAD_PATH, gascodeloc)

if nprocs > 1
    addprocs(nprocs - 1)  # Since we are *adding* processors.
end
# @everywhere using Pkg
# @everywhere Pkg.activate("/home/nick/Programming/git/gascode2/Chocobo/")

#@everywhere using Chocobo

# User input
input = Input()
read_input(input, controlf)
println(input)

# Get mesh size
println("Constructing mesh")
nnodes, ncells = get_mesh_size(meshf)

# Create mesh storage
mesh = Mesh{ncells, nnodes}()

read_mesh!(
    meshf,
    nnodes, ncells,
    mesh.x, mesh.y,
    mesh.nodelist, mesh.type,
    mesh.region, mesh.regioncelliterator,
    mesh.material, mesh.elelconn,
    mesh.nodnodconn,
    input.material_list
)

set_initial_conditions(eosf, mesh)

println("Initialisation")
# Set initial areas and volumes
calculate_area_volume!(
    mesh.x, mesh.y, mesh.nodelist,
    mesh.volume, mesh.area
)

# Set initial mass
calculate_mass!(mesh.volume, mesh.ρ, mesh.mass)

# Set initial pressure from ideal gas EoS
calculate_energy!(
    mesh.energy, mesh.pressure, mesh.material, mesh.ρ, input, mesh.gamma
)

println("Main loop")
mainloop(mesh, problem, input)
end

function julia_main()::Cint
    try
        println("ChocoboJ")
        run_chocobo()
    catch
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
    return 0 # if things finished successfully
end

end # module
