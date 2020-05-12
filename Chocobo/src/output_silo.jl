using Libdl

# File parameters
global const SILO_DB_CLOBBER = 0
global const SILO_DB_LOCAL = 0

# Data types
global const SILO_DB_NULL = -99
global const SILO_DB_DOUBLE = 20
global const SILO_DB_INTEGER = 16

# Mesh options
global const SILO_DBOPT_CYCLE = 263
global const SILO_DBOPT_DTIME = 280



# Mesh shapes
global const SILO_DB_ZONETYPE_QUAD = 24

# Variable centring
global const SILO_DB_NODECENT = 110
global const SILO_DB_ZONECENT = 111

# Keep track of number of files produced
silo_count = 0


function output_silo(mesh::Mesh, step::Int32, time::Float64, filestub::String)
global silo_count += 1
name = "$filestub" * fmt("04d",silo_count) * ".silo"

silo = Libdl.dlopen("/usr/local/prod/silo-4.10.2-so/lib/libsiloh5")
dbID = create_silo_file(silo, name)

create_silo_mesh(
    silo,
    "Mesh",
    mesh.nodelist, mesh.x, mesh.y,
    step, time,
    dbID
)

nodeids = NodeQuant{Int32}(mesh.nnodes)
cellids = NodeQuant{Int32}(mesh.ncells)
@inbounds for i in 1:mesh.nnodes
    nodeids[i] = i
end

@inbounds for i in 1:mesh.ncells
    cellids[i] = i
end

# nodeids = SharedArray([Int32(i) for i in 1:mesh.nnodes])
# cellids = SharedArray([Int32(i) for i in 1:mesh.ncells])

put_silo_material(silo, "Material", mesh.material, dbID, "Mesh")

put_silo_variable(silo, "Density", mesh.ρ, SILO_DB_ZONECENT, dbID, "Mesh")
put_silo_variable(silo, "Pressure", mesh.pressure, SILO_DB_ZONECENT, dbID, "Mesh")
put_silo_variable(silo, "Energy", mesh.energy, SILO_DB_ZONECENT, dbID, "Mesh")
put_silo_variable(silo, "Region", mesh.region, SILO_DB_ZONECENT, dbID, "Mesh")
put_silo_variable(silo, "NodeType", mesh.type, SILO_DB_NODECENT, dbID, "Mesh")
put_silo_variable(silo, "NodeID", nodeids, SILO_DB_NODECENT, dbID, "Mesh")
put_silo_variable(silo, "CellID", cellids, SILO_DB_ZONECENT, dbID, "Mesh")

end


function create_silo_file(silo, name::String)
# To get the Fortran routine to correctly set the database ID, we need to use
# a one-element array...
dbID::Array{Int32,1} = [0]
l = Ref{Int32}(length(name))
a = Ref{Int32}(SILO_DB_CLOBBER)
b = Ref{Int32}(SILO_DB_LOCAL)
l2 = Ref{Int32}(0)
c = Ref{Int32}(7)

DBCreate = Libdl.dlsym(silo, :dbcreate_)
status = ccall(
    DBCreate,
    Int32,
    (Ptr{UInt8}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
    name, l, a, b, "", l2, c, dbID
)
return dbID
end


function create_silo_mesh(
    silo,
    name::String, nodelist::CellQuant, x::NodeQuant, y::NodeQuant,
    step::Int32, time::Float64,
    dbID::Array{Int32}
)
DBPutZL2 = Libdl.dlsym(silo, :dbputzl2_)
DBPutUM = Libdl.dlsym(silo, :dbputum_)
DBMkOptList = Libdl.dlsym(silo, :dbmkoptlist_)
DBAddDOpt = Libdl.dlsym(silo, :dbadddopt_)
DBAddIOpt = Libdl.dlsym(silo, :dbaddiopt_)

ncells = size(nodelist, 2)
nnodes = size(x, 1)

# Create a 1D connectivity
connectivity::Array{Int32, 1} = Array{Int32}(undef, ncells*4)

k = 1
for i in 1:ncells
    for j in 1:4
        connectivity[k] = nodelist[j,i]
        k += 1
    end
end

# Set mesh options
optID::Array{Int32,1} = [0]

# dbmkoptlist(2, mesh_optlistID)
status = ccall(
    DBMkOptList,
    Int32,
    (Ref{Int32}, Ref{Int32}),
    2, optID
)

# SILO has a garbage interface here (treating Fortran arguments as pointer,
# which is wrong). Works "fine," except for int arguments, which it thinks
# are the value of the pointer. We can get around this by passing an array.
stepPtr::Array{Int32,1} = [step]
status = ccall(
    DBAddIOpt,
    Int32,
    (Ref{Int32}, Ref{Int32}, Ptr{Cint}),
    optID, SILO_DBOPT_CYCLE, stepPtr
)

status = ccall(
    DBAddDOpt,
    Int32,
    (Ref{Int32}, Ref{Int32}, Ref{Float64}),
    optID, SILO_DBOPT_DTIME, time
)

# Create a zonelist
shapetype::Array{Int32,1} = [SILO_DB_ZONETYPE_QUAD]
shapesize::Array{Int32,1} = [4]
shapecount::Array{Int32,1} = [ncells]
status2 = 0
status = ccall(
    DBPutZL2,
    Int32,
    (Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Int32},
     Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
     Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
    dbID, "ZL", 2, ncells, 2, connectivity, 4*ncells, 1, 0, 0,
    shapetype, shapesize, shapecount, 1, SILO_DB_NULL, status2
)

status = ccall(
    DBPutUM,
    Int32,
    (
        Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
        Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8}, Ref{Int32},
        Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8},
        Ref{Int32}, Ref{Int32}, Ref{Int32}
    ),
    dbID, name, length(name), 2, Array(x.d), Array(y.d), 0, "X", 1, "Y", 1, "", 0,
    SILO_DB_DOUBLE, nnodes, ncells, "ZL", 2, "", 0,
    optID, status2
)

end

function put_silo_material(
    silo,
    name::String, material::CellQuant,
    dbID::Array{Int32}, meshname::String,
)
DBPutMat = Libdl.dlsym(silo, :dbputmat_)
nmat::Int32 = maximum(material)
matnos = [Int32(i) for i in 1:nmat]
ndims::Int32 = 1
dims = [Int32(length(material)),]
status2::Int32 = 0

status = ccall(
    DBPutMat,
    Int32,
    (
        Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8}, Ref{Int32},
        Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
        Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
        Ref{Int32}, Ref{Int32}
    ),
    dbID, name, length(name), meshname, length(meshname),
    nmat, matnos, Array(material.d), dims, ndims, SILO_DB_NULL,
    SILO_DB_NULL, SILO_DB_NULL, SILO_DB_NULL, 0, SILO_DB_NULL, SILO_DB_NULL,
    status2
)

end


function put_silo_variable(
    silo,
    name::String, data::Quant, centring::Int64, dbID::Array{Int32},
    meshname::String
)
    # Work out data type
    if isa(data, Quant{Float64})
        put_silo_variable_f64(silo, name, data, centring, dbID, meshname)
    elseif isa(data, Quant{Int32})
        put_silo_variable_i32(silo, name, data, centring, dbID, meshname)
    else
        error("Unknown data type")
    end
end


function put_silo_variable_f64(
    silo,
    name::String, data::Quant, centring::Int64, dbID::Array{Int32},
    meshname::String
)
    DBPutUV1 = Libdl.dlsym(silo, :dbputuv1_)
    status2::Int32 = 0
    status = ccall(
        DBPutUV1,
        Int32,
        (
            Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8}, Ref{Int32},
            Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32}
        ),
        dbID, name, length(name), meshname, length(meshname),
        Array(data.d), length(data), 0, 0, SILO_DB_DOUBLE,
        centring, SILO_DB_NULL, status2
    )
end

function put_silo_variable_i32(
    silo,
    name::String, data::Quant, centring::Int64, dbID::Array{Int32},
    meshname::String
)
    DBPutUV1 = Libdl.dlsym(silo, :dbputuv1_)
    status2::Int32 = 0
    status = ccall(
        DBPutUV1,
        Int32,
        (
            Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ptr{UInt8}, Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32}
        ),
        dbID, name, length(name), meshname, length(meshname),
        Array(data.d), length(data), 0, 0, SILO_DB_INTEGER,
        centring, SILO_DB_NULL, status2
    )
end
