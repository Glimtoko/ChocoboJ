# N.b. Arrays are always mutable, and the data in them is the only part of these
# structures which should ever change, so the structs themselves do not need to
# be mutable.

using Parameters
using SharedArrays

abstract type Quant{T<:Number} end

struct NodeQuant{T} <: Quant{T}
    dims::Int32
    n::Int32
    n1::Int32
    n2::Int32
    d::SharedArray{T}
    NodeQuant{T}(n) where {T<:Number} = new(1, n, n, n, zeros(n))
    NodeQuant{T}(n1, n2) where {T<:Number} = new(2, n1*n2, n1, n2, zeros((n1, n2)))
end

struct CellQuant{T} <: Quant{T}
    dims::Int32
    n::Int32
    n1::Int32
    n2::Int32
    d::SharedArray{T}
    CellQuant{T}(n) where {T<:Number} = new(1, n, n, 1, zeros(n))
    CellQuant{T}(n1, n2) where {T<:Number} = new(2, n1*n2, n1, n2, zeros((n1, n2)))
end

function Base.setindex!(Q::Quant, v::Number, i::Number)
     1 <= i <= Q.n || throw(BoundsError(Q, i))
     Q.d[i] = v
end

function Base.setindex!(Q::Quant, v::Number, i1::Number, i2::Number)
     1 <= i1 <= Q.n1 || throw(BoundsError(Q, i1))
     1 <= i2 <= Q.n2 || throw(BoundsError(Q, i2))
     Q.d[i1,i2] = v
end

Base.getindex(Q::Quant, i::Number) = Q.d[i]
Base.getindex(Q::Quant, i1::Number, i2::Number) = Q.d[i1,i2]

Base.iterate(Q::Quant, state=1) = state > Q.n ? nothing : (Q.d[state], state+1)
Base.keys(Q::Quant) = Base.LinearIndices(Q.d)
Base.length(Q::Quant) = Q.n
Base.size(Q::Quant) = (Q.n1, Q.n2)
Base.size(Q::Quant, I::Integer) = ([Q.n1, Q.n2][I])

function copyquant!(q1::Quant, q2::Quant)
    @inbounds @sync @distributed for i in 1:q1.n
        q2[i] = q1[i]
    end
end

@with_kw mutable struct Mesh{NCells, NNodes}
    # Mesh dimensions
    ncells = NCells
    nnodes = NNodes

    # Nodal coordinates, velocities and node type
    x = NodeQuant{Float64}(NNodes)
    y = NodeQuant{Float64}(NNodes)
    u = NodeQuant{Float64}(NNodes)
    v = NodeQuant{Float64}(NNodes)
    type = NodeQuant{Int32}(NNodes)

    # Cell quants - Pressure, density, energy, etc.
    ρ = CellQuant{Float64}(NCells)
    pressure = CellQuant{Float64}(NCells)
    energy = CellQuant{Float64}(NCells)
    volume = CellQuant{Float64}(NCells)
    soundspeed = CellQuant{Float64}(NCells)
    area = CellQuant{Float64}(NCells)
    mass = CellQuant{Float64}(NCells)
    q = CellQuant{Float64}(NCells)

    # Half timestep quantities
    x½ = NodeQuant{Float64}(NNodes)
    y½ = NodeQuant{Float64}(NNodes)
    ρ½ = CellQuant{Float64}(NCells)
    pressure½ = CellQuant{Float64}(NCells)
    energy½ = CellQuant{Float64}(NCells)
    volume½ = CellQuant{Float64}(NCells)

    # Arrays for FEM - Shape Functions etc.
    ∫N = CellQuant{Float64}(4, NCells)
    ∫∂N∂x = CellQuant{Float64}(4, NCells)
    ∫∂N∂y = CellQuant{Float64}(4, NCells)
    ∂N∂x = CellQuant{Float64}(4, NCells)
    ∂N∂y = CellQuant{Float64}(4, NCells)
    element_weightc = CellQuant{Float64}(4, NCells)

    ▽●v = CellQuant{Float64}(NCells)
    ∫▽●V = CellQuant{Float64}(NCells)

    # Material and region lists
    region = CellQuant{Int32}(NCells)
    material = CellQuant{Int32}(NCells)

    # Connectivity
    nodelist = CellQuant{Int32}(4, NCells)
    regioncelliterator = Dict()
    elelconn = CellQuant{Int32}(4, NCells)
    nodnodconn = NodeQuant{Int32}(4, NNodes)

    # Stored (old) data
    volume_store = CellQuant{Float64}(NCells)
    u_store = NodeQuant{Float64}(NNodes)
    v_store = NodeQuant{Float64}(NNodes)

    # Averaged velocities
    ubar = NodeQuant{Float64}(NNodes)
    vbar = NodeQuant{Float64}(NNodes)

    # Work arrays
    mass_scatter_to_nodes = NodeQuant{Float64}(NNodes)
    force_scatter_to_nodes_x = NodeQuant{Float64}(NNodes)
    force_scatter_to_nodes_y = NodeQuant{Float64}(NNodes)

    # Materials
    gamma = []
end
