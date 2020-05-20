function get_mesh_size(filename::String)
    """
        get_mesh_size(filename)

        Gets the mesh size (numbers of nodes and cells) from the mesh
        file.
    """
    key = open(filename) do file
        readlines(file)
    end


    ncells = 0
    nnodes = 0

    innodelist = false
    incelllist = false
    for entry in key
        # First of all, assume that any "*" line ends the current section
        if occursin("*", entry)
            innodelist = false
            incelllist = false
        end

        # Parse nodelist lines
        if innodelist
            nnodes += 1
        end

        # Parse celllist lines
        if incelllist
            ncells += 1
        end

        # Detect start of nodelist
        if occursin("*NODE", entry)
            innodelist = true
        end

        # Detect start of region
        if occursin("PhysicalSurface", entry)
            incelllist = true
        end
    end

    return nnodes, ncells
end


function read_mesh!(
    filename::String,
    nnodes::Integer, ncells::Integer,
    x::NodeQuant, y::NodeQuant,
    nodelist::CellQuant, type::NodeQuant,
    region::CellQuant, regioncelliterator::Dict,
    material::CellQuant, elelconn::CellQuant,
    nodnodconn::NodeQuant,
    material_list::Dict
)
    boundaries = Dict(
    "XLOW" => -1,
    "XHIGH" => -1,
    "YLOW" => -2,
    "YHIGH" => -2,
)

    key = open(filename) do file
        readlines(file)
    end

    innodelist = false
    incelllist = false
    inbcs = false

    nl1 = Array{Int32, 1}()
    nl2 = Array{Int32, 1}()
    nl3 = Array{Int32, 1}()
    nl4 = Array{Int32, 1}()

    # type =
    reg = nothing
    bc = nothing
    ignore_next = false
    index_cell = 0
    index_node = 0
    for entry in key
        # First of all, assume that any "*" line ends the current section
        if occursin("*", entry)
            innodelist = false
            incelllist = false
        end

        # Parse nodelist lines
        if innodelist
            index_node += 1
            data = split(entry, ",")

            xp = parse(Float64, data[2])
            yp = parse(Float64, data[3])
            x[index_node] = xp
            y[index_node] = yp
        end

        # Parse celllist lines
        if incelllist
            index_cell += 1
            data = split(entry, ",")

            cell = parse(Int32, data[1])
            n1 = parse(Int32, data[3])
            n2 = parse(Int32, data[4])
            n3 = parse(Int32, data[5])
            n4 = parse(Int32, data[6])

            nodelist[1, index_cell] = n1
            nodelist[2, index_cell] = n2
            nodelist[3, index_cell] = n3
            nodelist[4, index_cell] = n4
            region[index_cell] = reg
        end

        # Boundary conditions
        if inbcs
            if ignore_next
                ignore_next = false
            else
                if occursin("#", entry)
                    bc = split(entry, "\$# ")[2]
                    ignore_next = true
                elseif !occursin("*", entry)
                    nodes = [parse(Int32, n) for n in split(entry, ",")]
                    for node in nodes
                        type[node] += boundaries[bc]
                    end
                end
            end
        end

        # Detect start of nodelist
        if occursin("*NODE", entry)
            println("Node list found")
            innodelist = true
        end

        # Detect start of region
        if occursin("PhysicalSurface", entry)
            reg = parse(Int32, split(entry, "PhysicalSurface")[2])÷11
            println("Found Region $reg")
            incelllist = true
        end

        # Detect node lists
        if occursin("*SET_NODE_LIST", entry)
            inbcs = true
        end
    end

    # Determine region/cell iterator forms
    nreg = maximum(region)
    println("Mesh contains $nreg regions")
    canuserange = true  # Assume we can use simple range interators
    for reg in 1:nreg
        println("Checking region $reg")
        first = findfirst(isequal(reg), region.d)
        last = findlast(isequal(reg), region.d)
        println("Lowest index: $first")
        println("Highest index: $last")
        @inbounds for i in first:last
            if region[i] != reg
                println("Error: Non-consecutive region cell iteration not coded")
                canuserange = false
            end
        end
        if !canuserange
            exit
        else
            regioncelliterator[reg] = first:last
        end
    end

    # Set the material list
    for i in 1:length(region)
        material[i] = material_list[region[i]]
    end

    # Calculate element-element connectivity
    # Sides
    sides = Dict(
        1 => ((1, 4), (2, 3)),
        2 => ((2, 1), (3, 4)),
        3 => ((4, 1), (3, 2)),
        4 => ((1, 2), (4, 3))
    )

    @inbounds for cell in 1:ncells
        @views n1 = nodelist.d[:, cell]
        @inbounds for cell2 in 1:ncells
            @views n2 = nodelist.d[:, cell2]

            for side in keys(sides)
                p11 = sides[side][1][1]
                p12 = sides[side][1][2]
                p21 = sides[side][2][1]
                p22 = sides[side][2][2]

                if x[n1[p11]] == x[n2[p12]] && x[n1[p21]] == x[n2[p22]]
                    if y[n1[p11]] == y[n2[p12]] && y[n1[p21]] == y[n2[p22]]
                        elelconn[side,cell] = cell2
                    end
                end
            end
        end
    end


    # Calculate node-node conneectivity
    for c0 in 1:ncells
        n0_tr = nodelist[3, c0]
        n1 = nodelist[2, c0]
        nodnodconn[1, n0_tr] = n1
        nodnodconn[3, n1] = n0_tr

        cr = elelconn[2, c0]
        if cr > 0
            n2 = nodelist[3, cr]
            nodnodconn[2, n0_tr] = n2
            nodnodconn[4, n2] = n0_tr
        else
            nodnodconn[2, n0_tr] = 0
        end

        cu = elelconn[3, c0]
        if cu > 0
            n3 = nodelist[3, cu]
            nodnodconn[3, n0_tr] = n3
            nodnodconn[1, n3] = n0_tr
        else
            nodnodconn[3, n0_tr] = 0
        end

        n4 = nodelist[4, c0]
        nodnodconn[4, n0_tr] = n4
        nodnodconn[2, n4] = n0_tr

        n0_bl = nodelist[1, c0]

        cd = elelconn[1, c0]
        if cd > 0
            n1 = nodelist[1, cd]
            nodnodconn[1, n0_bl] = n1
            nodnodconn[3, n1] = n0_bl
        else
            nodnodconn[1, n0_bl] = 0
        end

        n2 = nodelist[2, c0]
        nodnodconn[2, n0_bl] = n2
        nodnodconn[4, n2] = n0_bl

        n3 = nodelist[4, c0]
        nodnodconn[3, n0_bl] = n3
        nodnodconn[1, n3] = n0_bl

        cl = elelconn[4, c0]
        if cl > 0
            n4 = nodelist[1, cl]
            nodnodconn[4, n0_bl] = n4
            nodnodconn[2, n4] = n0_bl
        else
            nodnodconn[4, n0_bl] = 0
        end


    end

    # Check node-node connectivity
    for node in 1:nnodes
        nzero = sum(x->x==0, nodnodconn.d[:, node])
        if nzero == 1
            if type[node] != -1 && type[node] != -2
                error("Invalid nodnodconn for node $node")
            end
        elseif nzero == 2
            if type[node] != -3
                error("Invalid nodnodconn for node $node")
            end
        elseif nzero >= 3
            error("Invalid nodnodconn for node $node")
        end
    end
end
