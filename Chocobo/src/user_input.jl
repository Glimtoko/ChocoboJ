using Parameters

# Custom types
@enum ANTI_HG_TYPES begin
    ANTI_HG_NONE = 1
    ANTI_HG_DYNA_1 = 2
    ANTI_HG_DYNA_2 = 3
    ANTI_HG_BF = 4
end

@with_kw mutable struct Input
    boundary_values = Dict(
        "XLOW" => -1,
        "XHIGH" => -1,
        "YLOW" => -2,
        "YHIGH" => -2,
    )

    material_list = Dict(
        1 => 1,
        2 => 2
    )

    # Artificial Viscosity
    qform = 2
    CL = 0.5
    CQ = 0.75

    # Robustness
    anti_hg = [ANTI_HG_BF, ANTI_HG_BF]
    κ = 0.7

    # Mesh Extents
    x0 = 0.0
    x1 = 1.0
    y0 = 0.0
    y1 = 1.0

    # Timestep, etc.
    t0 = 0.0
    tend = 0.205
    dtinit = 0.0001
    dtmax = 0.01
    growth = 1.02
    dtmin_hg = 1.0e-6

    # Output frequency
    dtsilo = 0.01

    # Profiling/Debug settings
    debug_step_count = 0   # Set to zero to disable debugging
    profile = false

    # Cutoffs
    ρcutoff = 1.0e-6
end


function read_input(input::Input, filename::String)
    data = readlines(filename)
    fields = map((x) -> string(x), fieldnames(Input))
    for line in data
        line = split(line, "#")[1]   # Remove comments
        if occursin("=", line)
            # Setting something
            if occursin("{", line)
                # Dictionary
                key, value = split(line, "=")
                key = strip(key)

                if key in fields
                    ks = Symbol(key)
                    d = getproperty(input, ks)
                    value = replace(value, "{" => "")
                    value = replace(value, "}" => "")
                    items = split(value, ",")
                    for item in items
                        dkey, dvalue = split(item, ":")
                        dkey = try
                            parse(Int32, dkey)
                        catch
                            uppercase(strip(dkey))
                        end
                        if haskey(d, dkey)
                            t = typeof(d[dkey])
                            dvalue = parse(t, dvalue)
                            d[dkey] = dvalue
                        else
                            error("Invalid key $dkey provided to $key")
                        end
                    end
                else
                    error("Unexpected input item: $key")
                end
            elseif occursin("[", line)
                # Array - will contain enum values
            else
                # Regular input
                key, value = split(line, "=")
                key = strip(key)
                value = strip(value)
                if key in fields
                    # Get type of data currently held in input
                    ks = Symbol(key)
                    t = typeof(getproperty(input, ks))

                    # Convert data appropriately
                    value = parse(t, value)

                    # Set new value
                    setproperty!(input, ks, value)
                else
                    error("Unexpected input item: $key")
                end

            end
        end
    end

end
