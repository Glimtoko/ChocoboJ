function set_initial_conditions(eos_file::String, mesh::Mesh)
    # Read the EoS file
    include(eos_file)

    nmat = length(eos_data)
    println("\nFound data for $nmat material")

    global gamma = zeros(nmat)
    for mat in 1:nmat
        println("Material $mat")
        println("  Initial ρ: $(eos_data[mat].ρ0)")
        println("  Initial P: $(eos_data[mat].P0)")
        println("  Perfect Gas γ: $(eos_data[mat].γ)")
        gamma[mat] = eos_data[mat].γ

        @inbounds for cell in 1:mesh.ncells
            if mesh.material[cell] == mat
                mesh.ρ[cell] = eos_data[mat].ρ0
                mesh.pressure[cell] = eos_data[mat].P0
            end
        end
    end
end 
