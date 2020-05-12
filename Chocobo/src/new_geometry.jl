function calculate_area_volume!(
    x::NodeQuant, y::NodeQuant, nodelist::CellQuant,
    volume::CellQuant, area::CellQuant
)
    ncells = size(volume, 1)

    @inbounds @sync @distributed for i in 1:ncells
        n = nodelist.d[:, i]

        a1 = (-x[n[1]] + x[n[2]] + x[n[3]] - x[n[4]])/4.0
        a3 = (-x[n[1]] - x[n[2]] + x[n[3]] + x[n[4]])/4.0

        b1 = (-y[n[1]] + y[n[2]] + y[n[3]] - y[n[4]])/4.0
        b3 = (-y[n[1]] - y[n[2]] + y[n[3]] + y[n[4]])/4.0

        volume[i] = 4.0(a1*b3 - a3*b1)
        area[i] = volume[i]
    end
end


function calculate_mass!(volume::CellQuant, ρ::CellQuant, mass::CellQuant)
    ncells = size(volume, 1)  # All these arrays are cell-based
    @inbounds @sync @distributed for i in 1:ncells
        mass[i] = volume[i]*ρ[i]
    end
end

function calculate_energy!(
    energy::CellQuant, pressure::CellQuant,
    material::CellQuant, ρ::CellQuant,
    input::Input, gamma
)
    ncells = size(energy, 1)

    @inbounds @sync @distributed for i in 1:ncells
        energy[i] = pressure[i]/(gamma[material[i]] - 1.0)ρ[i]
    end
end
