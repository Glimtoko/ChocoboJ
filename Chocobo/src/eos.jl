function perfect_gas!(
    energy::CellQuant, ρ::CellQuant, pressure::CellQuant,
    material::CellQuant, gamma
)
    ncells = size(energy, 1)

    @inbounds @sync @distributed for i in 1:ncells
        pressure[i] = (gamma[material[i]] - 1.0)*ρ[i]*energy[i]
    end
end
