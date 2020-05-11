function get_finite_elements!(
    x::NodeQuant, y::NodeQuant, nodelist::CellQuant,
    ∫N::CellQuant, ∫∂N∂x::CellQuant, ∫∂N∂y::CellQuant,
    ∂N∂x::CellQuant, ∂N∂y::CellQuant, element_weightc::CellQuant
)
    ncells = size(∫N, 2)

    @inbounds @distributed for i in 1:ncells
        @views n = nodelist.d[:, i]

        # Calculate coefficients for shape function integeral
        a1 = (-x[n[1]] + x[n[2]] + x[n[3]] - x[n[4]])/4.0
        a2 = (x[n[1]] - x[n[2]] + x[n[3]] - x[n[4]])/4.0
        a3 = (-x[n[1]] - x[n[2]] + x[n[3]] + x[n[4]])/4.0

        b1 = (-y[n[1]] + y[n[2]] + y[n[3]] - y[n[4]])/4.0
        b2 = (y[n[1]] - y[n[2]] + y[n[3]] - y[n[4]])/4.0
        b3 = (-y[n[1]] - y[n[2]] + y[n[3]] + y[n[4]])/4.0

        # Shape function integrals
        ∫N[1,i] = ((3b3-b2)*(3a1-a2)-(3a3-a2)*(3b1-b2))/9.0
        ∫N[2,i] = ((3b3+b2)*(3a1-a2)-(3a3+a2)*(3b1-b2))/9.0
        ∫N[3,i] = ((3b3+b2)*(3a1+a2)-(3a3+a2)*(3b1+b2))/9.0
        ∫N[4,i] = ((3b3-b2)*(3a1+a2)-(3a3-a2)*(3b1+b2))/9.0

        # Partial derivative integral terms
        ∫∂N∂x[1,i] = -b3+b1
        ∫∂N∂x[2,i] = b3+b1
        ∫∂N∂x[3,i] = b3-b1
        ∫∂N∂x[4,i] = -b3-b1

        ∫∂N∂y[1,i] = a3-a1
        ∫∂N∂y[2,i] = -a3-a1
        ∫∂N∂y[3,i] = -a3+a1
        ∫∂N∂y[4,i] = a3+a1

        # Partial derivatives
        J = a1*b3 - a3*b1
        @inbounds for n in 1:4
            ∂N∂x[n,i] = 0.25*∫∂N∂x[n,i]/J
            ∂N∂y[n,i] = 0.25*∫∂N∂y[n,i]/J
        end

        # Element weights
        @inbounds for n in 1:4
            element_weightc[n, i] = ∫N[n, i]
        end
    end
end

function get_▽●V!(
    u::NodeQuant, v::NodeQuant, ∂N∂x::CellQuant,
    ∂N∂y::CellQuant, nodelist::CellQuant, ▽●v::CellQuant
)
    ncells = size(∂N∂x,2)
    nnodes = size(u, 1)
    @inbounds @distributed for i in 1:nnodes
        ▽●v[i] = 0.0
    end

    # N.b. This is roughly 68 times slower than using a loop, even before
    # we add in the @distributed
    #▽●v.d .= 0.0

    @inbounds @distributed for i in 1:ncells
        @inbounds for j in 1:4
            node = nodelist[j,i]
            ▽●v[i] += u[node]*∂N∂x[j,i] + v[node]*∂N∂y[j,i]
        end
    end
end

function set_soundspeed!(
    pressure::CellQuant, ρ::CellQuant,
    material::CellQuant, soundspeed::CellQuant
)
    ncells = size(pressure, 1)

    @inbounds @distributed for i in 1:ncells
        soundspeed[i] = sqrt(gamma[material[i]]*pressure[i]/ρ[i])
    end
end
