const TINY = 1.0e-10  # Used to avoid zeros in divisions

function get_finite_elements!(
    x::NodeQuant, y::NodeQuant, nodelist::CellQuant,
    ∫N::CellQuant, ∫∂N∂x::CellQuant, ∫∂N∂y::CellQuant,
    ∂N∂x::CellQuant, ∂N∂y::CellQuant, element_weightc::CellQuant
)
    ncells = size(∫N, 2)

    @inbounds @sync @distributed for i in 1:ncells
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
    @inbounds @sync @distributed for i in 1:ncells
        ▽●v[i] = 0.0
    end

    # N.b. This is roughly 68 times slower than using a loop, even before
    # we add in the @sync @distributed
    #▽●v.d .= 0.0

    @inbounds @sync @distributed for i in 1:ncells
        @inbounds for j in 1:4
            node = nodelist[j,i]
            ▽●v[i] += u[node]*∂N∂x[j,i] + v[node]*∂N∂y[j,i]
        end
    end
end

function set_soundspeed!(
    pressure::CellQuant, ρ::CellQuant,
    material::CellQuant, soundspeed::CellQuant,
    gamma
)
    ncells = size(pressure, 1)

    @inbounds @sync @distributed for i in 1:ncells
        soundspeed[i] = sqrt(gamma[material[i]]*pressure[i]/ρ[i])
    end
end

function get_q!(
    ρ::CellQuant, soundspeed::CellQuant, ▽●v::CellQuant,
    x::NodeQuant, y::NodeQuant,
    u::NodeQuant, v::NodeQuant,
    area::CellQuant, nodelist::CellQuant, elelconn::CellQuant,
    boundary::NodeQuant,
    q::CellQuant,
    input::Input
)
    ncells = size(ρ, 1)

    if input.qform == 1
        @inbounds for i in 1:ncells
            if ▽●v[i] < 0.0
                ∂u∂x = sqrt(area[i])*▽●v[i]
                q[i] = input.cq*ρ[i]*∂u∂x^2 + input.cl*ρ[i]*soundspeed[i]*abs(∂u∂x)
            else
                q[i] = 0.0
            end
        end
    elseif input.qform == 2
        ΔuΔx_b = zeros(Float64, ncells)
        ΔuΔx_t = zeros(Float64, ncells)
        ΔuΔx_r = zeros(Float64, ncells)
        ΔuΔx_l = zeros(Float64, ncells)
        Δxtb = zeros(Float64, ncells)
        Δxlr = zeros(Float64, ncells)

        @inbounds @sync @distributed for i in 1:ncells
            @views n = nodelist.d[:, i]

            # Calculate leg lengths
            Lhor_y = -(x[n[4]] + x[n[3]] - x[n[2]] - x[n[1]])
            Lhor_x = y[n[4]] + y[n[3]] - y[n[2]] - y[n[1]]
            Lver_y = x[n[3]] + x[n[2]] - x[n[4]] - x[n[1]]
            Lver_x = -(y[n[3]] + y[n[2]] - y[n[4]] - y[n[1]])

            # Area-weighted length scales used in side qs
            d1 = sqrt(Lhor_y^2 + Lhor_x^2) + TINY
            d2 = sqrt(Lver_y^2 + Lver_x^2) + TINY
            Δxtb[i] = area[i]/d1
            Δxlr[i] = area[i]/d2

            # Velocity gradients
            ΔuΔx_b[i] = (Lhor_x*(u[n[2]] - u[n[1]]) + Lhor_y*(v[n[2]] - v[n[1]]))/area[i]
            ΔuΔx_l[i] = (Lhor_x*(u[n[4]] - u[n[1]]) + Lhor_y*(v[n[4]] - v[n[1]]))/area[i]
            ΔuΔx_t[i] = (Lhor_x*(u[n[3]] - u[n[4]]) + Lhor_y*(v[n[3]] - v[n[4]]))/area[i]
            ΔuΔx_r[i] = (Lhor_x*(u[n[3]] - u[n[2]]) + Lhor_y*(v[n[3]] - v[n[2]]))/area[i]
        end

        # Christensen monotonic limit to the velocity gradients
        Φb = zeros(Float64, ncells)
        Φr = zeros(Float64, ncells)
        Φt = zeros(Float64, ncells)
        Φl = zeros(Float64, ncells)

        @inbounds @sync @distributed for i in 1:ncells
            @views n = nodelist.d[:, i]
            @views e = elelconn.d[:, i]

            rlb = 0.0        # edge=bottom ratio=left
            rrb = 0.0        # edge=bottom ratio=right
            rbr = 0.0        # edge=right ratio=bottom
            rtr = 0.0        # edge=right ratio=top
            rlt = 0.0        # edge=top ratio=left
            rrt = 0.0        # edge=top ratio=right
            rbl = 0.0        # edge=left ratio=bottom
            rtl = 0.0        # edge=left ratio=top

            # For all the following cases, we assume that e[x] = 0 only for boundary
            # nodes, as we don't support disjoints.

            # Bottom edge
            if e[1] != 0
                rlb = nudge(ΔuΔx_l[e[1]]) / nudge(ΔuΔx_l[i])
                rrb = nudge(ΔuΔx_r[e[1]]) / nudge(ΔuΔx_r[i])
            else
                if (boundary[n[1]] < 0 || boundary[n[2]] < 0) && e[3] != 0
                    # External boundary
                    rlb = rrb = 1.0
                end
            end

            # Right edge
            if e[2] != 0
                rbr = nudge(ΔuΔx_l[e[2]]) / nudge(ΔuΔx_l[i])
                rtr = nudge(ΔuΔx_r[e[2]]) / nudge(ΔuΔx_r[i])
            else
                if (boundary[n[2]] < 0 || boundary[n[3]] < 0) && e[4] != 0
                    # External boundary
                    rbr = rtr = 1.0
                end
            end

            # Top edge
            if e[3] != 0
                rlt = nudge(ΔuΔx_l[e[3]]) / nudge(ΔuΔx_l[i])
                rrt = nudge(ΔuΔx_r[e[3]]) / nudge(ΔuΔx_r[i])
            else
                if (boundary[n[3]] < 0 || boundary[n[4]] < 0) && e[1] != 0
                    # External boundary
                    rlt = rrt = 1.0
                end
            end

            # Left edge
            if e[4] != 0
                rbl = nudge(ΔuΔx_l[e[4]]) / nudge(ΔuΔx_l[i])
                rtl = nudge(ΔuΔx_r[e[4]]) / nudge(ΔuΔx_r[i])
            else
                if (boundary[n[1]] < 0 || boundary[n[4]] < 0) && e[2] != 0
                    # External boundary
                    rbl = rtl = 1.0
                end
            end

            Φb[i] = max(0.0, min(0.5(rbl+rbr), 2rbl, 2rbr, 1.0))
            Φt[i] = max(0.0, min(0.5(rtl+rtr), 2rtl, 2rtr, 1.0))
            Φl[i] = max(0.0, min(0.5(rlt+rlb), 2rlt, 2rlb, 1.0))
            Φr[i] = max(0.0, min(0.5(rrt+rrb), 2rrt, 2rrb, 1.0))
        end

        # Calculate Qs for the element sides
        @inbounds @sync @distributed for i in 1:ncells
            # If ΔuΔx > 0, then this direction is expanding, so set to zero
            if ▽●v[i] < 0.0
                ΔuΔx_l[i] >= -TINY && (ΔuΔx_l[i] = 0.0)
                ΔuΔx_r[i] >= -TINY && (ΔuΔx_r[i] = 0.0)
                ΔuΔx_t[i] >= -TINY && (ΔuΔx_t[i] = 0.0)
                ΔuΔx_b[i] >= -TINY && (ΔuΔx_b[i] = 0.0)

                qb = input.cq*ρ[i]*(ΔuΔx_b[i]*Δxtb[i])^2*(1.0 - Φb[i]^2) +
                     input.cl*ρ[i]*soundspeed[i]*abs(ΔuΔx_b[i]*Δxtb[i])*(1.0 - Φb[i])

                qt = input.cq*ρ[i]*(ΔuΔx_t[i]*Δxtb[i])^2*(1.0 - Φt[i]^2) +
                     input.cl*ρ[i]*soundspeed[i]*abs(ΔuΔx_t[i]*Δxtb[i])*(1.0 - Φt[i])

                ql = input.cq*ρ[i]*(ΔuΔx_l[i]*Δxlr[i])^2*(1.0 - Φl[i]^2) +
                     input.cl*ρ[i]*soundspeed[i]*abs(ΔuΔx_l[i]*Δxlr[i])*(1.0 - Φl[i])

                qr = input.cq*ρ[i]*(ΔuΔx_r[i]*Δxlr[i])^2*(1.0 - Φr[i]^2) +
                     input.cl*ρ[i]*soundspeed[i]*abs(ΔuΔx_r[i]*Δxlr[i])*(1.0 - Φr[i])

                q[i] = 0.5(qb + qt + qr + ql)
            else
                q[i] = 0.0
            end
        end
    end
end

function nudge(a::Number)
    """
        nudge(a)

    Returns the value of a, nudged slightly away from zero
    """
    fsign(abs(a)+TINY, a)
end

function fsign(a::Number, b::Number)
    """
        fsign(a, b)

    Mimics the Fortran version of the sign function - returns a, but with the
    sign of b. I.e. sign(b)*abs(a) if b != 0, abs(a) if b == 0
    """
    if b != 0
        sign(b)*abs(a)
    else
        abs(a)
    end
end

function get_dt(
    area::CellQuant, soundspeed::CellQuant, ρ::CellQuant,
    q::CellQuant, dtold::Float64, time::Float64,
    input::Input
)
    ncells = size(area, 1)

    δt::Float64 = 0.0

    dtmin::Float64 = 1.0
    dt::Float64 = 0.0
    control::Int32 = 0
    ρcutoff::Float64 = 1.0e-6

    # WARNING: This cannot be a distributed loop
    @inbounds for i in 1:ncells
        if area[i] <= 0.0
            println("Negative area in cell $i. Area = $(area[i])")
            δt = 9999.9
        else
            δt = sqrt(area[i]/max(ρcutoff, soundspeed[i]^2 + 2q[i]/ρ[i]))/2.0
        end
        if δt < dtmin
            dtmin = δt
            control = i
        end
    end

    if time <= input.t0
        dt = min(dtmin, input.dtmax, input.dtinit)
    else
        dt = min(dtmin, input.dtmax, dtold*input.growth)
        if dt == dtold*input.growth
            control = -1
        end
    end

    return dt, control
end

function move_nodes!(
    dt::Float64, x::NodeQuant, y::NodeQuant,
    u::NodeQuant, v::NodeQuant,
    xnew::NodeQuant, ynew::NodeQuant
)

    nnodes = size(x, 1)

    @inbounds @sync @distributed for i in 1:nnodes
        xnew.d[i] = x[i] + dt*u[i]
        ynew.d[i] = y[i] + dt*v[i]
    end

end

function get_ρ!(
    mass::CellQuant, volume::CellQuant, ρ::CellQuant
)
    ncells = size(mass, 1)

    @inbounds @sync @distributed for i in 1:ncells
        ρ[i] = mass[i]/volume[i]
    end

end


function get_∫▽●V!(
    u::NodeQuant, v::NodeQuant, nodelist::CellQuant,
    ∫∂N∂x::CellQuant, ∫∂N∂y::CellQuant,
    ∫▽●V::CellQuant
)
"""
Calculates the summed quantity in equation 3.51 of Andy's thesis
"""
    ncells = size(nodelist, 2)

    @inbounds @sync @distributed for i in 1:ncells
        ∫▽●V[i] = 0.0
    end

    @inbounds @sync @distributed for i in 1:ncells
        @inbounds for j in 1:4
            node = nodelist[j,i]
            ∫▽●V[i] += u[node]*∫∂N∂x[j,i] + v[node]*∫∂N∂y[j,i]
        end
    end

end

function get_energy!(
    dt::T, pressure::CellQuant, q::CellQuant, mass::CellQuant,
    energy0::CellQuant, ∫▽●V::CellQuant, energy::CellQuant
) where {T <: Real}
    ncells = size(pressure, 1)

    @inbounds @sync @distributed for i in 1:ncells
        energy[i] = energy0[i] - dt*(pressure[i] + q[i])*∫▽●V[i]/mass[i]
    end

end

function perform_momentum_update!(
    dt::Float64, u::NodeQuant, v::NodeQuant,
    ρ::CellQuant, pressure::CellQuant, q::CellQuant,
    nodelist::CellQuant, ∫N::CellQuant, ∫∂N∂x::CellQuant,
    ∫∂N∂y::CellQuant, boundary::NodeQuant,
    area::CellQuant, soundspeed::CellQuant,
    x05::NodeQuant, y05::NodeQuant,
    uout::NodeQuant, vout::NodeQuant,
    region::CellQuant, regioncelliterator::Dict,
    mass_scatter_to_nodes::NodeQuant,
    force_scatter_to_nodes_x::NodeQuant,
    force_scatter_to_nodes_y::NodeQuant,
    input::Input
)

    nnodes = size(u, 1)
    ncells = size(ρ, 1)
    nreg = maximum(region)

    @inbounds @sync @distributed for i = 1:nnodes
        mass_scatter_to_nodes[i] = 0.0
        force_scatter_to_nodes_x[i] = 0.0
        force_scatter_to_nodes_y[i] = 0.0
    end

    # Scatter masses and forces to nodes.
    @inbounds @sync @distributed for i in 1:ncells
        @inbounds for j in 1:4
            node = nodelist[j,i]
            mass_scatter_to_nodes[node] += ρ[i]*∫N[j,i]
            force_scatter_to_nodes_x[node] += (pressure[i] + q[i])*∫∂N∂x[j,i]
            force_scatter_to_nodes_y[node] += (pressure[i] + q[i])*∫∂N∂y[j,i]
        end
    end

    # Anti-hourglass - Including this within the force update function is
    # Significantly faster than a nested function call...
    for reg in 1:nreg
        if input.anti_hg[reg] == ANTI_HG_DYNA_1 || input.anti_hg[reg] == ANTI_HG_DYNA_2
            if input.anti_hg[reg] == ANTI_HG_DYNA_1
                f = (c) -> -1.0*input.κ*ρ[c]*abs(area[c])/max(dt, dtmin_hg)
            else
                f = (c) -> -1.0*input.κ*ρ[c]*sqrt(area[c])*soundspeed[c]
            end

            @inbounds @sync @distributed for cell in regioncelliterator[reg]
                @views n = nodelist.d[:, cell]
                uγ = u[n[1]] - u[n[2]] + u[n[3]] - u[n[4]]
                vγ = v[n[1]] - v[n[2]] + v[n[3]] - v[n[4]]

                factor = f(cell)

                force_scatter_to_nodes_x[n[1,cell]] += factor*uγ
                force_scatter_to_nodes_x[n[2,cell]] -= factor*uγ
                force_scatter_to_nodes_x[n[3,cell]] += factor*uγ
                force_scatter_to_nodes_x[n[4,cell]] -= factor*uγ

                force_scatter_to_nodes_y[n[1,cell]] += factor*vγ
                force_scatter_to_nodes_y[n[2,cell]] -= factor*vγ
                force_scatter_to_nodes_y[n[3,cell]] += factor*vγ
                force_scatter_to_nodes_y[n[4,cell]] -= factor*vγ
            end
        elseif input.anti_hg[reg] == ANTI_HG_BF
            @inbounds @sync @distributed for cell in regioncelliterator[reg]
                @views n = nodelist.d[:, cell]
                a1 = 0.25(-x05[n[1]] + x05[n[2]] + x05[n[3]] - x05[n[4]])
                a2 = 0.25(x05[n[1]] - x05[n[2]] + x05[n[3]] - x05[n[4]])
                a3 = 0.25(-x05[n[1]] - x05[n[2]] + x05[n[3]] + x05[n[4]])
                b1 = 0.25(-y05[n[1]] + y05[n[2]] + y05[n[3]] - y05[n[4]])
                b2 = 0.25(y05[n[1]] - y05[n[2]] + y05[n[3]] - y05[n[4]])
                b3 = 0.25(-y05[n[1]] - y05[n[2]] + y05[n[3]] + y05[n[4]])

                bb = 4(b3^2 + b1^2 + a3^2 + a1^2)

                Δx = x05[n[1]] - x05[n[2]] + x05[n[3]] - x05[n[4]]
                Δy = y05[n[1]] - y05[n[2]] + y05[n[3]] - y05[n[4]]
                γ1 = 0.5 - 0.5(∫∂N∂x[1,cell]*Δx + ∫∂N∂y[1,cell]*Δy)/area[cell]
                γ2 = -0.5 - 0.5(∫∂N∂x[2,cell]*Δx + ∫∂N∂y[2,cell]*Δy)/area[cell]
                γ3 = 0.5 - 0.5(∫∂N∂x[3,cell]*Δx + ∫∂N∂y[3,cell]*Δy)/area[cell]
                γ4 = -0.5 - 0.5(∫∂N∂x[4,cell]*Δx + ∫∂N∂y[4,cell]*Δy)/area[cell]

                # qx, qy different forms according to whether artificial damping or
                # stiffness are required
                q = -input.κ*ρ[cell]*soundspeed[cell]*sqrt(bb)/10.0   # damping
                q -= 0.5*input.κ*ρ[cell]*dt*soundspeed[cell]^2*bb/area[cell] # stiffness

                qy = q*(γ1*v[n[1]] + γ2*v[n[2]] + γ3*v[n[3]] + γ4*v[n[4]])
                qx = q*(γ1*u[n[1]] + γ2*u[n[2]] + γ3*u[n[3]] + γ4*u[n[4]])

                force_scatter_to_nodes_x[n[1]] += γ1*qx
                force_scatter_to_nodes_y[n[1]] += γ1*qy
                force_scatter_to_nodes_x[n[2]] += γ2*qx
                force_scatter_to_nodes_y[n[2]] += γ2*qy
                force_scatter_to_nodes_x[n[3]] += γ3*qx
                force_scatter_to_nodes_y[n[3]] += γ3*qy
                force_scatter_to_nodes_x[n[4]] += γ4*qx
                force_scatter_to_nodes_y[n[4]] += γ4*qy
            end
        end
    end

    # F = ma => v = Fδt/m
    @inbounds @sync @distributed for node in 1:nnodes
        uout[node] = u[node] + dt*force_scatter_to_nodes_x[node]/mass_scatter_to_nodes[node]
        vout[node] = v[node] + dt*force_scatter_to_nodes_y[node]/mass_scatter_to_nodes[node]

        # Apply boundary conditions
        if boundary[node] == -1 || boundary[node] == -3
            uout[node] = u[node]
        end
        if boundary[node] == -2 || boundary[node] == -3
            vout[node] = v[node]
        end
    end

end
