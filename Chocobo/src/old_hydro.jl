const TINY = 1.0e-10  # Used to avoid zeros in divisions

function get_finite_elements_OLD!(
    x::NodeQuant, y::NodeQuant, nodelist::CellQuant,
    ∫N::CellQuant, ∫∂N∂x::CellQuant, ∫∂N∂y::CellQuant,
    ∂N∂x::CellQuant, ∂N∂y::CellQuant, element_weightc::CellQuant
)
    ncells = size(∫N, 2)

    @inbounds for i in 1:ncells
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

function get_▽●V_OLD!(
    u::NodeQuant, v::NodeQuant, ∂N∂x::CellQuant,
    ∂N∂y::CellQuant, nodelist::CellQuant, ▽●v::CellQuant
)
    ncells = size(∂N∂x,2)
    ▽●v.d .= 0.0

    @inbounds for i in 1:ncells
        @inbounds for j in 1:4
            node = nodelist[j,i]
            ▽●v[i] += u[node]*∂N∂x[j,i] + v[node]*∂N∂y[j,i]
        end
    end
end

function set_soundspeed_OLD!(
    pressure::CellQuant, ρ::CellQuant,
    material::CellQuant, soundspeed::CellQuant
)
    ncells = size(pressure, 1)

    @inbounds for i in 1:ncells
        soundspeed[i] = sqrt(gamma[material[i]]*pressure[i]/ρ[i])
    end
end

function get_q_OLD!(
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
                q[i] = input.CQ*ρ[i]*∂u∂x^2 + input.CL*ρ[i]*soundspeed[i]*abs(∂u∂x)
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

        @inbounds for i in 1:ncells
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

        @inbounds for i in 1:ncells
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
        for i in 1:ncells
            # If ΔuΔx > 0, then this direction is expanding, so set to zero
            if ▽●v[i] < 0.0
                ΔuΔx_l[i] >= -TINY && (ΔuΔx_l[i] = 0.0)
                ΔuΔx_r[i] >= -TINY && (ΔuΔx_r[i] = 0.0)
                ΔuΔx_t[i] >= -TINY && (ΔuΔx_t[i] = 0.0)
                ΔuΔx_b[i] >= -TINY && (ΔuΔx_b[i] = 0.0)

                qb = input.CQ*ρ[i]*(ΔuΔx_b[i]*Δxtb[i])^2*(1.0 - Φb[i]^2) +
                     input.CL*ρ[i]*soundspeed[i]*abs(ΔuΔx_b[i]*Δxtb[i])*(1.0 - Φb[i])

                qt = input.CQ*ρ[i]*(ΔuΔx_t[i]*Δxtb[i])^2*(1.0 - Φt[i]^2) +
                     input.CL*ρ[i]*soundspeed[i]*abs(ΔuΔx_t[i]*Δxtb[i])*(1.0 - Φt[i])

                ql = input.CQ*ρ[i]*(ΔuΔx_l[i]*Δxlr[i])^2*(1.0 - Φl[i]^2) +
                     input.CL*ρ[i]*soundspeed[i]*abs(ΔuΔx_l[i]*Δxlr[i])*(1.0 - Φl[i])

                qr = input.CQ*ρ[i]*(ΔuΔx_r[i]*Δxlr[i])^2*(1.0 - Φr[i]^2) +
                     input.CL*ρ[i]*soundspeed[i]*abs(ΔuΔx_r[i]*Δxlr[i])*(1.0 - Φr[i])

                q[i] = 0.5(qb + qt + qr + ql)
            else
                q[i] = 0.0
            end
        end
    end
end

function get_dt_OLD(
    area::CellQuant, soundspeed::CellQuant, ρ::CellQuant,
    q::CellQuant, dtold::Float64, time::Float64,
    input::Input
)
    ncells = size(area, 1)

    δt::Float64 = 0.0

    dtmin::Float64 = 1.0
    dt::Float64 = 0.0
    control::Int32 = 0

    # WARNING: This cannot be a distributed loop
    @inbounds for i in 1:ncells
        if area[i] <= 0.0
            println("Negative area in cell $i. Area = $(area[i])")
            δt = 9999.9
        else
            δt = sqrt(area[i]/max(input.ρcutoff, soundspeed[i]^2 + 2q[i]/ρ[i]))/2.0
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

function move_nodes_OLD!(
    dt::Float64, x::NodeQuant, y::NodeQuant,
    u::NodeQuant, v::NodeQuant,
    xnew::NodeQuant, ynew::NodeQuant
)

    nnodes = size(x, 1)

    @inbounds for i in 1:nnodes
        xnew[i] = x[i] + dt*u[i]
        ynew[i] = y[i] + dt*v[i]
    end

end
