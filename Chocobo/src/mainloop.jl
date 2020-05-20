using Printf

function mainloop(mesh::Mesh, problem::String, input::Input)
    step = 0
    cend = 0

    time::Float64 = input.t0
    lastsilo::Float64 = input.t0
    dt::Float64 = input.dtinit
    step::Int32 = 1
    output_dumped = false

    # Create output directory
    try
        mkdir("$problem/results")
    catch
        nothing
    end

    @time while time <= input.tend
        # Get FEM elements
        get_finite_elements!(
            mesh.x, mesh.y, mesh.nodelist,
            mesh.∫N, mesh.∫∂N∂x, mesh.∫∂N∂y, mesh.∂N∂x, mesh.∂N∂y,
            mesh.element_weightc
        )

        # Get divergence
        get_▽●V!(
            mesh.u, mesh.v, mesh.∂N∂x, mesh.∂N∂y,
            mesh.nodelist, mesh.▽●v
        )

        # Set soundspeed
        set_soundspeed!(mesh.pressure, mesh.ρ, mesh.material, mesh.soundspeed, mesh.gamma)

        # Calculate artificial viscosity
        get_q!(
            mesh.ρ, mesh.soundspeed, mesh.▽●v,
            mesh.x, mesh.y,
            mesh.u, mesh.v,
            mesh.area, mesh.nodelist, mesh.elelconn,
            mesh.type,
            mesh.q,
            input
        )

        # Calculate timestep
        dt, control = get_dt(
            mesh.area, mesh.soundspeed, mesh.ρ, mesh.q, dt, time, input
        )
        time += dt
        @printf("%4d  %10.7f  %15.12f  %4d\n", step, time, dt, control)

        # Half timestep positions
        dt½ = dt/2
        move_nodes!(
            dt½, mesh.x, mesh.y,
            mesh.u, mesh.v,
            mesh.x½, mesh.y½
        )

        # Get half-timestep FEM elements
        get_finite_elements!(
            mesh.x½, mesh.y½, mesh.nodelist,
            mesh.∫N, mesh.∫∂N∂x, mesh.∫∂N∂y, mesh.∂N∂x, mesh.∂N∂y,
            mesh.element_weightc
        )

        # Half-timestep volume
        copyquant!(mesh.volume, mesh.volume_store)

        calculate_area_volume!(
            mesh.x½, mesh.y½, mesh.nodelist,
            mesh.volume½, mesh.area
        )

        # Half-timestep density
        get_ρ!(mesh.mass, mesh.volume½, mesh.ρ½)

        # Integrate divergence of v
        get_∫▽●V!(
            mesh.u, mesh.v, mesh.nodelist,
            mesh.∫∂N∂x, mesh.∫∂N∂y,
            mesh.∫▽●V
        )

        # Half-timestep energy
        get_energy!(
            dt½, mesh.pressure, mesh.q, mesh.mass, mesh.energy,
            mesh.∫▽●V, mesh.energy½
        )

        # Use EoS to calculate half-timestep pressure
        perfect_gas!(mesh.energy½, mesh.ρ½, mesh.pressure½, mesh.material, mesh.gamma)

        copyquant!(mesh.u, mesh.u_store)
        copyquant!(mesh.v, mesh.v_store)

        # Momentum update
        perform_momentum_update!(
            dt, mesh.u_store,  mesh.v_store, mesh.ρ½, mesh.pressure½, mesh.q,
            mesh.nodelist, mesh.∫N, mesh.∫∂N∂x, mesh.∫∂N∂y,
            mesh.type,
            mesh.area, mesh.soundspeed,
            mesh.x½, mesh.y½,
            mesh.u, mesh.v,
            mesh.region, mesh.regioncelliterator,
            mesh.mass_scatter_to_nodes,
            mesh.force_scatter_to_nodes_x,
            mesh.force_scatter_to_nodes_y,
            input
        )

        # Time-averaged velocities
        @inbounds @sync @distributed for i in 1:mesh.nnodes
            mesh.ubar[i] = (mesh.u_store[i] + mesh.u[i])/2
            mesh.vbar[i] = (mesh.v_store[i] + mesh.v[i])/2
        end

        move_nodes!(
            dt, mesh.x, mesh.y,
            mesh.ubar, mesh.vbar,
            mesh.x, mesh.y
        )

        # Get FEM elements, again
        get_finite_elements!(
            mesh.x, mesh.y, mesh.nodelist,
            mesh.∫N, mesh.∫∂N∂x, mesh.∫∂N∂y, mesh.∂N∂x, mesh.∂N∂y,
            mesh.element_weightc
        )

        # Full timestep volume
        calculate_area_volume!(
            mesh.x, mesh.y, mesh.nodelist,
            mesh.volume, mesh.area
        )

        # Full-timestep density
        get_ρ!(mesh.mass, mesh.volume, mesh.ρ)

        # Integrate divergence of v
        get_∫▽●V!(
            mesh.ubar, mesh.vbar, mesh.nodelist,
            mesh.∫∂N∂x, mesh.∫∂N∂y,
            mesh.∫▽●V
        )

        # Full-timestep energy
        get_energy!(
            dt, mesh.pressure½, mesh.q, mesh.mass, mesh.energy,
            mesh.∫▽●V, mesh.energy
        )

        # Use EoS to calculate full-timestep pressure
        perfect_gas!(mesh.energy, mesh.ρ, mesh.pressure, mesh.material, mesh.gamma)

        if !input.profile
            if (time >= 0.2 && !output_dumped) || input.debug_step_count > 0
                println("Text output file generated at time = $time")
                output_text(mesh, time, "$problem/results")
                output_dumped = true
            end

            if (time - lastsilo) >= input.dtsilo
                println("SILO Output file generated at time = $time")
                lastsilo = lastsilo + input.dtsilo
                output_silo(mesh, step, time, "$problem/results/visdump")
            end
        end
        if step == input.debug_step_count
            break
        end

        step += 1
    end
end
