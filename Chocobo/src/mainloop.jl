function mainloop(mesh::Mesh, problem::String)
    count = 0
    cend = 0
    while count <= cend
        # Get FEM elements
        @time get_finite_elements!(
            mesh.x, mesh.y, mesh.nodelist,
            mesh.∫N, mesh.∫∂N∂x, mesh.∫∂N∂y, mesh.∂N∂x, mesh.∂N∂y,
            mesh.element_weightc
        )

        # Get divergence
        @time get_▽●V!(
            mesh.u, mesh.v, mesh.∂N∂x, mesh.∂N∂y,
            mesh.nodelist, mesh.▽●v
        )

        # Set soundspeed
        @time set_soundspeed!(mesh.pressure, mesh.ρ, mesh.material, mesh.soundspeed)

        # Calculate artificial viscosity
        # @time get_q_OLD!(
        #     mesh.ρ, mesh.soundspeed, mesh.▽●v,
        #     mesh.x, mesh.y,
        #     mesh.u, mesh.v,
        #     mesh.area, mesh.nodelist, mesh.elelconn,
        #     mesh.type,
        #     mesh.q
        # )

        count += 1
    end
end
