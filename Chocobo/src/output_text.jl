using Formatting

count = -1

function output_text(mesh::Mesh, time, loc)

global count += 1

# Physical Quantities
name = "$loc/mypre_no" * fmt("03d",count) * ".txt"
open(name, "w") do io
    for i in  1:mesh.ncells
        line = format(
            "{1:6d}  {2:17.12f}  {3:17.12f}  {4:19.12f}\n",
            i, mesh.pressure[i], mesh.ρ[i], mesh.energy[i]
        )
        write(io, line)
    end
end

name = "$loc/myvel_no" * fmt("03d",count) * ".txt"
open(name, "w") do io
    for i in  1:mesh.nnodes
        line = format(
            "{1:17.12f}  {2:17.12f}  {3:17.12f}\n",
            time, mesh.u[i], mesh.v[i]
        )
        write(io, line)
    end
end


name = "$loc/myx_no" * fmt("03d",count) * ".txt"
open(name, "w") do io
    for i in  1:mesh.nnodes
        line = format(
            "{1:17.12f}  {2:17.12f}\n",
            mesh.x[i], mesh.y[i]
        )
        write(io, line)
    end
end

name = "$loc/nodelist_no" * fmt("03d",count) * ".txt"
open(name, "w") do io
    for i in  1:mesh.ncells
        line = format(
            "{1:6d}  {2:6d}  {3:6d}  {4:6d}\n",
            mesh.nodelist[1,i], mesh.nodelist[2,i], mesh.nodelist[3,i], mesh.nodelist[4,i],
        )
        write(io, line)
    end
end

end
