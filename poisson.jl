using SparseArrays, PyPlot
# UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>
# Arun Poisson Solver with Deformed Surface
# Poisson solver will evaluate the potential with the Robin boundary condition
#form: uxx = -f

function assemblePoissonA(n, nx, ny)
    N = (n+1)^2
    umap = reshape(1:N, n+1, n+1) # Index mapping from 2D grid to vector
    A = Tuple{Int64,Int64,Float64}[] # Array of matrix elements (row,col,value)
    # Main loop, insert stencil in matrix for each node point
    for j = 1:n+1
        for i = 1:n+1
        row = umap[i,j]
            if boundary[i, j] == 1
                #enforce +/- v as the potential on the right/left boundary
                if j == n+1
                    push!(A, (row, row, 1.0))
                else
                    push!(A, (row, row, 1.0))
                end

            elseif boundary_label[i, j] == 10
                #vertical right periodic t
                push!(A, (row, row, -(delta * epsilon * nx[i, j] + dx)))
                push!(A, (row, umap[i, j+1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[i+1, j], -delta/2 * epsilon * ny[i, j]))
                push!(A, (row, umap[n+1, j], delta/2 * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 11
                #vertical left
                push!(A, (row, row, -delta * epsilon * nx[i, j] + dx))
                push!(A, (row, umap[i, j-1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[i+1, j], delta/2 * epsilon * ny[i, j]))
                push!(A, (row, umap[n+1, j], -delta/2 * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 12
                #vertical right at bottom
                 push!(A, (row, row, -(delta * epsilon * nx[i, j] + dx)))
                 push!(A, (row, umap[i, j+1], delta * epsilon * nx[i, j]))
                 push!(A, (row, umap[1, j], -delta/2 * epsilon * ny[i, j]))
                 push!(A, (row, umap[i-1, j], delta/2 * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 13
                #vertical left
                push!(A, (row, row, -delta * epsilon * nx[i, j] + dx))
                push!(A, (row, umap[i, j-1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[1, j], delta/2 * epsilon * ny[i, j]))
                push!(A, (row, umap[i-1, j], -delta/2 * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 2
                #periodic condition for top boundary
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[n+1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))

            elseif boundary_label[i, j] == 3
                #periodic condition for bottom boundary
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))


            elseif boundary_label[i, j] == 4
                #vertical right
                 push!(A, (row, row, -(delta * epsilon * nx[i, j] + dx)))
                 push!(A, (row, umap[i, j+1], delta * epsilon * nx[i, j]))
                 push!(A, (row, umap[i+1, j], -delta/2 * epsilon * ny[i, j]))
                 push!(A, (row, umap[i-1, j], delta/2 * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 5
                #vertical left
                push!(A, (row, row, -delta * epsilon * nx[i, j] + dx))
                push!(A, (row, umap[i, j-1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[i+1, j], delta/2 * epsilon * ny[i, j]))
                push!(A, (row, umap[i-1, j], -delta/2 * epsilon * ny[i, j]))
                                
            elseif boundary_label[i, j] == 6
                #robin type condition for diagonal down
                push!(A, (row, row, delta * epsilon * ny[i, j] - delta * epsilon * nx[i, j] - dx))
                push!(A, (row, umap[i, j+1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[i+1, j], -delta * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 7
                #robin type condition for diagonal up
                push!(A, (row, row, -(delta * epsilon * nx[i, j] + delta * epsilon * ny[i, j] + dx)))
                push!(A, (row, umap[i, j+1], delta * epsilon * nx[i, j]))
                push!(A, (row, umap[i-1, j], delta * epsilon * ny[i, j]))

            elseif boundary_label[i, j] == 8
                #robin type condition for horizontal up
                push!(A, (row, row, -dx - delta * epsilon * ny[i, j]))
                push!(A, (row, umap[i, j+1], delta/2 * epsilon * nx[i, j]))
                push!(A, (row, umap[i-1, j], delta * epsilon * ny[i, j]))
                push!(A, (row, umap[i, j-1], -delta/2 * epsilon * nx[i, j]))

            elseif boundary_label[i, j] == 9 
                #robin type condition for horizontal down
                push!(A, (row, row, -dx + delta * epsilon * ny[i, j]))
                push!(A, (row, umap[i, j+1], delta/2 * epsilon * nx[i, j]))
                push!(A, (row, umap[i+1, j], -delta * epsilon * ny[i, j]))
                push!(A, (row, umap[i, j-1], -delta/2 * epsilon * nx[i, j]))

            else
                # Interior nodes, 5-point stencil
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))

            end
        end
    end
    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N, N)
    return A
end

function assemblePoissonb(n, f)
    N = (n+1)^2
    umap = reshape(1:N, n+1, n+1) # Index mapping from 2D grid to vector
    b = zeros(N)
    # Main loop, insert stencil in matrix for each node point
    for j = 1:n+1
        for i = 1:n+1
        row = umap[i,j]
            if boundary[i, j] == 1
                #enforce +/- v as the potential on the right/left boundary
                if j == n+1
                    b[row] = v
                else
                    b[row] = -v
                end

            elseif boundary_label[i, j] == 0 || boundary_label[i, j] == 2 || boundary_label[i, j] == 3
                #periodic condition
                b[row] = f[i, j] * dx^2 / epsilon^2

            else
                #everything else
                b[row] = v * dx
            end
        end
    end
    return b
end
