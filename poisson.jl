using SparseArrays, PyPlot
# UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>
# Arun Poisson Solver with Deformed Surface
# Poisson solver will evaluate the potential with the Robin boundary condition
#form: uxx = -f

function assemblePoisson(n, f)
    N = (n+1)^2
    umap = reshape(1:N, n+1, n+1) # Index mapping from 2D grid to vector
    A = Tuple{Int64,Int64,Float64}[] # Array of matrix elements (row,col,value)
    b = zeros(N)
    # Main loop, insert stencil in matrix for each node point
    for j = 1:n+1
        for i = 1:n+1
        row = umap[i,j]
            if grid[i, j] == 1
                #enforce +/- v as the potential on the right/left boundary
                if j == n+1
                    push!(A, (row, row, 1.0))
                    b[row] = v
                else
                    push!(A, (row, row, 1.0))
                    b[row] = -v
                end

            elseif grid[i, j] == 2
                #periodic condition for top boundary
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[n+1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] * dx^2 / epsilon^2

            elseif grid[i, j] == 3
                #periodic condition for bottom boundary
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] * dx^2 / epsilon^2

            elseif grid[i, j] == 4
                #vertical right
                 push!(A, (row, row, -(delta * epsilon + dx)))
                 push!(A, (row, umap[i, j+1], delta * epsilon))
                 b[row] = v * dx

            elseif grid[i, j] == 5
                #vertical left
                push!(A, (row, row, delta * epsilon + dx))
                push!(A, (row, umap[i, j-1], -delta * epsilon))
                b[row] = v * dx
                                
            elseif grid[i, j] == 6
                #robin type condition for diagonal down
                push!(A, (row, row, -2*delta * epsilon - sqrt(2)*dx))
                push!(A, (row, umap[i, j+1], delta * epsilon))
                push!(A, (row, umap[i+1, j], delta * epsilon))
                b[row] = sqrt(2) * v * dx

            elseif grid[i, j] == 7
                #robin type condition for diagonal up
                push!(A, (row, row, -2*delta * epsilon - sqrt(2)*dx))
                push!(A, (row, umap[i, j+1], delta * epsilon))
                push!(A, (row, umap[i-1, j], delta * epsilon))
                b[row] = sqrt(2) * v * dx

            elseif grid[i, j] == 8
                #robin type condition for horizontal up
                push!(A, (row, row, -(delta * epsilon + dx)))
                push!(A, (row, umap[i-1, j], delta * epsilon))
                b[row] = v * dx

            elseif grid[i, j] == 9 
                #robin type condition for horizontal down
                push!(A, (row, row, -(delta * epsilon + dx)))
                push!(A, (row, umap[i+1, j], delta * epsilon))
                b[row] = v * dx

            else
                # Interior nodes, 5-point stencil
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] * dx^2 / epsilon^2
            end
        end
    end
    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N, N)
    return A, b
end