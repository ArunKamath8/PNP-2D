function parseGrid(grid)
    n = num_spaces
    new = zeros(Int8, n+1, n+1)
    boundary_indices = []
    for j = 1:n+1
        for i = 1:n+1
            if grid[i, j]
                continue

            elseif i == 1 && j != 2 && j != n
                #top boundary but not next to electrode
                new[i, j] = 2

            elseif i == n+1 && j != 2 && j != n
                #bottom boundary but not next to electrode
                new[i, j] = 3

            elseif i == 1
                if j == 2
                    new[i, j] = 10
                    push!(boundary_indices, [i, j])
                else
                    new[i, j] = 11
                    push!(boundary_indices, [i, j])
                end
                
            elseif i == n+1
                if j == 2
                    new[i, j] = 12
                    push!(boundary_indices, [i, j])
                else
                    new[i, j] = 13
                    push!(boundary_indices, [i, j])
                end 
                
            elseif grid[i-1, j] && !grid[i, j+1] && !grid[i+1, j] && grid[i, j-1]
                #diagonal down
                new[i, j] = 6
                push!(boundary_indices, [i, j])

            elseif !grid[i-1, j] && !grid[i, j+1] && grid[i+1, j] && grid[i, j-1]
                #diagonal up
                new[i, j] = 7
                push!(boundary_indices, [i, j])

            elseif !grid[i-1, j] && !grid[i, j+1] && !grid[i+1, j] && grid[i, j-1]
                #vertical right
                new[i, j] = 4
                push!(boundary_indices, [i, j])

            elseif !grid[i-1, j] && grid[i, j+1] && !grid[i+1, j] && !grid[i, j-1]
                #vertial left
                new[i, j] = 5
                push!(boundary_indices, [i, j])

            elseif !grid[i-1, j] && !grid[i, j+1] && grid[i+1, j] && !grid[i, j-1]
                #horizontal up
                new[i, j] = 8
                push!(boundary_indices, [i, j])

            elseif grid[i-1, j] && !grid[i, j+1] && !grid[i+1, j] && !grid[i, j-1]
                #horizontal down
                new[i, j] = 9
                push!(boundary_indices, [i, j])
                
            end
        end
    end
    return new, boundary_indices
end