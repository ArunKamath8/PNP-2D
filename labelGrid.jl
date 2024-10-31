#grid parser -> given a grid, will label points that are either periodic, (top and bottom boundary) or next to a boundary
function parseGrid(grid)
    n = num_spaces
    for j = 1:n+1
        for i = 1:n+1
            if grid[i, j] == 1
                continue

            elseif i == 1 && j != 2 && j != n
                #top boundary but not next to electrode
                grid[i, j] = 2

            elseif i == n+1 && j != 2 && j != n
                #bottom boundary but not next to electrode
                grid[i, j] = 3

            elseif i == 1
                if j == 2
                    grid[i, j] = 4
                else
                    grid[i, j] = 5
                end
                
            elseif i == n+1
                if j == 2
                    grid[i, j] = 4
                else
                    grid[i, j] = 5
                end 
                
            elseif grid[i-1, j] == 1 && grid[i, j+1] != 1 && grid[i+1, j] != 1 && grid[i, j-1] == 1 
                #diagonal down
                grid[i, j] = 6

            elseif grid[i-1, j] != 1 && grid[i, j+1] != 1 && grid[i+1, j] == 1 && grid[i, j-1] == 1 
                #diagonal up
                grid[i, j] = 7

            elseif grid[i-1, j] != 1 && grid[i, j+1] != 1 && grid[i+1, j] != 1 && grid[i, j-1] == 1 
                #vertical right
                grid[i, j] = 4

            elseif grid[i-1, j] != 1 && grid[i, j+1] == 1 && grid[i+1, j] != 1 && grid[i, j-1] != 1 
                #vertial left
                grid[i, j] = 5

            elseif grid[i-1, j] != 1 && grid[i, j+1] != 1 && grid[i+1, j] == 1 && grid[i, j-1] != 1
                #horizontal up
                grid[i, j] = 8 

            elseif grid[i-1, j] == 1 && grid[i, j+1] != 1 && grid[i+1, j] != 1 && grid[i, j-1] != 1
                #horizontal down
                grid[i, j] = 9
                
            end
        end
    end
    return grid
end