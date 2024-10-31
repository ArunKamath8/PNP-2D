#defining differential operators for a grid
function grad2(var)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 1:num_spaces+1
            #ignore boundary
            if grid[i, j] == 1
                continue
            #check if normal
            elseif grid[i, j] == 0
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[i+1, j] - 4 * var[i, j])
            #check if top row
            elseif grid[i, j] == 2
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[num_spaces+1, j] + var[i+1, j] - 4 * var[i, j])
            #check if bottom row
            elseif grid[i, j] == 3
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[1, j] - 4 * var[i, j])
            end
        end
    end
    return x
end


function grad(var)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 1:num_spaces+1
            #ignore boundary
            if grid[i, j] == 1
                continue

            elseif grid[i, j] == 0
                #standard 4 - point gradient
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] - var[i+1, j] + var[i-1, j])

            elseif grid[i, j] == 2
                #periodic definition of gradient
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] - var[i+1, j] + var[num_spaces+1, j])

            elseif grid[i, j] == 3
                #periodic definition of gradient
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] - var[1, j] + var[i-1, j])

            end
        end
    end
    return x
end

function graddotgrad(var1, var2)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 1:num_spaces+1
            #ignore boundary
            if grid[i, j] == 1
                continue

            elseif grid[i, j] == 0
                #normal gradient dot gradient
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[i-1, j] - var1[i+1, j]) * (var2[i-1, j] - var2[i+1, j]))

            elseif grid[i, j] == 2
                #periodic definition
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[num_spaces+1, j] - var1[i+1, j]) * (var2[num_spaces+1, j] - var2[i+1, j]))

            elseif grid[i, j] == 3
                #periodic definition
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[i-1, j] - var1[1, j]) * (var2[i-1, j] - var2[1, j]))

            end
        end
    end
    return x
end

