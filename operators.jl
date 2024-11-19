function norms(levelset, boundary_indices)   
    normx = zeros(Float64, num_spaces+1, num_spaces+1)
    normy = zeros(Float64, num_spaces+1, num_spaces+1)
    for index in boundary_indices
        i = index[1]
        j = index[2]
        if i == 1
            absx = (levelset[i, j+1] - levelset[i, j-1]) / (2*dx)
            absy = (levelset[num_spaces+1, j] - levelset[i+1, j]) / (2*dx)
        elseif i == num_spaces+1
            absx = (levelset[i, j+1] - levelset[i, j-1]) / (2*dx)
            absy = (levelset[i-1, j] - levelset[1, j]) / (2*dx)
        else
            absx = (levelset[i, j+1] - levelset[i, j-1]) / (2*dx)
            absy = (levelset[i-1, j] - levelset[i+1, j]) / (2*dx)
        end
        normx[i, j] = absx / sqrt(absx^2 + absy^2)
        normy[i, j] = absy / sqrt(absx^2 + absy^2)
    end
    return normx, normy
end

#defining differential operators for a grid
function grad2(var)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 1:num_spaces+1
            #ignore boundary
            if boundary[i, j]
                continue
            #check if normal
            elseif boundary_label[i, j] == 0
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[i+1, j] - 4 * var[i, j])
            #check if top row
            elseif boundary_label[i, j] == 2
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[num_spaces+1, j] + var[i+1, j] - 4 * var[i, j])
            #check if bottom row
            elseif boundary_label[i, j] == 3
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[1, j] - 4 * var[i, j])
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
            if boundary[i, j]
                continue

            elseif boundary_label[i, j] == 0
                #normal gradient dot gradient
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[i-1, j] - var1[i+1, j]) * (var2[i-1, j] - var2[i+1, j]))

            elseif boundary_label[i, j] == 2
                #periodic definition
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[num_spaces+1, j] - var1[i+1, j]) * (var2[num_spaces+1, j] - var2[i+1, j]))

            elseif boundary_label[i, j] == 3
                #periodic definition
                x[i, j] = (1/(4*dx^2)) * ((var1[i, j+1] - var1[i, j-1]) * (var2[i, j+1] - var2[i, j-1])+(var1[i-1, j] - var1[1, j]) * (var2[i-1, j] - var2[1, j]))

            end
        end
    end
    return x
end
