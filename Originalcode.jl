#This is my original code I used to solve a 2D Blocking electrode problem
#contains all necessary parts
#2D grid, blocking electrodes.
using LinearAlgebra
using PyPlot
using Dates
fig, ax = subplots(figsize=(7,7)) 

#grid creation
num_spaces = 100
grid = zeros(Int8, num_spaces+1, num_spaces+1)
dx = 2 / num_spaces

#parameters
epsilon = 0.05
delta = 0.1
dt = 10^-5
v = 0.1

function grad2(var)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 2:num_spaces
            if i == 1
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[num_spaces+1, j] + var[i+1, j] - 4 * var[i, j])
            elseif i == num_spaces+1
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[1, j] - 4 * var[i, j])
            else
                x[i, j] = (1/dx^2)*(var[i, j+1] + var[i, j-1] + var[i-1, j] + var[i+1, j] - 4 * var[i, j])
            end
        end
    end
    return x
end

function grad(var)   
    x = zeros(Float64, num_spaces+1, num_spaces+1)
    for i in 1:num_spaces+1
        for j in 2:num_spaces
            if i == 1
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] + var[i+1, j] - var[num_spaces+1, j])
            elseif i == num_spaces+1
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] + var[1, j] - var[i-1, j])
            else
                x[i, j] = (1/(2*dx))*(var[i, j+1] - var[i, j-1] + var[i+1, j] - var[i-1, j])
            end
        end
    end
    return x
end

#initialize
c = ones(Float64, num_spaces+1, num_spaces+1)
rho = zeros(Float64, num_spaces+1, num_spaces+1)
#phi will get updated according to poissons equation
phi = zeros(Float64, num_spaces+1, num_spaces+1)



for time in 1:10^4
##enforce boundary conditions
    #solve poissons equation
    A, b = assemblePoisson(num_spaces, rho)
    phi = reshape(A \ b, num_spaces+1, num_spaces+1)

    #enforcing no flux condition at boundary: j = 1
    for i = 1:num_spaces+1
        deltaphi = phi[i, 2] - phi[i, 1]
        rho[i, 1] = (c[i, 2] * deltaphi + rho[i, 2])/(1 - deltaphi^2)
        c[i, 1] = -1 * (rho[i, 2]-rho[i, 1]) / deltaphi
    end

    #enforcing no flux condition at boundary: j = n+1
    for i = 1:num_spaces+1
        deltaphi = phi[i, num_spaces+1] - phi[i, num_spaces]
        c[i, num_spaces+1] = (c[i, num_spaces] - rho[i, num_spaces] * deltaphi)/(1 - deltaphi^2)
        rho[i, num_spaces+1] = (c[i, num_spaces] - c[i, num_spaces+1])/(deltaphi)
    end

    #update c and rho; boundary is taken care of at the top of the loop
    gradc = grad(c)
    gradrho = grad(rho)
    gradphi = grad(phi)
    grad2c = grad2(c)
    grad2rho = grad2(rho)
    grad2phi = grad2(phi)

    dcdt = epsilon .* (grad2c .+ gradrho .* gradphi .+ rho .* grad2phi)
    drhodt = epsilon .* (grad2rho .+ gradc .* gradphi .+ c .* grad2phi)
    c[:,2:num_spaces] .+= (dcdt[:,2:num_spaces] .* dt)
    rho[:,2:num_spaces] .+= (drhodt[:,2:num_spaces] .* dt)
end

using SparseArrays, PyPlot
# UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

function assemblePoisson(n, f)
    h = dx
    N = (n+1)^2
    umap = reshape(1:N, n+1, n+1) # Index mapping from 2D grid to vector
    A = Tuple{Int64,Int64,Float64}[] # Array of matrix elements (row,col,value)
    b = zeros(N)
    # Main loop, insert stencil in matrix for each node point
    for j = 1:n+1
        for i = 1:n+1
        row = umap[i,j]
            if i == 1 && j != 1 && j != n+1
                #periodic condition
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[n+1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] * dx^2 / epsilon^2
            elseif i == n+1 && j != 1 && j != n+1
                #periodic condition
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] * dx^2 / epsilon^2
            elseif j == 1
                # Robin boundary condition
                push!(A, (row, row, delta * epsilon + dx))
                push!(A, (row, umap[i, j+1], -delta * epsilon))
                b[row] = -v * dx
            elseif j == n+1
                # Robin boundary condition
                push!(A, (row, row, delta * epsilon + dx))
                push!(A, (row, umap[i, j-1], -delta * epsilon))
                b[row] = v * dx
            else
                # Interior nodes, 5-point stencil
                push!(A, (row, row, 4.0))
                push!(A, (row, umap[i+1,j], -1.0))
                push!(A, (row, umap[i-1,j], -1.0))
                push!(A, (row, umap[i,j+1], -1.0))
                push!(A, (row, umap[i,j-1], -1.0))
                b[row] = f[i, j] / epsilon^2 * dx^2
            end
        end
    end
    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N, N)
    return A, b
end