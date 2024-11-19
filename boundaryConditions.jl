function enforceNoFlux(c, rho, phi, nx, ny, boundary_indices)
    c_copy = copy(c)
    rho_copy = copy(rho)

    for index in boundary_indices
        i = index[1]
        j = index[2]

        if boundary_label[i, j] == 10
            #vertical right
            ngradphi = (phi[num_spaces+1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [nx[i, j] -ngradphi; -ngradphi nx[i, j]]
            soln = [(c_copy[num_spaces+1, j] - c_copy[i+1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j]; (rho_copy[num_spaces+1, j] - rho_copy[i+1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 11
            #vertical left
            ngradphi = (phi[num_spaces+1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [-nx[i, j] -ngradphi; -ngradphi -nx[i, j]]
            soln = [(c_copy[num_spaces+1, j] - c_copy[i+1, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j]; (rho_copy[num_spaces+1, j] - rho_copy[i+1, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 12
            #vertical right
            ngradphi = (phi[i-1, j] - phi[1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [nx[i, j] -ngradphi; -ngradphi nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j]; (rho_copy[i-1, j] - rho_copy[1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 13
            #vertical left
            ngradphi = (phi[i-1, j] - phi[1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [-nx[i, j] -ngradphi; -ngradphi -nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j]; (rho_copy[i-1, j] - rho_copy[i, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 4
             #vertical right
            ngradphi = (phi[i-1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [nx[i, j] -ngradphi; -ngradphi nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i+1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j]; (rho_copy[i-1, j] - rho_copy[i+1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 5
            #vertical left
            ngradphi = (phi[i-1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [-nx[i, j] -ngradphi; -ngradphi -nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i+1, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j]; (rho_copy[i-1, j] - rho_copy[i+1, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 6
            #diagonal down
            ngradphi = (phi[i, j] - phi[i+1, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [(nx[i, j]-ny[i, j]) -ngradphi; -ngradphi (nx[i, j]-ny[i, j])]
            soln = [-c_copy[i+1, j] * ny[i, j] + c_copy[i, j+1] * nx[i, j]; -rho_copy[i+1, j] * ny[i, j] + rho_copy[i, j+1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 7
            #diagonal up
            ngradphi = (phi[i-1, j] - phi[i, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [(nx[i, j]+ny[i, j]) -ngradphi; -ngradphi (nx[i, j]+ny[i, j])]
            soln = [c_copy[i-1, j] * ny[i, j] + c_copy[i, j+1] * nx[i, j]; (rho_copy[i-1, j]) * ny[i, j] + rho_copy[i, j+1] * nx[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 8
            #horizup 
            ngradphi = (phi[i-1, j] - phi[i, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j-1])/2 * nx[i, j]
            matrix = [ny[i, j] -ngradphi; -ngradphi ny[i, j]]
            soln = [(c_copy[i, j+1] - c_copy[i, j-1])/2 * nx[i, j] + c_copy[i-1, j] * ny[i, j]; (rho_copy[i, j+1] - rho_copy[i, j-1])/2 * nx[i, j] + rho_copy[i-1, j] * ny[i, j]]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 9
            #horizdown
            ngradphi = (phi[i, j] - phi[i+1, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j-1])/2 * nx[i, j]
            matrix = [-ny[i, j] -ngradphi; -ngradphi -ny[i, j]]
            soln = [(c_copy[i, j+1] - c_copy[i, j-1])/2 * nx[i, j] - c_copy[i+1, j] * ny[i, j]; (rho_copy[i, j+1] - rho_copy[i, j-1])/2 * nx[i, j] - rho_copy[i+1, j] * ny[i, j]]
            c[i, j], rho[i, j] = matrix \ soln
        end
    end
    return c, rho
end

function enforceBV(c, rho, phi, nx, ny, boundary_indices)
    c_copy = copy(c)
    rho_copy = copy(rho)
    Rs = []

    for index in boundary_indices
        i = index[1]
        j = index[2]

        if boundary_label[i, j] == 10
            #vertical right
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[num_spaces+1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [nx[i, j]+A*dx A*dx-ngradphi; A*dx-ngradphi nx[i, j]+A*dx]
            soln = [(c_copy[num_spaces+1, j] - c_copy[i+1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j] + B*dx; (rho_copy[num_spaces+1, j] - rho_copy[i+1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 11
            #vertical left
            A = exp(-(v-phi[i, j])/2) * kf
            B = exp((v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[num_spaces+1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [-nx[i, j]+A*dx A*dx-ngradphi; A*dx-ngradphi -nx[i, j]+A*dx]
            soln = [(c_copy[num_spaces+1, j] - c_copy[i+1, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j] + B*dx; (rho_copy[num_spaces+1, j] - rho_copy[i+1, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 12
            #vertical right
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [A*dx+nx[i, j] A*dx-ngradphi; A*dx-ngradphi A*dx+nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j] + B*dx; (rho_copy[i-1, j] - rho_copy[1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 13
            #vertical left
            A = exp(-(v-phi[i, j])/2) * kf
            B = exp((v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [A*dx-nx[i, j] A*dx-ngradphi; A*dx-ngradphi A*dx-nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j] + B*dx; (rho_copy[i-1, j] - rho_copy[i, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 4
            #vertical right
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [A*dx+nx[i, j] A*dx+-ngradphi; A*dx+-ngradphi A*dx+nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i+1, j])/2 * ny[i, j] + c_copy[i, j+1] * nx[i, j] + B*dx; (rho_copy[i-1, j] - rho_copy[i+1, j])/2 * ny[i, j] + rho_copy[i, j+1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 5
            #vertical left
            A = exp(-(v-phi[i, j])/2) * kf
            B = exp((v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[i+1, j])/2 * ny[i, j] + (phi[i, j] - phi[i, j-1]) * nx[i, j]
            matrix = [A*dx-nx[i, j] A*dx-ngradphi; A*dx-ngradphi A*dx-nx[i, j]]
            soln = [(c_copy[i-1, j] - c_copy[i+1, j])/2 * ny[i, j] - c_copy[i, j-1] * nx[i, j] + B*dx; (rho_copy[i-1, j] - rho_copy[i+1, j])/2 * ny[i, j] - rho_copy[i, j-1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 6
            #diagonal down
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i, j] - phi[i+1, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [A*dx+(nx[i, j]-ny[i, j]) A*dx+-ngradphi; A*dx+-ngradphi A*dx+(nx[i, j]-ny[i, j])]
            soln = [-c_copy[i+1, j] * ny[i, j] + c_copy[i, j+1] * nx[i, j] + B*dx; -rho_copy[i+1, j] * ny[i, j] + rho_copy[i, j+1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 7
            #diagonal up
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[i, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j]) * nx[i, j]
            matrix = [A*dx+(nx[i, j]+ny[i, j]) A*dx+-ngradphi; A*dx+-ngradphi A*dx+(nx[i, j]+ny[i, j])]
            soln = [c_copy[i-1, j] * ny[i, j] + c_copy[i, j+1] * nx[i, j] + B*dx; (rho_copy[i-1, j]) * ny[i, j] + rho_copy[i, j+1] * nx[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 8
            #horizup 
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i-1, j] - phi[i, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j-1])/2 * nx[i, j]
            matrix = [A*dx+ny[i, j] A*dx+-ngradphi; A*dx+-ngradphi A*dx+ny[i, j]]
            soln = [(c_copy[i, j+1] - c_copy[i, j-1])/2 * nx[i, j] + c_copy[i-1, j] * ny[i, j] + B*dx; (rho_copy[i, j+1] - rho_copy[i, j-1])/2 * nx[i, j] + rho_copy[i-1, j] * ny[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln

        elseif boundary_label[i, j] == 9
            #horizdown
            A = exp(-(-v-phi[i, j])/2) * kf
            B = exp((-v-phi[i, j])/2) * kb / Cb
            ngradphi = (phi[i, j] - phi[i+1, j]) * ny[i, j] + (phi[i, j+1] - phi[i, j-1])/2 * nx[i, j]
            matrix = [A*dx-ny[i, j] A*dx-ngradphi; A*dx-ngradphi A*dx-ny[i, j]]
            soln = [(c_copy[i, j+1] - c_copy[i, j-1])/2 * nx[i, j] - c_copy[i+1, j] * ny[i, j] + B*dx; (rho_copy[i, j+1] - rho_copy[i, j-1])/2 * nx[i, j] - rho_copy[i+1, j] * ny[i, j] + B*dx]
            c[i, j], rho[i, j] = matrix \ soln
        end
        push!(Rs, A*(rho[i, j] + c[i, j]) - B)
    end
    return c, rho, Rs
end