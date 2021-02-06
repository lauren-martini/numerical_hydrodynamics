function flux_step(Q, istart, iend)
    flux = zero(Q)
    for i in range(istart, iend, step=1)
        # Velocity across the boundary
        u_boundary = 0.5((Q[2, i]/Q[1, i]) + (Q[2, i-1]/Q[1, i-1]))
        # u_boundary = calc_u_boundary(Q, i)

        if u_boundary >= 0
            flux[:, i] = Q[:, i-1] * u_boundary
        else
            flux[:, i] = Q[:, i] * u_boundary
        end
    end
    return flux
end
