function to_flow_animation(grid, res, init; frame_interval=10)

    num_res = size(res)[2]
    println()
    anim = @animate for i in 1:num_res
        plot(grid, res[:, i], ylims=(0.95, 1.35), color="grey", label="Density", legend=:topleft)
        #plot!(size=(900, 400))
        plot!(grid, init[1, :], label="Initial Condition", line=[:dash], color="blue")
        xlabel!("x")
        ylabel!("Q(x, t)")
        title!("Frame $i out of $num_res")

        #println()
        print("\rAnimating frame $i out of $num_res")
    end every frame_interval

    return anim
end

function to_shock_animation(grid, res, init; frame_interval=10)

    num_res = size(res)[2]
    println()
    anim = @animate for i in 1:num_res
        plot(grid, res[:, i], ylims=(0.95, 2.1), color="grey", label="Density", legend=:topright)
        #plot!(size=(900, 400))
        plot!(grid, init[1, :], label="Initial Condition", line=[:dash], color="blue")
        xlabel!("x")
        ylabel!("Q(x, t)")
        title!("Frame $i out of $num_res")

        print("\rAnimating frame $i out of $num_res")
    end every frame_interval

    return anim
end
