function grapfOutput!(data::alpheusData)
    x = data.preprocessor.x
    H = data.output.H
    U = data.output.U
    J = data.input.J
    upstream = data.options.upstream
    dowstream = data.options.downstream
    solver = data.options.solver
    z = data.output.z


    plot(x, H ,title="Altura",xlabel="X (m)", ylabel="H (m)", label="Altura (m)",legend=:topleft)
    scatter!(x[z.==1.0], H[z.==1.0], color = "green", label = "Turbina")
    savefig("altura_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    plot(x, U ,title="Velocidade",xlabel="X (m) ", ylabel="U (m/s)", label="Velocidade (m/s)")
    scatter!(x[z.==1.0], U[z.==1.0], color = "green", label = "Turbina")
    savefig("velocidade_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    nothing
end

function textOutput!(data::alpheusData)
    x = data.preprocessor.x
    H = data.output.H
    U = data.output.U
    W = data.output.W
    J = data.input.J
    upstream = data.options.upstream
    dowstream = data.options.downstream
    solver = data.options.solver
    z = data.output.z

    open("output_$(J)","w") do file
        write(file, "$(J) $(data.preprocessor.x'z) $(H'z) $(U'z) $(W'z)")
    end


    nothing
end

function grapfOutputEnum!(data::alpheusData)
    x = data.preprocessor.x
    W = data.output.W_enum
    J = data.input.J
    upstream = data.options.upstream
    dowstream = data.options.downstream
    solver = data.options.solver

    W = W*g*density*data.input.Q

    plot(x, W ,title="Potência",xlabel="X (m)", ylabel="Y (W)", label="Potência(W)",legend=:topright)
    savefig("energiaEnum_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    nothing
end