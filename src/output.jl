function grapfOutput!(data::alpheusData)
    x = data.preprocessor.x
    H = data.output.H
    U = data.output.U
    J = data.input.J
    upstream = data.options.upstream
    dowstream = data.options.downstream
    solver = data.options.solver
    z = data.output.z


    plot(x, H ,title="Altura",xlabel="X", ylabel="H", label="Altura(m)")
    scatter!(x[z.==1.0], H[z.==1.0], color = "green", label = "Turbina")
    savefig("altura_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    plot(x, U ,title="Velocidade",xlabel="X", ylabel="U", label="Velocidade(m/s)")
    scatter!(x[z.==1.0], U[z.==1.0], color = "green", label = "Turbina")
    savefig("velocidade_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    nothing
end

function grapfOutputEnum!(data::alpheusData)
    x = data.preprocessor.x
    W = data.output.W_enum
    J = data.input.J
    upstream = data.options.upstream
    dowstream = data.options.downstream
    solver = data.options.solver

    plot(x, W ,title="Energia",xlabel="X", ylabel="W", label="Energia(m)")
    savefig("energiaEnum_$(solver)_$(upstream)_$(dowstream)_$(J).png")

    nothing
end