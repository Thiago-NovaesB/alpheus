function createStructures()
    options = alpheusOptions()
    input = alpheusInput()
    preprocessor = alpheusPreprocessor()
    output = alpheusOutput()

    data = alpheusData(options,input,preprocessor,output)

    return data
end


function fillInput!(data::alpheusData)
    input = data.input

    input.J = J
    input.B = B
    input.θ = θ
    input.M = M
    input.Gt = Gt
    input.L = L
    input.n0 = n0
    input.nt = nt
    input.Q = Q
    input.ηu = ηu
    input.ηd = ηd
    input.f = f
    input.Hu = Hu
    input.Hd = Hd
    input.Vmin = Vmin
    input.Vmax = Vmax

    nothing
end

function fillPreprocessor!(data::alpheusData)
    preprocessor = data.preprocessor
    input = data.input

    preprocessor.I = tan(input.θ);
    preprocessor.Δx = input.L/(input.J-1);
    preprocessor.x = LinRange(0,input.L,input.J);
    preprocessor.E = (input.L.-preprocessor.x)*preprocessor.I;

    nothing
end

function fillOptions!(data::alpheusData)
    options = data.options

    options.upstream = upstream
    options.downstream = downstream
    options.enumerate = enumerate
    options.solver = solver
    options.parameters = parameters

    nothing
end