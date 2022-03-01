function singleModel!(data::alpheusData)
    options = data.options
    input = data.input
    preprocessor = data.preprocessor

    ### initialize model ###
    if options.solver == bonmin
        model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe));
    elseif options.solver == couenne
        model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe));
    elseif options.solver == ipopt
        model = Model(Ipopt.Optimizer)
    end

    ### set parameters ###
    for (key, value) in options.parameters
        set_optimizer_attribute(model, key, value)
    end

    ### create variables ###
    @variable(model, 0.0 <= H[1:J] <= deepest, start = 1.0);
    @variable(model, 0.0 <= V[1:J] <= g*deepest, start = 1.0);
    @variable(model, z[1:J], Bin);
    @variable(model, 0.0 <= W[1:J], start = 1.0 );

    ### create expressions ###
    #=5.9=#@NLexpression(model, A[j in 1:J], input.B*H[j]);
    #=5.10=#@NLexpression(model, R[j in 1:J], A[j]/(input.B+2*H[j]));
    #=5.11=#@NLexpression(model, S[j in 1:J], V[j]*(input.n0^2+input.nt^2*z[j]+2*input.n0*input.nt*z[j])/(R[j]^(4/3)));
    #=5.12=#@NLexpression(model, SM[j in 1:J-1],(S[j] + S[j+1])/2);

    ### create constraint ###
    if !options.upstream
        @constraint(model, H[1] <= input.Hu);
    end
    if !options.downstream
        @constraint(model, H[J] >= input.Hd);
    end
    #=5.6=#@constraint(model, sum(z) == 1);

    #=5.8=#@constraint(model, z[1] == 0);
    #=5.8=#@constraint(model, z[J] == 0);
    #=5.14=#@constraint(model, [j in 1:J-1], -(1-z[j])*input.M*H[j]-z[j]*deepest <= H[j+1] - H[j]);
    #=5.14=#@constraint(model, [j in 1:J-1], H[j+1] - H[j] <= (1-z[j])*input.M*H[j]+z[j]*deepest);
    #=5.7=#@constraint(model, [j in 1:J], V[j] <= g*H[j]);
    #=5.13=#@NLconstraint(model,[j in 1:J], V[j]*A[j]^2 == input.Q^2);
    #=5.5=#@NLconstraint(model, [j in 1:J-1], V[j]/(2*g) + H[j] + preprocessor.E[j] == V[j+1]/(2*g) + H[j+1] + preprocessor.E[j+1] + SM[j]*preprocessor.Δx + W[j]*z[j]);

    ### create objective ###

    @NLobjective(model, Max,sum(W[j]*z[j] for j in 1:J) -(input.Hu-H[1])*input.ηu*options.upstream + (H[J]-input.Hd)*input.ηd*options.downstream)
    
    optimizeModel!(model)
    fillOutput!(model,data)
    nothing
end

function optimizeModel!(model::Model)
    @time status = optimize!(model)
    nothing
end

function fillOutput!(model::Model,data::alpheusData)
    output = data.output

    output.W = value.(model[:W]);
    output.H = value.(model[:H]);
    output.z = value.(model[:z]);
    V = value.(model[:V]);
    output.U = sqrt.(V)

    nothing
end

