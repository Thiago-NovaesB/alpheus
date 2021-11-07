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

    ### create expressions ###
    @NLexpression(model, W[j in 1:J], z[j]*5*(V[j]/2.8)^(3/2)/8);
    @NLexpression(model, A[j in 1:J], input.B*H[j]);
    @NLexpression(model, R[j in 1:J], A[j]/(input.B+2*H[j]));
    @NLexpression(model, S[j in 1:J], V[j]*(input.n0^2+input.nt^2*z[j]+2*input.n0*input.nt*z[j])/(R[j]^(4/3)));
    @NLexpression(model,SM[j in 1:J-1],(S[j] + S[j+1])/2);

    ### create constraint ###
    if !options.upstream
        @constraint(model, H[1] <= input.Hu);
    end
    if !options.downstream
        @constraint(model, H[J] >= input.Hd);
    end
    @constraint(model, sum(z) <= input.f);

    @constraint(model, [j in 1:J], V[j] >= input.Vmin*z[j]);
    @constraint(model, [j in 1:J], V[j] <= (input.Vmax - g*deepest)*z[j] + g*deepest);

    @constraint(model, z[1] == 0);
    @constraint(model, z[J] == 0);
    @constraint(model, [j in 1:J], H[j] >= input.Gt*z[j]);
    @constraint(model, [j in 1:J-1], -input.M*H[j] <= H[j+1] - H[j]);
    @constraint(model, [j in 1:J-1], H[j+1] - H[j] <= input.M*H[j]);
    @constraint(model, [j in 1:J], V[j] <= g*H[j]);
    @NLconstraint(model,[j in 1:J], V[j]*A[j]^2 == input.Q^2);
    @NLconstraint(model, [j in 1:J-1], V[j]/(2*g) + H[j] + preprocessor.E[j] == V[j+1]/(2*g) + H[j+1] + preprocessor.E[j+1] + SM[j]*preprocessor.Δx + W[j]*z[j]);

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

