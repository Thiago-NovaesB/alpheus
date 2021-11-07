function explicityEnumeration!(data::alpheusData)

    options = data.options
    input = data.input
    preprocessor = data.preprocessor
    initialiseEnumerateOutput!(data)

    for k in 2:input.J-1
        ### enumerate ###
        z = zeros(J)
        z[k] = 1

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

        @constraint(model, [j in 1:J], V[j] >= input.Vmin*z[j]);
        @constraint(model, [j in 1:J], V[j] <= (input.Vmax - g*deepest)*z[j] + g*deepest);

        @constraint(model, [j in 1:J], H[j] >= input.Gt*z[j]);
        @constraint(model, [j in 1:J-1], -input.M*H[j] <= H[j+1] - H[j]);
        @constraint(model, [j in 1:J-1], H[j+1] - H[j] <= input.M*H[j]);
        @constraint(model, [j in 1:J], V[j] <= g*H[j]);
        @NLconstraint(model,[j in 1:J], V[j]*A[j]^2 == input.Q^2);
        @NLconstraint(model, [j in 1:J-1], V[j]/(2*g) + H[j] + preprocessor.E[j] == V[j+1]/(2*g) + H[j+1] + preprocessor.E[j+1] + SM[j]*preprocessor.Δx + W[j]*z[j]);

        ### create objective ###

        @NLobjective(model, Max,sum(W[j]*z[j] for j in 1:J) -(input.Hu-H[1])*input.ηu*options.upstream + (H[J]-input.Hd)*input.ηd*options.downstream)
        
        optimize!(model)
        status = termination_status(model)
        @show status
        if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            fillOutputEnumeration!(model,data,k)
        end
        
    end
    nothing
end

function initialiseEnumerateOutput!(data::alpheusData)
    output = data.output
    input = data.input

    output.W_enum = zeros(input.J)
    output.H_enum = zeros(input.J,input.J)
    output.U_enum = zeros(input.J,input.J)
    output.z_enum = zeros(input.J,input.J)
    nothing
end

function fillOutputEnumeration!(model::Model,data::alpheusData,k::Int)
    output = data.output

    output.W_enum[k] = value.(model[:W])[k];
    output.H_enum[k,:] = value.(model[:H]);
    V = value.(model[:V]);
    output.U_enum[k,:] = sqrt.(V)
    output.z_enum[k,k] = 1
    nothing
end

