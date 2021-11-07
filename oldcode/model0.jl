using JuMP, AmplNLWriter, Bonmin_jll, Couenne_jll
using Plots

J = 50;
B = 10.0;
θ = 0.0024;
M = zeros(J).+0.1;
Gt = 1.0;
g = 9.81;
L = 1000.0;
n0 = 0.033;
nt = 0.0133;
Q = 20;
ηu = 0.0;
ηd = 0.0;
f = 1;

I = tan(θ);
Δx = L/(J-1);
x = LinRange(0,L,J);
E = (L.-x)*I;

model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe));
# model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe));

# MOI.set(model, MOI.RawParameter("bonmin.algorithm"), "B-BB")
# MOI.set(model, MOI.RawParameter("bonmin.num_resolve_at_root"), 5)
MOI.set(model, MOI.RawParameter("max_iter"), 5000)



@variable(model, 0.0 <= H[1:J] <= Q, start = 2.0);
@variable(model, 0.0 <= A[1:J] <= Q*Q, start = 20.0);
@variable(model, 0.0 <= SM[1:J] <= Q, start = 0.01);
@variable(model, 0.0 <= V[1:J] <= Q, start = 1.0);
# @variable(model, 0.0 <= W[1:J] <= Q, start = 0.01);
@variable(model, z[1:J], Bin);
@NLexpression(model,W[j in 1:J], 5*(V[j]/2.8)^(3/2)/8);


@NLconstraint(model, [j in 1:J-1], V[j]/(2*g) + H[j] + E[j] == V[j+1]/(2*g) + H[j+1] + E[j+1] + SM[j]*Δx + W[j]*z[j]);
@constraint(model, sum(z) <= f);
@constraint(model, z[1] == 0);
@constraint(model, z[J] == 0);
@constraint(model,[j in 1:J], A[j] == B*H[j]);
@NLconstraint(model,[j in 1:J], V[j]*A[j]^2 == Q^2);
@NLexpression(model, R[j in 1:J], A[j]/(B+2*H[j]));
@constraint(model, [j in 1:J], V[j] <= g*H[j]);
@NLexpression(model, S[j in 1:J], V[j]*(n0^2+nt^2*z[j]+2*n0*nt*z[j])/(R[j]^(4/3)));
@NLconstraint(model,[j in 1:J-1], SM[j] == (S[j] + S[j+1])/2);
@constraint(model, [j in 1:J], H[j] >= Gt*z[j]);
@constraint(model, [j in 1:J-1], -M[j]*H[j] <= H[j+1] - H[j]);
@constraint(model, [j in 1:J-1], H[j+1] - H[j] <= M[j]*H[j]);

@NLobjective(model, Max, -H[1]*ηu + H[J]*ηd + sum(W[j]*z[j] for j in 1:J) )

@time status = optimize!(model)

W = value.(model[:W]);
H = value.(model[:H]);
V = value.(model[:V]);
A = value.(model[:A]);
z = value.(model[:z]);
S = value.(model[:S]);
SM = value.(model[:SM]);
U = sqrt.(V)



plot(x,[H,U] ,title=["Altura" "Velocidade"],xlabel="X", ylabel=["H" "U"],layout=(2,1))
# scatter!(x,H)
# scatter!(x,U)
# annotate!(x[4],H[4],"Turbine")
# annotate!(x[4],U[4],"Turbine")
# for i in 1:J
#     if z[i] == 1.0
#         @show i
#         annotate!(x[i],H[i],"Turbine")
#         annotate!(x[i],U[i],"Turbine")
#     end
# end