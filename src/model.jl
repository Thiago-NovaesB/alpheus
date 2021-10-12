using JuMP, AmplNLWriter, Bonmin_jll

J = 1001;
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

# MOI.set(model, MOI.RawParameter("bonmin.algorithm"), "B-BB")
# MOI.set(model, MOI.RawParameter("bonmin.num_resolve_at_root"), 5)
MOI.set(model, MOI.RawParameter("max_iter"), 50000)



@variable(model, 0.0 <= H[1:J] <= L, start = 2.0);
@variable(model, 0.0 <= A[1:J] <= L, start = 20.0);
@variable(model, 0.0 <= SM[1:J] <= L, start = 0.0);
@variable(model, 0.0 <= V[1:J] <= L, start = 1.0);
@variable(model, 0.0 <= W[1:J] <= L, start = 0.0);
@variable(model, z[1:J], Bin);

@constraint(model, [j in 1:J-1], V[j]/(2*g) + H[j] + E[j] == V[j+1]/(2*g) + H[j+1] + E[j+1] + SM[j]*Δx + W[j]*z[j]);
@constraint(model, sum(z) <= f);
@constraint(model, z[1] == 0);
@constraint(model, z[J] == 0);
@constraint(model, H[J] == 2);
@constraint(model,[j in 1:J], A[j] == B*H[j]);
@NLexpression(model, R[j in 1:J], A[j]/(B+2*H[j]));
@constraint(model, [j in 1:J], V[j] <= g*H[j]);
@NLexpression(model, S[j in 1:J], V[j]*(n0^2+nt^2*z[j]+2*n0*nt*z[j])/(R[j]^(4/3)));
@NLconstraint(model,[j in 1:J-1], SM[j] == (S[j] + S[j+1])/2);
@constraint(model, [j in 1:J], H[j] >= Gt*z[j]);
@constraint(model, [j in 1:J-1], -M[j]*H[j] <= H[j+1] - H[j]);
@constraint(model, [j in 1:J-1], H[j+1] - H[j] <= M[j]*H[j]);

@NLobjective(model, Max, -H[1]*ηu + H[J]*ηd + sum(W[j]*z[j] for j in 1:J) - sum((V[j]*A[j]^2 - Q^2)^2 for j in 1:J)/J)

 @time status = optimize!(model)

W = value.(model[:W]);
H = value.(model[:H]);
V = value.(model[:V]);
A = value.(model[:A]);
z = value.(model[:z]);
S = value.(model[:S]);
SM = value.(model[:SM]);

