using CPLEX
using JuMP
using Plots
using StatsPlots
using Random

T = 12
M = 4
E = 3
Random.seed!(0)
d = [rand(20:70) for _ in 1:T]
f = [10,30,60,90]
e = [8,6,4,2]
R = 2

model = Model(CPLEX.Optimizer)

@variable(model, x[1:T, 1:M] ≥ 0)
@variable(model, y[1:T, 1:M], Bin)
@variable(model, s[1:T] ≥ 0)

@objective(model, Min, sum(f[m]*y[t,m] for m in 1:M, t in 1:T) + sum(s))


@constraint(model, [t=2:T], sum(x[t,m] for m in 1:M) - s[t] + s[t-1] == d[t])
@constraint(model, sum(x[1,m] for m in 1:M) - s[1]  == d[1] ) #t=1

@constraint(model, [t=1:T, m=1:M], x[t,m] ≤ (sum(d[t_] for t_ in t:T)*y[t,m]) )

#Contrainte carbone glissante
@constraint(model, [t=1:T-R], sum( (e[m]-E)*x[t_,m] for m in 1:M, t_ in t:t+R) ≤ 0 )

optimize!(model)

groupedbar(value.(x))
