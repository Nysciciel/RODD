using GLPK
using JuMP

model = Model(GLPK.Optimizer)

N = 10
K = 6

rare = [species <= 3 for species in 1:K]
alpha = [rare[species] ? 0.7 : 0.5 for species in 1:K]
p = [[0.7 for i in 1:N] for k in 1:K]
c = [1 for i in 1:N]

function neighbours(site::Int)
    return [site-1, site, site+1]
end

@variable(model, x[1:N], Bin) # Le site est protégé
@variable(model, y[1:N], Bin) # Le site est central

@constraint(model,[k in 1:K], sum((rare[k] ? y : x)[i]*log(1-p[k][i]) for i in 1:N) <= log(1-alpha[k]))

@constraint(model, [i in 1:N, j in neighbours(i)], y[i] <= x[j])

@objective(model, Min, sum(c[i]*x[i] for i in 1:N))

optimize!(model)
if has_values(model)
    print(solution_summary(model, verbose=true))
else
    print("No solution")
end