using CPLEX
using JuMP

include("ExplForet_opl.dat")
# include("inst2.dat")
model = Model(CPLEX.Optimizer)

@variable(model, d[1:m,1:n] >= 0)
@variable(model, x[1:m,1:n], Bin)

function neighbours(i::Int, j::Int)
    L = [(k,l) for k in [i+1,i,i-1] for l in Dict(i+1 => [j], i => [j-1,j+1], i-1 => [j])[k]]
    r = copy(L)
    for (i,j) in L
        if i==0 || i==n+1 || j==0 || j== m+1
            filter!(e->e!=(i,j),r)
        end
    end
    return r
end

@constraint(model, [i=1:m,j=1:n], d[i,j] >= sum(x[k,l] for (k,l) in neighbours(i,j)) - length(neighbours(i,j))*(1-x[i,j]) )
@constraint(model, sum(x) >= 60)

@objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:m, j in 1:n) + w2*g*L*sum(4*x-d) )

optimize!(model)
println("objective value= ", objective_value(model))
println("nombre de parcelles non coupées = ", Int(sum(value.(x))))
println("Nombre de lisières : ", Int.(sum(4*value.(x)-value.(d))))
round.(Int,value.(x))