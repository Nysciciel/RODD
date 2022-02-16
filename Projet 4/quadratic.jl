using CPLEX
using JuMP

include("ExplForet_opl.dat")

model = Model(CPLEX.Optimizer)

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

@variable(model, x[i=1:n,j=1:m], Bin)
@variable(model, y[i=1:n,j=1:m,k=1:n,l=1:m] >=0)

#@objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:n,j in 1:m) + w2*g*L* sum( x[i,j]*(5-length(neighbours(i,j))) - y[i,j,k,l] for i in 1:n,j in 1:m for (k,l) in neighbours(i,j) ) )
#@objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:n,j in 1:m) + w2*g*L* sum( x[i,j]*(sum(1-x[k,l] for (k,l) in neighbours(i,j))+4-length(neighbours(i,j))) for i in 1:n,j in 1:m  ) )
@objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:n,j in 1:m) + w2*g*L* ( sum(x[i,j]*(4-length(neighbours(i,j))) for i in 1:n,j in 1:m) + sum(x[i,j]-y[i,j,k,l] for i in 1:n, j in 1:m for (k,l) in neighbours(i,j) )  ) )

@constraint(model, [i=1:n,j=1:m,k=1:n,l=1:m;i!=k || j!=l], y[i,j,k,l] >= x[i,j] + x[k,l] - 1 )


optimize!(model)
if has_values(model)
    print(prod((e->get(Dict(1=>"■ ",0=>"□ ",-1=>"\n"),e," ")).(Int.(round.(hcat(value.(x),(e->-1).(Matrix{Int}(undef, size(value.(x))[1],1)))')))))
end