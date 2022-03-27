using CPLEX
using JuMP
include("../utilitaire.jl")
include("instance.jl")


function runPL(inst::Instance)
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    m = inst.m
    n = inst.n
    w1 = inst.w1
    w2 = inst.w2
    L = inst.L
    g = inst.g
    t = inst.t
    min_x = inst.min_x


    @variable(model, d[1:m,1:n] >= 0)
    @variable(model, x[1:m,1:n], Bin)

    function neighbours(i::Int, j::Int)
        L_ = [(k,l) for k in [i+1,i,i-1] for l in Dict(i+1 => [j], i => [j-1,j+1], i-1 => [j])[k]]
        r = copy(L_)
        for (i,j) in L_
            if i==0 || i==n+1 || j==0 || j== m+1
                filter!(e->e!=(i,j),r)
            end
        end
        return r
    end

    @constraint(model, [i=1:m,j=1:n], d[i,j] >= sum(x[k,l] for (k,l) in neighbours(i,j)) - length(neighbours(i,j))*(1-x[i,j]) )
    @constraint(model, sum(x) >= min_x)
    @objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:m, j in 1:n) + w2*g*L*sum(4*x-d))

    optimize!(model)
    return Solution(round(solve_time(model), digits = 6), node_count(model), value.(x))
end
# if has_values(model)
#     print(prod((e->get(Dict(1=>"■ ",0=>"□ ",-1=>"\n"),e," ")).(Int.(round.(hcat(value.(x),(e->-1).(Matrix{Int}(undef, size(value.(x))[1],1)))')))))
# end