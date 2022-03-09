using CPLEX
using JuMP
using Plots
using StatsPlots
using Random



function rolling(;
    returns::String="cost", 
    verbose::Bool=false, 
    T::Int = 12, 
    M::Int = 4, 
    E::Int = 3, 
    f::Vector{Int64} = [10,30,60,90], 
    e::Vector{Int64} = [8,6,4,2], 
    R::Int = 2, 
    seed::Int = 0)


    Random.seed!(seed)
    d = [rand(20:70) for _ in 1:T]

    model = Model(CPLEX.Optimizer)


    set_optimizer_attribute(model, MOI.Silent(), !verbose)

    @variable(model, x[1:T, 1:M] ≥ 0)
    @variable(model, y[1:T, 1:M], Bin)
    @variable(model, s[1:T] ≥ 0)

    @objective(model, Min, sum(f[m]*y[t,m] for m in 1:M, t in 1:T) + sum(s))


    @constraint(model, [t=2:T], sum(x[t,m] for m in 1:M) - s[t] + s[t-1] == d[t])
    @constraint(model, sum(x[1,m] for m in 1:M) - s[1]  == d[1] ) #t=1

    @constraint(model, [t=1:T, m=1:M], x[t,m] ≤ (sum(d[t_] for t_ in t:T)*y[t,m]) )

    #Contrainte carbone glissante
    @constraint(model, [t=1:T-R+1], sum( (e[m]-E)*x[t_,m] for m in 1:M, t_ in t:t+R-1) ≤ 0 )

    optimize!(model)


    if verbose
        modes = ["Mode "*string(m) for m in 1:M]

        p = groupedbar( value.(x) , bar_position=:stack, label = reshape(modes,(1,M)) )

        ys=[]
        xs=[]
        list_of_xs = 0.5:1:T+0.5
        for i in 1:length(d)
            append!(ys,d[i])
            append!(ys,d[i])
            append!(ys,NaN)
            append!(xs,list_of_xs[i])
            append!(xs,list_of_xs[i+1])
            append!(xs,NaN)

        end

        scatter(p, xs,  ys, linetype=:steppost, label = "Demande", linewidth = 3)
    else

        if returns == "cost"
            return objective_value(model)
        elseif returns == "carbon"
            return sum(e[m]*value(x[t,m]) for m in 1:M, t in 1:T)/T
        else
            return
        end
    end

end


p1 = plot(1:12, [rolling(R=R) for R in 1:12], label="Cost")
p2 = plot(1:12, [rolling(R=R, returns="carbon") for R in 1:12], label="Carbon average")

plot(p1,p2,layout=2)