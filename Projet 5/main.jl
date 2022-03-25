using CPLEX
using JuMP
using Plots
using StatsPlots
using Random



function rolling(;
    verbose::Bool=false, 
    T::Int = 12, 
    M::Int = 4, 
    E::Int = 3, 
    f::Vector{Float64} = [10.,30.,60.,90.], 
    e::Vector{Float64} = [8.,6.,4.,2.], 
    R::Int = 2, 
    seed::Union{UInt64, Nothing} = UInt(0))


    seed = seed==nothing ? rand(UInt) : seed
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
        p = scatter(p, xs,  ys, linetype=:steppost, label = "Demande", linewidth = 3, xlabel="Time", ylabel="Production")
        savefig(p, "plot.pdf")
        display(p)
    end
    return [objective_value(model), sum(e[m]*value(x[t,m]) for m in 1:M, t in 1:T)/T]
end


function graphs(;
    T::Int = 12,
    seed::Union{UInt64, Nothing} = UInt(0))
    seed = seed==nothing ? rand(UInt64) : seed

    arguments = Dict{Symbol, Any}(:seed=>seed, :T=>T)

    results = Array{Float64}(undef, 2, T)
    
    results = cat(results, hcat([rolling(;arguments..., R=R) for R in 1:T]...), dims=3)

    arguments[:E] = 4
    results = cat(results, hcat([rolling(;arguments..., R=R) for R in 1:T]...), dims=3)

    arguments[:E] = 3
    arguments[:M] = 5
    arguments[:f] = [10.,30.,50.,60.,90.]
    arguments[:e] = [8.,6.,5.,4.,2.]
    results = cat(results, hcat([rolling(;arguments..., R=R) for R in 1:T]...), dims=3)

    arguments[:M] = 6
    arguments[:f] = [10.,30.,50.,60.,90.,95.]
    arguments[:e] = [8.,6.,5.,4.,2.,1.]
    results = cat(results, hcat([rolling(;arguments..., R=R) for R in 1:T]...), dims=3)

    arguments[:M] = 7
    arguments[:f] = [1.,10.,30.,50.,60.,90.,95.]
    arguments[:e] = [9.,8.,6.,5.,4.,2.,1.]
    results = cat(results, hcat([rolling(;arguments..., R=R) for R in 1:T]...), dims=3)


    results = results[:,:,2:end]
    p1 = plot(1:T, results[1,:,:], linewidth=2, title="Cost", labels=["Base" "E=4" "M=5" "M=6" "M=7"], xticks=1:T, xlabel="R")
    p2 = plot(1:T, results[2,:,:], linewidth=2, title="Carbon average", legend=false, xticks=1:T, xlabel="R")
    plot!(size=(80*T,round(Int, 80*T/sqrt(2))))
    p = plot(p1, p2, layout=2)
    savefig(p, "plot.pdf")
    display(p)
end
graphs(T=20)
# rolling()