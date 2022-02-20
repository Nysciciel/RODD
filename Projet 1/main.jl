include("io.jl")
include("../utilitaire.jl")

using CPLEX
using JuMP
using Random



function solve_model(m,n,proba::Vector{Vector{Vector{Float64}}}, c::Vector{Vector{Int}}, alpha::Vector{Float64}, p::Int, K::Int, verbose::Bool=false)
    rare = [species <= p for species in 1:K] # 1 si espèce rare
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    x = @variable(model, x[1:n,1:m], Bin) # Le site est protégé
    y = @variable(model, y[1:n,1:m], Bin) # Le site est central

    @objective(model, Min, sum(c[i][j]*x[i,j] for i in 1:n, j in 1:m))

    probabilities = @constraint(model,  [k in 1:K], base_name = "probabilities",
        sum((rare[k] ? y[i,j] : x[i,j])*log(1-proba[k][i][j]) for i in 1:n, j in 1:m) <= log(1-alpha[k]))

    #Si un voisin de (i,j) n'est pas protégé, alors (i,j) n'est pas centrale
    voisins_1 = @constraint(model, [i in 1:n, j in 1:m, (k,l) in neighbours(i,j,n,m)], 
                base_name = "voisins_1", y[i,j] <= x[k,l])
    #Si tous les voisins de (i,j) sont protégés, alors (i,j) est centrale
    voisins_2 = @constraint(model, [i in 1:n, j in 1:m], base_name="voisins_2", 
                y[i,j] >= sum(x[k,l] for (k,l) in neighbours(i,j,n,m)) - 8)#length(neighbours(i,j))+1)
    

    # Les bords ne peuvent pas être centraux
    @constraint(model,[j in 1:m], y[1,j]==0)
    @constraint(model,[j in 1:m], y[n,j]==0)
    @constraint(model,[i in 1:n], y[i,1]==0)
    @constraint(model,[i in 1:n], y[i,m]==0)

    set_silent(model)
    if verbose
        unset_silent(model)
    end

    optimize!(model)

    resolution_time = round(solve_time(model), digits=2)
    noeuds = 0#node_count(model)

    if verbose
        println("resolution_time = ", resolution_time)
        println("Noeuds : ", noeuds )
    end

    if has_values(model)
        x = value.(x)
        y = value.(y)
        
        cout = round(Int,objective_value(model))

        #Probabilités de survie
        surv_proba = [ round(1-prod(1-proba[k][i][j]*value( rare[k] ? y[i,j] : x[i,j]) for i in 1:n, j in 1:m ),digits=2 ) for k in 1:K ]
        
        if verbose
            println("Fonction objetif : ", cout)
            affichage(x, y)

            println("aplha = \t", alpha)
            println("proba de survie : ", surv_proba)
        end
    else #On affiche les contraintes en conflit

        if verbose
            println("No solution")
            conflict_constraint_list = ConstraintRef[]
            println(compute_conflict!(model))
            for (F, S) in list_of_constraint_types(model)
                for con in all_constraints(model, F, S)
                    if MOI.get(model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
                        push!(conflict_constraint_list, con)
                        println(con)
                    end
                end
            end
        end
        x = nothing
        y = nothing
        cout = Inf
        
        surv_proba = zeros(K)
    end
    
    return x,y,cout,resolution_time,noeuds, surv_proba 
end


function instance_1_4(verbose::Bool=false)
    #Fichier de données : m, n, p, q, alpha, proba
    include("ReserveNaturelles_opl.dat")
    K = p + q #Nombre d'espèces à protéger 

    alphas = [  [0.5,0.5,0.5,0.5,0.5,0.5],
                [0.9,0.9,0.9,0.5,0.5,0.5],
                [0.5,0.5,0.5,0.9,0.9,0.9],
                [0.8,0.8,0.8,0.6,0.6,0.6]]
    for i in 1:length(alphas)
        if verbose
            println("------------------- Instance ", i, " ------------------")
        end
        x,y,cout,resolution_time,noeuds, surv_proba = solve_model(m,n,proba, c, alphas[i], p, K, verbose)
        write_solution(x,y,cout,resolution_time,noeuds,surv_proba, 1, m)
    end
end

function generation_instances(m,n,p,K::Int)
    proba =  [ [ [0.0 for _ in 1:n] for _ in 1:m] for _ in 1:K]
    for k in 1:K
    #Pour chaque espèce rare(commune), on sélectionne 5-10%(10-15%)
    # de cases sur lesquelles elles ont une proba non nulles
        borne_inf = round(Int, (k<=p ? 0.05 : 0.1)*m*n)
        borne_sup = round(Int, (k<=p ? 0.1 : 0.15)*m*n)
        for _ in 1:rand(borne_inf:borne_sup)
            proba[k][rand(1:n)][rand(1:m)] = rand(0.1:0.1:0.5)
        end
    end
    c = [ [rand(1:n) for i in 1:n] for j in 1:m]
    return proba, c
end


function comportement(alphas::Vector{Vector{Float64}}, p::Int, K::Int, verbose::Bool=false)
    for α in 1:length(alphas)
        result_folder = "res/alpha_"*string(α)
        if verbose
            println("--------------------------- ", α, " ---------------------------" )
        end
        for m in 10:50
            if verbose
                println("---------------------- ", m, " -----------------------" )
            end
            proba, c = generation_instances(m,m,p,K)
            x,y,cout,resolution_time,noeuds, surv_proba = solve_model(m,m,proba,c,alphas[α],p,K,verbose)
            write_solution(x,y,cout,resolution_time,noeuds,surv_proba, α,m)
        end
        
        filename = "res/results_alpha_" * string(α) * "_instances_10_50"
        results_tex(result_folder, filename, alphas[α])
    end
    
end

p = 3
K = 6
alphas = [  [0.5,0.5,0.5,0.5,0.5,0.5],
            [0.9,0.9,0.9,0.5,0.5,0.5],
            [0.5,0.5,0.5,0.9,0.9,0.9],
            [0.8,0.8,0.8,0.6,0.6,0.6]]

comportement(alphas, p, K, true)
