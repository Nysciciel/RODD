include("../io.jl")

using CPLEX
using JuMP
using Random

function neighbours(i::Int, j::Int)
    """
    Renvoie les cases voisines de (i,j) elle même incluse
    """
    i_m = max(1, i-1)
    j_m = max(1, j-1)
    i_p = min(n, i+1)
    j_p = min(m, j+1)
    return [(k,l) for k in i_m:i_p for l in j_m:j_p]
end

function solve_model(proba, c, alpha, p, K)
    rare = [species <= p for species in 1:K] # 1 si espèce rare
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    x = @variable(model, x[1:n,1:m], Bin) # Le site est protégé
    y = @variable(model, y[1:n,1:m], Bin) # Le site est central

    @objective(model, Min, sum(c[i][j]*x[i,j] for i in 1:n, j in 1:m))

    probabilities = @constraint(model,  [k in 1:K], base_name = "probabilities",
        sum((rare[k] ? y[i,j] : x[i,j])*log(1-proba[k][i][j]) for i in 1:n, j in 1:m) <= log(1-alpha[k]))

    #Si un voisin de (i,j) n'est pas protégé, alors (i,j) n'est pas centrale
    voisins_1 = @constraint(model, [i in 1:n, j in 1:m, (k,l) in neighbours(i,j)], 
                base_name = "voisins_1", y[i,j] <= x[k,l])
    #Si tous les voisins de (i,j) sont protégés, alors (i,j) est centrale
    voisins_2 = @constraint(model, [i in 1:n, j in 1:m], base_name="voisins_2", 
                y[i,j] >= sum(x[k,l] for (k,l) in neighbours(i,j)) - 8)#length(neighbours(i,j))+1)
    

    # Les bords ne peuvent pas être centraux
    @constraint(model,[j in 1:m], y[1,j]==0)
    @constraint(model,[j in 1:m], y[n,j]==0)
    @constraint(model,[i in 1:n], y[i,1]==0)
    @constraint(model,[i in 1:n], y[i,m]==0)

    optimize!(model)
    time = round(solve_time(model), digits=2)
    println("solve_time = ", time)
    noeuds = 0#node_count(model)
    println("Noeuds : ", noeuds )

    if has_values(model)
        x = value.(x)
        y = value.(y)
        affichage(value.(x), value.(y) )
        cout = round(objective_value(model),digits=2)
        println("Fonction objetif : ", cout)
        
        #Probabilités de survie
        surv_proba = [ round(1-prod(1-proba[k][i][j]*value( rare[k] ? y[i,j] : x[i,j]) for i in 1:n, j in 1:m ),digits=2 ) for k in 1:K ]
        println("aplha = \t", alpha)
        println("proba de survie : ", surv_proba)
    else #On affiche les contraintes en conflit
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
        println("No solution")
        x = nothing
        y = nothing
        cout = Inf
        surv_proba = zeros(K)
    end
    return x,y,cout,time,noeuds, surv_proba 
end


function instance_1_4()
    #Fichier de données : m, n, p, q, alpha, proba
    include("ReserveNaturelles_opl.dat")
    K = p + q #Nombre d'espèces à protéger 

    #alpha_k = 0.5 pour tout k
    println("------------------- Instance 1 ------------------")
    x,y,cout,time,noeuds, surv_proba = solve_model(proba, c, alpha, p, K)
    write_solution(alpha,x,y,cout,time,noeuds,surv_proba, "instance_1")

    # #2ème instance
    println("------------------- Instance 2 ------------------")
    alpha = [0.9,0.9,0.9,0.5,0.5,0.5]
    x,y,cout,time,noeuds, surv_proba = solve_model(proba, c, alpha, p, K)
    write_solution(alpha,x,y,cout,time,noeuds,surv_proba, "instance_2")

    # #3ème instance
    println("------------------- Instance 3 ------------------")
    alpha = [0.5,0.5,0.5,0.9,0.9,0.9]
    x,y,cout,time,noeuds, surv_proba = solve_model(proba, c, alpha, p, K)
    write_solution(alpha,x,y,cout,time,noeuds,surv_proba, "instance_3")

    # #4ème instance
    println("------------------- Instance 4 ------------------")
    alpha = [0.8,0.8,0.8,0.6,0.6,0.6]
    x,y,cout,time,noeuds, surv_proba = solve_model(proba, c, alpha, p, K)
    write_solution(alpha,x,y,cout,time,noeuds,surv_proba, "instance_4")

    results_tex()
end

function generation_instances(m,n, K::Int)
    p = rand(1:K)
    alpha_rare = round(rand(0.5:0.1:0.9),digits=2)
    alpha_commun = round(rand(0.5:0.1:0.9),digits=2)
    alpha = [   alpha_rare*ones(p,1)' alpha_commun*ones(K-p,1)' ]
    proba =  [ [ [0.0 for i in 1:n] for j in 1:m] for k in 1:K]
    for k in 1:K
    #Pour chaque espèce rare(commune), on sélectionne 5-10%(10-15%)
    # de cases sur lesquelles elles ont une proba non nulles
        borne_inf = round(Int, (k<=p ? 0.05 : 0.1)*m*n)
        borne_sup = round(Int, (k<=p ? 0.1 : 0.15)*m*n)
        for i in 1:rand(borne_inf:borne_sup)
            proba[k][rand(1:n)][rand(1:m)] = rand(0.1:0.1:0.5)
        end
    end
    c = [ [rand(1:n) for i in 1:n] for j in 1:m]
    return proba, c, alpha, p
end

m = 10
n = 10
K = 6
# Random.seed!(0002)
proba, c, alpha, p = generation_instances(m,n,K)
x,y,cout,time,noeuds, surv_proba = solve_model(proba,c,alpha,p,K)