using CPLEX
using JuMP
using Random
using Distributions
include("../utilitaire.jl")

struct Instance
    N::Int
    C::Int
    G::Int
    A::Int
    h::Int
    B::Int
    init::Float64
    individu
    α::Float64
    function Instance(B::Int=8, α::Float64=0.,N::Int=8, h::Int=10, G::Int=5, alea::Bool=true)
        Random.seed!(0)
        if !alea
            individu =  [[[2, 1], [1, 1], [2, 1], [2, 1], [1, 2]]],
                        [[[2, 2], [2, 1], [1, 2], [1, 2], [1, 1]]],
                        [[[2, 2], [2, 1], [2, 1], [1, 2], [1, 2]]],
                        [[[2, 2], [1, 1], [1, 2], [2, 2], [2, 2]]],
                        [[[2, 1], [1, 1], [1, 1], [2, 2], [1, 2]]],
                        [[[1, 2], [1, 1], [2, 2], [2, 1], [2, 1]]],
                        [[[1, 1], [1, 1], [2, 2], [1, 1], [1, 1]]],
                        [[[2, 2], [1, 1], [2, 2], [2, 1], [2, 1]]]
        else
            #α entre 0 et 1: traduit la rareté de l'allèle 1 des gènes. La rareté diminue avec le nombre de colonnes (ie le nombre de gènes)
            function p(g)
                m = 1 * G * 2 #C*G*A
                Bernoulli( (g+α*m)/(2*m*(1+α)) )
            end
            individu = [[[[1+rand( p(g) ), 1+rand( p(g) )] for g in 1:G]] for _ in 1:N]
        end
        return new(N, 1, G, 2, h, B, 0.001, individu, α)
    end
end

struct Solution
    time::Float64
    nodes::Int
    probs
    expec::Float64
    LB::Float64
    nb_parents::Int
end

""" Affichage de la Solution du model """
function affichage(model, instance, probas)
    x = model[:x]
    p = model[:p]
    N = instance.N
    Nm = round(N/2)
    A = instance.A
    m = instance.C * instance.G * A

    if has_values(model)
        #X
        for i in 1:N
            print(name(x[i])*"\t")
            if i == Nm
                print(" |\t")
            end
        end
        println("\n")
        for i in 1:N
            print(round(Int,value(x[i])),"\t")
            if i == Nm
                print(" |\t")
            end
        end
        println()
        #PROBAS
        for k in 1:m
            print(name(p[k])*"\t")
            if rem(k,2) == 0
                print("\t")
            end
        end
        println("\n")
        for k in 1:m
            print(round(value.(p[k]), digits = 5), "\t")
            if rem(k,2) == 0
                print("\t")
            end
        end
        println("\nBorne inférieure du problème : ", round(objective_value(model), digits=5), "\n")

        #PROBAS REELLES
        for k in 1:m
            print(name(p[k])*"\t")
            if rem(k,2) == 0
                print("\t")
            end
        end
        println("\n")
        for k in 1:m
            print( round( probas[k]  , digits = 5), "\t")
            if rem(k,2) == 0
                print("\t")
            end
        end
        println("\nEspérance du nombre d'allèles perdus : ", sum(probas))
    else
        println("No Solution")
    end
end

function solve_PL(inst::Instance, verbose::Bool=false)
    start = time()
    N = inst.N
    C = inst.C
    G = inst.G
    A = inst.A
    B = inst.B
    individu = inst.individu

    Nm = round(Int, N/2)
    m = C*G*A #nombre d'allèles
    #Nombre d'allèles k pour l'individu i
    n = [ [length(findall(x->x==1+rem(k+1,A),individu[i][1][ceil(Int,k/A)])) for k in 1:m] for i in 1:N]

    I_2 = [ [i for i in 1:N if n[i][k] == 2] for k in 1:m]  #I_2[k] = individus qui ont 2 fois l'allèle k 
    I_1 = [ [i for i in 1:N if n[i][k] == 1] for k in 1:m]  #I_1[k] = individus qui ont 1 fois l'allèle k

    #THETA
    h = inst.h
    θ = zeros(h)
    θ[1] = inst.init
    for r in 2:h
        θ[r] = θ[1]^((h-r)/(h-1))
    end

    #MODEL
    model = Model(CPLEX.Optimizer)
    set_silent(model)

    @variable(model, 0 ≤ x[1:N] ≤ B, Int)   #nb d'enfants de l'individu i
    @variable(model, p[1:m] ≥ 0)    #Proba de disparition de l'allèle k
    @variable(model, t[1:m] ≥ 0)    #

    @constraint(model, sum(x[1:Nm]) == N)   #Les mâles font N enfants
    @constraint(model, sum(x[Nm+1:N]) == N) #Les femelles font N enfants
    @constraint(model, [k=1:m], p[k] >= t[k] - sum(x[i] for i in I_2[k]))
    @constraint(model, [k in 1:m, r in 1:h], log(θ[r]) + (t[k] - θ[r])/θ[r] >= log(0.5)*sum(x[i] for i in I_1[k])) #Approx

    @objective(model, Min, sum(p))

    optimize!(model)

    #Calcul des probas de disparition sur la solution optimale
    probas = Float64[]
    for k in 1:m
        proba = 1
        for i in 1:N
            if (i ∈ I_2[k]) && (value(x[i])>0)
                proba = 0
                break
            elseif i ∈ I_1[k]
                proba *= 1/(2^value(x[i]))
            end
        end
        push!(probas, proba)
    end
    if verbose
        affichage(model, inst, probas)
    end

    return Solution(round(time()-start, digits=2), node_count(model), probas, round(sum(probas), digits=6), round(objective_value(model), digits=6), sum(value.(x).≥1))
end

function instances1_2()
	instances = Instance[Instance(3), Instance(2)]
	
	rows = Vector{Vector{String}}(undef, length(instances))
	for i in 1:length(instances)
		inst = instances[i]
		println(inst)
		sol = solve_PL(inst)
		rows[i] = [string(inst.N), string(inst.h) , string(sol.time), string(sol.nodes), string(sol.LB), string(sol.expec)]
	end
	titles = ["Taille", "Morceaux", "Temps(s)", "Noeuds", "Borne inférieure", "Espérance"]
	write_table_tex("instances1-2", "Résultats obtenus sur les instances générées aléatoirement", titles, rows)
end
@enum Param alpha hh taille_pop genes
function instances_alea(param, verbose::Bool=false)
	instances_1 = Instance[]
    h = 9
    N = 8
    α = 0.5
    G = 500
    # d = Dict(:α => 0:0.1:1, :h=>1:10, :N=>8:2:20)

    if param == alpha
        liste = 0:0.1:1
        for par in liste
            push!(instances_1, Instance(N, par, N, h, 5, true))
        end
    elseif param == hh
        liste = 1:10
        for par in liste
            push!(instances_1, Instance(N, α, N, par, 5, true))
        end
    elseif param == taille_pop
        liste = 8:5:100
        for par in liste
            push!(instances_1, Instance(N, α, par, h, G, true))
        end
    elseif param == genes
        liste = 10:10:60
        for par in liste
            push!(instances_1, Instance(N, α, N, h, par, true))
        end
    else
        error("Le paramètre n'est pas reconnu")
    end
    rows_1 = Vector{Vector{String}}(undef, length(instances_1))
    solve_PL(instances_1[1])
	for i in 1:length(instances_1)
		global inst = instances_1[i]
        if verbose
    		println(inst)
        end
        sol = solve_PL(inst, verbose)
		rows_1[i] = [string(inst.N), string(inst.G), string(inst.h), string(inst.α), string(sol.time), string(sol.nodes), string(sol.LB), string(sol.expec), string(sol.nb_parents)]
	end
	titles = ["Taille", "Gènes", "Morceaux", "\$\\alpha\$", "Temps(s)", "Noeuds", "Borne inférieure", "Allèles disparus", "Parents"]
	write_table_tex("instances_alea", "Résultats obtenus pour ", titles, rows_1)


    
end

instances_alea(taille_pop)