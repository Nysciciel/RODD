using CPLEX
using JuMP

include("DivGenetique_opl.dat")

m = C*G*A #nombre d'allèles
#Nombre d'allèles k pour l'individu i
n = [ [length(findall(x->x==1+rem(k+1,A),individu[i][1][ceil(Int,k/A)])) for k in 1:m] for i in 1:N]

I_2 = [ [i for i in 1:N if n[i][k] == 2] for k in 1:m]  #I_2[k] = individus qui ont 2 fois l'allèle k 
I_1 = [ [i for i in 1:N if n[i][k] == 1] for k in 1:m]  #I_1[k] = individus qui ont 1 fois l'allèle k

h=50
θ = zeros(h)
θ[1] = 0.001
for r in 2:h
    θ[r] = θ[1]^((h-r)/(h-1))
end
model = Model(CPLEX.Optimizer)
set_silent(model)

@variable(model, 0 ≤ x[1:N] ≤ 3, Int)   #nb d'enfants de l'individu i
@variable(model, p[1:m] ≥ 0)    #Proba de disparition de l'allèle k
@variable(model, t[1:m] ≥ 0)    #

@constraint(model, sum(x[1:Nm]) == N)   #Les mâles font N enfants
@constraint(model, sum(x[Nm+1:N]) == N) #Les femelles font N enfants
@constraint(model, [k=1:m], p[k] >= t[k] - sum(x[i] for i in I_2[k]))
@constraint(model, [k in 1:m, r in 1:h], log(θ[r]) + (t[k] - θ[r])/θ[r] >= log(0.5)*sum(x[i] for i in I_1[k])) #Approx

@objective(model, Min, sum(p))

optimize!(model)
if has_values(model)
    # print(solution_summary(model, verbose=true))
    #X
    for i in 1:N
        print(name(x[i])*"\t")
        if i == Nm
           print(" |\t")
       end
    end
    println("\n")#,"-"^(8*N+3))
    for i in 1:N
        print(round(Int,value(x[i])),"\t")
        if i == Nm
           print(" |\t")
       end
    end
    println("\nNombre d'enfants = ", round(Int,value(sum(x))/2), "\n")
    #PROBAS
    for k in 1:m
        print(name(p[k])*"\t")
        if rem(k,2) == 0
           print("\t")
       end
    end
    println("\n")#,"-"^(8*N+3))
    for k in 1:m
        print(round(value.(p[k]), digits = 6), "\t")
        if rem(k,2) == 0
           print("\t")
       end
    end
    println("\nBorne inférieure du problème : ", round(objective_value(model), digits=6), "\n")

    #PROBAS REELLES
    for k in 1:m
        print(name(p[k])*"\t")
        if rem(k,2) == 0
            print("\t")
        end
    end
    println("\n")#,"-"^(8*N+3))
    probas = Float64[]
    for k in 1:m
#         p = prod( ((i ∈ I_2[k]) && (value(x[i])>0) ) ? 0 : (i in I_1[k] ? (1/(2^value(x[i])) : 1)) for i in 1:N)
        proba = 1
        for i in 1:N
            if (i ∈ I_2[k]) && (value(x[i])>0)
                proba = 0
            else
                proba *= 1/(2^value(x[i]))
            end
        end
        push!(probas, proba)
        print( round( probas[k]  , digits = 6), "\t")
        if rem(k,2) == 0
            print("\t")
        end
    end
    println("\nEspérance du nombre d'allèles perdus : ", sum(probas))

else
    print("No solution")
end