#=
Main:
- Julia version: 
- Author: Vincent
- Date: 2022-02-09
=#
using CPLEX
using JuMP

include("DivGenetique_opl.dat")

model = Model(CPLEX.Optimizer)
m = C*G*A #nombre d'allèles
@variable(model, 3>=x[1:N]>=0, Int)
@variable(model, 0<= p[1:m])
@variable(model, t[1:m])

@constraint(model, sum(x[1:Nm]) == N)
@constraint(model, sum(x[Nm+1:N]) == N)
#Nombre d'allèles k pour l'individu i
n = [ [length(findall(x->x==1+rem(k+1,A),individu[i][1][ceil(Int,k/A)])) for k in 1:m] for i in 1:N]
I_2 = [ [i for i in 1:N if n[i][k] == 2] for k in 1:m]
@constraint(model, [k=1:m], p[k] >= t[k] - sum(x[i] for i in I_2[k]))
I_1 = [ [i for i in 1:N if n[i][k] == 1] for k in 1:m]
#Approx
h=100
θ = Array(1/h:1/h:1)
@constraint(model, [k in 1:m, θr in θ], log(θr) + (t[k] - θr)/θr >= sum(x[i]*log(0.5) for i in I_1[k]))



@objective(model, Min, sum(p))

optimize!(model)
if has_values(model)
    print(solution_summary(model, verbose=true))
    #X
    for i in 1:N
        print(name(x[i])*"\t")
        if i == Nm
           print(" |\t")
       end
    end
    println("\n","-"^(8*N+3))
    for i in 1:N
        print(round(Int,value(x[i])),"\t")
        if i == Nm
           print(" |\t")
       end
    end
    println("\nNombre d'enfants = ", round(Int,value(sum(x))/2))
    #PROBAS
    for k in 1:m
        print(name(p[k])*"\t")
        if rem(k,2) == 0
           print(" |\t")
       end
    end
    println("\n","-"^(8*N+3))
    for k in 1:m
        print(round(value.(p[k]), digits = 5), "\t")
        if rem(k,2) == 0
           print(" |\t")
       end
    end
    println("\nBorne inférieure du problème : ", objective_value(model))

    #PROBAS REELLES
    for k in 1:m
        print(name(p[k])*"\t")
        if rem(k,2) == 0
           print(" |\t")
       end
    end
    println("\n","-"^(8*N+3))
    for k in 1:m
#         p = prod( ((i ∈ I_2[k]) && (value(x[i])>0) ) ? 0 : (i in I_1[k] ? (1/(2^value(x[i])) : 1)) for i in 1:N)
        p = 1
        for i in 1:N
            if (i ∈ I_2[k]) && (value(x[i])>0)
                p = 0
                break
            p *= 1/(2^value(x[i]))
            end
        end
        print( round( p  , digits = 10), "\t")
        if rem(k,2) == 0
           print(" |\t")
       end
    end

else
    print("No solution")
end