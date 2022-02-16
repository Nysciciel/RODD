using CPLEX
using JuMP

#Fichier de données : m, n, p, q, alpha, proba
include("ReserveNaturelles_opl.dat")

model = Model(CPLEX.Optimizer)

K = p + q #Nombre d'espèces à protéger 
rare = [species <= p for species in 1:K] # indices des espèces rares

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

x = @variable(model, x[1:n,1:m], Bin) # Le site est protégé
y = @variable(model, y[1:n,1:m], Bin) # Le site est central

survival_probabilities = @constraint(model,[k in 1:K],
 sum((rare[k] ? y : x)[i,j]*log(1-proba[k][i][j]) for i in 1:n, j in 1:m) <= log(1-alpha[k]))

@constraint(model, [i in 1:n, j in 1:m, (k,l) in neighbours(i,j)], y[i,j] <= x[k,l])
# Les bords ne peuvent pas être centraux
@constraint(model,[j in 1:m], y[1,j]==0)
@constraint(model,[j in 1:m], y[n,j]==0)
@constraint(model,[i in 1:n], y[i,1]==0)
@constraint(model,[i in 1:n], y[i,m]==0)

@objective(model, Min, sum(c[i][j]*x[i,j] for i in 1:n, j in 1:m))

optimize!(model)
if has_values(model)
    print(solution_summary(model, verbose=true))
    for i in 1:n
        for j in 1:m
            if round(value(y[i, j]))==1
                print("O")
            elseif round(value(x[i, j]))==1
                print("o")
            else
                print("_")
            end
        end
        println()
    end
else
    print("No solution")
end