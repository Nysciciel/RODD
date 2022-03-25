using CPLEX
using JuMP
include("../utilitaire.jl")
include("MinFragmentation_opl.dat")

mutable struct Instance
	Amin::Int
	Amax::Int
	B::Int
end
mutable struct Solution
	model::Model
	res_time::Float64
	nb_nodes::Int
	nb_it::Int
	DMPPV::Float64
	parcelles::Matrix{Int}
end
function f(x,y)
	return sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i in 1:n, k in 1:n, j in 1:m, l in 1:m if (i!=k || j!=l))
end
function g(x,y)
	return sum(x)
end


function solve_PL(lambda::Float64, instance::Instance)
	model = Model(CPLEX.Optimizer)
	set_silent(model)
	@variable(model, x[1:n, 1:m], Bin )
	@variable(model, y[1:n, 1:m, 1:n, 1:m], Bin )

	@constraint(model, sum(x) <= instance.Amax)
	@constraint(model, sum(x) >= instance.Amin)
	@constraint(model, sum(x[i,j]*(10*c[i][j]) for i in 1:n, j in 1:m) <= instance.B)
	@constraint(model, [i=1:n,j=1:m,k=1:n,l=1:m; i!=k || j!=l], y[i,j, k ,l] <= x[k,l])
	@constraint(model, [i=1:n,j=1:m], sum(y[i,j,k,l] for k in 1:n, l in 1:m if (i!=k || j!=l)) == x[i,j] )

	@objective(model, Min, f(x,y) - lambda*g(x,y) )

	optimize!(model)
	return model
end

function Dinkelbach(lambda_0::Float64, instance::Instance)
	start = time()
	lambda = lambda_0
	nb_it = 0
	nb_nodes = 0
	while true
		nb_it += 1
		model = solve_PL(lambda, instance)
		nb_nodes += node_count(model)
		if has_values(model)
			x_lambda = value.(model[:x])
			y_lambda = value.(model[:y])
			lambda = f(x_lambda,y_lambda)/g(x_lambda,y_lambda)
			println("objective = ", round(objective_value(model),digits=2))
			println("λ = ", round(lambda,digits=2))
		else
			println("No solution")
			return Solution(model, round(time()-start, digits=2), nb_nodes, nb_it, round(lambda,digits=2), zeros(Int, (n,m)) )
		end
			
		if objective_value(model) >= -1e-10
			return Solution(model, round(time()-start, digits=2), nb_nodes, nb_it, round(lambda,digits=2), round.(value.(model[:x])) )
		end
		
	end

end

function run(instance::Instance)
	
	sol = Dinkelbach(50., instance)
	
	if has_values(sol.model)
		# print(solution_summary(model, verbose=false))
		println("time = ", round(sol.res_time, digits=2))
		println("nodes : ", sol.nb_nodes)
		println("nb_it = ", sol.nb_it)
		println("DMPPV = ", round(sol.DMPPV, digits=2))
		println("nb parcelles : ", sum(sol.parcelles))
		affichage(sol.parcelles)
	else
		print("No solution")
	end
	return sol
end

instances = Instance[Instance(30,35,920), Instance(20,21,520), Instance(70,75,3500)]
rows = Vector{Vector{String}}(undef, length(instances))
for i in 1:length(instances)
	inst = instances[i]
	println(inst)
	sol = run(inst)
	rows[i] = [string(inst.Amin) * ", " * string(inst.Amax) * ", " * string(inst.B), string(sol.res_time), string(sol.nb_nodes), string(sol.nb_it), string(sol.DMPPV)]
end
titles = ["\$A_{min}, A_{max}, B\$", "Temps(s)", "Noeuds", "Itérations", "DMPPV"]
write_table_tex("instances1-3", "Résultats obtenus sur les 3 instances", titles, rows)