using CPLEX
using JuMP
include("../utilitaire.jl")
include("MinFragmentation_opl.dat")

mutable struct Instance
	Amin::Int
	Amax::Int
	B::Int
	n::Int
	m::Int
	c
	function Instance(Amin::Int, Amax::Int, B::Int, n::Int=10, m::Int=10; alea::Bool=false)
		if alea
			c = rand(1:round(Int, sqrt(n*m)), (n,m))
		else
			@assert(n==m==10)
			c = [7 3 10 10 2 8 6 4 5 5;
				7 7 10 5 2 8 6 3 9 9;
				7 3 4 6 3 2 4 9 7 8;
				6 2 7 6 4 7 5 10 7 8;
				2 4 3 4 9 6 4 9 8 4;
				7 5 2 9 8 9 5 6 10 10;
				5 2 3 7 9 9 4 9 6 3;
				5 2 9 4 2 8 6 9 3 4;
				9 6 5 4 5 6 8 9 6 6;
				8 8 7 7 3 5 8 3 9 9]
		end
		new(Amin, Amax, B, n, m, 10*c)
	end
end
mutable struct Solution
	#model::Model
	res_time::Float64
	nb_nodes::Int
	nb_it::Int
	DMPPV::Float64
	parcelles::Matrix{Int}
end
function f(x,y)
	n,m = size(x)
	return sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i in 1:n, k in 1:n, j in 1:m, l in 1:m if (i!=k || j!=l))
end
function g(x,y)
	return sum(x)
end


function solve_PL(lambda::Float64, instance::Instance)
	model = Model(CPLEX.Optimizer)
	set_silent(model)
	m = instance.m
	n = instance.n
	c = instance.c
	@variable(model, x[1:n, 1:m], Bin )
	@variable(model, y[1:n, 1:m, 1:n, 1:m], Bin )

	@constraint(model, sum(x) <= instance.Amax)
	@constraint(model, sum(x) >= instance.Amin)
	@constraint(model, sum(x[i,j]*(c[i,j]) for i in 1:n, j in 1:m) <= instance.B)
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
		model = solve_PL(lambda, instance)
		nb_nodes += node_count(model)
		if has_values(model)
			x_lambda = value.(model[:x])
			y_lambda = value.(model[:y])
			lambda = f(x_lambda,y_lambda)/g(x_lambda,y_lambda)
		else
			println("No solution")
			return Solution(round(time()-start, digits=2), nb_nodes, nb_it, round(lambda,digits=2), zeros(Int, (n,m)) ), false
		end
			
		if abs(objective_value(model)) <= 1e-10
			println("objective = ", round(objective_value(model),digits=2))
			println("λ = ", round(lambda,digits=2))
			return Solution(round(time()-start, digits=2), nb_nodes, nb_it, round(lambda,digits=2), round.(value.(model[:x])) ), true
		end
		println("objective = ", round(objective_value(model),digits=2))
		println("λ = ", round(lambda,digits=2))
		nb_it += 1
		
	end

end

function run(instance::Instance)
	
	sol, has_val = Dinkelbach(10000., instance)
	
	if has_val
		# print(solution_summary(model, verbose=false))
		println("time = ", round(sol.res_time, digits=2))
		println("nodes : ", sol.nb_nodes)
		println("nb_it = ", sol.nb_it)
		println("DMPPV = ", round(sol.DMPPV, digits=2))
		println("nb parcelles : ", sum(sol.parcelles))
		affichage(sol.parcelles)
	else
		println("No solution")
	end
	return sol
end

function instances1_3()


	instances = Instance[Instance(30,35,920), Instance(20,21,520), Instance(70,75,3500)]
	rows = Vector{Vector{String}}(undef, length(instances))
	for i in 1:length(instances)
		inst = instances[i]
		println(inst)
		sol = run(inst)
		rows[i] = [string(inst.n) * "\$\\times\$" * string(inst.m), string(inst.Amin) * ", " * string(inst.Amax) * ", " * string(inst.B), string(sol.res_time), string(sol.nb_nodes), string(sol.nb_it), string(sol.DMPPV)]
	end
	titles = ["Taille", "\$A_{min}, A_{max}, B\$", "Temps(s)", "Noeuds", "Itérations", "DMPPV"]
	write_table_tex("instances1-3", "Résultats obtenus sur les 3 instances", titles, rows)
end

# run(Instance(20,21,520))


function instances_alea(n::Int=15)
	instances = Instance[]
	for dim in 10:n
		for _ in 1:5
			A = rand(1:100)
			Amin = round(Int, dim^2*rand(1:A)/100)
			Amax = round(Int, dim^2*rand(A:100)/100)
			B = round(Int, dim^2*20*A/100)
			push!(instances, Instance(Amin, Amax, B, dim, dim, alea=true))
		end
	end
	
	rows = Vector{Vector{String}}(undef, length(instances))
	for i in 1:length(instances)
		inst = instances[i]
		println(inst)
		sol = run(inst)
		rows[i] = [string(inst.n) * "\$\\times\$" * string(inst.m), string(inst.Amin) * ", " * string(inst.Amax) * ", " * string(inst.B), string(sol.res_time), string(sol.nb_nodes), string(sol.nb_it), string(sol.DMPPV)]
	end
	titles = ["Taille", "\$A_{min}, A_{max}, B\$", "Temps(s)", "Noeuds", "Itérations", "DMPPV"]
	write_table_tex("instancesAlea", "Résultats obtenus sur les instances générées aléatoirement", titles, rows)
end


instances_alea(16)	
