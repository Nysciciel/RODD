include("MinFragmentation_opl.dat")

#2. calculer v(lambda) = min{f(x,y)-lambda*g(x,y)}: x \in X}
function f(x,y)
	return sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i in 1:n, k in 1:n, j in 1:m, l in 1:m if (i!=k || j!=l))
end
function g(x,y)
	return sum(x)
end

function solve_PL(lambda)
	model = Model(CPLEX.Optimizer)
	@variable(model, x[1:n, 1:m], Bin )
	@variable(model, y[1:n, 1:m, 1:n, 1:m], Bin )

	@constraint(model, sum(x) <= Amax)
	@constraint(model, sum(x) >= Amin)
	@constraint(model, sum(map(*,x,reduce(vcat,transpose.(10*c)))) <= B)
	@constraint(model, [i=1:n,j=1:m,k=1:n,l=1:m; i!=k || j!=l], y[i,j, k ,l] <= x[k,l])
	@constraint(model, [i=1:n,j=1:m], sum(y[i,j,k,l] for k in 1:n, l in 1:m if (i!=k || j!=l)) == x[i,j] )
	#@constraint(model, [i=1:n,j=1:m], y[i,j,i,j] == 0)

	@objective(model, Min, f(x,y) - lambda*g(x,y) )

	optimize!(model)
	
	println("Objectif : ", v_lambda)
	return model
end

function Dinkelbach(lambda_0)
	lambda = lambda_0
	
	while true
		model = solve_PL(lambda)
		if objective_value(model) == 0
			return model
		end
		lambda = f(x_lambda,y_lambda)
		println("lambda = ", lambda)
	end
end

model = Dinkelbach(50)

if has_values(model)
    print(solution_summary(model, verbose=false))
    for i in 1:n
        for j in 1:m
            if round(value(x[i, j]))==1
                print("O")
            else
                print("_")
            end
        end
        println()
    end
else
    print("No solution")
end
