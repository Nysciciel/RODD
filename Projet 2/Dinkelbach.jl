#1. lambda = lambda_0 #on va faire decroitre lambda

#2. calculer v(lambda) = min{f(x,y)-lambda*g(x,y)}: x \in X}
model = Model(CPLEX.Optimizer)
@variable(model, x[1:n, 1:m], Bin )
@variable(model, y[1:n, 1:m, 1:n, 1:m], Bin )

@constraint(model, sum(x) <= Amax)
@constraint(model, sum(x) >= Amin)
@constraint(model, sum(map(*,x,reduce(vcat,transpose.(c)))) <= B)
@constraint(model, [i=1:n,j=1:m,k=1:n,l=1:m], y[i,j, k ,l] <= x[k,l])
@constraint(model, [i=1:n,j=1:m], sum(y[i,j,:,:]) == x[i,j] )

@objective(model, Min, sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i,k in 1:n, k in 1:n, j in 1:m, l in 1:m) )

optimize!(model)
if has_values(model)
    print(solution_summary(model, verbose=true))
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