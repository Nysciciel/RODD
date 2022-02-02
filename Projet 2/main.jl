#cd("D:\\M2\\RODD\\RODD\\1")
using CPLEX
using JuMP

model = Model(CPLEX.Optimizer)

io = open("temp", "w")
write(io, replace(read(open("MinFragmentation_opl.dat"), String), ";"=>""))
close(io)
include("temp")
rm("temp")

# set_time_limit_sec(model, 30)

@variable(model, x[1:n,1:m], Bin)
@variable(model, y[1:n,1:m])

@constraint(model, sum(x) <= Amax)
@constraint(model, sum(x) >= Amin)

@constraint(model, sum(map(*,x,reduce(vcat,transpose.(c)))) <= B)

D = sqrt(m^2 + n^2)
@constraint(model, [i=1:n,j=1:m], y[i,j] <= D)
@constraint(model, [i=1:n,j=1:m], y[i,j] <= D*x[i,j])
@constraint(model, [i=1:n,j=1:m,k=1:n,l=1:m;i!=k || j!=l],
    y[i,j] <= sqrt((i-k)^2+(j-l)^2)*x[k,l] + (1-x[k,l])*D)

@objective(model, Max, sum(y))

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