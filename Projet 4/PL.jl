using CPLEX
using JuMP

include("ExplForet_opl.dat")

model = Model(CPLEX.Optimizer)

@variable(model, d[1:m,1:n] >= 0)
@variable(model, x[1:m,1:n], Bin)

function neighbours(i::Int, j::Int)
	if i == 1
		line = [2]
	elseif i == n
		line = [n-1]
	else
		line = [i-1, i+1]
	end
	
    if j == 1
        col = [2]
    elseif j == m
        col = [m-1]
	else
		col = [j-1, j+1]
    end
    return [(k,l) for k in line for l in col]
end

@constraint(model, [i=1:m,j=1:n], d[i,j] >= sum(x[k,l] for (k,l) in neighbours(i,j)))

@objective(model, Max, w1*sum(t[i][j]*(1-x[i,j]) for i in 1:m, j in 1:n) + w2*g*L*sum(4*x-d) )

