struct Instance
    P
    m
    n
    w1
    w2
    L
    g
    t
    min_x
end


struct Solution
    time::Float64
    nodes::Int
    x
end

function gen_random(sizes, bound)
    instances = []
    P=1
    L=3
    g=1.26157
    for m in sizes
        n=m

        w1=rand()
        w2=rand()
        t= [[rand(round(Int,m*n*0.7):m*n) for _ in 1:n] for _ in 1:m]
        inst = Instance(P,m,n,w1,w2,L,g,t,bound)

        push!(instances, inst)

    end
    return instances
end