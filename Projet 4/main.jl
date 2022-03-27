include("quadratic.jl")
include("PL.jl")

include("ExplForet_opl.dat")
inst1 = Instance(P,m,n,w1,w2,L,g,t,0)
inst1_bounded = Instance(P,m,n,w1,w2,L,g,t,60)

include("inst2.dat")
inst2 = Instance(P,m,n,w1,w2,L,g,t,0)


instances = Instance[inst1, inst2]
instances_bounded = Instance[inst1_bounded]
instances_alea_unbounded = gen_random(10:5:50, 0)
instances_alea_bounded = gen_random(10:5:50, 60)


function test_instances(instances, filename::String, label::String)
	rows = Vector{Vector{String}}(undef, length(instances))
	for i in 1:length(instances)
		inst = instances[i]
		solPL = runPL(inst)
		solQr = runQuad(inst)
		rows[i] = [string(inst.m)*"\$\\times\$"*string(inst.n), string(solQr.time), string(solQr.nodes), string(solPL.time), string(solPL.nodes)]
	end
	titles = ["Taille", "Quadratique", "PLNE"]
    subtitles = ["", "Temps(s)", "Noeuds", "Temps(s)", "Noeuds"]
	write_table_tex(filename, label, titles, rows, num_col_titles = [1,2,2], subtitles = subtitles)
end

test_instances(instances, "Instances", "Résultats obtenus sur les instances données")
test_instances(instances_bounded, "Instances_bounded", "Résultats obtenus avec contrainte")
test_instances(instances_alea_unbounded, "random_instances", "Résultats obtenus sur les instances aléatoires sans contrainte")
test_instances(instances_alea_bounded, "random_instances_unbounded", "Résultats obtenus sur les instances aléatoires avec contrainte")