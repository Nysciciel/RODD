using CPLEX
using JuMP
include("../utilitaire.jl")

mutable struct Sommet
    l::Int #âge jachère
    a::Int #âge culture
    j::Int #culture ou jachère en cours
end
""" Surcharge de l'opérateur == """
function Base.:(==)(x::Sommet, y::Sommet)
    x.l == y.l && x.a == y.a && x.j == y.j
end

mutable struct Arc
    i::Sommet #initial
    f::Sommet #final
    rend::Int #rendement
end
""" Surcharge de l'opérateur == """
function Base.:(==)(x::Arc, y::Arc)
    x.i == y.i && x.f == y.f && x.rend == y.rend
end

include("inst") #import de l'instance

function arcs_entrants(sommet::Sommet)
    return unique([k for k in 1:length(Arcs) if sommet == Arcs[k].i])
end
function arcs_sortants(sommet::Sommet)
    return unique([k for k in 1:length(Arcs) if sommet == Arcs[k].f])
end

########## Paramètres de l'instance (en plus de ceux dans le fichier inst)
P = NbParcelles::Int #nb parcelles
D = Demande
Cultures = append!(C...)
Demande = zeros(Int, (length(Cultures), T)) #Demande[Culture,t]
initArcs = [1, 2]
##############################################################

########## MODEL #############################################
model = Model(CPLEX.Optimizer)
@variable(model, x[1:length(Arcs), 1:T, 1:P], Bin)

init_flot = @constraint(model, [p in 1:P], sum(x[u, 1, p] for u in initArcs) ≤ 1)
arcs_interdits = @constraint(model, [p in 1:P], sum(x[u, 1, p] for u in 1:length(Arcs) if u ∉ initArcs) == 0)
conserv = @constraint(model, [sommet in Sommets, t = 1:T-1, p = 1:P],
    sum(x[u, t, p] for u in arcs_sortants(sommet)) == sum(x[u, t+1, p] for u in arcs_entrants(sommet)))
demande = @constraint(model, [j in Cultures, t = 1:T], sum(Arcs[u].rend * x[u, t, p] for u in 1:length(Arcs), p in 1:P if Arcs[u].f.j == j) ≥ D[j][t])

@objective(model, Min, sum(x[u, 1, p] for u in initArcs, p in 1:P))
##############################################################

########## SOLVE #############################################
optimize!(model)
##############################################################

########## SAVE THE RESULTS ##################################
row = [string(round(solve_time(model), digits=6)), string(node_count(model)), string(round(Int, sum(value(model[:x][arc, 1, p]) for p in 1:P, arc in initArcs)))]
titles = ["Temps(s)", "Noeuds", "Parcelles"]
write_table_tex("res", "Résultats pour l'instance Test 2 cultures", titles, [row])
##############################################################
