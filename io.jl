#print( prod( (e->get(Dict(1=>"■ ",0=>"□ ",-1=>"\n"),e," ") ).([Int.(round.(value.(x))) .-ones(Int,size(x,1),1) ]') ) )

function affichage(x,y=zeros(Int,size(x)))

    matrice_x = [Int.(round.(x)) .-ones(Int,size(x,1),1) ]'
    dic = Dict(1=>"■ ",0=>"□ ",-1=>"\n")
    matrice_x = (e->get(dic,e," ") ).(matrice_x)
    matrice = copy(matrice_x)
    #fusion des matrices
    for i in 1:size(x,1)
        for j in 1:size(x,2)
            matrice[j,i] = round(y[i,j])==1 ? "x " : matrice_x[j,i]
        end
    end
    print( prod( matrice) )
end

function write_solution(x,y,cout,time,noeuds, surv_proba, α::Int, taille::Int)
    filename = "res/alpha_" * string(α) * "/instance" * "_" * string(taille) * ".txt"
    fout = open(filename, "w")
    println(fout, "x = ", x)
    println(fout, "y = ", y)
    println(fout, "cout = ", cout)      
    println(fout, "resolution_time = ", time)
    println(fout, "nb_noeuds = ", noeuds)
    println(fout, "surv_proba =", surv_proba)
    close(fout)
end

"""
Create a latex file which contains an array with the results
"""
function results_tex(α::Int, taille_min::Int, taille_max::Int)
    println("--------- RESULTS TEX  ---------")

    # Open the latex output file
    filename = "res/results_alpha_" * string(α) * "_instances_" * string(taille_min) * "_" * string(taille_max)
    fout = open( filename * ".tex", "w")

    #Create the alpha result_folder
    result_folder = "res/alpha_"*string(α)
    if !isdir(result_folder)
        mkdir(result_folder)
    end

    # Print the latex file output
    println(fout, raw"""\documentclass[main.tex]{subfiles}
    
	\begin{document}""")

    header = raw"""
\begin{table}
\centering
\caption{Résultats obtenus pour $\alpha = """
alphas = [  [0.5,0.5,0.5,0.5,0.5,0.5],
            [0.9,0.9,0.9,0.5,0.5,0.5],
            [0.5,0.5,0.5,0.9,0.9,0.9],
            [0.8,0.8,0.8,0.6,0.6,0.6]]
    header *= string(alphas[α]) * "\$"
    header *= raw"""}
\begin{tabular}{|c||c|c|c|c|c|c||c|c|c|}
\hline
\textbf{Taille instance}&\multicolumn{6}{c}{\textbf{Proba de survie (par espèce)}}&\textbf{Temps(s)}&\textbf{Noeuds}&\textbf{Coût}\\
&1&2&3&4&5&6& & &\\
\hline
"""

    footer = raw"""
\hline
\end{tabular}
\end{table}
"""
    println(fout, header)

    maxInstancePerPage = 50
    id = 1
    
    for file in readdir(result_folder)
        path = result_folder * "/" * file
        if isfile(path)
            println(path)
            include(path)
            
            #Taille instance
            taille = split(file[1:end-4], "_")[2]
            print(fout, taille, "\$\\times\$ ", taille, " &")
            
            #Proba survie
            if surv_proba == zeros(length(surv_proba))
                print(fout, " - & - & - & - & - & - &")
            else     
                for k in 1:length(surv_proba)
                    print(fout, surv_proba[k], " &")
                end
            end

            #Resolution_time, nodes, value
            print(fout, resolution_time, " &", nb_noeuds, " &")
            if cout == Inf
                println(fout, "-\\\\")
            else
                println(fout, cout, "\\\\")
            end

            #If we need to start a new page
            if rem(id, maxInstancePerPage) == 0
                println(fout, footer, "\\newpage")
                println(fout, header)
            end
            id += 1
        end
    end


    # Print the end of the latex file
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end
