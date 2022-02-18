#print( prod( (e->get(Dict(1=>"■ ",0=>"□ ",-1=>"\n"),e," ") ).([Int.(round.(value.(x))) .-ones(Int,size(x,1),1) ]') ) )

function affichage(x,y=zeros(Int,size(x)))
    matrice_x = [Int.(round.(x)) .-ones(Int,size(x,1),1) ]'
    dic = Dict(1=>"■ ",0=>"□ ",-1=>"\n")
    matrice_x = (e->get(dic,e," ") ).(matrice_x)
    matrice = copy(matrice_x)
    #fusion des matrices
    for i in 1:size(x,1)
        for j in 1:size(x,2)
            matrice[j,i] = y[i,j]==1 ? "x " : matrice_x[j,i]
        end
    end
    print( prod( matrice) )
end

function write_solution(alpha, x,y,cout,time,noeuds, surv_proba, instance::String)
    fout = open("res/" * instance * ".txt", "w")
    println(fout, "alpha = ", alpha)
    println(fout, "x = ", x)
    println(fout, "y = ", y)
    println(fout, "cout = ", cout)      
    println(fout, "time = ", time)
    println(fout, "nb_noeuds = ", noeuds)
    println(fout, "surv_proba =", surv_proba)
    close(fout)
end

"""
Create a latex file which contains an array with the results
"""
function results_tex()
    println("--------- RESULTS TEX  ---------")

    # Open the latex output file
    fout = open("instances_1_4.tex", "w")

    # Print the latex file output
    println(fout, raw"""\documentclass[main.tex]{subfiles}
	\begin{document}""")

    header = raw"""
	\begin{center}
	\begin{tabular}{|l||l|c|c|c|c|c|c||c|c|c|}
    \hline
    \textbf{Instance}&\textbf{Espèce}&1&2&3&4&5&6&\textbf{Temps(s)}&\textbf{Noeuds}&\textbf{Coût}\\\\"""

    footer = raw"""\hline
    \end{tabular}
	\end{center}
	"""
    println(fout, header)

    # For each file in the result folder
    for file in readdir("res")
        path = "res/"*file
        if isfile(path)
            println(path)
            include(path)

            print(fout, "\n\\hline\\hline\n", replace(file, "_" => "\\_"), " &\$\\alpha\$&")
            for α in alpha
                print(fout, α, " &")
            end
            print(fout, time, " &", nb_noeuds, " &", cout, "\\\\\n")
            print(fout, "&Probabilités de survie &")
            for k in 1:length(surv_proba)
                print(fout, surv_proba[k], " &")
            end
            println(fout, "& &\\\\")
        end
    end


    # Print the end of the latex file
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end
