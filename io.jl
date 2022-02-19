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
function results_tex(filename::String)
    println("--------- RESULTS TEX  ---------")

    # Open the latex output file
    fout = open(filename*".tex", "w")

    # Print the latex file output
    println(fout, raw"""\documentclass[main.tex]{subfiles}
    
	\begin{document}""")

    header = raw"""
\begin{center}
\begin{tabular}{|c||r|c|c|c|c|c|c||c|c|c|}
\hline
\textbf{Taille instance}&\textbf{Espèce}&1&2&3&4&5&6&\textbf{Temps(s)}&\textbf{Noeuds}&\textbf{Coût}\\
\hline
"""

    footer = raw"""
\end{tabular}
\end{center}
"""
    println(fout, header)

    maxInstancePerPage = 6
    id = 1

    instances = ["" for i in 10:50]
    # For each file in the result folder
    for instance in 10:50
        
        println(fout, "\\hline")
        print(fout, "\\multirow{8}{*}{", instance, "\$\\times\$ ", instance, "} &")

        for i in 1:4
            path = "res/instance_"*string(instance)*"_"*string(i)*".txt"
            if isfile(path)
                println(path)
                include(path)
                # taille = split(file, "_")[2]
                # alpha  = split(file, "_")[3]
                if i > 1
                    print(fout, " &")
                end
                print(fout, " \$\\alpha\$ &")
                for α in alpha
                    print(fout, α, " &")
                end
                print(fout, "\\multirow{2}{*}{",time, "} &\\multirow{2}{*}{", nb_noeuds, "} & \\multirow{2}{*}{")
                if cout == Inf
                    println(fout, "-}\\\\")
                else
                    println(fout, cout, "}\\\\")
                end

                print(fout, " &Proba survie &")
                if surv_proba == zeros(length(surv_proba))
                    print(fout, " - & - & - & - & - & - &")
                else     
                    for k in 1:length(surv_proba)
                        print(fout, surv_proba[k], " &")
                    end
                end
            end
            println(fout, " & & &\n\\cline{2-11}")
        end
        println(fout, "\\hline")
        
        if rem(id, maxInstancePerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end
        id += 1
        
    end


    # Print the end of the latex file
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end
