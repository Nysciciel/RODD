#print( prod( (e->get(Dict(1=>"■ ",0=>"□ ",-1=>"\n"),e," ") ).([Int.(round.(value.(x))) .-ones(Int,size(x,1),1) ]') ) )

function affichage(x,y=zeros(Int,size(x)))
    matrice = [Int.(round.(x)) .-ones(Int,size(x,1),1) ]'
    dic = Dict(1=>"■ ",0=>"□ ",-1=>"\n")
    matrice = (e->dic[e] ).(matrice)
    #fusion des matrices
    matrice[[CartesianIndex(j,i) for i in 1:10, j in 1:10 if round(x[i,j]*y[i,j])==1]] .= "x "
    print( prod( matrice) )
end

function write_solution(x,y,cout,time,noeuds, surv_proba, α::Int, taille::Int)
    result_folder = "res/alpha_" * string(α)
    #Create the alpha result_folder
    if !isdir(result_folder)
        mkdir(result_folder)
    end
    filename = result_folder * "/instance" * "_" * string(taille) * ".txt"
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
function results_tex(result_folder::String, filename::String, alpha::Vector{Float64})
    println("--------- RESULTS TEX  ---------")

    raws = Vector{String}[]
    for file in readdir(result_folder)
        path = result_folder * "/" * file
        if isfile(path)
            println(path)
            include(path)
            
            raw = []

            #Taille instance
            taille = split(file[1:end-4], "_")[2]
            push!(raw, string(taille) * "\$\\times\$ " * string(taille))
            
            #Proba survie
            if surv_proba == zeros(length(surv_proba))
                push!(raw," - ", " - "," - "," - "," - "," - ")
            else     
                for k in 1:length(surv_proba)
                    push!(raw, string(surv_proba[k]) )
                end
            end

            #Resolution_time, nodes, value
            push!(raw, string(resolution_time), string(nb_noeuds) )
            if cout == Inf
                push!(raw, " - ")
            else
                push!(raw, string(cout) )
            end
            push!(raws, raw)
        end
    end

    caption = "Résultats obtenus pour \$\\alpha = " * string(alpha) * "\$"
    titres = ["Taille instance", "Proba de survie (par espèce)", "Temps(s)", "Noeuds", "Coût"]
    subtitles = ["", "1", "2", "3", "4", "5", "6", "", "", ""]
    num_col_titles = [1,6,1,1,1]
    write_table_tex(filename, caption, titres, raws, subtitles, num_col_titles)
end



function write_table_tex(output::String, caption::String, titles::Array{String}, raws::Vector{Vector{String}},
    subtitles::Array{String}=String[], num_col_titles::Array{Int}=ones(Int,length(titles)), num_col_sub::Array{Int}=ones(Int,length(subtitles)) )
    fout = open(output * ".tex", "w")

    println(fout, raw"""\documentclass[main.tex]{subfiles}

\begin{document}
""")

    #HEADER OF TABLE
    header = raw"""
\begin{table}
    \centering
    \caption{"""
    header *= caption
    header *= raw"""
}
    \begin{tabular}{|"""

    header *= "c|"^sum(num_col_titles) * "}\n\t\\hline\n\t"

    for i in 1:length(titles)
        if num_col_titles[i] > 1
            header *= "\\multicolumn{" * string(num_col_titles[i]) * "}{c}{"
        end
        header *= "\\textbf{" * titles[i] * "}"
        if num_col_titles[i] > 1
            header *= "}"
        end
        if i < length(titles)
            header *= " &"
        end
    end
    header *= "\\\\\n\t\\hline"

    #SUBHEADERS
    subheader = "\n\t"
    if length(subtitles) > 0
        for i in 1:length(subtitles)
            if num_col_sub[i] > 1
                subheader *= "\\multicolumn{" * string(num_col_sub[i]) * "}{c}{"
            end
            subheader *= subtitles[i]
            if num_col_sub[i] > 1
                subheader *= "}"
            end
            if i < sum(num_col_titles)
                subheader *= " &"
            end
        end
        subheader *= "\\\\\n\t\\hline\n"
    end

    #FOOTER OF TABLES
    footer = raw"""
    \hline
    \end{tabular}
\end{table}
"""

    print(fout, header)
    println(fout, subheader)


    maxRawsPerPage = 50
    id = 1

    #CONTENT
    for raw in raws
        for i in 1:length(raw)
            print(fout, raw[i])
            if i < length(raw)
                print(fout, " &")
            end
        end
        println(fout, "\\\\")

        
        #If we need to start a new page
        if rem(id, maxRawsPerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end
        id += 1
    end
    
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end
