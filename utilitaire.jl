function neighbours(i::Int, j::Int, n::Int, m::Int)
    """
    Renvoie les cases voisines de (i,j) elle même incluse
    """
    i_m = max(1, i-1)
    j_m = max(1, j-1)
    i_p = min(n, i+1)
    j_p = min(m, j+1)
    return [(k,l) for k in i_m:i_p for l in j_m:j_p]
end

function affichage(x,y=zeros(Int,size(x)))
    matrice = [Int.(round.(x)) .-ones(Int,size(x,1),1) ]'
    dic = Dict(1=>"■ ",0=>"□ ",-1=>"\n")
    matrice = (e->dic[e] ).(matrice)
    #fusion des matrices
    matrice[[CartesianIndex(j,i) for i in 1:10, j in 1:10 if round(x[i,j]*y[i,j])==1]] .= "x "
    println( prod( matrice) )
end

"""
Write a table in a .tex file
Input:
    - output = filename of the output file
    - caption = caption of the table
    - titles = header of the table
        - num_col_titles = number of cols for each title
    - subtitles = subheader of the table
        - num_col_sub = number of cols for each subtitle
"""
function write_table_tex(output::String, caption::String, titles::Array{String}, rows::Vector{Vector{String}};
    subtitles::Array{String}=String[], subsubtitles::Array{String}=String[], 
    num_col_titles::Array{Int}=ones(Int,length(titles)), num_col_sub::Array{Int}=ones(Int,length(subtitles)),
    alignment::String="c"^sum(num_col_titles), lines::Array{String}=fill("",length(rows)), maxRawsPerPage::Int=50 )

    fout = open(output * ".tex", "w")

    println(fout, raw"""\documentclass[main.tex]{subfiles}
\newmargin{2cm}{2cm}
\setlength{\voffset}{-1.5cm}
\begin{document}
\thispagestyle{empty}
""")

    #HEADER OF TABLE
    header = raw"""
\begin{table}
    \centering
    \caption{"""
    header *= caption
    header *= raw"""
}
    \begin{tabular}{
    """

    header *= alignment * "}\n\\hline\t\n\t"

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
    header *= "\\\\\n\t\\hline\n"

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
            if i < length(subtitles)
                subheader *= " &"
            end
        end
        subheader *= "\\\\\n\t"
    end

    #SUBSUBHEADERS
    subsubheader = "\n\t"
    if length(subsubtitles) > 0
        for i in 1:length(subsubtitles)
            subsubheader *= subsubtitles[i]
            
            if i < length(subsubtitles)
                subsubheader *= " &"
            end
        end
        subsubheader *= "\\\\\n\t\\hline\n"
    end

    #FOOTER OF TABLES
    footer = raw"""
    \end{tabular}
\end{table}
"""

    print(fout, header)
    println(fout, subheader)
    println(fout, subsubheader)

    id = 1

    #CONTENT
    for j in 1:length(rows)
        for i in 1:length(rows[j])
            print(fout, rows[j][i])
            if i < length(rows[j])
                print(fout, " &")
            end
        end
        
        println(fout, "\\\\" * lines[j])
        
        
        #If we need to start a new page
        if rem(id, maxRawsPerPage) == 0
            println(fout, footer, "\\newpage\n\\thispagestyle{empty}")
            println(fout, header)
            println(fout, subheader)
            println(fout, subsubheader)
        end
        id += 1
    end
    
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end