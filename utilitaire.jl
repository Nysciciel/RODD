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