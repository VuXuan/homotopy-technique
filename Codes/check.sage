def test1(v,q,G): ### test for start matrix after lifting 
    nr = G.nrows()
    nd = G.ncols()
    
    G1 = G.substitute(x2 = v[0,0], x3 = v[1,0])
    
    G2 = Matrix(parent(q), nr,nd)
    for i in range(nr):
    	for j in range(nd):
	    G2[i,j] = G1[i,j]

    G3 = Matrix(parent(q), nr, nd)
    for i in range(nr):
    	for j in range(nd):
	    G3[i,j] = G1[i,j].mod(q)

    G4 = G3.substitute(T=0)
    return G4
