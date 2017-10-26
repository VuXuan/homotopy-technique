
#### Check Chinese Remainder Theorem  ####


def chekcrv(w1,qq1,w12, prec):
    ww1 = w1.mod(qq1)
    cww1 = ww1.coefficients(sparse = False)

    cw12 = w12.coefficients(sparse = False)

    dv1 = []
    for i in range(len(cww1)):
	dv1.append(Matrix(k,1,1))
    for i in range(len(cww1)):
	dv1[i][0,0] = cww1[i]

    dc = []
    for i in range(len(cw12)):
	dc.append(Matrix(k,1,1))

    for i in range(len(cw12)):
	dc[i][0,0] = cw12[i]

    dvv = []
    for i in range(len(dc)):
	dvv.append(dv1[i][0,0] - dc[i][0,0])

    dv = [dvv[i].mod(T^prec) for i in range(len(dvv))]

    return dv


def chekcrq(q1,q2, q3,q12,prec):
    qq = q1*q2*q3
    cqq = qq.coefficients(sparse = False)
    
    cq12 = q12.coefficients(sparse = False)
    
    dq = []
    for i in range(len(cqq)):
	dq.append(Matrix(k,1,1))
    for i in range(len(cqq)):
	dq[i][0,0] = cqq[i]

    dcq = []
    for i in range(len(cq12)):
	dcq.append(Matrix(k,1,1))
    for i in range(len(cq12)):
	dcq[i][0,0] = cq12[i]

    dv = []
    for i in range(len(dcq)):
	dv.append(dq[i][0,0] - dcq[i][0,0])

    d = []
    for i in range(len(dv)):
	d.append(dv[i].mod(T^prec))  

    return d

#######################################################################


#### Check Chinese Remainder Theorem  ####


def chekcrv(w1,qq1,w12, prec):
    ww1 = w1.mod(qq1)
    cww1 = ww1.coefficients(sparse = False)

    cw12 = w12.coefficients(sparse = False)

    dv1 = []
    for i in range(len(cww1)):
	dv1.append(Matrix(k,1,1))
    for i in range(len(cww1)):
	dv1[i][0,0] = cww1[i]

    dc = []
    for i in range(len(cw12)):
	dc.append(Matrix(k,1,1))

    for i in range(len(cw12)):
	dc[i][0,0] = cw12[i]

    dvv = []
    for i in range(len(dc)):
	dvv.append(dv1[i][0,0] - dc[i][0,0])

    dv = [dvv[i].mod(T^prec) for i in range(len(dvv))]

    return dv


def chekcrq(q1,q2, q3,q12,prec):
    qq = q1*q2*q3
    cqq = qq.coefficients(sparse = False)
    
    cq12 = q12.coefficients(sparse = False)
    
    dq = []
    for i in range(len(cqq)):
	dq.append(Matrix(k,1,1))
    for i in range(len(cqq)):
	dq[i][0,0] = cqq[i]

    dcq = []
    for i in range(len(cq12)):
	dcq.append(Matrix(k,1,1))
    for i in range(len(cq12)):
	dcq[i][0,0] = cq12[i]

    dv = []
    for i in range(len(dcq)):
	dv.append(dq[i][0,0] - dcq[i][0,0])

    d = []
    for i in range(len(dv)):
	d.append(dv[i].mod(T^prec))  

    return d

#######################################################################



def teststart(v,w,q,G): ### test for start matrix after lifting 
    nr = G.nrows()
    nd = G.ncols()
    
    G1 = G.substitute(x2 = v, x3 = w)
    
    G2 = Matrix(parent(q), nr,nd)
    for i in range(nr):
    	for j in range(nd):
	    G2[i,j] = G1[i,j]

    G3 = Matrix(parent(q), nr, nd)
    for i in range(nr):
    	for j in range(nd):
	    G3[i,j] = G2[i,j].mod(q)

    G4 = G3.substitute(T=0)
    return G4


####################################################################

def teststartT(v1,w1,q1,G): ### test for start matrix after lifting 
    nr = G.nrows()
    nd = G.ncols()

    v = v1.substitute(T=0)
    w = w1.substitute(T=0)
    q = q1.substitute(T=0)
    
    G1 = G.substitute(x2 = v, x3 = w)
    
    G2 = Matrix(parent(q), nr,nd)
    for i in range(nr):
    	for j in range(nd):
	    G2[i,j] = G1[i,j]

    G3 = Matrix(parent(q), nr, nd)
    for i in range(nr):
    	for j in range(nd):
	    G3[i,j] = G2[i,j].mod(q)

    #G4 = G3.substitute(T=0)
    return G3


#####################################################################


def testwithoutdeno(E1, v, w, nq): ### test the rank of the matrix after substitution
    nr = G.nrows()
    nd = G.ncols()
    F1 = E1.substitute(x2 = v, x3 = w)
    F2 = Matrix(parent(nq), nr, nd)
    
    return F2.minors(nr)


def testwithdono(E1,v,w,nq):
    F1 = E1.substitute(x2 = v, x3 = w)

    listminor = F1.minors(2)

    listnum = []
    for i in range(len(listminor)):
    	listnum.append(listminor[i].numerator())

    listdeno = []
    for i in range(len(listminor)):
    	listdeno.append(listminor[i].denominator())
    
    listnumA = []
    for i in range(len(listnum)):
    	listnumA.append(Matrix(AA,1,1))
    for i in range(len(listnum)):
    	listnumA[i][0,0] = listnum[i]

    listnumAA = []
    for i in range(len(listnumA)):
    	listnumAA.append(listnumA[i][0,0])
    
    listdenoA = []
    for i in range(len(listdeno)):
    	listdenoA.append(Matrix(AA,1,1))
    for i in range(len(listdeno)):
    	listdenoA[i][0,0] = listdeno[i]

    listdenoAA = []
    for i in range(len(listdenoA)):
    	listdenoAA.append(listdenoA[i][0,0])

    return listnumAA, listdenoAA

#num1 = test2(E1, vfin, wfin, qfin)[0]
#deno1 = test2(E1, vfin, wfin, qfin)[1]

qA = Matrix(AA,1,1)
qA[0,0] = qfin

qAA = qA[0,0]


