k.<t>=PolynomialRing(GF(9001))
K.<T,x>=PolynomialRing(GF(9001),2)
bF = GF(9001)

FF.<T,x0,x1,x2,x3> = PolynomialRing(GF(9001),5)
kk.<t,x0,x1,x2,x3> = PolynomialRing(GF(9001),5)
R.<x0,x1,x2,x3>=PolynomialRing(GF(9001),4)

#H.<t,T,x0,x1,x2,x3> = PolynomialRing(GF(9001),6)

z = PowerSeriesRing(bF,'z').gen()



D = [2,1]
p = 2
q = 5
n = q-p+1



def randpoly(d): 
   return add([add([GF(9001).random_element()*R({tuple(a):1}) for a in WeightedIntegerVectors(i,[1 for i in range(R.ngens())])]) for i in range(d+1)])



def nonrand(bF): ## random element in Q - nonzero
    r = bF.random_element()
    while r == 0:   
        r = bF.random_element()
    return r



def randlin(R): ## linear form with random coefficients nonzero
    n = R.ngens()
    Q = R.base_ring()
    return sum(nonrand(bF)* R.gens()[i] for i in range(n)) + nonrand(bF)



def multi(R,d): ## product of d linear form with generic coefficients
    a = randlin(R)
    for i in range(1,d):
        a = a*randlin(R)
    return a



def inputmat(D): ## Input matrix F - random
    mat = Matrix(R, q, p, [[randpoly(D[i]) for i in range(p)] for j in range(q)])
    C = transpose(mat)
    return C

 

def listrandom(R, D, n, p, q): # list of all linear random for the start matrix 
    d = n*(sum(D[k] for k in range(len(D))))
    ld = []
    for i in range(d):
        ld.append(randlin(R))
    return ld 

ld = listrandom(R, D, n, p, q)



def startmat(R, D,ld, p, q):  ## start matrix - p∗q polynomial matrix with the entry (i,j) is product of Di linear form with generic coefficients, first block is a diag matrix
    p = p
    q = q
    d = sum(D[i] for i in range(p))  
    A = Matrix(R, p, q)
    A[0,0] = prod(ld[k] for k in (0..(D[0] - 1)))
    for j in (1..q-1):
	if j > p-1:
	    A[0,j] = prod(ld[k] for k in ((j-p+1)*D[0]..((j-p+2)*D[0]-1)))
	else: 
	    A[0,j] = 0
    for i in (1..p-1):
	a = p*sum(D[k] for k in range(i))
	for j in (1..q-1):
	    if j == i:
		A[i,j] = prod(ld[k] for k in (a..(a+D[j]-1)))
	    elif j > p-1:
		A[i,j] = prod(ld[t] for t in ((a+(j-p+1)*D[i])..(a+(j-p+2)*D[i]-1)))
	    else: 
		A[i,j] = 0
    return A


A = startmat(R,D,ld,2,5)
C = inputmat(D)
B = (1-FF.gens()[0])*A + FF.gens()[0]*C ## Homotopy matrix



def startsol(A): #Zero-dimensional parametrizations for start matrix
    soll1 = sol1(A)
    soll2 = sol2(A)
    rsol = ratrec(soll2,255)
    denolist = deno(rsol)
    lcmlist = ss(denolist)
    seccase = deldeno(rsol,lcmlist)
    return soll1+seccase



def sol1(A):
    s11 = []  #g_00 = g_0,p = ... = g_0,q-1 = 0 
    for i in range(2):
	for t in (2..3):
	    for l in (4..5):
		for m in (6..7):
		    s11.append(ideal(ld[i],ld[t],ld[l],ld[m]).variety())
    s12 = [] #g_11 = g_1,p = ... = g_1,q-1 = 0
    s12.append(ideal(ld[8],ld[9],ld[10],ld[11]).variety())
    s = s11+s12
    return condict2(s)



def condict2(s): #convert to dictionary for sol1(A)
    J.<t,x0,x1,x2,x3>=PolynomialRing(GF(9001),5)
    k.<t>=PolynomialRing(GF(9001))
    l = []
    for i in range(len(s)):
	ll = []
	for j in range(4):
	    ll.append(k(s[i][0][J.gen(j+1)]))
	l.append(ll)
    return l



def sol2(A): #g_00 = g_11 = 0
    return sol21(A)+sol22(A)



def ratrec(sol2,e): # rational reconstruction for sol2(A) at degree(e,e)
    rsol1 = []
    for i in range(len(sol2)):
	rsol2 = []
	for j in range(len(sol2[i])):
	    rsol2.append(sol2[i][j].substitute(t = z))
	rsol1.append(rsol2)
    rsol = []
    for i in range(len(rsol1)):
	rsol3 = []
	for j in range(len(rsol1[i])):
	    rsol3.append(rsol1[i][j].pade(e,e))
	rsol.append(rsol3)
    return rsol




def deno(rsol): # denominators list
    deno1 = []
    for i in range(len(rsol)): 
        deno2 = []
	for j in range(len(rsol[i])): 
	    deno2.append(denominator(rsol[i][j]))
        deno1.append(deno2) 
    return deno1



def ss(deno): #lcm list
    a = []
    for i in range(len(deno)):
	a.append(lcm(deno[i][j] for j in range(len(deno[i])))) 
    return a



def deldeno(a,b): # delete denominatiors
    a1 = []
    for i in range(len(a)): 
	a2 = []
	for j in range(len(a[i])):
	    a2.append(a[i][j]*b[i])
        a1.append(a2)
    return a1



def sol21(A): #g_00 = g_11 = 0 (1st case)
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(GF(9001))
    J.<t,x0,x1,x2,x3>=PolynomialRing(GF(9001))
    G = copy(A)
    s01 = []
    m = x0 - 1/ld[0].coefficients()[0]*ld[0]
    ld[8] = ld[8].substitute(x0 = m)
    x12 = x1 - 1/ld[8].coefficients()[0]*ld[8]
    x02 = m.substitute(x1 = x12)
    for j in (p..q-1):
	for i in range(p):
	    G[i,j] = G[i,j].substitute(x0 = x02, x1 = x12)
    E = G[0:2,2:5]	
    P = sub(D,2,3,E)
    F = []
    for i in range(len(P)):
	F1 = [x02.substitute(x2 = P[i][0],x3 = P[i][1]),x12.substitute(x2 = P[i][0],x3 = P[i][1]),P[i][0],P[i][1]]
    	F.append(F1)
    return F



def sol22(A): #g_00 = g_11 = 0 (2nd case) 
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(GF(9001))
    J.<t,x0,x1,x2,x3>=PolynomialRing(GF(9001),5)
    G = copy(A)
    s01 = []
    m = x0 - 1/ld[1].coefficients()[0]*ld[1]
    ld[8] = ld[8].substitute(x0 = m)
    x12 = x1 - 1/ld[8].coefficients()[0]*ld[8]
    x02 = m.substitute(x1 = x12)
    for j in (p..q-1):
	for i in range(p):
	    G[i,j] = G[i,j].substitute(x0 = x02, x1 = x12)
    E = G[0:2,2:5]	
    P = sub(D,2,3,E)
    F = []
    for i in range(len(P)):
	F1 = [x02.substitute(x2 = P[i][0],x3 = P[i][1]),x12.substitute(x2 = P[i][0],x3 = P[i][1]),P[i][0],P[i][1]]
    	F.append(F1)
    return F



def sub(D,p,q,l2): #2*3 submatrix
    S.<x2,x3>=PolynomialRing(GF(9001),2)
    randlin(S)
    ld = listrandom(S, D, 2, 2, 3)
    F = l2
    G = startmat(S,D,ld,2,3)
    H = (1-T)*G+T*F
    s11 = []
    for m in range(D[0]):
	for t in (D[0]..2*D[0]-1):
	    s11.append(ideal(ld[m], ld[t]).variety()) 

    H11 = [H[0,0]*H[1,1],-H[1,1]*H[0,2]]
    l11 = condict(s11)
    lift(H11,l11,64)
    s12 = [] 
    s12.append(ideal(ld[4],ld[5]).variety())
    H12 = [H[0,0]*H[1,1],H[0,0]*H[1,2]]
    l12 = condict(s12)
    lift(H12,l12,64)
    s13 = []
    for m in range(D[0]):
	s13.append(ideal(ld[m], ld[4]).variety())
    H13 = [H[0,0]*H[1,2],-H[1,1]*H[0,2]]
    l13 = condict(s13)
    lift(H13,l12,64)
    listlift = lift(H11,l11,64)+lift(H12,l12,64)+lift(H13,l13,64)
    rsol = ratrec(listlift,31)
    denolist = deno(rsol)
    lcmlist = ss(denolist)
    seccase = deldeno(rsol,lcmlist)
    return seccase #lift(H11,l11,64)+lift(H12,l12,64)+lift(H13,l13,64)



def condict(a): #convert to dictionary
    S.<x2,x3>=PolynomialRing(GF(9001),2)
    k.<t>=PolynomialRing(GF(9001))
    l11 = []
    for i in range(len(a)):
	l21 = []
	for j in range(2):
	    l21.append(k(a[0][0][S.gen(j)]))
	l11.append(l21)
    return l11



def lift(H,ll,prec):
    RS = []
    for i in range(len(ll)):
        RS.append(lift1(H,ll[i],prec))
    return RS


def lift1(H,sol,prec): #lift solutions for submatrix of dimension 2*3 for start matrix
    if prec == 1:
        return sol
    else:
        sol_half = lift1(H,sol,prec//2)
        I = quotient(parent(sol[0]),t**prec)
        m = jacobian(H, [FF.gen(0),FF.gen(3),FF.gen(4)])
        dic = dict()
        dic[FF.gen(0)] = I.0
	dic[FF.gen(3)] = sol_half[0]
	dic[FF.gen(4)] = sol_half[1]
        n1 = Matrix(I,len(H),len(H))
        for i in range(len(H)):
            for j in range(len(H)):
                n1[i,j] = m[i,j].subs(dic)
        n = inver(n1, sol, prec)
        a = matrix(len(H),1,[sol_half[i] for i in range(len(sol_half))])
        v = Matrix(I,len(H),1)
        for i in range(1,len(H)):
            v[i,0] = H[i].subs(dic)
        p = a - n*v
        sol_new = []
        for i in range(len(H)):
            sol_new.append(p[i][0])
        result = []
        for i in range(len(sol_new)):
            result.append(sol_new[i].lift())    
        return result 



def inver(N,sol,deg):
    k.<t>=PolynomialRing(GF(9001))
    I = quotient(parent(sol[0]),t**(deg))
    Nbar = []
    P = []
    for k in range(deg):
        Nbar.append(matrix([[N[i][j][k] for j in range(len(sol))] for i in range(len(sol))]))       
    P.append(Nbar[0].inverse())
    for i in range(1,deg):
       P.append(Nbar[0].inverse()*(-sum(Nbar[i-j]*P[j] for j in range(i))))
    return sum(P[i]*(I.0**i) for i in range(0,len(P)))
        
     






















