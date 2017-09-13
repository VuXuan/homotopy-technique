R.<x0,x1,x2,x3>=PolynomialRing(GF(9001),4)
k.<t>=PolynomialRing(GF(9001))
FF.<T,x0,x1,x2,x3,t> = PolynomialRing(GF(9001),6)
bF = GF(9001)
z = PowerSeriesRing(bF,'z').gen()



D = [2,1]
p = 2
q = 5
n = q-p+1



def randpoly(d): 
   return add([add([bF.random_element()*R({tuple(a):1}) for a in WeightedIntegerVectors(i,[1 for i in range(R.ngens())])]) for i in range(d+1)])




def nonrand(bF): ## random element in bF - nonzero
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



def condict2(s): ##convert to dictionary for sol1(A)
    J.<t,x0,x1,x2,x3>=PolynomialRing(bF,5)
    k.<t>=PolynomialRing(bF)
    l = []
    for i in range(len(s)):
	ll = []
	for j in range(4):
	    ll.append(k(s[i][0][J.gen(j+1)]))
	l.append(ll)
    return l




def sol21(A): ##g_00 = g_11 = 0 (1st case) 
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(bF)
    J.<t,x0,x1,x2,x3>=PolynomialRing(bF,5)
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
    return E


E = sol21(A)



S.<x2,x3>=PolynomialRing(bF,2)
ld1 = listrandom(S, D, 2, 2, 3)
G = startmat(S,D,ld1,2,3)


def condict(a): ##convert to dictionary
    S.<x2,x3>=PolynomialRing(bF,2)
    l11 = []
    for i in range(len(a)):
	l21 = []
	for j in range(2):
	    l21.append((a[i][0][S.gen(j)]))	    
	    #l21.append(k(a[0][0][S.gen(j)]))
	l11.append(l21)
    return l11




def start231(D,2,3,G,ld1):     ## zero-dimensional parametrization for s1 for substart_00 = substart_02 = 0
    #H = (1-T)*G+T*E   
    #H1 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,1]*H[1,2]-H[1,1]*H[0,2]]
    s11 = []
    for m in range(D[0]):
	for t in (D[0]..2*D[0]-1):
	    s11.append(ideal(ld1[m], ld1[t]).variety()) 
    s1 = condict(s11)
    f = []
    for m in range(len(s1)):
	f.append(s1[m][1]*(prod((x2-s1[l][0])/(s1[m][0] - s1[l][0]) for l in range(len(s1)) if l !=m)))
    f1 = sum(f[m] for m in range(len(f)))
    z1 = prod(x2 - s1[m][0] for m in range(len(s1)))
    l1 = [z1,f1]
    q = z1.substitute(x2 = FF.gen(5))
    w1 = x2.substitute(x2 = FF.gen(5))
    v1 = f1.substitute(x2 = FF.gen(5))
    l = [q,w1,v1] 
    return l


l1 = start231(D,2,3,G,ld1)


def start232(D,2,3,G,ld1):     ## zero-dimensional parametrization for s1 for substart_11 = substart_12 = 0
    #H = (1-T)*G+T*E   
    #H2 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,0]*H[1,2] - H[1,0]*H[0,2]]
    s11 = []
    s11.append(ideal(ld1[4], ld1[5]).variety()) 
    s1 = condict(s11)
    f = []
    for m in range(len(s1)):
	f.append(s1[m][1]*(prod((x2-s1[l][0])/(s1[m][0] - s1[l][0]) for l in range(len(s1)) if l !=m)))
    f1 = sum(f[m] for m in range(len(f)))
    z1 = prod(x2 - s1[m][0] for m in range(len(s1)))
    l1 = [z1,f1]
    q = z1.substitute(x2 = FF.gen(5))
    w1 = x2.substitute(x2 = FF.gen(5))
    v1 = f1.substitute(x2 = FF.gen(5))
    l = [q,w1,v1] 
    return l


l2 = start232(D,2,3,G,ld1)



def start233(D,2,3,G,ld1): ## zero-dimensional parametrization for s1 for substart_00 = substart_11 = 0
    #H = (1-T)*G+T*E   
    #H3 = [H[0,0]*H[1,2] - H[1,0]*H[0,2],H[0,0]*H[1,2]-H[1,1]*H[0,2]]
    s11 = []
    for m in range(D[0]):
	for t in (2*D[0]..2*D[0]+D[1]-1):
	    s11.append(ideal(ld1[m], ld1[t]).variety()) 
    s1 = condict(s11)
    f = []
    for m in range(len(s1)):
	f.append(s1[m][1]*(prod((x2-s1[l][0])/(s1[m][0] - s1[l][0]) for l in range(len(s1)) if l !=m)))
    f1 = sum(f[m] for m in range(len(f)))
    z1 = prod(x2 - s1[m][0] for m in range(len(s1)))
    l1 = [z1,f1]
    q = z1.substitute(x2 = FF.gen(5))
    w1 = x2.substitute(x2 = FF.gen(5))
    v1 = f1.substitute(x2 = FF.gen(5))
    l = [q,w1,v1] 
    return l


l3 = start233(D,2,3,G,ld1)



H = (1-T)*G+T*E
H1 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,1]*H[1,2]-H[1,1]*H[0,2]] #square system for substart_00 = substart_02 = 0
H2 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,0]*H[1,2] - H[1,0]*H[0,2]] #square system for substart_11 = substart_12 = 0
H3 = [H[0,0]*H[1,2] - H[1,0]*H[0,2],H[0,0]*H[1,2]-H[1,1]*H[0,2]] #square system for substart_00 = substart_11 = 0


def lift(H1,l1,prec): #lift from zero-dimensional parametrization for s1 for substart_00 = substart_02 = 0 to zero-dimensional parametrization for 2*3 dense matrix 

#prec = 64???

###TODO: finish

    sol = [l1[1],l1[2]]
    FF.<T,x0,x1,x2,x3,t> = PolynomialRing(bF,6)
    if prec == 1:
	return sol
    else:
	sol_half = lift(H1,l1,prec//2)
	m = jacobian(H1,[FF.gen(3),FF.gen(4)])
    	n = m.substitute(x2 = sol[0],x3 = sol[1])
	g = n.inverse() 
	v = Matrix(FF,len(H1),1)
        for i in range(len(H1)):
            v[i,0] = H1[i].substitute(x2 = sol[0],x3 = sol[1])
	w = Matrix(FF,len(H1),1)
	for i in range(len(H1)):
	    w[i,0] = sol[i]
    return n


n = lift(H1,l1,2)




P.<T> = PolynomialRing(bF)
I.<t> = PowerSeriesRing(bF,1)
J.<T> = PolynomialRing(I)
#S.<t> = PolynomialRing(bF,1)

fR = FractionField(J)

def randpow(d): 
    return sum(I.random_element(3)*((J.gen())^i) for i in range(d+1))


d = [1,2]

def mat(J,p,d):
    M = Matrix(J, p, p)
    for i in range(p):	
	for j in range(p):
	    M[i,j] = randpow(d[i])
    return M
    
M = mat(J,2,d)


def inver(fR,M,p):
    a = M.det()
    B = Matrix(fR, p, p)
    for i in range(p):
	for j in range(p):
	    N = M.delete_columns([j])
	    K = N.delete_rows([i])
	    B[i,j] = (-1)^(i+j)*(K.det())/a
    return B.transpose()





def sol22(A): ##g_00 = g_11 = 0 (2nd case) 
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(bF)
    J.<t,x0,x1,x2,x3>=PolynomialRing(bF,5)
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
    return E


E1 = sol22(A)


### Use the same G as the 1st case
HH = (1-T)*G+T*E1
HH1 = [HH[0,0]*HH[1,1] - HH[1,0]*HH[0,1],HH[0,1]*HH[1,2]-HH[1,1]*HH[0,2]] #square system for substart_00 = substart_02 = 0
HH2 = [HH[0,0]*HH[1,1] - HH[1,0]*HH[0,1],HH[0,0]*HH[1,2] - HH[1,0]*HH[0,2]] #square system for substart_11 = substart_12 = 0
HH3 = [HH[0,0]*HH[1,2] - HH[1,0]*HH[0,2],HH[0,0]*HH[1,2]-HH[1,1]*HH[0,2]] #square system for substart_00 = substart_11 = 0



def ratrec(sol2,e): ## rational reconstruction for sol2(A) at degree(e,e)
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




def deno(rsol): ## denominators list
    deno1 = []
    for i in range(len(rsol)): 
        deno2 = []
	for j in range(len(rsol[i])): 
	    deno2.append(denominator(rsol[i][j]))
        deno1.append(deno2) 
    return deno1



def ss(deno): ## lcm list
    a = []
    for i in range(len(deno)):
	a.append(lcm(deno[i][j] for j in range(len(deno[i])))) 
    return a



def deldeno(a,b): ## delete denominatiors
    a1 = []
    for i in range(len(a)): 
	a2 = []
	for j in range(len(a[i])):
	    a2.append(a[i][j]*b[i])
        a1.append(a2)
    return a1


