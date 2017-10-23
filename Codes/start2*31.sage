
bF = GF(9007)

R.<x0,x1,x2,x3>=PolynomialRing(bF,4)

AA.<t> = PolynomialRing(bF)
BB.<T> = PolynomialRing(AA)

FF.<x0,x1,x2,x3> = PolynomialRing(BB,4)

fR = FractionField(BB)
JJ.<T,t> = PolynomialRing(bF,2)

S.<x2,x3> = PolynomialRing(bF,2)

R1.<x2,x3,t,T> = PolynomialRing(bF,4)


u = x2

R2.<x0,x1,x2,x3,t,T> = PolynomialRing(bF,6)

k.<T> = PolynomialRing(bF)
k1.<t> = PolynomialRing(k)
K = FractionField(k);
A1.<t> = PolynomialRing(K) 



D = [2,1]
p = 2
q = 5
n = q-p+1


###################################################################

def star1(S,D):
    ld1 = listrandom(S, D, 2, 2, 3)

    I = []
    for m in range(D[0]):
	for t in (D[0]..2*D[0]-1):
	    I.append(ideal(ld1[m], ld1[t]))         

    for i in range(2*D[0]):
        if I[i].dimension() != 0: 
	    ld1 = listrandom(S, D, 2, 2, 3)

    if ideal(ld1[4], ld1[5]).dimension() !=0:
	ld1 = listrandom(S, D, 2, 2, 3)

    G = startmat(S,D,ld1,2,3)

    J = []
    for m in range(D[0]):
	for t in (2*D[0]..2*D[0]+D[1]-1):
	    J.append(ideal(ld1[m], ld1[t])) 

    for i in range(len(J)):
        if J[i].dimension() != 0: 
	    ld1 = listrandom(S, D, 2, 2, 3)

    return ld1,G


ld1,G = star1(S,D)

 #### zero-dim para components for 2*3 start matrix G ####



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


##################################################################
## zero-dimensional para for first case ###

def ch(AA,S,D,ld1,G):  #### substart_00 = substart_02 = 0
    s11 = []   
    for m in range(D[0]):
	for t in (D[0]..2*D[0]-1):
	    s11.append(ideal(ld1[m], ld1[t]).variety()) 
    s1 = condict(s11)
    
    v = AA.lagrange_polynomial(s1)
    
    z = prod(x2 - s1[m][0] for m in range(len(s1)))
    w = x2.substitute(x2 = AA.gen(0))
    q = z.substitute(x2 = AA.gen(0))

    polelim11 = Matrix(R1,1,1)
    polelim11[0,0] = q

    polelim1 = polelim11[0,0]

    param1 = Matrix(R1,2,1)
    param1[0,0] = w
    param1[1,0] = v

    return param1, polelim1


#############z zero-dim para for second case ######


def ch2(AA,S,D,ld1,G):  ###substart_11 = substart_12 = 0
    s12 = []
    s12.append(ideal(ld1[4], ld1[5]).variety())

    s2 = condict(s12)

    v2 = AA.lagrange_polynomial(s2);
    
    z2 = prod(x2 - s2[m][0] for m in range(len(s2)))
    w2 = x2.substitute(x2 = AA.gen(0))
    q2 = z2.substitute(x2 = AA.gen(0))

    polelim22 = Matrix(R1,1,1)
    polelim22[0,0] = q2

    polelim2 = polelim22[0,0]

    param2 = Matrix(R1,2,1)
    param2[0,0] = w2
    param2[1,0] = v2
    
    return param2,polelim2



		###############################

##### zero-dim para for the 3rd case #####



def ch3(AA,S,D,ld1,G): ###substart_00 = substart_11 = 0
    s13 = []
    for m in range(D[0]):
	for t in (2*D[0]..2*D[0]+D[1]-1):
	    s13.append(ideal(ld1[m], ld1[t]).variety()) 

    s3 = condict(s13)

    v3 = AA.lagrange_polynomial(s3);
    
    z3 = prod(x2 - s3[m][0] for m in range(len(s3)))
    w3 = x2.substitute(x2 = AA.gen(0))
    q3 = z3.substitute(x2 = AA.gen(0))

    polelim33 = Matrix(R1,1,1)
    polelim33[0,0] = q3

    polelim3 = polelim33[0,0]

    param3 = Matrix(R1,2,1)
    param3[0,0] = w3
    param3[1,0] = v3
    
    return param3, polelim3
    
#####################################################################

### Find square system to lifting ####


def chek1(H,GG,EE,a):

    p1 = GG.nrows()
    q1 = GG.ncols()

    #G = GG[0:p1,p1:q1]

    l = a[0]
    q = a[1]

    J = GG.substitute(x2 = l[0,0], x3 = l[1,0])
    
    p2 = J.nrows()
    q2 = J.ncols()

    J1 = Matrix(R1,p2,q2)
    for i in range(p2):
    	for j in range(q2):
	    J1[i,j] = J[i,j].mod(q)

    E = J1.rref()
    
    Q = []
    for i in range(p2):
    	for j in range(q2):
            if E[i, j] != 0:
               Q.append(j)
	       break
    A1 = range(q1)

    
    N = diff(A1,Q)
    
    B = []
    
    for i in range(len(N)):
    	B.append(sorted(Q+[N[i]]))

    W = []
    for i in range(len(B)):
    	W.append(H.matrix_from_columns(B[i]).det())

    WG = []
    for i in range(len(B)):
    	WG.append(GG.matrix_from_columns(B[i]).det())


    W1 = []
    for i in range(len(B)):
    	W1.append(EE.matrix_from_columns(B[i]).det())
    return W,WG,W1



def diff(first, second): ## compute list difference
        second = set(second)
        return [item for item in first if item not in second]
	
#####################################################################

GG = Matrix(R1,2,3)
for i in  range(2):
    for j in range(3):
	GG[i,j] = G[i,j]



EE = E1.substitute(x2 = R1.gen(0), x3 = R1.gen(1))
H = (1-R1.gen(3))*GG+R1.gen(3)*EE

####################################################################

a1 = ch(AA,S,D,ld1,G)
a2 = ch2(AA,S,D,ld1,G)
a3 = ch3(AA,S,D,ld1,G)


H1 = chek1(H,GG,EE,a1)[0]

H2 = chek1(H,GG,EE,a2)[0]

H3 = chek1(H,GG,EE,a3)[0]
