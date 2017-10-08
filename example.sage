bF = GF(9001)

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
    m = Matrix(R, p, q)
    for i in range(p): 
	for j in range(q):
	    m[i,j] = randpoly(D[i])
    #C = transpose(m)
    return m




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




AM = startmat(R,D,ld,2,5)

C = inputmat(D)

Hmat = (1-R2.gen(5))*AM+R2.gen(5)*C

def sol21(A): ##g_00 = g_11 = 0 (1st case) 
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(bF)
    J.<t,x0,x1,x2,x3>=PolynomialRing(bF,5)
    	
    G5 = copy(A)
    s01 = []
    m = x0 - 1/ld[0].coefficients()[0]*ld[0]
    ld[8] = ld[8].substitute(x0 = m)
    x12 = x1 - 1/ld[8].coefficients()[0]*ld[8]
    x02 = m.substitute(x1 = x12)
    for j in (p..q-1):
	for i in range(p):
	    G5[i,j] = G5[i,j].substitute(x0 = x02, x1 = x12)
    E1 = G5[0:2,2:5]

    return E1



E1 = sol21(AM)



def sol22(AM): ##am_00 = am_11 = 0 (2nd case) 
    (x0, x1, x2, x3) = R.gens()

    G = copy(AM)
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



E2 = sol22(AM)





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

def ch(AA,S,D,ld1,G):  #### substart_00 = substart_02 = 0
    s11 = []   
    for m in range(D[0]):
	for t in (D[0]..2*D[0]-1):
	    s11.append(ideal(ld1[m], ld1[t]).variety()) 
    s1 = condict(s11)
    
    v = AA.lagrange_polynomial(s1);
    
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

#l1 = ch(AA,S,D,ld1,G)


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
    #return q2,w2,v2

#l2 = ch2(AA,S,D,ld1,G)

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



#l3 = ch3(AA,S,D,ld1,G)


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


GG = Matrix(R1,2,3)
for i in  range(2):
    for j in range(3):
	GG[i,j] = G[i,j]



EE = E1.substitute(x2 = R1.gen(0), x3 = R1.gen(1))
H = (1-R1.gen(3))*GG+R1.gen(3)*EE




EE2 = E2.substitute(x2 = R1.gen(0), x3 = R1.gen(1))
HH = (1-R1.gen(3))*GG+R1.gen(3)*EE2




def soll(H,prec): ##ld[0] = ld[8] = 0
    H1 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,1]*H[1,2]-H[1,1]*H[0,2]]

    (param1, polelim1) = ch(AA,S,D,ld1,G)

    l1 = HenselLift(H1, param1, polelim1, u, prec) 

    ll1 = []
    for i in range(len(l1)):
	ll1.append(l1[i].substitute(T = 1))

    H2 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,0]*H[1,2]-H[1,0]*H[0,2]]

    (param2, polelim2) = ch2(AA,S,D,ld1,G)
  
    l2 = HenselLift(H2, param2, polelim2, u, prec)
    
    ll2 = []
    for i in range(len(l2)):
	ll2.append(l2[i].substitute(T = 1))
	
    H3 = [H[0,0]*H[1,2] - H[1,0]*H[0,2],H[0,1]*H[1,2]-H[1,1]*H[0,2]]

    (param3, polelim3) = ch3(AA,S,D,ld1,G)

    l3 = HenselLift(H3, param3, polelim3, u, prec)

    ll3 = []
    for i in range(len(l3)):
	ll3.append(l3[i].substitute(T = 1))

    l = Matrix(R1,2,1) ## TODO: Check l[0,0] (i.e. x2)
    l[0,0] = ll1[0][0,0]
    l[1,0] = ll1[0][1,0]*ll2[0][1,0]*ll3[0][1,0]

    q = ll1[1]*ll2[1]*ll3[1]
    return ll1,ll2,ll3#l,q



def subb(AM,H,prec): ## FINAL SOLUTION  am_00 = am_11 = 0 (1st case) 
    (x0, x1, x2, x3) = R.gens()
    k.<t>=PolynomialRing(bF)
    J.<t,x0,x1,x2,x3>=PolynomialRing(bF,5)
    	
    G5 = copy(AM)
    s01 = []
    m = x0 - 1/ld[0].coefficients()[0]*ld[0]
    ld[8] = ld[8].substitute(x0 = m)
    x12 = x1 - 1/ld[8].coefficients()[0]*ld[8]
    x02 = m.substitute(x1 = x12)

#    (ll1,ll2,ll3) = soll(H,prec)

    (ll1,ll2,ll3) = soll(H,prec)

    

    v01 = x02.substitute(x2 = ll1[0][0,0], x3 = ll1[0][1,0])
    v11 = x12.substitute(x2 = ll1[0][0,0], x3 = ll1[0][1,0])

    a1 = [(v01,v11,ll1[0][0,0],ll1[0][1,0]),ll1[1]]

    v02 = x02.substitute(x2 = ll2[0][0,0], x3 = ll2[0][1,0])
    v12 = x12.substitute(x2 = ll2[0][0,0], x3 = ll2[0][1,0])

    a2 = [(v02,v12,ll2[0][0,0],ll2[0][1,0]),ll2[1]]

    v03 = x02.substitute(x2 = ll3[0][0,0], x3 = ll3[0][1,0])
    v13 = x12.substitute(x2 = ll3[0][0,0], x3 = ll3[0][1,0])

    a3 = [(v03,v13,ll3[0][0,0],ll3[0][1,0]),ll3[1]]
    return [a1,a2,a3]


def sol2(HH,prec): ##ld[1] = ld[8]
    HH1 = [HH[0,0]*HH[1,1] - HH[1,0]*HH[0,1],HH[0,1]*HH[1,2]-H[1,1]*HH[0,2]]

    (param1, polelim1) = ch(AA,S,D,ld1,G)

    l12 = HenselLift(HH1, param1, polelim1, u, prec) 

    HH2 = [HH[0,0]*HH[1,1] - HH[1,0]*HH[0,1],HH[0,0]*HH[1,2]-HH[1,0]*HH[0,2]]

    (param2, polelim2) = ch2(AA,S,D,ld1,G)
  
    l22 = HenselLift(HH2, param2, polelim2, u, prec)


    HH3 = [HH[0,0]*HH[1,2] - HH[1,0]*HH[0,2],HH[0,1]*HH[1,2]-HH[1,1]*HH[0,2]]

    (param3, polelim3) = ch3(AA,S,D,ld1,G)

    l32 = HenselLift(HH3, param3, polelim3, u, prec)


    lll2 = Matrix(R1,2,1)
    lll2[0,0] = l12[0][0,0]
    lll2[1,0] = l12[0][1,0]*l22[0][1,0]*l32[0][1,0]

    q2 = l12[1]*l22[1]*l32[1]
    return lll2,q2


def subb2(AM,HH,prec): ##FINAL SOLUTION am_00 = am_11 = 0 (2nd case) 
    (x0, x1, x2, x3) = R.gens()

    G = copy(AM)
    s01 = []
    m = x0 - 1/ld[1].coefficients()[0]*ld[1]
    ld[8] = ld[8].substitute(x0 = m)
    x12 = x1 - 1/ld[8].coefficients()[0]*ld[8]
    x02 = m.substitute(x1 = x12)
   
    (l2,q2) = sol2(HH,prec)
    
    v00 = x02.substitute(x2 = l[0,0], x3 = l[1,0])
    v11 = x12.substitute(x2 = l[0,0], x3 = l[1,0])
    return [(v00,v11,l2[0],l2[1]), q2]


def soll1(H,prec): 
    H1 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,1]*H[1,2]-H[1,1]*H[0,2]]
    (param1, polelim1) = ch(AA,S,D,ld1,G)
    l1 = HenselLift(H1, param1, polelim1, u, prec) 

    return l1


def soll2(H,prec): ### TODO: Check the solution

    H2 = [H[0,0]*H[1,1] - H[1,0]*H[0,1],H[0,0]*H[1,2]-H[1,0]*H[0,2]]

    (param2, polelim2) = ch2(AA,S,D,ld1,G)

    l2 = HenselLift(H2, param2, polelim2, u, prec)

    return l2



def soll3(H,prec):
    H3 = [H[0,0]*H[1,2] - H[1,0]*H[0,2],H[0,1]*H[1,2]-H[1,1]*H[0,2]]
    (param3, polelim3) = ch3(AA,S,D,ld1,G)
    l3 = HenselLift(H3, param3, polelim3, u, prec)
    return l3



def HenselLift(F, param, polelim, u, myprec): 
    
    prec = 1
    V = param
   # print "V = ", V, "\n"

    qq = Matrix(A1,1,1)
    qq[0,0] = polelim
    q = qq[0,0]
    
    #print "q = ", q, "\n"

    while prec <= myprec:

    	#print "prec = ", prec, "\n"
	
        J = jacobian(F,(R1.gen(0),R1.gen(1)))
        JJ = J.substitute(x2 = V[0,0], x3 = V[1,0])
        d = JJ.det()

	dd1 = Matrix(A1,1,1)
        dd1[0,0] = d
        dd0 = dd1[0,0]
	
	#print "dd0 = ", dd0, "\n"

	s = gc(dd0,q, 2*prec)
	
	#print "s = ", s, "\n"
	
	sJ = JJ.inverse()
 	
	#print "sJ = ", sJ, "\n"
	
	sF = Matrix(R1,1,2)
        sF[0,0] =F[0].substitute(x2 = V[0,0], x3 = V[1,0])
        sF[0,1] =F[1].substitute(x2 = V[0,0], x3 = V[1,0])

	#print "sF = ", sF, "\n"
	
	NEW = sJ*(sF.transpose())

	#print "NEW = ", NEW, "\n"

	dd2 = NEW[0,0].denominator()

	#print "dd2 = ", dd2, "\n"

        dd3 = Matrix(A1,1,1)
        dd3[0,0] = dd2
        dd = dd3[0,0]

	s0 = gc(dd,q,2*prec)
	
	s01 = Matrix(R1,1,1)
	s01[0,0] = s0
	s = s01[0,0]
       
	#print "s = ", s, "\n"

	ddd2 = NEW[1,0].denominator()

	#print "ddd2 = ", ddd2, "\n"

        ddd3 = Matrix(A1,1,1)
        ddd3[0,0] = ddd2
        ddd = ddd3[0,0]
	
	s1 = gc(ddd,q,2*prec)

	s32 = Matrix(R1,1,1)
	s32[0,0] = s1
	s3 = s32[0,0]
	
	#print "s3 = ", s3, "\n"

	NEWV1 = [ V[0,0] - s*NEW[0,0].numerator(), V[1,0] - s3*NEW[1,0].numerator() ] 

	#print "NEWV1 = ", NEWV1, "\n"

	q12 = Matrix(k1,1,1)
	q12[0,0] = q
	q100 = q12[0,0]

	N0 = Matrix(k1,1,1)
	N0[0,0] = NEWV1[0]
	N00 = N0[0,0]
	
	r0 = N00.mod(q100)

	
	#print "r0 = ", r0, "\n"
	

	r011 = r0//(T^(2*prec))
	r01 = r0-r011*(T^(2*prec))

	#print "r01 = ", r01, "\n"

	N1 = Matrix(k1,1,1)
	N1[0,0] = NEWV1[1]
	N10 = N1[0,0]
	
	r1 = N10.mod(q100)

	#print "r1 = ", r1, "\n"
	
	r111 = r1//(T^(2*prec))
	r11 = r1-r111*(T^(2*prec))

	#print "r11 = ", r11, "\n"

	NEWW = Matrix(R1,1,2)
	NEWW[0,0] = r01
	NEWW[0,1] = r11

	r01 =  NEWW[0,0]
	r11 = NEWW[0,1]

	NEWV = [r01,r11]

	#print "NEWV = ", NEWV, "\n"

	B = Matrix(R1,1,1)
        B[0,0] = t

	delta = u.substitute(x2 = NEWV[0], x3 = NEWV[1]) - B[0,0]
	#print "delta = ", delta, "\n"

	Y01 = delta*derivative(NEWV[0],R1.gen(2))

	Y0 = Matrix(k1,1,1)
	Y0[0,0] = Y01
	Y01 = Y0[0,0]
	
	Y0 = Y01.mod(q100)

	Y000 = Y0//(T^(2*prec))
	Y00 = Y0-Y000*(T^(2*prec))

	#print "Y00 = ", Y00, "\n"

	X01 = Matrix(R1,1,1)
	X01[0,0] = Y00
	Y00 = X01[0,0]

	Y0 =  NEWV[0] - Y00
	#print "Y0 = ", Y0, "\n"


	Y11 = delta*derivative(NEWV[1],R1.gen(2))

	Y1 = Matrix(k1,1,1)
	Y1[0,0] = Y11
	Y11 = Y1[0,0]
	
	Y1 = Y11.mod(q100)

	
	Y111 = Y1//(T^(2*prec))
	Y11 = Y1-Y111*(T^(2*prec))

	#print "Y11 = ", Y11, "\n"
	
	X11 = Matrix(R1,1,1)
	X11[0,0] = Y11
	Y11 = X11[0,0]

	Y1 =  NEWV[1] - Y11

	#print "Y1 = ", Y1, "\n"

	VERYNEWV = [Y0,Y1]

	#print "VERYNEWV = ", VERYNEWV, "\n"
        
	del1 = Matrix(k1,1,1)
	del1[0,0] = delta
	delta1 = del1[0,0]

	#print "delta1 = ", delta1, "\n"

	q01 = q100 - (delta1*derivative(q100)).mod(q100)

	q11 = q01//(T^(2*prec))
	NEWq = q01-q11*(T^(2*prec))

	V[0,0] = VERYNEWV[0]
        V[1,0] = VERYNEWV[1]

	NEWqq = Matrix(A1,1,1)
	NEWqq[0,0] = NEWq

        q = NEWqq[0,0]

        prec = 2*prec

	#print "V = ", V, "\n"
	#print "q = ", q, "\n"
	
    return V,q



def gc(p5,q5,prec):
    RR.<T> = PowerSeriesRing(bF)
    (g,s,w) = xgcd(p5,q5)
    d = s.degree()
    L = s.coefficients(sparse=False)
    a = []
    for i in range(len(L)): 
        a.append((L[i].subs(T=T).O(prec)).truncate(prec))

    m = sum(a[i]*t^(i) for i in range(len(a)))

    return m






ll1 = soll1(H,2)

l1 = []
for i in range(len(ll1)):
    l1.append(ll1[i].substitute(T = 1))


ll2 = soll2(H,2)

l2 = []

for i in range(len(ll2)):
    l2.append(ll2[i].substitute(T = 1))

ll3 = soll3(H,2)

l3 = []
for i in range(len(ll3)):
    l3.append(ll3[i].substitute(T = 1))


def uni(p1,p2):
    q1 = p1[1]
    q2 = p2[1]
    d = gcd(q1,q2)
    if d.degree() == 0:
	q = q1*q2
    else:
	q = 1
	a1 = list(q1.factor())
	a2 = list(q2.factor())
    return q


#prec = 2

#aa1 = subb(AM,H,prec)

#print "1st component zero-dim para for start matrix  is aa1 =  ", aa1, "\n"

#aa2 = subb(AM,HH,prec)

#print "2nd component zero-dim para for start matrix  is aa2 =  ", aa2, "\n"




a100 = subb(AM,H,2)

uu = x2




def HenselLift4(F, param, polelim, uu, myprec): 
    
    prec = 1
    V = param
    print "V = ", V, "\n"

    qq = Matrix(A,1,1)
    qq[0,0] = polelim
    q = qq[0,0]
    
    print "q = ", q, "\n"

    while prec <= myprec:

        print "prec = ", prec, "\n"
	
	J = jacobian(F,(x0,x1,x2,x3))
	JJ = J.substitute(x0 = V[0,0], x1 = V[1,0], x2 = V[2,0], x3 = V[3,0])
	d = JJ.det()

	dd1 = Matrix(A,1,1)
	dd1[0,0] = d
	dd0 = dd1[0,0]
	
	print "dd0 = ", dd0, "\n"

	s = gc(dd0,q, 2*prec)
	
	print "s = ", s, "\n"
	
	sJ = JJ.inverse()
 	
	print "sJ = ", sJ, "\n"
	
	sF = Matrix(R2,1,4)
    
	for i in range(4):
	    sF[0,i] = F[i].substitute(x0 = V[0,0], x1 = V[1,0], x2 = V[2,0], x3 = V[3,0])
   
	print "sF = ", sF, "\n"
	
	NEW = sJ*(sF.transpose())

	print "NEW = ", NEW, "\n"

	deno1 = []
	for i in range(4):
	    deno1.append(NEW[i,0].denominator())
    

	deno2 = []
	for i in range(len(deno1)):
	    deno2.append(Matrix(A,1,1))
	for i in range(len(deno1)):
	    deno2[i][0,0] = deno1[i]

        deno = []
	for i in range(len(deno1)):
	    deno.append(deno2[i][0,0])    
       

	print "deno = ", deno, "\n"


	nume = [] ##numerator of NEW
	for i in range(4):
	    nume.append(NEW[i,0].numerator())

        print "nume = ",nume, "\n"


	lis1 = [] ##list of gcd between deno.NEW and q
	for i in range(4):
	    lis1.append(gc(deno[i],q,prec))


	lis = []
	for i in range(4):
	    lis.append(lis1[i].substitute(t = R2.gen(4), T = R2.gen(5)))

	print "lis = ",lis, "\n"

	NEWV1 = []
	for i in range(4): 
	    NEWV1.append(V[i,0] - lis[i]*nume[i])

	print "NEWV1 = ", NEWV1, "\n"

    
	nn1 = []
	for i in range(len(NEWV1)):
	    nn1.append(Matrix(k1,1,1))

	for i in range(len(NEWV1)):
	    nn1[i][0,0] = NEWV1[i]

	print "nn1 = ", nn1, "\n"

	nn2 = []
	for i in range(len(NEWV1)):
	    nn2.append(nn1[i][0,0])

	q12 = Matrix(k1,1,1)
	q12[0,0] = q
	q100 = q12[0,0]

	print "q100 = ", q100, "\n"

	nn3 = []
	for i in range(len(nn1)):
	    nn3.append(nn2[i].mod(q100))

	print "nn3 = ", nn3, "\n"    

### take prec T^(2*prec) for any nn3[i]

	nn4 = []
	for i in range(len(nn3)): 
	    nn4.append(nn3[i]//(T^(2*prec)))

	nn = []
	for i in range(len(nn4)):
	    nn.append(nn3[i] - nn4[i]*(T^(2*prec)))

	print "nn = ", nn, "\n"

	pNEWV = []
	for i in range(len(nn)):
	    pNEWV.append(Matrix(R2,1,1))

	for i in range(len(nn)): 
	    pNEWV[i][0,0] = nn[i]

	NEWV = []
	for i in range(len(nn)):
	    NEWV.append(pNEWV[i][0,0])

	print "NEWV = ", NEWV, "\n"

	B = Matrix(R2,1,1)
        B[0,0] = t

        delta = u.substitute(x0 = NEWV[0], x1 = NEWV[1], x2 = NEWV[2], x3 = NEWV[3]) - B[0,0]
    
        print "delta = ", delta, "\n"

        der = []
        for i in range(len(NEWV)): 
	    der.append(delta*(derivative(NEWV[i],R2.gen(4))))
   
        print "der = ", der, "\n"

        Y0 = []
        for i in range(len(der)): 
	    Y0.append(Matrix(k1,1,1))
        for i in range(len(der)): 
	    Y0[i][0,0] = der[i]

        print "Y0 = ", Y0, "\n"

        Y1 = []
        for i in range(len(Y0)): 
	    Y1.append(Y0[i][0,0].mod(q100))

        print "Y1 = ", Y1, "\n"

        Y000 = []
        for i in range(len(Y1)): 
	    Y000.append(Y1[i]//(T^(2*prec)))
   
        Y00 = []
        for i in range(len(Y1)):
	    Y00.append(Y1[i] - Y000[i]*T^(2*prec)) 

        print "Y00 = ", Y00, "\n"

        X01 = []
        for i in range(len(Y00)): 
	    X01.append(Matrix(R2,1,1))
        for i in range(len(Y00)):
	    X01[i][0,0] = Y00[i]

        Y01 = []
        for i in range(len(X01)):
	    Y01.append(X01[i][0,0])

        VERYNEWV  = []
        for i in range(len(Y01)): 
	    VERYNEWV .append(NEWV[i] - Y01[i])

        print "VERYNEWV  = ", VERYNEWV , "\n"


        del1 = Matrix(k1,1,1)
        del1[0,0] = delta
        delta1 = del1[0,0]

        print "delta1 = ", delta1, "\n"

        q01 = q100 - (delta1*derivative(q100)).mod(q100)

        q11 = q01//(T^(2*prec))
        NEWq = q01-q11*(T^(2*prec))


        NEWqq = Matrix(A,1,1)
        NEWqq[0,0] = NEWq

        q = NEWqq[0,0]
    
        print "q = ", q, "\n"
  
        V = Matrix(R2,4,1)
        for i in range(4): 
	    V[i,0] = VERYNEWV[i]
	
	prec = 2*prec
    	
	print "V = ", V, "\n"

    return V,q


def qua(Hmat,AM,a100):
    H = copy(Hmat)
    G = copy(AM)
    ll = copy(a100)

    l = []
    for i in range(len(ll)):
	l.append(ll[i][0])


    J = []      ## List of all p*q matrices over K after evaluating in G = AM
    for i in range(len(l)): 
    	J.append(G.substitute(x0 = l[i][0], x1 = l[i][1], x2 = l[i][2], x3 = l[i][3]))
	

    E = []      ## row echelon form
    for i in range(len(J)):
    	E.append(J[i].rref())



    """
    For x∗ be a solution of IG, we  construct the system of n equations from H

    The list Q: with Q[i] is the set of index, namely C, for H¯
    for each x∗=l[i]

    """



    Q = [] 
    for k in range(len(E)):
    	A = []
    	for i in range(p):
            for j in range(q):
            	if E[k][i, j] != 0:
                    A.append(j)
                    break
        
        Q.append(A)
    
    A1 = range(q)


    
    N = []
    for i in range(len(Q)):
    	N.append(diff(A1,Q[i])) 
    

    """
    For each k, the construction of C[k] contains all of the indices of p×p-minors (equations in the systems)

    C[k] for the solution l[k] of I_G
    """

    C = []
    for k in range(len(Q)): 
    	B = []
    	for i in range(len(N[k])):
            B.append(Q[k]+[N[k][i]])
        C.append(B)

    """
    Z[k] is the system of n equations for the solution l[k]
    """
    #aa = H.matrix_from_columns(C[0][i]).det()

    
    return C,N,Q,E



    Z = []
    for k in range(len(C)):
    	W = []
    	for i in range(len(C[k])):
    	    W.append(H.matrix_from_columns(C[k][i]).det())
	Z.append(W)
	    
    return C



def diff(first, second): ## compute list difference
        second = set(second)
        return [item for item in first if item not in second]












