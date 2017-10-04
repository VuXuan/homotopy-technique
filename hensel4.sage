bF = GF(9001)

S.<t> = PolynomialRing(bF)
k.<T> = PolynomialRing(bF);
K = FractionField(k);
A.<t> = PolynomialRing(K)


#q = t^2-1
#p = T^2+T*t-1


R2.<x0,x1,x2,x3,t,T> = PolynomialRing(bF,6)

f1 = (x0-1)^2 + (x1-1)^2+ (x2-1)^2+(x3-1)^2-4-T-T^2

f2 = (x0+1)^2 + (x1+1)^2+ (x2+1)^2 + (x3+1)^2-4-T

f3 = (2*x0 - 3)^3 + (2*x1 - 3)^3+ (2*x2-3)^3 + (2*x3)^3-4-T/2

f4 = (3*x0+2)^3 + (3*x1+2)^3+ (3*x2+2)^3 + (3*x3+2)^3 - 4-T-T/3^2

F = [f1,f2,f3,f4]
polelim = t^2-1


k1.<t> = PolynomialRing(k)
K = FractionField(k);
A.<t> = PolynomialRing(K)

RR.<T> = PowerSeriesRing(bF)


param = Matrix(R2,4,1)

param[0,0] = t+1
param[1,0] = t-1
param[2,0] = t
param[3,0] = -t

u = x2

def HenselLift4(F, param, polelim, u, myprec): 
    
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














def gc(p,q,prec):
    RR.<T> = PowerSeriesRing(bF)
    (g,s,w) = xgcd(p,q)
    d = s.degree()
    L = s.coefficients(sparse=False)
    a = []
    for i in range(len(L)): 
        a.append((L[i].subs(T=T).O(prec)).truncate(prec))

    m = sum(a[i]*t^(i) for i in range(len(a)))

    return m

