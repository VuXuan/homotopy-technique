
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

####################################################################


def HenselLift(F, param, polelim, u, myprec): 
    
    prec = 1
    V = param
    #print "V = ", V, "\n"

    qq = Matrix(A1,1,1)
    qq[0,0] = polelim
    q = qq[0,0]
    
    #print "q = ", q, "\n"

    while prec <= myprec:

    #	print "prec = ", prec, "\n"
	
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

#	print "V = ", V, "\n"
#	print "q = ", q, "\n"
	
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
