bF = QQ

S.<T> = PolynomialRing(bF)
k.<t> = PolynomialRing(bF);
K = FractionField(k);
A.<T> = PolynomialRing(K)
#RR.<t> = PowerSeriesRing(bF)


q = T^2-1
p = t^2+T*t-1



#bF = QQ
R.<x1,x2,T,t> = PolynomialRing(bF,4)

f1 = (x1-1)^2+(x2-1)^2-4-t-t^2
f2 = (x1+1)^2 + (x2+1)^2-4-t

F = [f1,f2]
polelim = T^2-1


S.<T> = PolynomialRing(bF)
k.<t> = PolynomialRing(bF);
K = FractionField(k);
A.<T> = PolynomialRing(K)
RR.<t> = PowerSeriesRing(bF)



param = Matrix(R,2,1)
param[0,0] = T
param[1,0] = -T

 

u = x1

def HenselLift(F, param, polelim, u, myprec): 
    
    prec = 1
    V = param
   # print prec


    while prec <= myprec:

    	print "prec = ", prec, "\n"
	
        J = jacobian(F,(x1,x2))
        JJ = J.substitute(x1 = V[0,0], x2 = V[1,0])
        d = JJ.det()

	#print d

        dd1 = Matrix(A,1,1)
        dd1[0,0] = d
        dd = dd1[0,0]
	
        
        qq = Matrix(A,1,1)
        qq[0,0] = polelim
        q = qq[0,0]
	
  
        sJ = JJ.inverse()
 	
	print "sJ = ", sJ, "\n"

        #sF = Matrix(R,2,1)
        #sF[0,0] =F[0].substitute(x1 = V[0,0], x2 = V[1,0])
        #sF[1,0] =F[1].substitute(x1 = V[0,0], x2 = V[1,0])

	sF = Matrix(R,1,2)
        sF[0,0] =F[0].substitute(x1 = V[0,0], x2 = V[1,0])
        sF[0,1] =F[1].substitute(x1 = V[0,0], x2 = V[1,0])

	print "sF = ", sF, "\n"
	
        NEW = sJ*(sF.transpose())

	print "NEW = ", NEW, "\n"

 
        dd2 = NEW[0,0].denominator()

        dd3 = Matrix(A,1,1)
        dd3[0,0] = dd2
        dd = dd3[0,0]
	
	s0 = gc(dd,q,2*prec)
	
	s01 = Matrix(R,1,1)
	s01[0,0] = s0
	s = s01[0,0]
       
	print "s = ", s, "\n"

	ddd2 = NEW[1,0].denominator()

        ddd3 = Matrix(A,1,1)
        ddd3[0,0] = ddd2
        ddd = ddd3[0,0]
	
	s1 = gc(ddd,q,2*prec)

	s32 = Matrix(R,1,1)
	s32[0,0] = s1
	s3 = s32[0,0]
	
	print "s3 = ", s3, "\n"

        NEWV1 = [ V[0,0] - s*NEW[0,0].numerator(), V[1,0] - s3*NEW[1,0].numerator() ] 

	print "NEWV1 = ", NEWV1, "\n"

        NEWV = [ (NEWV1[0].mod(q)).mod(t^(2*prec)), (NEWV1[1].mod(q)).mod(t^(2*prec)) ] 
	
	print "NEWV = ", NEWV, "\n"
	#print prec	

        B = Matrix(R,1,1)
        B[0,0] = T
        delta = u.substitute(x1 = NEWV[0], x2 = NEWV[1]) - B[0,0]
	print "delta = ", delta, "\n"
        
	VERYNEWV = [ (NEWV[0] - (delta*derivative(V[0,0],R.gen(2))).mod(q)).mod(t^(2*prec)), (NEWV[1] - (delta*derivative(V[1,0],R.gen(2))).mod(q)).mod(t^(2*prec)) ]
	
	print "VERYNEWV = ", VERYNEWV, "\n"
        q1 = Matrix(R,1,1)
        q1[0,0] = q
        q = q1[0,0]

	NEWq = (q - (delta*derivative(q,R.gen(2))).mod(q)).mod(t^(2*prec))

        V[0,0] = VERYNEWV[0]
        V[1,0] = VERYNEWV[1]
   
        q = NEWq
        prec = 2*prec

	print "V = ", V, "\n"
	print "q = ", q, "\n"

    return V,q


#c = HenselLift(F, param, polelim, u, 8)
def gc(p,q,prec):
    RR.<t> = PowerSeriesRing(bF)
    (g,s,w) = xgcd(p,q)
    d = s.degree()
    L = s.coefficients()
    a = []
    for i in range(len(L)): 
        a.append((L[i].subs(t=t).O(prec)).truncate(prec))

    m = sum(a[i]*T^(d-i) for i in range(len(a)))

    return m











