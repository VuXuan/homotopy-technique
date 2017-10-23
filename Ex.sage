bF = QQ
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
    while prec <= myprec:
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
	

        (g,s,w) = xgcd(dd,q) ## s*dd+w*qq = g
	#c = s.lc()
	#return c

	de = s.degree()
	c1 = diff( s, T, de)/(de.factorial())
	a = Matrix(K,1,1)
	a[0,0] = c1
	c = a[0,0]

        ss = c.subs(t = t).O(2*prec)
        s = ss.truncate(2*prec)*T^de

	#print s
  
        sJ = JJ.inverse()
 	
	#print sJ

        sF = Matrix(R,2,1)
        sF[0,0] =F[0].substitute(x1 = V[0,0], x2 = V[1,0])
        sF[1,0] =F[1].substitute(x1 = V[0,0], x2 = V[1,0])
	
        NEW = sJ*sF

	#return NEW

 
        dd2 = NEW[0,0].denominator()

        dd3 = Matrix(A,1,1)
        dd3[0,0] = dd2
        dd = dd3[0,0]

        (g,s,w) = xgcd(dd,q)

	#print s

 	de = s.degree()
	c1 = diff( s, T, de)/(de.factorial())
	a = Matrix(K,1,1)
	a[0,0] = c1
	c = a[0,0]

        ss = c.subs(t = t).O(2*prec)
        s0 = ss.truncate(2*prec)*T^de
	
	s01 = Matrix(R,1,1)
	s01[0,0] = s0
	s = s01[0,0]
       

	ddd2 = NEW[1,0].denominator()

        ddd3 = Matrix(A,1,1)
        ddd3[0,0] = ddd2
        ddd = ddd3[0,0]

        (g,s1,w) = xgcd(ddd,q)

	de1 = s1.degree()
	c11 = diff( s1, T, de1)/(de1.factorial())
	a1 = Matrix(K,1,1)
	a1[0,0] = c11
	cc = a1[0,0]

        sss = cc.subs(t = t).O(2*prec)
        s33 = sss.truncate(2*prec)*T^de1

	s32 = Matrix(R,1,1)
	s32[0,0] = s33
	s3 = s32[0,0]
 

        NEWV1 = [ V[0,0] - s*NEW[0,0].numerator(), V[1,0] - s3*NEW[1,0].numerator() ] 

	#print NEWV

        NEWV = [ (NEWV1[0].mod(q)).mod(t^(2*prec)), (NEWV1[1].mod(q)).mod(t^(2*prec)) ] 
	

	#print NEWV	

        B = Matrix(R,1,1)
        B[0,0] = T
        delta = u.substitute(x1 = NEWV[0], x2 = NEWV[1]) - B[0,0]

        VERYNEWV = [ (NEWV[0] - (delta*derivative(V[0,0],R.gen(2))).mod(q)).mod(t^(2*prec)), (NEWV[1] - (delta*derivative(V[1,0],R.gen(2))).mod(q)).mod(t^(2*prec)) ]
	
        q1 = Matrix(R,1,1)
        q1[0,0] = q
        q = q1[0,0]

	NEWq = (q - (delta*derivative(q,R.gen(2))).mod(q)).mod(t^(2*prec))

        V[0,0] = VERYNEWV[0]
        V[1,0] = VERYNEWV[1]
   
        q = NEWq
        prec = 2*prec

    return NEWV





