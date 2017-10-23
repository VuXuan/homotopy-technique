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

## rational reconstruction for coefficients of q and v1

def ratrec(solv,e): ## for q
    P.<T> = PowerSeriesRing(bF, default_prec = 64)

    soll = Matrix(A1,1,1)
    soll[0,0] = solv

    sol = soll[0,0]

    m = sol.coefficients()

    m1 = []
    for i in range(len(m)):
    	m1.append(Matrix(k,1,1))

    for i in range(len(m1)):
    	m1[i][0,0] = m[i]

    m2 = []
    for i in range(len(m1)):
    	m2.append(m1[i][0,0])

    m3 = []
    for i in range(len(m2)):
    	m3.append( m2[i].substitute(T = P.gen()))

    m4 = []
    for i in range(len(m3) - 1):
    	m4.append(m3[i].pade(e,e))

    a = Matrix(parent(m4[0]),1,1)
    a[0,0] = 1
    
    m4.append(a[0,0])
    return m4


def ratrecvv(solv,e): ## for v1
    P.<T> = PowerSeriesRing(bF,default_prec = 64)

    soll = Matrix(A1,1,1)
    soll[0,0] = solv

    sol = soll[0,0]
    
    m = sol.coefficients()
    
    m1 = []
    for i in range(len(m)):
    	m1.append(Matrix(k,1,1))

    for i in range(len(m1)):
    	m1[i][0,0] = m[i]

    m2 = []
    for i in range(len(m1)):
    	m2.append(m1[i][0,0])

    m3 = []
    for i in range(len(m2)):
    	m3.append( m2[i].substitute(T = P.gen()))
    
    m4 = []
    for i in range(len(m3)):
    	m4.append(m3[i].pade(e,e))
   
    return m4



##################################################################





def CRT(a1,H1,a2,H2,a3,H3,prec):
    v1,q1 = HenselLift(H1,a1[0],a1[1],u,prec)
    v2,q2 = HenselLift(H2,a2[0], a2[1], u, prec)
    v3,q3 = HenselLift(H3,a3[0],a3[1], u,prec)


    vv1 = emph(v1,q1)
    vv2 = emph(v2,q2)
    vv3 = emph(v3,q3)


    aa,nq = com1(vv1,q1,vv2,q2,vv3,q3,prec)
    
    return aa,nq



def emph(v1,q1):
    v11 = Matrix(parent(q1),2,1)
    for i in range(2):
    	v11[i,0] = v1[i,0]	
    return v11




def c1(v1,q1,v2,q2, prec):
    P.<T> = PowerSeriesRing(bF,default_prec = 64)
    n1,n2 = gmm(q1,q2,prec)
 
    c = v2*n1*q1 + v1*n2*q2

    dc = c.degree()
    lc1 = c.coefficients(sparse = False)

    ac = []
    for i in range(len(lc1)):
	ac.append(Matrix(k,1,1))
    for i in range(len(lc1)):
	ac[i][0,0] = lc1[i]

    ac1 = []
    for i in range(len(ac)):
	ac1.append(ac[i][0,0].substitute(T = P.gen()))

    ac2 = []
    for i in range(len(ac1)):
	ac2.append((ac1[i].O(prec)).truncate(prec))

    nc = sum(ac2[i]*t^(i) for i in range(len(ac2)))
    return nc



def crth(vv1,q1,vv2,q2,prec):
    RR.<T> = PowerSeriesRing(bF, default_prec = 64)
    
    vv = Matrix(parent(q1), 2,1)
    vv[0,0] = c1(vv1[0,0],q1,vv2[0,0],q2, prec)
    vv[1,0] = c1(vv1[1,0],q1,vv2[1,0],q2, prec)

    q12 = q1*q2
    

    dq = q12.degree()
    Lq = q12.coefficients(sparse = False)

    aq = []
    for i in range(len(Lq)):
    	aq.append((Lq[i].subs(T=T).O(prec)).truncate(prec))

    nq = sum(aq[i]*t^(i) for i in range(len(aq)))
    return  vv,nq

    

def com1(vv1,q1,vv2,q2,vv3,q3,prec):
    c12,q12 =  crth(vv1,q1,vv2,q2,prec)
    c,qq = crth(c12,q12,vv3,q3,prec)
    return c,qq


def gmm(q1,q2,prec):
    RR.<T> = PowerSeriesRing(bF, default_prec = 64)
    (g,m1,m2) = xgcd(q1,q2)
    
    d1 = m1.degree()
    L1 = m1.coefficients(sparse=False)
    a1 = []
    for i in range(len(L1)): 
        a1.append((L1[i].subs(T=T).O(prec)).truncate(prec))

    n1 = sum(a1[i]*t^(i) for i in range(len(a1)))

    d2 = m2.degree()
    L2 = m2.coefficients(sparse=False)
    a2 = []
    for i in range(len(L2)):
    	a2.append((L2[i].subs(T=T).O(prec)).truncate(prec))

    n2 = sum(a2[i]*t^(i) for i in range(len(a2)))
    
    return n1,n2


#####################################################################
### reconstruct a zero-dim para S with coeffs in K(T), deg(num) <= e, deg(denum) <=e ###

def PaS(a1,H1,a2,H2,a3,H3,prec):
    aa,nq = CRT(a1,H1,a2,H2,a3,H3,prec)
    mq = ratrec(nq,prec//2)

    dq = nq.degree()

    sq = sum(mq[i]*t^(i) for i in (0..dq))
    
    mv0 = ratrecvv(aa[0,0],prec//2)

    dv0 = aa[0,0].degree()

    sv0 = sum(mv0[i]*t^(i) for i in (0..dv0))

    mv1 = ratrecvv(aa[1,0], prec//2)

    dv1 = aa[1,0].degree()

    sv1 = sum(mv1[i]*t^(i) for i in (0..dv1))
    
    mv = Matrix(parent(sq),2,1)
    mv[0,0] = sv0
    mv[1,0] = sv1
    
    return mv,sq




#######################################################################

#### Chekck###


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

#def crth(vv1,q1,vv2,q2,prec):
 #   RR.<T> = PowerSeriesRing(bF, default_prec = 64)
  #  n1,n2 = gmm(q1,q2,prec)
    
  #  c = Matrix(parent(q1),2,1)
   # c[0,0] = vv1[0,0]*n2*q2 + vv2[0,0]*n1*q1
   # c[1,0] = vv1[1,0]*n2*q2 + vv2[1,0]*n1*q1

   # dd0 = c[0,0].degree()
   # LL0 = c[0,0].coefficients(sparse = False)

   # aa0 = []
   # for i in range(len(LL0)):
    #	aa0.append((LL0[i].subs(T=T).O(prec)).truncate(prec))

   # nn0 = sum(aa0[i]*t^(i) for i in range(len(aa0)))

   # dd1 = c[1,0].degree()
   # LL1 = c[1,0].coefficients(sparse = False)

   # aa1 = []
   # for i in range(len(LL1)):
   # 	aa1.append((LL1[i].subs(T=T).O(prec)).truncate(prec))

 #   nn1 = sum(aa1[i]*t^(i) for i in range(len(aa1)))

  #  aa = Matrix(parent(q1),2,1)
   # aa[0,0] = nn0
  #  aa[1,0] = nn1

    #q12 = q1*q2
    

   # dq = q12.degree()
   # Lq = q12.coefficients(sparse = False)

#    aq = []
 #   for i in range(len(Lq)):
  #  	aq.append((Lq[i].subs(T=T).O(prec)).truncate(prec))

  #  nq = sum(aq[i]*t^(i) for i in range(len(aq)))
  #  return  aa,nq
