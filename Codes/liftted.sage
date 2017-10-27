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


##################################################################





def CRT(H, prec, prec2):
    a1 = ch(AA,S,D,ld1,G)
    H1 = chek1(H,GG,EE,a1)[0]
    
    v1,q1 = HenselLift(H1,a1[0],a1[1],u,prec)
    
    a2 = ch2(AA,S,D,ld1,G)
    H2 = chek1(H,GG,EE,a2)[0]

    v2,q2 = HenselLift(H2,a2[0], a2[1], u, prec)

    a3 = ch3(AA,S,D,ld1,G)
    H3 = chek1(H,GG,EE,a3)[0]
    

    v3,q3 = HenselLift(H3,a3[0],a3[1], u,prec)


    vv1 = emph(v1,q1)
    vv2 = emph(v2,q2)
    vv3 = emph(v3,q3)


    aa,nq = com1(vv1,q1,vv2,q2,vv3,q3,prec2)
    
    return aa,nq



def emph(v1,q1):
    v11 = Matrix(parent(q1),2,1)
    for i in range(2):
    	v11[i,0] = v1[i,0]	
    return v11



########## Chinese Remainder Theorem
## Input v1, v2
## Output v s.t v = v1 mod q1, v = v2 = q1


def crmth2(v1,q1,v2,q2):
    (s,n1,n2) = xgcd(q1,q2)
    v = (v1*n2*q2 + v2*n1*q1).mod(q1*q2)
    return v


### Input: vv1, vv2, vv3, q1, q1, q3
### Output: the union of those components

def CRT2(H, prec):
    a1 = ch(AA,S,D,ld1,G)
    H1 = chek1(H,GG,EE,a1)[0]

    v1,q1 = HenselLift(H1,a1[0],a1[1],u,prec)
    
    a2 = ch2(AA,S,D,ld1,G)
    H2 = chek1(H,GG,EE,a2)[0]

    v2,q2 = HenselLift(H2,a2[0], a2[1], u, prec)

    a3 = ch3(AA,S,D,ld1,G)
    H3 = chek1(H,GG,EE,a3)[0]
    

    v3,q3 = HenselLift(H3,a3[0],a3[1], u,prec)


    vv1 = emph(v1,q1)
    vv2 = emph(v2,q2)
    vv3 = emph(v3,q3)
    d1 = vv1[0,0]
    d2 = vv2[0,0]

    dd = crmth2(d1, q1, d2, q2)

    e1 = vv1[1,0]
    e2 = vv2[1,0]

    ee = crmth2(e1, q1, e2, q2)


    d3 = vv3[0,0]

    e3 = vv3[1,0]

    d123 = crmth2(dd,zz,d3,q3)
    
    e123 = crmth2(ee, zz, e3, q3)

    q123 = q1*q2*q3
    return d123, e123, q123
###############

######################################################################################"""
##### Solutions x2, x3, q have coefficients in polynomials of T of degrees at most 2e

def rat(sol, prec1): ##### should take prec1 = 40
    P.<T> = PowerSeriesRing(bF,default_prec = 64)
    csol = sol.coefficients(sparse = False)
    dsol = sol.degree()

    lsol = []
    for i in range(len(csol)):
    	lsol.append((csol[i].subs(T=T).O(prec1)).truncate(prec1))

    nsol = sum(lsol[i]*t^(i) for i in range(len(lsol)))
    return nsol


#ratv123 = rat(v123, 40)
#ratw123 = rat(w123, 40)
#ratq123 = rat(q123, 40)

###############################
### Rational reconstruction at prec = e, e = 19


def rat_recons1(nsol, e):
    P.<T> = PowerSeriesRing(bF,default_prec = 64)

    soll = Matrix(A1,1,1)
    soll[0,0] = nsol

    sol = soll[0,0]

    m = sol.coefficients(sparse = False)

    m1 = []
    for i in range(len(m)):
    	m1.append(Matrix(k,1,1))

    for i in range(len(m1)):
    	m1[i][0,0] = m[i]

    m2 = []
    for i in range(len(m1)):
    	m2.append(m1[i][0,0])
    return m2

def rat_recons2(m2, e):
    m4 = []
    for i in range(len(m2) - 1):
    	m4.append(m2[i].rational_reconstruct(T^(2*e+1),e,e))
 
    m51 = []
    for i in range(len(m4)):
	m51.append(m4[i][0])

    m52 = []
    for i in range(len(m4)):
	m52.append(m4[i][1])

    m5 = []
    for i in range(len(m52)):
    	m5.append(m51[i]/m52[i])
 
    nc = sum(m5[i]*t^(i) for i in range(len(m5)))
    return nc


def rat_recons(sol, e):
    m = rat_recons1(sol, e)
    mm = rat_recons2(m, e)
    return mm

##############################


def c1(v1,q1,v2,q2,prec):
    P.<T> = PowerSeriesRing(bF,default_prec = 64)
    n1,n2 = gmm(q1,q2,prec)
    cc = v2*n1*q1 + v1*n2*q
    q12 = q1*q2
    c = cc.mod(q12)
    
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
    
    ww = Matrix(parent(q1), 2,1)
    ww[0,0] = c1(vv1[0,0],q1,vv2[0,0],q2, prec)
    ww[1,0] = c1(vv1[1,0],q1,vv2[1,0],q2, prec)

    q12 = q1*q2
    

    dq = q12.degree()
    Lq = q12.coefficients(sparse = False)

    aq = []
    for i in range(len(Lq)):
    	aq.append((Lq[i].subs(T=T).O(prec)).truncate(prec))

    nq = sum(aq[i]*t^(i) for i in range(len(aq)))
    return  ww,nq

    

def com1(vv1,q1,vv2,q2,vv3,q3,prec):
    c12,q12 =  crth(vv1,q1,vv2,q2,prec)
    c,qq = crth(c12,q12,vv3,q3,prec)
    return c,qq


def gmm(q1,q2,prec):
    RR.<T> = PowerSeriesRing(bF, default_prec = 64)
    (g,m1,m2) = xgcd(q1,q2)
    
    d1 = m1.degree()
    L1 = m1.coefficients(sparse=False)
    a11 = []
    for i in range(len(L1)): 
        a11.append((L1[i].subs(T=T).O(prec)).truncate(prec))

    n1 = sum(a11[i]*t^(i) for i in range(len(a11)))

    d2 = m2.degree()
    L2 = m2.coefficients(sparse=False)
    a2 = []
    for i in range(len(L2)):
    	a2.append((L2[i].subs(T=T).O(prec)).truncate(prec))

    n2 = sum(a2[i]*t^(i) for i in range(len(a2)))
    
    return n1,n2


## rational reconstruction for coefficients of q and v1

def ratrec1(solv,e): ## for q
    P.<T> = PowerSeriesRing(bF, default_prec = 64)

    soll = Matrix(A1,1,1)
    soll[0,0] = solv

    sol = soll[0,0]

    m = sol.coefficients(sparse = False)

    m1 = []
    for i in range(len(m)):
    	m1.append(Matrix(k,1,1))

    for i in range(len(m1)):
    	m1[i][0,0] = m[i]

    m2 = []
    for i in range(len(m1)):
    	m2.append(m1[i][0,0])

    return m2
    #m3 = []
    #for i in range(len(m2)):
    	#m3.append( m2[i].substitute(T = P.gen()))

def ratrec2(m2,e):
    m4 = []
    for i in range(len(m2) - 1):
    	m4.append(m2[i].rational_reconstruct(T^(2*e+1),e,e))
 
    a = Matrix(k,1,1)
    a[0,0] = 1
    
    m4.append(a[0,0].rational_reconstruct(T^(2*e+1),e,e))

    
    return m4


#def ratrecc(nq,)


def ratrec(nq,e):
    m2 = ratrec1(nq,e)
    m4 = ratrec2(m2,e)
    m51 = []
    for i in range(len(m4)):
	m51.append(m4[i][0].substitute(T=1))

    m52 = []
    for i in range(len(m4)):
	m52.append(m4[i][1].substitute(T=1))

    m5 = []
    for i in range(len(m4)):
	m5.append(m51[i]/m52[i])
    
    nc = sum(m5[i]*t^(i) for i in range(len(m5)))
    
    return nc
  


########################################################################


def ratrecvv(solv,e): ## for v1
    P.<T> = PowerSeriesRing(bF,default_prec = 64)

    soll = Matrix(A1,1,1)
    soll[0,0] = solv

    sol = soll[0,0]
    
    m = sol.coefficients(sparse = False)
    
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



def ratrecv1(solv,e):
    m4 = ratrecvv(solv,e)

    #mc4 = sum(m4[i]*t^(i) for i in range(len(m4)))

    #return mc4
    
    m51 = []
    for i in range(len(m4)):
	m51.append(m4[i].numerator())

    m52 = []
    for i in range(len(m4)):
	m52.append(m4[i].denominator())
	
    #m5 = []
   # for i in range(len(m5)):
    #	m5.append(m51[i]/m52[i])

    #mc = sum(m5[i]*t^(i) for i in range(len(m5)))
    #return mc
    
    m61 = []
    for i in range(len(m51)):
	m61.append(m51[i].substitute(T=1))

    m62 = []
    for i in range(len(m52)):
	m62.append(m52[i].substitute(T=1))

    m6 = []
    for i in range(len(m61)):
	m6.append(m61[i]/m62[i])

    nc = sum(m6[i]*t^(i) for i in range(len(m6)))
    return nc


#####################################################################
### reconstruct a zero-dim para S with coeffs in K(T), deg(num) <= e, deg(denum) <=e ###




def ratrecvv2(m4, e):
    mv = sum(m4[i]*t^i for i in range(len(m4)))
    return mv
    
aa, nq = CRT(H,32, 40)

def par(aa,nq):
    dq = nq.derivative()
    v = (aa[0,0]*dq).mod(nq)
    w = (aa[1,0]*dq).mod(nq)
    return v,w

vv, ww =  par(aa,nq)

qq = ratrec(nq, 19)
aa1 = ratrecv1(vv, 19)
aa2 = ratrecv1(ww, 19)

qdev = qq.derivative()

#Matt = Matrix(AA,3,1)
#Matt[0,0] = qq
#Matt[1,0] = aa1
#Matt[2,0] = aa2

#qfin = Matt[0,0]
#vfin = Matt[1,0]
#wfin = Matt[2,0]


qfin = qq
vfin = aa1/qdev
wfin = aa2/qdev





