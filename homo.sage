IntegerRing()
bF = QQ
var('n')
q = 4
p = 2
n = q - p + 1
R = PolynomialRing(bF, n, 'x')
S = PolynomialRing(bF, n+1, 'x')
D = [1,1]
AC = bF.algebraic_closure()
T = PolynomialRing(AC,1, "z")



def randpoly(d): ## Random polynomial in R
   return add([add([bF.random_element()*R({tuple(a):1}) for a in WeightedIntegerVectors(i,[1 for i in range(R.ngens())])]) for i in range(d+1)])

def randmat(D): ## Input matrix F
    mat = Matrix(R, q, p, [[randpoly(D[i]) for i in range(p)] for j in range(q)])
    F = transpose(mat)
    return F
    
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
    
def G1(R, D, q): ## $p*q$ polynomial matrix with the entry (i,j) is product of $D_i$ linear form with generic coefficients
    p = len(D)
    A = Matrix(R, p, q)
    for i in range(p):
        for j in range(q):
           A[i,j] = multi(R,D[i])
    return A

def G2(R, D, q):  ## p∗q polynomial matrix with the entry (i,j) is product of Di linear form with generic coefficients, first block is a diag matrix
    p = len(D)
    B = Matrix(R, p, q)
    for i in range(p):
        for j in range(q):
            if j == i:
                B[i,j] = multi(R,D[i]) 
            elif j > p-1:
                B[i,j] = multi(R,D[i]) 
            else:
                B[i,j] = 0
    return B
 
def diff(first, second): ## compute list difference
        second = set(second)
        return [item for item in first if item not in second]

G = G2(R,D,q)    
m=ideal(G.minors(2)) ## TODO: 
l = m.variety()      ## TODO:

def var1(A,R,p,q): ## g_11 = 0 or g_22 = 0
    x0 = R.gen(0)
    G = copy(A)
    ll = []
    for i in range(p):
        ll.append(x0 - 1/G[i,i].coefficients()[0]*G[i,i])
    for j in range(p,q):
        for i in range(p):
            G[i, j] = G[i,j].substitute(x0 = ll[i])
    l2 = []
    for i in range(p):
        l2.append(ideal(G[i,2], G[i,3], x0 - ll[i]).variety())
    return l2
  
list1 = var1(G,R,p,q)

def var2(A,R,p,q): ## g_11 = g_22 = 0
    (x0, x1, x2) = R.gens()
    z = T.gen()
    G = copy(A)
    m = x0 - 1/G[0,0].coefficients()[0]*G[0,0]
    G[1,1] = G[1,1].substitute(x0 = m)
    x12 = x1 - 1/G[1,1].coefficients()[0]*G[1,1]
    x02 = m.substitute(x1 = x12)
    for j in range(p,q):
        for i in range(p):
            G[i, j] = G[i,j].substitute(x0 = x02, x1 = x12)
    g = G.submatrix(0,2,2,2).det().substitute(x2  = z)
    t = ideal(g).variety()
    l4 = []
    l31 = [t[0][T.gen()], x12.substitute(x2 = t[0][T.gen()]),  x02.substitute(x2 =  t[0][T.gen()]),] 
    d1 = dict()
    for i in range(len(l31)):
        d1[R.gen(i)] = l31[i]
    l32 = [t[1][T.gen()], x12.substitute(x2 = t[1][T.gen()]),  x02.substitute(x2 =  t[1][T.gen()]),] 
    d2 = dict()
    for i in range(len(l32)):
        d2[R.gen(i)] = l32[i]
    return [[d2]] + [[d1]] 

list2 = var2(G,R,p,q)
    
F = randmat(D)
H = (1-S.gens()[n])*G + S.gens()[n]*F ## Homotopy


"""
Find the indices as in \bar{G}
"""

J = []      ## List of all p*q matrices over K after evaluating in G
for i in range(len(l)): 
    J.append(G.subs(l[i]))
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
    A
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

Z = []
for k in range(len(C)):
    W = []
    for i in range(len(C[k])):
        W.append(H.matrix_from_columns(C[k][i]).det())
    Z.append(W)   
































