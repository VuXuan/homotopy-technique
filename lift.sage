k.<t>=PolynomialRing(GF(9001))
K.<T,x>=PolynomialRing(GF(9001),2)


FF.<T,x0,x1,x2> = PolynomialRing(GF(9001),4)
R.<x0,x1,x2>=PolynomialRing(GF(9001),3)
g1 = (x0 + x1 + x2 + 1)*(x0 + 2*x1+ 2**2*x2 + 2**3)
g2 = (x0 + 3*x1 + 3**2*x2 + 3**3)*(x0 + 4*x1+ 4**2*x2 + 3**3)
g3 = (x0 + 5*x1 + 5**2*x2 + 5**3)
G = [g1,g2,g3]

    
def randpoly(d): 
   return add([add([GF(9001).random_element()*R({tuple(a):1}) for a in WeightedIntegerVectors(i,[1 for i in range(R.ngens())])]) for i in range(d+1)])

   
f1 = randpoly(2)
f2 = randpoly(2)
f3 = randpoly(1)
F = [f1,f2,f3]


def path(G,F):
    H = []
    for i in range(3):
        H.append((1-T)*G[i] + T*F[i])
    return H
        
H = path(G,F)


I1 = ideal(g1,g2,g3)
list = I1.variety()


ll = []
for i in range(len(list)-1):
    l2 = []
    for j in range(1,len(H)+1):
        l2.append(k(list[i][FF.gen(len(H)+1-j)]))
    ll.append(l2)  
  
      
def lift(H,ll,prec):
    RS = []
    for i in range(len(ll)):
        RS.append(lift1(H,ll[i],prec))
    return RS

def lift1(H,sol,prec):
    if prec == 1:
        return sol
    else:
        sol_half = lift1(H,sol,prec//2)
        I = quotient(parent(sol[0]),t**prec)
        m = jacobian(H, [FF.gen(i) for i in range(1,len(sol)+1)])
        dic = dict()
        dic[FF.gen(0)] = I.0
        for i in range(1, len(ll[0])+1):
            dic[FF.gen(i)] = sol_half[i-1]
        n1 = Matrix(I,len(H),len(H))
        for i in range(len(H)):
            for j in range(len(H)):
                n1[i,j] = m[i,j].subs(dic)
        n = inver(n1, sol, prec)
        a = matrix(len(H),1,[sol_half[i] for i in range(len(sol_half))])
        v = Matrix(I,len(H),1)
        for i in range(1,len(H)):
            v[i,0] = H[i].subs(dic)
        p = a - n*v
        sol_new = []
        for i in range(len(H)):
            sol_new.append(p[i][0])
        result = []
        for i in range(len(sol_new)):
            result.append(sol_new[i].lift())    
        return result 


def inver(N,sol,deg):
    k.<t>=PolynomialRing(GF(9001))
    I = quotient(parent(sol[0]),t**(deg))
    Nbar = []
    P = []
    for k in range(deg):
        Nbar.append(matrix([[N[i][j][k] for j in range(len(sol))] for i in range(len(sol))]))       
    P.append(Nbar[0].inverse())
    for i in range(1,deg):
       P.append(Nbar[0].inverse()*(-sum(Nbar[i-j]*P[j] for j in range(i))))
    return sum(P[i]*(I.0**i) for i in range(0,len(P)))
        
        
       
#n1 = matrix([f(I.0,sol_half[0],sol_half[1],sol_half[2]) for f in m])

#v = matrix(len(H),1,[H[i](I.0, sol_half[0],  sol_half[1],  sol_half[2]) for i in range(len(H))])

       
#h1 = (1-T)*g1 + T*f1
#h2 = (1-T)*g2 + T*f2
#h3 = (1-T)*g3 + T*f3
       
#ll = [[k(list[0][x2])]+[k(list[0][x1])] + [k(list[0][x0])], [k(list[1][x2])]+[k(list[1][x1])] + [k(list[1][x0])], [k(list[2][x2])]+[k(list[2][x1])] + [k(list[2][x0])]]
        
#[sol_new[0].lift(),sol_new[1].lift(),sol_new[2].lift()]     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
