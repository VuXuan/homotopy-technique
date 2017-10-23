F:=convert([(x1-1)^2+(x2-1)^2-4-t-t^2,(x1+1)^2+(x2+1)^2-4-t],Matrix);
V:=convert([T, -T],Vector); q:=T^2-1;
with(linalg):
with(LinearAlgebra):
    
HenselLift:=proc(F, param, polelim, u, myprec)
local prec, V;
prec:=1: V:=param: q:=polelim:
while prec<=myprec do
  J:=convert(jacobian(convert(F,list),[x1,x2]), Matrix);

  ls:={x1=V[1], x2=V[2]};
  dd:=det(subs(ls,evalm(J)));
  gcdex(dd,q,T,'s','w');
  lprint(s);
  
  s:=convert(series(s, t, 2*prec),polynom);
  sJ:=convert(MatrixInverse(subs(ls,J)),Matrix); sF:=subs(ls, F);
  #sJ:=MatrixScalarMultiply(sJ, s):

 # lprint(sJ);
 # lprint(Transpose(sF));
  NEW:=MatrixVectorMultiply(sJ, Transpose(sF));
  dd:=denom(NEW[1,1]):
  #(s,w):=gcdex(dd,q,T);
  gcdex(dd,q,T,'s','w');
  lprint(s);
  s:=convert(series(s, t, 2*prec),polynom);

NEWV:=[ V[1]-s*numer(NEW[1,1]), V[2]-s*numer(NEW[2,1]) ];
NEWV:=[ rem(rem(NEWV[1],q, T), t^(2*prec), t), rem(rem(NEWV[2],q,T), t^(2*prec), t) ];
#NEWV:=map(numer, NEWV);
#NEWV:=[s*NEWV[1], s*NEWV[2]];

#lprint(NEWV);

den:=denom(NEWV[1]);
#NEWV:=[rem(expand(NEWV[1]), q, T), rem(expand(NEWV[2]), q, T)];
#NEWV:=[rem(NEWV[1], t^(2*prec), t), rem(NEWV[2], t^(2*prec), t)];
delta:=subs({x1=NEWV[1], x2=NEWV[2]}, u)-T;

    
VERYNEWV:=[
    rem(NEWV[1]-rem(expand((delta*(diff(V[1], T)))), q, T),t^(2*prec),t),
    rem(NEWV[2]-rem(expand((delta*(diff(V[2], T)))), q, T),t^(2*prec),t)
    ];
NEWq:=rem(q-(rem(delta*diff(q, T),q,T)), t^(2*prec), t);


V:=convert(VERYNEWV, Vector); q:=NEWq: prec:=2*prec:
print(rem(subs({x1=V[1],x2=V[2]}, F)[1,1],q,T));
od;

return V, q;
end proc;

debug(HenselLift);
HenselLift(F, V, q, x1,2);


#qq := -1/8*T^2*t^2 - 1/2*T^3 + 1/8*T*t^2 + 1/4*T*t + 3/2*T;
#p := T^2-1;

#a:=rem(qq,p,T); 
#aa:=rem(a,t^2,t);

#q100 := -1/8*T^2*t^2 + 1/2*T^3 - 1/8*T*t^2 - 1/4*T*t - 3/2*T;

#a100 := rem(q100,p,T);
#aa100:=rem(a100,t^2,t);
