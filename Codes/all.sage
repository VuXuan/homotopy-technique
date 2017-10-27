load ("data.sage")
load ("2*31.sage")
load ("start2*31.sage")
load ("Hensel2.sage")
load ("liftted.sage")
load ("check.sage")


### Zero-dim parametrization after lifting
##### Test CRT

v1,q1 = HenselLift(H1,a1[0],a1[1],u,32)
v2,q2 = HenselLift(H2,a2[0], a2[1], u, 32)
v3,q3 = HenselLift(H3,a3[0],a3[1], u,32)


vv1 = emph(v1,q1)
vv2 = emph(v2,q2)
vv3 = emph(v3,q3)

d1 = vv1[0,0]
d2 = vv2[0,0]

dd = crmth2(d1, q1, d2, q2)

e1 = vv1[1,0]
e2 = vv2[1,0]

ee = crmth2(e1, q1, e2, q2)

zz = q1*q2

te1 = teststart(d1, e1, q1, G)


print "te1 = ", te1.minors(2), "\n"

te2 = teststart(d2, e2, q2, G)


te2.minors(2)

print "te2 = ", te2.minors(2), "\n"

#testq = teststart(dd, ee, zz, G)

testq1 = teststart(dd, ee, q1, G)


testq2 = teststart(dd, ee, q2, G)

#print "testq = ", testq.minors(2), "\n"

print "testq1 = ", testq1.minors(2), "\n"

print "testq2 = ", testq2.minors(2), "\n"




##################################################################
#### Test CRT

###

d3 = vv3[0,0]

e3 = vv3[1,0]

te3 = teststart(d3, e3, q3, G)

print "te3 = ", te3.minors(2), "\n"

d123 = crmth2(dd,zz,d3,q3) ### parametrization for x2

e123 = crmth2(ee, zz, e3, q3) ### parametrization for x3


res1 = teststart(d123, e123, q1, G)

print "res1 = ", res1.minors(2), "\n"

res2 = teststart(d123, e123, q2, G)

print "res2 = ", res2.minors(2), "\n"

res3 = teststart(d123, e123, q3, G)

print "res3 = ", res3.minors(2), "\n"

#######################################################

##### Test CRT2


v123, w123, q123 =  CRT2(H, 32)

relq1 = teststart(v123, w123, q1, G)

print "q1.testCRT2(1,2,3) = ", relq1.minors(2), "\n"

relq2 = teststart(v123, w123, q2, G)

print "q2.testCRT2(1,2,3) = ", relq2.minors(2), "\n"

relq3 = teststart(v123, w123, q3, G)

print "q3.testCRT2(1,2,3) = ", relq3.minors(2), "\n"

#####################################################################################""

### Test the result for the fuction rat(sol, prec1)


ratv123 = rat(v123, 40)
ratw123 = rat(w123, 40)
ratq123 = rat(q123, 40)

testratq1 = teststart(ratv123, ratw123, q1, G)
print "testratq1 = ", testratq1.minors(2), "\n"

testratq2 = teststart(ratv123, ratw123, q2, G)
print "testratq2 = ", testratq2.minors(2), "\n"

testratq3 = teststart(ratv123, ratw123, q3, G)
print "testratq3 = ", testratq3.minors(2), "\n"