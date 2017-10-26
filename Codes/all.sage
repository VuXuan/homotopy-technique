load ("data.sage")
load ("2*31.sage")
load ("start2*31.sage")
load ("Hensel2.sage")
load ("liftted.sage")
load ("check.sage")



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

testq = teststart(dd, ee, zz, G)

testq1 = teststart(dd, ee, q1, G)


testq2 = teststart(dd, ee, q2, G)

print "testq = ", testq.minors(2), "\n"

print "testq1 = ", testq1.minors(2), "\n"

print "testq2 = ", testq2.minors(2), "\n"




##################################################################
#### Test CRT

###

d3 = vv3[0,0]

e3 = vv3[1,0]

te3 = teststart(d3, e3, q3, G)

print "te3 = ", te3.minors(2), "\n"

d123 = crmth2(dd,zz,d3,q3)

e123 = crmth2(ee, zz, e3, q3)


res1 = teststart(d123, e123, q1, G)

print "res1 = ", res1.minors(2), "\n"

res2 = teststart(d123, e123, q2, G)

print "res2 = ", res2.minors(2), "\n"

res3 = teststart(d123, e123, q3, G)

print "res3 = ", res3.minors(2), "\n"