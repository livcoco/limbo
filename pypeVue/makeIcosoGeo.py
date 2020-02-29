#!/usr/bin/env python3

from pypeVue import Point, Layout, sssq, addEdges
from math import sqrt

def genTriangleK (k, v0, v1, v2, pn):
    def tellstuff(way, b):
        pass
        #print (f'Add{way} {pn:2} {b:2}    ro {ro:2}   rbp {rbp:2}   rb {rb:2}   re {re:2}')
    def genPoint(p, q, r):
        x = (p*v0.x + q*v1.x + r*v2.x)/k
        y = (p*v0.y + q*v1.y + r*v2.y)/k
        z = (p*v0.z + q*v1.z + r*v2.z)/k
        t = sssq(x,y,z)         # t = distance to origin
        LO.posts.append(Point(x/t,y/t,z/t))
        if pn-ro <= re:  addEdges(pn, pn-ro,   LO); tellstuff('/',  pn-ro)
        
        if pn-ro > rbp:  addEdges(pn, pn-ro-1, LO); tellstuff('\\', pn-ro-1)
        
        if pn    >  rb:  addEdges(pn, pn-1,    LO); tellstuff('-',  pn-1)
        return pn+1
    
    rbp = rb = 2; re = -1;  ro=0
    a,b,c = k, 0, 0;     pn=genPoint(a,b,c)
    rbp = re = 0; rb = 1
    for ro in range(1,1+k):
        a -= 1;  b += 1; pn=genPoint(a,b,c)
        e,f = b,c
        for co in range(ro):
            e -= 1;  f += 1; pn=genPoint(a,e,f)
        rbp = rb; re = pn-1; rb = pn

phi = (1+sqrt(5))/2
v0 = Point(0,   1,   phi)
v1 = Point(1,   phi, 0  )
v2 = Point(phi, 0,   1  )
LO = Layout()                   # Init an empty layout
k = 4

print ('k = {k}')
genTriangleK (k, v0, v1, v2, 0);
for p in sorted(LO.posts):
    print (f'{p.x:11.8f} {p.y:11.8f} {p.z:11.8f}')
print (LO.edgeList)

k=5
print ('k = {k}')
LO = Layout(posts=[], cyls=[],  edgeList={})    # Init an empty layout
genTriangleK (k, v0, v1, v2, 0);
for p in sorted(LO.posts):
    print (f'{p.x:11.8f} {p.y:11.8f} {p.z:11.8f}')
print (LO.edgeList)