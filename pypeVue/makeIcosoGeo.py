#!/usr/bin/env python3

from pypeVue import Point, Layout, ssq, sssq, addEdges
from math import sqrt

def genTriangleK (layout, k, v0, v1, v2, pn):
    def genPoint(p, q, r):
        x = (p*v0.x + q*v1.x + r*v2.x)/k
        y = (p*v0.y + q*v1.y + r*v2.y)/k
        z = (p*v0.z + q*v1.z + r*v2.z)/k
        t = sssq(x,y,z)         # t = distance to origin
        layout.posts.append(Point(x/t,y/t,z/t))
        if pn-ro <= re:  addEdges(pn, pn-ro,   layout)
        
        if pn-ro > rbp:  addEdges(pn, pn-ro-1, layout)
        
        if pn    >  rb:  addEdges(pn, pn-1,    layout)
        return pn+1
    
    rbp = rb = pn+2; re = -1;  ro=0
    a,b,c = k, 0, 0;     pn=genPoint(a,b,c)
    rbp = re = pn; rb = pn+1
    for ro in range(1,1+k):
        a -= 1;  b += 1; pn=genPoint(a,b,c)
        e,f = b,c
        for co in range(ro):
            e -= 1;  f += 1; pn=genPoint(a,e,f)
        rbp = rb; re = pn-1; rb = pn
    
def dedup(layi, layo):
    '''Given list of points via layi, return de-duplicated list etc'''
    L  = layi.posts
    pprev = Point(9e9, 8e8, 7e7);  eps = 0.001
    for n, p in enumerate(L): p.dex = n
    
    transi = {}
    for p in sorted(L):
        if ssq(*(p.diff(pprev))) > eps:
            layo.posts.append(p)
        transi[p.dex] = len(layo.posts)-1
        pprev = p

    print ('edges in:', sorted([(k,layi.edgeList[k]) for k in layi.edgeList.keys()]))
    print ('transi: ', sorted([(k,transi[k]) for k in transi.keys()]))
    el = layi.edgeList
    layo.edgeList = {}          # debug!!
    for i in el.keys():
        v = transi[i]
        print (f'{i:4} -> {v} for {el[i]}')
        for j in el[i]:
            w = transi[j]
            addEdges(v, w, layo)
    print ('edges out:', sorted([(len(layo.edgeList[k]), k, layo.edgeList[k]) for k in layo.edgeList.keys()]))

def genIcosahedron(layin, Vfreq, zMin):
    '''Generate points and edges for icosahedral faces having corners
    above the plane z=zMin, with faces triangulated at frequency Vfreq '''
    phi = (1+sqrt(5))/2
    icoCornerNote = 'oip ojp ojq oiq  poi qoi qoj poj  ipo jpo jqo iqo'
    facesNote = 'aij ajf afb abe aei bfk bkl ble cdh chl clk ckg cgd dgj dji dih dji dih elh ehi fjg fgk'
    corr1 = {'o':0, 'i':1, 'j':-1, 'p':phi, 'q':-phi}
    corners = [Point(corr1[i], corr1[j], corr1[k]) for i,j,k in icoCornerNote.split()]
    oa = ord('a')
    # Init an empty layout for local use (ie before deduplication)    
    laylo = Layout(posts=[], cyls=[],  edgeList={})
    for i,j,k in facesNote.split():
        p,q,r = corners[ord(i)-oa], corners[ord(j)-oa], corners[ord(k)-oa]
        pn = len(laylo.posts)
        if min(p.z, q.z, r.z) >= zMin:
            genTriangleK (laylo, Vfreq, p, q, r, pn)
            print (f'Have {len(laylo.posts):3} posts after face {i}{j}{k}')
    dedup(laylo, layin)     # Dedup laylo and copy points into layin

for Vfreq in (4,):
    zMin = 0.5
    print (f'Vfreq {Vfreq},  zMin {zMin}')
    LO = Layout(posts=[], cyls=[],  edgeList={})    # Init an empty layout
    genIcosahedron(LO, Vfreq, zMin)
    for p in LO.posts:
        print (f'{p.x:11.8f} {p.y:11.8f} {p.z:11.8f}  {p.dex}')
    print (LO.edgeList)
    
phi = (1+sqrt(5))/2
v0 = Point(0,   1,   phi)
v1 = Point(1,   phi, 0  )
v2 = Point(phi, 0,   1  )
for k in (4,5):
    if k>2: continue
    LO = Layout(posts=[], cyls=[],  edgeList={})    # Init an empty layout
    print (f'\nk = {k}')
    genTriangleK (k, v0, v1, v2, 0);
    for p in sorted(LO.posts):
        print (f'{p.x:11.8f} {p.y:11.8f} {p.z:11.8f}')
    print (LO.edgeList)