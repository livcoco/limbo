#!/usr/bin/env python3
# jiw - 29 Feb 2020

# This is rough code for testing generation of geodesic dome point
# coordinates and edges.  A z limit parameter allows rough control
# over how much of the sphere is generated.  The orientation of
# the icosahedron needs to have an orientation control but doesn't.
# At the moment, frequency and z limit are hardcoded, at about line
# 72. To run an example, try eg:
#       ./makeIcosoGeo.py > t1-v; ./pypeVue.py f=t1-v

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
    rbp = re = pn; rb = pn
    for ro in range(1,1+k):
        a -= 1;  b += 1; pn=genPoint(a,b,c)
        e,f = b,c
        for co in range(ro):
            e -= 1;  f += 1; pn=genPoint(a,e,f)
        rbp = rb; re = pn-1; rb = pn
    
def dedupClip(layi, layo, clip1, clip2):
    '''Given list of points via layi, return de-duplicated and clipped
    list etc'''
    xlo, xhi = min(clip1.x, clip2.x), max(clip1.x, clip2.x)
    ylo, yhi = min(clip1.y, clip2.y), max(clip1.y, clip2.y)
    zlo, zhi = min(clip1.z, clip2.z), max(clip1.z, clip2.z)
    L  = layi.posts
    pprev = Point(9e9, 8e8, 7e7);  eps = 0.001
    for n, p in enumerate(L): p.dex = n
    
    transi = {} # Make node-number translation table for merged points
    for p in sorted(L):
        me = p.dex; del p.dex
        if xlo<=p.x<=xhi and ylo<=p.y<=yhi and zlo<=p.z<=zhi:
            if ssq(*(p.diff(pprev))) > eps:
                layo.posts.append(p)
            transi[me] = len(layo.posts)-1
        else:          # p is out-of-box; remove its edge evidence
            nbrs = layi.edgeList[me]
            for nbr in nbrs:    # Get rid of all refs to me
                layi.edgeList[nbr].remove(me)
            del layi.edgeList[me] # Get rid of me
        pprev = p

    # Translate all edge numbers in the edgeList to merge points
    el = layi.edgeList
    for i in el.keys():
        v = transi[i]
        for j in el[i]:
            w = transi[j]
            addEdges(v, w, layo)

def genIcosahedron(layin, Vfreq, zMin):
    '''Generate points and edges for icosahedral faces having corners
    above the plane z=zMin, with faces triangulated at frequency Vfreq.
    Ref: "Geodesic Domes", by Tom Davis - a pdf file - pp. 5-10 '''
    phi = (1+sqrt(5))/2
    cornerNote = 'oip ojp ojq oiq  poi qoi qoj poj  ipo jpo jqo iqo'
    facesNote = 'aij ajf afb abe aei bfk bkl ble cdh chl clk ckg cgd dgj dji dih dji dih elh ehi fjg fgk'
    corr1 = {'o':0, 'i':1, 'j':-1, 'p':phi, 'q':-phi}
    corners = [Point(corr1[i], corr1[j], corr1[k]) for i,j,k in cornerNote.split()]
    oa = ord('a')
    # Init an empty layout for local use (ie before deduplication)    
    laylo = Layout(posts=[], cyls=[],  edgeList={})
    for i,j,k in facesNote.split():
        p,q,r = corners[ord(i)-oa], corners[ord(j)-oa], corners[ord(k)-oa]
        pn = len(laylo.posts)
        if min(p.z, q.z, r.z) >= zMin:
            genTriangleK (laylo, Vfreq, p, q, r, pn)
            print (f'=   {len(laylo.posts):3} posts after face {i}{j}{k}')
    # Dedup laylo and copy points into layin
    dedupClip(laylo, layin, Point(-.23,-1,-.23), Point(1,1,1))
    print (f'=   {len(layin.posts):3} posts after dedup')

for Vfreq in (5,):
    zMin = 0.5
    zMin = -1.1
    print (f'=  Vfreq {Vfreq},  zMin {zMin}')
    LO = Layout(posts=[], cyls=[],  edgeList={})    # Init an empty layout
    genIcosahedron(LO, Vfreq, zMin)

    print (f'=  Writing {len(LO.posts)} post coordinates')
    print ('=P  endGap=0 postAxial=f postLabel=f  pDiam=.01  endGap=0  postHi=.02 postDiam=.01 ')
    print ('=L O 0,0,0; C ', end='')
    for p in LO.posts:
        print (f'  {p.x:0.4f},{p.y:0.4f},{p.z:0.4}', end='')
    print (';\n=C  Mpaa')
    print ("=A gg['endGap']=0")
    out = 0
    for j in LO.edgeList.keys():
        for k in LO.edgeList[j]:
            print (f' {j:2} {k:2};', end='')
            out += 1
            if out%11 == 0: print()
    print()
    print (f'=  Wrote {out} cylinders')
