#!/usr/bin/env python3
# jiw - 29 Feb 2020

# This is rough code for testing generation of geodesic dome point
# coordinates and edges.  Clipping-box parameters allow some control
# over how much of the sphere is generated.  The orientation of the
# icosahedron has two controls: z angle and y angle rotation.  At
# present, frequency, limits, and angles are hardcoded, near line 114.
# To see an example, try eg:

#       ./makeIcosaGeo.py > t1-v; ./pypeVue.py f=t1-v

from pypeVue import Point, Layout, ssq, sssq, addEdges
from math import sqrt, pi, asin, sin, cos, atan2

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

    # Triangulate the face with corners v0,v1,v2; stepping row by row
    # from point v0 toward v1 (hence decreasing a, increasing b), and
    # in each row stepping across from the v0-v1 side toward the v0-v2
    # side (hence decreasing e, increasing f).  Neighbor tracking
    # looks to the left and looks up both left and right.  rb tracks
    # beginning of current row.  rbp and re track beginning and end of
    # previous row.
    rbp = rb = pn+2; re = -1;  ro=0
    a,b,c = k, 0, 0;     pn=genPoint(a,b,c)
    rbp = re = pn; rb = pn
    for ro in range(1,1+k):
        a -= 1;  b += 1; pn=genPoint(a,b,c)
        e,f = b,c
        for co in range(ro):
            e -= 1;  f += 1; pn=genPoint(a,e,f)
        rbp = rb; re = pn-1; rb = pn

# See if point p is in box with corners given by elements of clip1, clip2
def pointInBox (p, clip1, clip2):
    xlo, xhi = min(clip1.x, clip2.x), max(clip1.x, clip2.x)
    ylo, yhi = min(clip1.y, clip2.y), max(clip1.y, clip2.y)
    zlo, zhi = min(clip1.z, clip2.z), max(clip1.z, clip2.z)
    return xlo <= p.x <= xhi and ylo <= p.y <= yhi and zlo <= p.z <=zhi

def dedupClip(layi, layo, clip1, clip2):
    '''Given list of points via layi, return de-duplicated and clipped
    list etc'''
    def CCW(p):
        # Sort points by descending z and clockwise about z
        return -round(p.z*1000)-atan2(p.y,p.x)/(2*pi)
    L  = layi.posts
    eps = 0.00001   # Good enough for freq 36, where .001 is too big
    pprev = Point(9e9, 8e8, 7e7)
    for n, p in enumerate(L): p.dex = n
    L.sort(key=CCW)
    transi = {} # Make node-number translation table for merged points
    for p in L:
        me = p.dex; del p.dex
        if pointInBox (p, clip1, clip2):
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

def genIcosahedron(layin, Vfreq, clip1, clip2, rotay, rotaz):
    '''Generate points and edges for triangulated icosahedral faces.
    Rotate basic icosahedron faces about y by ry degrees, and about z
    by rz degrees. Use genTriangleK to triangulate feasible faces at
    frequency Vfreq.  Use dedupClip to discard corners outside of box
    clip1, clip2.  Ref: "Geodesic Domes", by Tom Davis - a pdf file -
    pp. 5-10 ; and note, "The vertices of an icosahedron centered at
    the origin with an edge-length of 2 and a circumradius of sqrt(phi
    +2) ~ 1.9 are described by circular permutations of (0, ±1, ±ϕ)
    where ϕ = 1 + √5/2 is the golden ratio", from
    https://en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
    '''
    phi = (1+sqrt(5))/2
    cornerNote = 'oip ojp ojq oiq  poi qoi qoj poj  ipo jpo jqo iqo'
    facesNote = 'aij ajf afb abe aei bfk bkl ble cdh chl clk ckg cgd dgj dji dih dji dih elh ehi fjg fgk'
    corr1 = {'o':0, 'i':1, 'j':-1, 'p':phi, 'q':-phi}
    corners = [Point(corr1[i], corr1[j], corr1[k]) for i,j,k in cornerNote.split()]
    # Rotate corners by rz, ry degrees. See:
    # https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations
    DtoR = pi/180;   ry, rz = rotay*DtoR, rotaz*DtoR
    sa, ca, sb, cb = sin(rz), cos(rz), sin(ry), cos(ry)
    rox = Point(ca*cb,  -sa,  ca*sb)
    roy = Point(sa*cb,   ca,  sa*sb) # Set up x,y,z rows
    roz = Point(-sb,      0,  cb)    #   of Z,Y rotation matrix
    oa = ord('a')
    # Init an empty layout for local use (ie before deduplication)    
    laylo = Layout(posts=[], cyls=[],  edgeList={})
    # Now make faces, and triangulate those that look feasible
    for i,j,k in facesNote.split():
        pp, qq, rr = corners[ord(i)-oa], corners[ord(j)-oa], corners[ord(k)-oa]
        # Rotate each of the pp, qq, rr faces in space by Z and Y degrees
        p = Point(rox.inner(pp), roy.inner(pp), roz.inner(pp))
        q = Point(rox.inner(qq), roy.inner(qq), roz.inner(qq))
        r = Point(rox.inner(rr), roy.inner(rr), roz.inner(rr))
        # Now maybe triangulate this face if any of its corners are in
        # box.  (When box & face intersection is strictly inside the
        # face we mess up and don't process it.  Oh well.)
        if pointInBox(p,clip1,clip2) or pointInBox(q,clip1,clip2) or pointInBox(r,clip1,clip2):
            pn = len(laylo.posts)
            genTriangleK (laylo, Vfreq, p, q, r, pn)
            print (f'=   {len(laylo.posts):3} posts after face {i}{j}{k}')
        else:            
            print (f'=   {len(laylo.posts):3} posts after face {i}{j}{k} skipped')
    # Have done all faces.  Now dedup & clip laylo and copy points into layin
    print (f'=  {len(laylo.posts)} posts before dedup and clip')
    dedupClip(laylo, layin, clip1, clip2)
    print (f'=  {len(layin.posts)} posts after dedup and clip')

if __name__ == '__main__':
    # This is a test section for genIcosahedron
    phi = (1+sqrt(5))/2;  r = sqrt(2+phi)
    yAngle, zAngle = asin(phi/r)*180/pi, -18 # ~ 58.2825, -18
    for Vfreq in (27,):
        clipLo = Point(-2,-2,-2)
        clipLo = Point(-2,-2,-0.001)
        #clipLo = Point(-2,-2,-0.2)
        clipHi = Point(2,2,2)
        print (f'=  Vfreq {Vfreq},   yAngle {yAngle},  zAngle {zAngle}')
        print (f'=  Clip box corners = {clipLo} and {clipHi}')
        LO = Layout(posts=[], cyls=[],  edgeList={}) # Init an empty layout
        # At present, genIcosahedron reports about
        # posts per face and about dedup/clip stats
        genIcosahedron(LO, Vfreq, clipLo, clipHi, yAngle, zAngle)

        print (f'=  Writing {len(LO.posts)} post coordinates')
        print ('=P  endGap=0 postAxial=f postLabel=f  pDiam=.01  endGap=0  postHi=.02 postDiam=.01 ')
        print ('=L O 0,0,0;')
        np = 0
        for p in LO.posts:
            if np%3==0: print (' C', end='')
            print (f' {p.x:0.5f},{p.y:0.5f},{p.z:0.5f} ', end='')
            np += 1
            if np%3==0: print (';')
        if np%4 !=0: print (';')
        print (";\n=A gg['endGap']=0")
        print ('=C  Mpaa')
        out = 0
        for j in sorted(LO.edgeList.keys()):
            for k in sorted(LO.edgeList[j]):
                if j<k:             # Both of j,k and k,j are in the list
                    print (f' {j:2} {k:2};', end='')
                    out += 1
                    if out%11 == 0: print()
        print()
        print (f'=  Wrote {out} cylinders')
