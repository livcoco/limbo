#!/usr/bin/env python3
# jiw - 29 Feb 2020

# This is rough code for testing generation of geodesic dome point
# coordinates and edges.  Clipping-box parameters allow some control
# over how much of the sphere is generated.  The orientation of the
# icosahedron has two controls: z angle and y angle rotation.
# To see an example, try eg:

# from the pypevue directory within the source tree:
#       ./makeIcosaGeo.py > t1-v; ./pypeVue.py f=t1-v
# after installation anywher on the system:
#       makeIcosaGeo > t1-v; pypeVue f=t1-v

from math import sqrt, pi, asin, sin, cos, atan2, tan, radians, pi, acos, degrees
from pypevue.pypeVue import Point, Layout, ssq, sssq
from pypevue.pypePlugins import FunctionList
from pypevue.pypePlugins.baseFuncs import addEdges
    
Vfreq = 1     # Arbitrary initial value for global Vfreq

class IcosaGeoPoint(Point):
    facess = [(1,2,3,4,5), (6,7,8,9,10,11,12,13,14,15), (16,17,18,19,20)]
    def __init__(self, x, y, z, freq, rank = None, face = None, step = None, stepInRank = None, num=None, nnbrs = None, dupl = None):
        super().__init__(x,y,z)
        self.freq = freq
        self.rank = rank
        self.face = face
        self.step = step
        self.stepInRank = stepInRank
        self.num = num
        self.nnbrs = nnbrs #the number of struts connected to this node
        self.dupl = dupl

    @property
    def stepsInRank(self):
        if self.rank <= self.freq:
            steps = self.rank * 5
        elif self.rank <= 2 * self.freq:
            steps = self.freq * 5
        else:
            steps = (self.freq - (self.rank - self.freq * 2)) * 5
        return steps

    @property
    def topFaces(self):
        return self.facess[0]
    @property
    def midFaces(self):
        return self.facess[1]
    @property
    def bottomFaces(self):
        return self.facess[2]
    @property
    def botFaces(self):
        return self.facess[2]
    @property
    def radius(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    def nutation(self, q):
        ''' 
        Angle from the tangent plane of the sphere at this point, 
        to the vector from this point to q.
        https://www.superprof.co.uk/resources/academic/maths/analytical-geometry/distance/angle-between-line-and-plane.html 
        '''
        # use p for the point of this class
        p = self
        # define the plane tangent to the sphere at p
        pl = Point(2*p.x, 2*p.y, 2*p.z)
        # the vector corresponding to the slope from p to q
        m = Point(q.x-p.x, q.y-p.y, q.z-p.z)
        # find the angle between the line and the plane
        angle = self._angleLineSlopeToPlane(pl, m)

        # determine if q is below or above the plane so we know the sign(+-) of the angle
        # https://math.stackexchange.com/questions/7931/point-below-a-plane
        # If v is the vector that points 'up' and p0 is some point on your plane, and finally p is the point that might be below the plane, compute the dot product v⋅(p−p0). This projects the vector to p on the up-direction. This product is {−,0,+} if p is below, on, above the plane, respectively. 
        qpDiff = q.diff(p)
        qpDiff = Point(qpDiff[0], qpDiff[1], qpDiff[2])
        plInnerQP = pl.inner(qpDiff)
        if plInnerQP < 0:
            angle = -angle
        return degrees(angle)

    def _angleLineSlopeToPlane(self, plane, lineSlope):
        '''
        return the angle (in radians) between a plane and a 3D line defined by it slope
        https://www.superprof.co.uk/resources/academic/maths/analytical-geometry/distance/angle-between-line-and-plane.html 
        '''
        show = False
        if show: print(f'in angle() with plane ({plane}), line ({lineSlope})')
        pl = plane
        m = lineSlope # slope
        # the dot product is proportional to the cosine of the angle between two vectors
        # in this case the vector defining the plane which is normal to the plane
        plInnerM = abs(pl.inner(m)) # inner product is equivalent to dot product
        #magnitude of the vector normal to the plane
        plMag = sqrt(pl.x**2 + pl.y**2 + pl.z**2)
        # magnitude of the vector corresponding to the slope from p to q
        mMag = sqrt(m.x**2 + m.y**2 + m.z**2)
        # normalize the dot product to the magnitude of the two vectors
        # since the plane angle is normal to the plane, we want the arcsine
        # instead of the arccosine 
        if show: print(f'  plInnerM {plInnerM:1.2f},  plMag {plMag:1.2f}, mMag {mMag:1.2f}, plInnerM / (plMag * mMag) {plInnerM / (plMag * mMag)}')
        sinOfAngle = plInnerM / (plMag * mMag)
        sinOfAngle = min(1, max(0, sinOfAngle))
        angle = asin(sinOfAngle)
        if show: print(f'  returning angle {degrees(angle):1.2f}deg')
        return angle

    def distanceToPlane(self, pl):
        '''
        return the shortest distance from myself to the plane pl
        '''
        pln = pl.norm()
        plNorm = Point(pln[0], pln[1], pln[2])
        return self.inner(plNorm)

    def angle(self, q):
        '''
        angle between vector me and another vector q
        '''
        angleCos = self.inner(q) / (self.mag() * q.mag())
        return degrees(acos(angleCos))
    
    def precession(self, q):
        '''
        Angle from the plane created by the tangent line to the circle on the sphere at height z and the origin, 
        to the vector from this point to q.
        https://www.youtube.com/watch?v=Zbiclnu2IlU @ 14:50 min:sec
        '''
        show = False
        p = self
        if show: print(f'in precession() with p ({p}), q ({q})')

        # we need to remove the effects of the nutation angle on the precession angle
        # define the plane tangent to the sphere at p
        plts = Point(2*p.x, 2*p.y, 2*p.z) # plane of tangent sphere
        pltsn = plts.norm()
        pltsNorm = Point(pltsn[0], pltsn[1], pltsn[2])
        nut = radians(self.nutation(q))

        #find the component of q-p perpendicular to the tangent plane of the sphere at p
        qMp = q.diff(p)
        qMp = Point(qMp[0], qMp[1], qMp[2])
        qMpMag = qMp.mag()
        pltsMag = sin(nut)*qMpMag 
        pqPerp = Point(pltsNorm.x*pltsMag, pltsNorm.y*pltsMag, pltsNorm.z*pltsMag)
        aa = qMp.diff(pqPerp)
        m = qMpMpqPerp = Point(aa[0], aa[1], aa[2])
        aa = p.add(m)
        qProj = Point(aa[0], aa[1], aa[2])
        if show:
            aa = p.add(qMpMpqPerp)
            qProj = Point(aa[0], aa[1], aa[2])
            nutCheck = radians(self.nutation(qProj))
            print(f'  plts ({plts}), pltsNorm ({pltsNorm}), nutation {degrees(nut):1.2f}deg, qMp ({qMp}), qMpMag {qMpMag:1.2f}, pltsMag = {pltsMag:1.2f}, pqPerp = ({pqPerp}), qMpMpqPerp = ({qMpMpqPerp}), qProj ({qProj}),  nutCheck {degrees(nutCheck):1.2f}deg (should be 0)')

        dx = 0.1
        if (-dx < self.x < dx) and (-dx < self.y < dx):
            # we are at the top center, define the plane as the y-axis
            pltc = Point(0,1,0)
            tv = Point(1,0,p.z)
        else:
            # find the plane which is formed by the line tangent to the circle on
            # the sphere with z = p.z and the origin.
            tv = Point(-p.y, p.x, 0) #vector in direction of tangent line
            t = Point(p.x + tv.x, p.y + tv.y, p.z) # point on the tangent line
            pl = p.cross(t)
            pltc = Point(2*pl[0], 2*pl[1], 2*pl[2]) # plane of tangent to cicle
            if show: print(f'  tangent vector ({tv}), points on tangent line: p ({p}), t ({t}), pltc ({pltc})')

        angle = self._angleLineSlopeToPlane(pltc, m)
        aa = p.cross(pltc)
        pXpltc = Point(aa[0], aa[1], aa[2])
        if show: print(f'  raw angle {degrees(angle):1.2f}deg, pltc.inner(qProj) {pltc.inner(qProj)}, tv.inner(qProj) {tv.inner(qProj)}, pXpltc.inner(qProj) {pXpltc.inner(qProj)}')
        # adjust the positive acute angle based on what quadrant it is in
        if pltc.inner(qProj) >= 0:
            # 0-180 degrees
            if pXpltc.inner(qProj) > 0:
                # >90
                angle = pi - angle
        else:
            # 180-360
            if pXpltc.inner(qProj) < 0:
                # >270
                angle = pi - angle
            angle += pi
        if angle >= 2*pi:
            angle -= 2*pi
        if show: print(f'  returning angle {degrees(angle):1.2f}deg')
        return degrees(angle)

    def __str__(self):
        return f'num {self.num}, rank {self.rank}, face {self.face}, step {self.step}, nnbrs {self.nnbrs}, dupl {self.dupl}, coords ({self.x:1.2f}, {self.y:1.2f}, {self.z:1.2f})'
    
def genTriangleK (layout, v0, v1, v2, pn):
    def genPoint(p, q, r):
        # Add portions p, q, r of corner points v0, v1, v2
        x = (p*v0.x + q*v1.x + r*v2.x)/Vfreq
        y = (p*v0.y + q*v1.y + r*v2.y)/Vfreq
        z = (p*v0.z + q*v1.z + r*v2.z)/Vfreq
        t = sssq(x,y,z)         # t = distance to origin
        layout.posts.append(IcosaGeoPoint(x/t,y/t,z/t, Vfreq))
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
    a,b,c = Vfreq, 0, 0;     pn=genPoint(a,b,c)
    rbp = re = pn; rb = pn
    for ro in range(1,1+Vfreq):
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

def CCW(p):    # Sort points by descending z and clockwise about z
#    return -round(p.z*1000)-atan2(p.y,p.x)/8
    return -round(p.z*100000)-atan2(p.y,p.x)/8

# FRA2 -- Sort points by rank from 0 and clockwise about z always
# starting a new rank on an icosahedron edge.  Note: this only works
# with zAngle == 0, so should probably remove the rotation option
def FRA2(p):
    angle = atan2(p.y,p.x) #angle from -pi to +pi (-180deg to 180deg)
    # our vertical edge 0 point is at 180 degrees.  If a rounding error occurs, this could
    # become -180 degrees.  Adjust the value to make sure this does not happen
    tweak = (2*pi) / (Vfreq * 5) / 4
    angle -= tweak
    if angle < -pi:
        angle += 2*pi
    if Vfreq < p.rank <= Vfreq * 2: # our 0 edge is on a ~36deg angle
        # each segment makes an approximate angle of 360 / number of segments around the dome
        # each rank has a number of segments offset from the center, e.g. 0.5, 1, 1.5, 2, ...
        segAngle = 2*pi / (Vfreq * 5)
        segsOffset = ((p.rank - Vfreq) / 2)
        angle += segAngle * segsOffset
    elif p.rank > Vfreq * 2:
        # the center line is offset by half a pentagon segment (360/10)
        angle += (2*pi) / 10
    if angle > pi:
        angle -= 2*pi
    sortVal = p.rank - angle/8
    return sortVal

def dedupClip(phase, layi, layo, clip1, clip2):
    '''Given list of points via layi, return de-duplicated and clipped
    list etc'''
    L  = layi.posts
    # eps should be smaller than actual point to point distances,
    # but larger than possible floating point rounding error.
    # For example, .001 is too big to work ok at freq=36.
    eps = 0.00001
    pprev = IcosaGeoPoint(9e9, 8e8, 7e7, Vfreq)
    for n, p in enumerate(L): p.dex = n
    L.sort(key = CCW if phase==1 else FRA2)
    transi = {} # Make node-number translation table for merged points
    if phase==1:
        for p in L:  p.dupl = 1 # how many instances of point
    for p in L:
        me = p.dex; del p.dex
        if pointInBox (p, clip1, clip2):
            if ssq(*(p.diff(pprev))) > eps:
                layo.posts.append(p)
            elif phase==1: pprev.dupl += 1
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

def genIcosahedron(layin, VfreqPar, clip1, clip2, rotay, rotaz):
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
    global Vfreq; Vfreq = VfreqPar
    phi = (1+sqrt(5))/2
    cornerNote = 'oip ojp ojq oiq  poi qoi qoj poj  ipo jpo jqo iqo'
    facesNote = 'aij ajf afb abe aei bfk bkl ble cdh chl clk ckg cgd dgj dji dih elh ehi fjg fgk'
    corr1 = {'o':0, 'i':1, 'j':-1, 'p':phi, 'q':-phi}    
    corners = [IcosaGeoPoint(corr1[i], corr1[j], corr1[k], Vfreq) for i,j,k in cornerNote.split()]
    # Rotate corners by rz, ry degrees. See:
    # https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations
    DtoR = pi/180;   ry, rz = rotay*DtoR, rotaz*DtoR
    sa, ca, sb, cb = sin(rz), cos(rz), sin(ry), cos(ry)
    rox = IcosaGeoPoint(ca*cb,  -sa,  ca*sb, Vfreq)
    roy = IcosaGeoPoint(sa*cb,   ca,  sa*sb, Vfreq) # Set up x,y,z rows
    roz = IcosaGeoPoint(-sb,      0,  cb, Vfreq)    #   of Z,Y rotation matrix
    oa = ord('a')
    # Init empty layouts for local use (ie before deduplication)
    laylo1 = Layout(posts=[], cyls=[],  edgeList={})
    laylo2 = Layout(posts=[], cyls=[],  edgeList={})
    # Now make faces, and triangulate those that look feasible
    for i,j,k in facesNote.split():
        pp, qq, rr = corners[ord(i)-oa], corners[ord(j)-oa], corners[ord(k)-oa]
        # Rotate each of the pp, qq, rr faces in space by Z and Y degrees
        p = IcosaGeoPoint(rox.inner(pp), roy.inner(pp), roz.inner(pp), Vfreq)
        q = IcosaGeoPoint(rox.inner(qq), roy.inner(qq), roz.inner(qq), Vfreq)
        r = IcosaGeoPoint(rox.inner(rr), roy.inner(rr), roz.inner(rr), Vfreq)
        # Now maybe triangulate this face if any of its corners are in
        # box.  (When box & face intersection is strictly inside the
        # face we mess up and don't process it.  Oh well.)
        if pointInBox(p,clip1,clip2) or pointInBox(q,clip1,clip2) or pointInBox(r,clip1,clip2):
            pn = len(laylo1.posts)
            genTriangleK (laylo1, p, q, r, pn)
            #print (f'=   {len(laylo1.posts):3} posts after face {i}{j}{k}')
        else:
            #print (f'=   {len(laylo1.posts):3} posts after face {i}{j}{k} skipped')
            pass
    # Have done all faces.  Now dedup & clip laylo and copy points into layin
    dedupClip(1, laylo1, laylo2, clip1, clip2)
    #print (f'=  Made {len(laylo2.posts)} posts for geodesic with frequency {Vfreq}')

    # Find ranks, or number of rows down from rank-0 center point.
    po = laylo2.posts; elo = laylo2.edgeList;  infin = Vfreq*20
    for k,p in enumerate(po):
        p.num=k;   p.pa = p.pb = p.rank = infin  # Set stuff large
    if len(po): po[0].rank = 0
    for p in po:
        p.nnbrs = len(elo[p.num])
        for dq in elo[p.num]:
            q = po[dq]
            if 1+p.rank < q.rank: # Is p the new pop of q?
                q.pa, q.pb, q.rank = p.num, p.num, 1+p.rank
            elif p.rank < q.rank: # Is p the new mom of q?
                q.pb = p.num
    dedupClip(2, laylo2, layin, clip1, clip2)

    #add face, step, stepInRank
    po = layin.posts; rank = -1; stepsPerFace = [0,0];
    for p in po:
        if p.rank != rank:
            rank = p.rank; faceIdx = 0; step = 0; stepInRank = 0
            if rank == 0:
                faces = p.topFaces
                stepsPerFace[0] = 0
            elif rank == Vfreq * 2:
                faces = p.bottomFaces
                stepsPerFace[0] = Vfreq
            elif rank == Vfreq + 1:
                faces = p.midFaces
                stepsPerFace = [Vfreq-1,1]
            elif rank <= Vfreq:
                stepsPerFace[0] += 1
            elif rank < Vfreq * 2:
                stepsPerFace[0] -= 1
                stepsPerFace[1] += 1
            else:
                stepsPerFace[0] -= 1

        #print(f'TMPDEBUG faces {faces}, faceIdx {faceIdx}')
        p.face = faces[faceIdx]
        p.step = step
        p.stepInRank = stepInRank
        step += 1
        stepInRank += 1
        #print(p)
        if rank <= Vfreq or rank >= Vfreq * 2:
            if step == stepsPerFace[0]:
                faceIdx += 1; step = 0
        else:
            idx = faceIdx % 2
            if step == stepsPerFace[idx]:
                faceIdx += 1; step = 0


if __name__ == '__main__':
    # This is a test section for genIcosahedron
    phi = (1+sqrt(5))/2;  r = sqrt(2+phi)
    yAngle, zAngle = asin(phi/r)*180/pi, 0 # ~ 58.2825
    Vfreq = 30
    Vfreq = 2
    clipLo = Point(-2,-2,-2)
    clipLo = Point(-2,-2,-0.65)
    clipLo = Point(-2,-2,-0.01)
    clipHi = Point(2,2,2)
    print (f'=  Vfreq {Vfreq},   yAngle {yAngle},  zAngle {zAngle}')
    print (f'=  Clip box corners = {clipLo} and {clipHi}')
    LO = Layout(posts=[], cyls=[],  edgeList={}) # Init an empty layout

    FunctionList.registrar('') # can add different functions here (see docs)
    # genIcosahedron will report posts per face and dedup/clip stats
    genIcosahedron(LO, Vfreq, clipLo, clipHi, yAngle, zAngle)

    print (f'=  Writing {len(LO.posts)} post coordinates')
    print ('=P  endGap=0 postAxial=f postLabel=f  pDiam=.01  endGap=0  postHi=.02 postDiam=.01 ')
    print ('=L O 0,0,0;')
    np = 0; quota=3
    for p in LO.posts:
        if np%quota==0: print (' C', end='')
        print (f'  {p.x:8.5f},{p.y:8.5f},{p.z:8.5f} ', end='')
        np += 1
        if np%quota==0: print (';')
    if np%quota !=0: print (';')
    print (f'=  Wrote {np} point locations')
    lopo = LO.posts;  loel = LO.edgeList;  quota = 8
    # Generate sets of cylinders in various colors.
  
    for co in ('Y', 'B', 'R', 'C'):
        print (f'=C  {co}pff')
        out = 0
        for j in sorted(loel.keys()):
            for k in sorted(loel[j]):
                if j<k:   # Both of j,k and k,j are in the list
                    p, q = lopo[j], lopo[k]
                    oB = p.rank == q.rank
                    oY = p.nnbrs==5 or q.nnbrs==5
                    oC = p.dupl>1 and q.dupl>1 and not (oB or oY)
                    oR = not (oB or oY or oC)
                    # Note, some spokes may satisfy multiple conditions.
                    # In next line, one can attach 'and not' clauses to
                    # suppress extra cylinders if necessary.
                    if (co=='B' and oB and not oY) or (co=='Y' and oY) or (co=='R' and oR) or (co=='C' and oC):
                        print (f' {j:2} {k:2};', end='')
                        out += 1
                        if out%quota == 0: print()
        if out%quota != 0: print()
        print (f'=  Wrote {out} {co} cylinders')
