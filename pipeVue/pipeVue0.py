#!/usr/bin/env python3

# jiw 26 Dec 2018
'''A program that uses SolidPython to generate OpenSCAD code for tubes
along selected edges between `posts` in a plane.  This supports
visualization of arrangements of edges in geodesic dome structures.'''

#  This program processes a layout script and a cylinders script.  See
#  `spec-layout-script.odt` for examples and details of both kinds of
#  scripts.

#  A layout script tells where to locate posts.  It has entries with
#  type of pattern (polygon, rectangular grid, triangular grid),
#  numbers of corners, rows, or columns, and radius or spacing.

#  A cylinders script tells what post-to-post cylinders to make.  It
#  has entries with optional <color>, <diam>, <post>, and <level>
#  elements.  Each semicolon terminator invokes cylinder production.

# Optional params for __main__:
#    designNum/name, endGap, postHi, pDiam, qDiam, SF

# [parameter handling revision, 12 Feb: allow params in any order; but
# require keyword=value forms -- eg `pDiam=0.07` -- where keyword
# should exactly match the name of a variable in the program.]

#  The command-line parameters default to designNum=0, endGap=.03,
#  postHi=.24, pDiam=.06, qDiam=.02, and SF=100 respectively, which
#  when scaled represent design number 0; 3-unit end gaps; 24-unit
#  post heights; thick and thin diameters of 6 and 2 units; and
#  100-unit scale factor.

#  Note, if the first parameter is a valid integer k (or blank, ie
#  k=0) pipeVue uses the k'th built-in script set for layout and
#  cylinders.  Else, the first parameter should be the name of a file
#  containing a layout script and a cylinders script.

#  Note, an end gap is a small gap between a post and a cylinder end.
#  With endGap=3, a gap of about 6 units is drawn between the ends of
#  cylinders meeting at the same point.  With endGap=0, there'd be no
#  gap.

# When modifying this code:
# (a) At outset (ie once only), at command prompt say:
#           V=0;  SCF=pipeVue$V.scad;  PYF=${SCF%.scad}.py
#           STF=${SCF%.scad}.stl; exec-on-change $PYF "python3 $PYF" &
#           echo $PYF, $SCF, $STF;  openscad $SCF &
#     [To avoid a 'Text file busy' shell error message, instead of
#      just saying  ./$PYF  the commands above use python to run $PYF]
# (b) After changes that you want to see the effect of, save the
#     file.  The exec-on-change script will be informed of the file
#     change and will run $PYF.  Then [if openscad `Design -> Automatic
#     Reload and Preview` option is on] openscad will see that $SCF
#     changed*, and re-render its image.
# (c) In openscad, press F6 to render details, then Export, as STL.
# (d) Say `craftware $STF &` then slice it and save gcode

from solid import cylinder, translate, rotate, scad_render_to_file, color, text
from solid.utils import down, up, left
from sys import argv
from math import sqrt, pi, cos, sin, asin, atan2
from collections import namedtuple

Point  = namedtuple('Point',  'x,y,z')
#  Design elements cSides, nPosts, pLayout, cSpec are two integers and
#  two strings. cSides is the number of central sides (allowing
#  central area to be pentagonal, hexagonal, etc.).  nPosts is the
#  total number of posts.  pLayout is a script for arrangement of
#  posts.  cSpec is a script for cylinders between points on posts.
#  Re script contents, see `spec-layout-script.odt`.
Design = namedtuple('Design', 'pLayout, cSpec')
Layout = namedtuple('Layout', 'BP, posts')

def rotate2(a,b,theta):
    st = sin(theta)
    ct = cos(theta)
    return  a*ct-b*st, a*st+b*ct

def ssq(x,y,z):    return x*x + y*y + z*z
def sssq(x,y,z):   return sqrt(ssq(x,y,z))

def produceOut(code, numText, LO):
    BP, posts = LO.BP, LO.posts
    bx, by, bz = BP.x, BP.y, BP.z
    nn = len(numText)
    def getNums(j, k):
        nums = []
        try:
            if j > nn or nn > k :
                raise ValueError;
            for ns in numText:
                nums.append(float(ns))
        except ValueError:
            print (f'Anomaly: code {code}, {numText} has wrong count or format')
            return None
        return nums
    
    if code=='B':               # Set base point, BP
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums: return Layout(Point(*nums), posts)

    if code=='C':               # Create a collection of posts
        nums = getNums(3,33333) # Need at least 3 numbers
        if nums:
            while len(nums) >= 3:
                x, y, z = nums[0]+bx,  nums[1]+by, nums[2]+bz
                posts.append(Point(x,y,z))
                nums = nums[3:]
            if len(nums)>0:
                print (f'Anomaly: code {code}, {numText} has {nums} left over')
            return Layout(BP, posts)

    if code=='L':               # Create a line of posts
        nums = getNums(4,4)     # Need exactly 4 numbers
        if nums:
            n, dx, dy, dz = int(nums[0]), nums[1], nums[2], nums[3]
            x, y, z = bx, by, bz
            for k in range(n):
                x, y, z = x+dx, y+dy, z+dz
                posts.append(Point(x,y,z))
            return Layout(BP, posts)

    if code=='P':               # Create a polygon of posts        
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:
            n, r, a0 = int(nums[0]), nums[1], nums[2]
            theta = 2*pi/n
            x, y = rotate2(r, 0, a0*pi/180) # a0 in degrees
            for post in range(n):
                posts.append(Point(bx+x, by+y, bz))
                x, y = rotate2(x, y,theta)
    
    if code in 'RT':            # Create an array of posts
        nums = getNums(4,4)     # Need exactly 4 numbers
        if nums:
            r, c, dx, dy = int(nums[0]), int(nums[1]), nums[2], nums[3]
            y, z = by, bz
            for rr in range(r):
                x, roLen = bx, c
                # For odd rows of triangular arrays, offset the row
                if code=='T' and (rr&1)==1:
                    x, roLen = bx - dx/2, c+1
                for cc in range(roLen):
                    posts.append(Point(x,y,z))
                    x += dx
                y += dy
            return Layout(BP, posts)
    return LO                   # No change if we fail or fall thru

def isTrue(x):
    '''Return false if x is None, or False, or an empty string, or a
    string beginning with f, F, N, or n.  Else, return True.    '''
    return not str(x)[:1] in 'fFNn'

def levelLet(lev):
    deltaHi = SF*postHi/(len(levels)-1)
    return (ord(lev)-ord(levels[0]))*deltaHi

def thickLet(thix):
    if thix=='p':
        return SF*pDiam
    else: # diameters q, r, s, t... scale geometrically
        expo = max(0, ord(thix)-ord('q'))
        return SF * qDiam * pow(dRatio,expo)
        
def doLayout(dz):
    LO = Layout(Point(0,0,0), [])
    pc, code, numbers = '?', '?', []
    codes, digits = 'BCLPRT', '01234356789+-.'
    
    for cc in dz.pLayout:       # Process current character
        # Add character to number, or store a number?
        if cc in digits:
            num = num + cc if pc in digits else cc
        elif pc in digits:
            numbers.append(num) # Add number to list of numbers

        # Process a completed entry, or start a new entry?
        if cc==';':
            LO = produceOut(code, numbers, LO)
        elif cc in codes:
            pc, code, numbers = '?', cc, []
        pc = cc                 # Prep to get next character
        
    # Now LO has an unscaled points list.  Create and return CSG
    posts = LO.posts
    for k in range(len(posts)):
        p = posts[k]
        p = Point(SF*p.x, SF*p.y, SF*p.z)
        posts[k] = p
        if isTrue(postList):
            #print (f'Post {k:<3} ({p.x:8.2f}, {p.y:8.2f}, {p.z:8.2f} )')
            print (f'p{k:<2}=Point( {p.x:8.2f}, {p.y:8.2f}, {p.z:8.2f})')
    assembly = None
    for k, p in enumerate(posts):
        tube = cylinder(d=SF*postDiam, h=SF*postHi)
        cyli = translate([p.x, p.y, p.z])(tube)
        assembly = assembly + cyli if assembly else cyli
        if isTrue(postLabel):
            cName = colorSet['B']
            thik  = thickLet('t')
            zd    = levelLet('e')
            for cc in postLabel:
                if cc in colors: cName = colorSet[cc]
                if cc in thixx:  thik  = thickLet(cc)
                if cc in levels: zd = levelLet(cc)
            tx =  color(cName)(text(text=str(k),size=thik))
            tr = translate([p.x-thik*(1+len(str(k))), p.y, zd+p.z])(tx)
            assembly = assembly + tr

    return assembly, Layout(LO.BP, posts)

def doCylinders(dz, LO, assembly):
    def oneCyl(listIt):   # Return a cylinder & its end-post #'s
        m, n = int(post1)%nPosts, int(post2)%nPosts
        p, q = posts[m], posts[n]
        za1 = levelLet(level1)
        za2 = levelLet(level2)
        pz, qz = za1 + p.z, za2 + q.z
        # p, q are scaled, so dx,dy,dz & L are too.
        dz, dx, dy = qz-pz, q.x-p.x,  q.y-p.y
        L = max(0.1, sssq(dx,  dy,  dz))
        cName = colorSet[colorr]
        alpha = SF*endGap/L     # endGap needs scaling
        # Inputs are scaled, so cx, cy, cz are too.
        cx, cy, cz = p.x+alpha*dx, p.y+alpha*dy, pz+alpha*dz
        if isTrue(listIt):
            print (f'Make  {cName:8} {thix} {m:2}{level1} {n:2}{level2}   Length {L:2.2f}')
        yAxisAngle = (pi/2 - asin(dz/L)) * 180/pi
        zAxisAngle =  atan2(dy, dx)      * 180/pi
        diam = thickLet(thix)
        tube = cylinder(d=diam, h=L-SF*2*endGap)
        colo = color(cName)(tube)
        tilt = rotate([0,yAxisAngle,zAxisAngle])(colo)
        # Return a ready-to-use cylinder
        return translate([cx,cy,cz])(tilt), m, n, L

    def addEdge(v,w):
        if v in edgeList:
            if w not in edgeList[v]:
                edgeList[v].append(w)
        else:
            edgeList[v] = [w]
           

    specs, posts = dz.cSpec, LO.posts
    colorr='G'; thix='p'; pc = None
    post1, post2, level1, level2 = '0', '1', 'c','c'
    nPosts = len(posts)
    nonPost = True
    edgeList, Lmax = {}, 0
    for cc in specs:
        if cc in colors: colorr = cc
        elif cc in thixx: thix  = cc
        elif cc in levels:
            level1, level2 = level2, cc
        elif cc in digits:
            if pc in digits:  post2 = post2 + cc
            else:             post1, post2 = post2, cc
            nonPost = False
        elif cc=='/':
            level1, level2 = level2, level1
        elif cc==';':
            if nonPost:
                post1, post2 = str(1+int(post1)), str(1+int(post2))
            cyli, p1, p2, L = oneCyl(cylList)
            assembly = assembly + cyli if assembly else cyli
            Lmax = max(L, Lmax)
            addEdge(p1, p2)  # Add edge to edges list
            addEdge(p2, p1)
            nonPost = True
        pc = cc
    # Finished with specs; now see if we need to auto-add cylinders
    cutoff = Lmax + autoTol
    if cutoff > 0:   # See if no way for any more edges
        print (f'In auto-add, cutoff distance is {cutoff:7.2f} = Lmax + autoTol = {Lmax:0.2f} + {autoTol}')
        cutoff2 = cutoff*cutoff
        print (edgeList)
        for pn in range(nPosts):
            p = posts[pn]
            for qn in range(1+pn, nPosts):
                q = posts[qn]
                dx, dy = p.x-q.x, p.y-q.y
                #print (f'pn {pn}  qn {qn}   cutoff {cutoff:7.2f}  dx {dx:7.2f}  dy {dy:7.2f} ')
                if abs(dx) > cutoff or abs(dy) > cutoff:
                    continue
                d2 = ssq(dx, dy, p.z-q.z)
                #print (f'cutoff2 {cutoff2:7.2f}  d2 {d2:7.2f} ')
                if d2 > cutoff2: continue
                if pn not in edgeList or qn not in edgeList[pn]:
                    post1, post2 = str(pn), str(qn)
                    cyli, p1, p2, L = oneCyl(autoList)
                    assembly = assembly + cyli if assembly else cyli
    return assembly

def installParams(parTxt):
    '''Given a string like "var1=val1 var2=val2 var3=val3 ...", extract
    the variable names and the values, convert the values to numeric
    forms, and store the values in the globals() dict.  Note, do not
    use any white space within any of the var=val strings.    '''
    plist = parTxt.split()      # Split the list on white space
    flubs, glob = '', globals()
    for vev in plist:
        p = pq = vev.split('=')          # Split each equation on = sign
        q, ok = '', False
        if len(pq)==2:
            p, q = pq
            if p in glob.keys():
                t, v = type(glob[p]), q
                try:
                    if   t==int:   v = int(q);     ok=True
                    elif t==float: v = float(q);   ok=True
                    elif t==bool:  v = isTrue(q);  ok=True
                    elif t==str:   v = q;          ok=True
                except:  pass
        if ok:
            glob[p] = v
        else:  flubs += f' [ {p} {q} ] '
    if flubs: print (f'Parameter-setting fail: {flubs}')

def loadScriptFile(fiName):
    '''Read parameters, layout script, and cylinders script from file'''
    mode = 0;                   # Start out in comments mode
    pt = los = cs = ''          # Start with empty scripts
    with open(fiName) as fi:
        for line in fi:
            ll, l1, l2 = len(line), line[:1], line[:2]
            if l1 == '=':   # Detect section change vs comment ...
                if   l2=='=P': mode = 1 # Parameters
                elif l2=='=L': mode = 2 # Layout
                elif l2=='=C': mode = 3 # Cylinders
            else:
                if   mode==1: pt  = pt  + line
                elif mode==2: los = los + line
                elif mode==3: cs  = cs  + line
    installParams(pt)           # Install params, if any
    return Design(los, cs)      # Return Design

colors, levels = 'GYRBCMW',  'abcde'
thixx,  digits = 'pqrstuvw', '01234356789'
colorSet = dict({'G':'Green', 'Y':'Yellow', 'R':'Red', 'B':'Blue', 'C':'Cyan', 'M':'Magenta', 'W':'White'})   

if __name__ == '__main__':
    # Set initial values of main parameters
    pDiam,   qDiam,    dRatio   = 0.06, 0.02, sqrt(2)
    endGap,  postHi,   postDiam = 0.03, 0.16, qDiam
    f,       SF,    cylSegments = '', 100, 30
    version, paramTxt, postLabel= 0, '','Bte' # Blue, size u, level e
    scadFile = f'pipeVue{version}.scad' # Name of scad output file
    postList = cylList = False # Control printing of post and cyl data
    autoTol, autoList = -1e9, False
    for k in range(1,len(argv)):
        paramTxt = paramTxt + ' ' + argv[k]
    installParams(paramTxt)     # Set params from command line
    
    if f == '':
        dz = Design('C 0,0,0; P5,1,0;', 'Gpae 1,2;;;;1;')
    else:
        dz = loadScriptFile(f)  # May install params from file.
    installParams(paramTxt)     # Again, set params from command line.

    assembly, LO = doLayout(dz)
    assembly = doCylinders(dz, LO, assembly)
    scad_render_to_file(assembly, scadFile,
                        file_header = f'$fn = {cylSegments};',
                        include_orig_code=False)
    print (f'Wrote scad code to {scadFile}')
