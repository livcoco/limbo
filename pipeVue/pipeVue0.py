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

# Optional parameters for this program are described in a table in
# file `spec-layout-script.odt`.  That documentation tells how to
# specify parameter values endGap, postHi, pDiam, qDiam, SF, and
# numerous others.  Parameter settings can appear on the command line
# as well as in =P lines of a scripts file.  A table in the
# documentation shows parameter names, parameter usage, and default
# values.

#  Note, an end gap is a small gap between a post and a cylinder end.
#  With endGap=3 default value, a gap of about 6 units is drawn
#  between the ends of cylinders meeting at the same point.  With
#  endGap=0, there'd be no gap.

# When modifying this code:
# (a) At outset (ie once only), at command prompt say:
#           V=0;  SCF=pipeVue$V.scad;  PYF=${SCF%.scad}.py
#           STF=${SCF%.scad}.stl; exec-on-change $PYF "python3 $PYF" &
#           echo $PYF, $SCF, $STF;  openscad $SCF &
#     [To avoid a 'Text file busy' shell error message, instead of
#      just saying  ./$PYF  the commands above use python to run $PYF]

#     [Bash script exec-on-change runs a command upon changes.  See
#      https://github.com/ghjwp7/plastics/blob/master/exec-on-change .
#      Or just run the program manually as needed.  When you are
#      working on a particular script S, you can replace the e-o-c
#      shown above with:    exec-on-change S "./pipeVue0 f=S" &
#      which will run pipeVue0 on S whenever you change the file S.]

# (b) After program changes that you want to see the effect of, save
#     the file.  The exec-on-change script will be informed of the
#     file change and will run $PYF.  Then [if openscad's `Design ->
#     Automatic Reload and Preview` option is on] openscad will see
#     that $SCF changed*, and re-render its image.

from solid import cylinder, translate, rotate, scad_render_to_file, color, text
from solid.utils import down, up, left
from sys import argv
from math import sqrt, pi, cos, sin, asin, atan2

class Point:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    def scale(self, s):
        self.x = s*self.x
        self.y = s*self.y
        self.z = s*self.z
    def scalexy(self, s):
        self.x = s*self.x
        self.y = s*self.y

class Layout:
    def __init__(self, BP=Point(0,0,0), OP=Point(0,0,0), posts=[], sRadi=0):
        self.BP = BP  # Current basepoint value
        self.OP = OP  # Origin point of net
        self.posts = posts
        self.sRadi = sRadi       # Radius of interior sphere, if any

    def get3(self):
        return  self.BP, self.OP, self.sRadi

def ssq(x,y,z):    return x*x + y*y + z*z
def sssq(x,y,z):   return sqrt(ssq(x,y,z))
def dist(p,q):     return sqrt(ssq(p.x-q.x, p.y-q.y, p.z-q.z))

def rotate2(a,b,theta):
    st = sin(theta)
    ct = cos(theta)
    return  a*ct-b*st, a*st+b*ct

def produceOut(code, numText, LO):
    '''Modify layout LO according to provided code and numbers'''
    BP, OP, SR = LO.get3()
    bx, by, bz = BP.x, BP.y, BP.z
    nn = len(numText)
    
    def getNums(j, k): # Get list of values from list of number strings
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
        if nums:   LO.BP = Point(*nums);  return

    if code=='O':               # Set origin point, OP
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:   LO.OP = Point(*nums);  return

    if code=='C':               # Create a collection of posts
        nums = getNums(3,33333) # Need at least 3 numbers
        if nums:
            while len(nums) >= 3:
                x, y, z = nums[0]+bx,  nums[1]+by, nums[2]+bz
                LO.posts.append(Point(x,y,z))
                nums = nums[3:]
            if len(nums)>0:
                print (f'Anomaly: code {code}, {numText} has {nums} left over')
            return

    if code=='L':               # Create a line of posts
        nums = getNums(4,4)     # Need exactly 4 numbers
        if nums:
            n, dx, dy, dz = int(nums[0]), nums[1], nums[2], nums[3]
            x, y, z = bx, by, bz
            for k in range(n):
                x, y, z = x+dx, y+dy, z+dz
                LO.posts.append(Point(x,y,z))
            return

    if code=='P':               # Create a polygon of posts        
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:
            n, r, a0 = int(nums[0]), nums[1], nums[2]
            theta = 2*pi/n
            x, y = rotate2(r, 0, a0*pi/180) # a0 in degrees
            for post in range(n):
                LO.posts.append(Point(bx+x, by+y, bz))
                x, y = rotate2(x, y,theta)
            return
    
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
                    LO.posts.append(Point(x,y,z))
                    x += dx
                y += dy
            return
    return                      # We might fail or fall thru

def isTrue(x):
    '''Return false if x is None, or False, or an empty string, or a
    string beginning with f, F, N, or n.  Else, return True.    '''
    return not str(x)[:1] in 'fFNn'

# Compute coords of a letter-point on post based at pb, with top at pt
def levelAt(lev, pb, pt):
    a = (ord(lev)-ord(levels[0]))/(len(levels)-1) # Get portion-of-post
    b = 1-a                                       # Get unused portion
    return b*pb.x+a*pt.x, b*pb.y+a*pt.y, b*pb.z+a*pt.z

def thickLet(thix):
    if thix=='p':
        return SF*pDiam
    else: # diameters q, r, s, t... scale geometrically
        expo = max(0, ord(thix)-ord('q'))
        return SF * qDiam * pow(dRatio,expo)

def postTop(p, OP):   # Given post location p, return loc. of post top
    x, y, z = p.x, p.y, p.z
    ox, oy, oz = (x, y, z-99) if postAxial else (OP.x, OP.y, OP.z)
    u = SF*postHi               # Distance from p to post-top
    v = sssq(x-ox, y-oy, z-oz)  # Distance from p to origin point
    if v>0.01:
        a, b = (u+v)/v, -u/v    # Extrapolation ratios a + b = 1
        tx, ty, tz = a*x+b*ox, a*y+b*oy, a*z+b*oz
    else:
        tx, ty, tz = x, y, z+u  # Fallback if p ~ OP
    siny = min(1, max(-1, (tz-z)/u)) # Don't let rounding error shut us down
    yAxisAngle = (pi/2 - asin(siny)) * 180/pi
    zAxisAngle =  atan2(ty-y, tx-x)  * 180/pi
    #print (f'OP is {ox} {oy} {oz},  xyz {x:0.1f} {y:0.1f} {z:0.1f}   txyz {tx:0.1f} {ty:0.1f} {tz:0.1f}')
    return tx, ty, tz, yAxisAngle, zAxisAngle

def doLayout(layoutScript):
    LO = Layout()
    pc, code, numbers = '?', '?', []
    codes, digits = 'BCLOPRST', '01234356789+-.'
    assembly = None
    
    for cc in layoutScript:   # Process current character of script
        # Add character to number, or store a number?
        if cc in digits:
            num = num + cc if pc in digits else cc
        elif pc in digits:
            numbers.append(num) # Add number to list of numbers

        # Process a completed entry, or start a new entry?
        if cc==';':
            produceOut(code, numbers, LO)
        elif cc in codes:
            pc, code, numbers = '?', cc, []
        pc = cc                 # Prep to get next character

    LO.OP.scale(SF)
    # ----------------------
    # Scale the set of posts 
    for k in range(len(LO.posts)):
        p = LO.posts[k]
        if isTrue(zSpread):
            zrat = 2/(1+p.z/zSize) # assumes z centers at z==0
            p.scalexy(zrat)
        p.scale(SF)
        LO.posts[k] = p
        if isTrue(postList):
            print (f'p{k:<2}=Point( {p.x:8.2f}, {p.y:8.2f}, {p.z:8.2f})')

    # ----------------------
    # Draw the set of posts
    pHi = SF*postHi
    for k, p in enumerate(LO.posts):
        tx, ty, tz, yAxisAngle, zAxisAngle = postTop(p, LO.OP)
        tube = cylinder(d=SF*postDiam, h=SF*postHi)
        tilt = rotate([0,yAxisAngle,zAxisAngle])(tube)
        cyli = translate([p.x, p.y, p.z])(tilt)
        assembly = assembly + cyli if assembly else cyli
        if isTrue(postLabel):
            # -------------------------
            # Draw labels for the posts
            top = Point (tx, ty, tz)
            cName = colorSet['B']
            thik  = thickLet('t')
            x,y,z = levelAt('e', p, top)
            for cc in postLabel:
                if cc in colors: cName = colorSet[cc]
                if cc in thixx:  thik  = thickLet(cc)
                if cc in levels: x,y,z = levelAt(cc, p, top)
            txt = color(cName)(text(text=str(k),size=thik))
            txt = translate([x-thik*len(str(k)), y, z])(txt)
            assembly = assembly + txt

    return assembly, Layout(LO.BP, LO.OP, LO.posts)

def doCylinders(cylScript, LO, assembly):
    def topPoint(p):
        tx, ty, tz, yAxisAngle, zAxisAngle = postTop(p, LO.OP)
        return Point(tx, ty, tz)
    
    def oneCyl(listIt):   # Return a cylinder & its end-post #'s
        m, n = int(post1)%nPosts, int(post2)%nPosts
        p, q = LO.posts[m], LO.posts[n]
        pTop, qTop = topPoint(p), topPoint(q)
        px, py, pz = levelAt(level1, p, pTop)
        qx, qy, qz = levelAt(level2, q, qTop)
        dx, dy, dz = qx-px, qy-py, qz-pz
        L = max(0.1, sssq(dx,  dy,  dz))
        cName = colorSet[colorr]
        alpha = SF*endGap/L     # endGap needs scaling
        # Inputs are scaled, so cx, cy, cz are too.
        cx, cy, cz = px+alpha*dx, py+alpha*dy, pz+alpha*dz
        if isTrue(listIt):
            print (f'Make  {cName:8} {thix} {m:2}{level1} {n:2}{level2}   Length {L:2.2f}')
        yAxisAngle = (pi/2 - asin(dz/L)) * 180/pi
        zAxisAngle =  atan2(dy, dx)      * 180/pi
        diam = thickLet(thix)
        tube = cylinder(d=diam, h=L-SF*2*endGap)
        colo = color(cName)(tube)
        tilt = rotate([0,yAxisAngle,zAxisAngle])(colo)
        # Return a ready-to-use cylinder
        return translate([cx,cy,cz])(tilt), m, n

    def addEdge(v,w):
        if v in edgeList:
            if w not in edgeList[v]:
                edgeList[v].append(w)
        else:
            edgeList[v] = [w]
           
    colorr='G'; thix='p'; pc = None
    post1, post2, level1, level2 = '0', '1', 'c','c'
    nPosts = len(LO.posts)
    nonPost = True
    edgeList, Lmax = {}, 0
    for cc in cylScript:
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
            cyli, p1, p2 = oneCyl(cylList)
            assembly = assembly + cyli if assembly else cyli
            addEdge(p1, p2)  # Add edge to edges list
            addEdge(p2, p1)
            nonPost = True
        pc = cc
    # Finished with cylScript; now see if we need to auto-add cylinders
    cutoff = autoMax
    if cutoff > 0:   # See if any way for any more edges
        print (f'In auto-add, cutoff distance autoMax is {cutoff:7.3f}')
        cutoff2 = cutoff*cutoff
        #print (edgeList)
        for pn in range(nPosts):
            p = LO.posts[pn]
            for qn in range(1+pn, nPosts):
                q = LO.posts[qn]
                dx, dy = p.x-q.x, p.y-q.y
                #print (f'pn {pn}  qn {qn}   cutoff {cutoff:7.2f}  dx {dx:7.2f}  dy {dy:7.2f} ')
                if abs(dx) > cutoff or abs(dy) > cutoff:
                    continue
                d2 = ssq(dx, dy, p.z-q.z)
                #print (f'cutoff2 {cutoff2:7.2f}  d2 {d2:7.2f} ')
                if d2 > cutoff2: continue
                if pn not in edgeList or qn not in edgeList[pn]:
                    post1, post2 = str(pn), str(qn)
                    cyli, p1, p2 = oneCyl(autoList)
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
    parTxt=''; layoutTxt=''; cylTxt=''  # Start with empty scripts
    with open(fiName) as fi:
        for line in fi:
            ll, l1, l2 = len(line), line[:1], line[:2]
            if l1 == '=':   # Detect section change vs comment ...
                if   l2=='=P': mode = 1 # Parameters
                elif l2=='=L': mode = 2 # Layout
                elif l2=='=C': mode = 3 # Cylinders
            else:
                if   mode==1: parTxt    = parTxt    + line
                elif mode==2: layoutTxt = layoutTxt + line
                elif mode==3: cylTxt    = cylTxt    + line
    installParams(parTxt)       # Install params, if any
    return layoutTxt, cylTxt    # Return layout & cylinder scripts

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
    autoMax, autoList  = 0, True
    zSpread, zSize, postAxial = False, 1, True
    for k in range(1,len(argv)):
        paramTxt = paramTxt + ' ' + argv[k]
    installParams(paramTxt)     # Set params from command line
    
    if f == '':
        lScript, cScript = 'C 0,0,0; P5,1,0;', 'Gpae 1,2;;;;1;Rea 1,2;;;;1;'
    else:
        lScript, cScript = loadScriptFile(f) # May override some params
    installParams(paramTxt)     # Again, set params from command line.

    assembly, LO = doLayout(lScript)
    assembly = doCylinders(cScript, LO, assembly)
    scad_render_to_file(assembly, scadFile,
                        file_header = f'$fn = {cylSegments};',
                        include_orig_code=False)
    print (f'Wrote scad code to {scadFile}')
