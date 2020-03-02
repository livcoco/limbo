#!/usr/bin/env python3

# jiw 26 Dec 2018
'''A program that generates OpenSCAD code for tubes along selected
edges between `posts` in a plane.  This supports visualization of
arrangements of edges in geodesic dome structures.'''

# This program processes a layout script and a cylinders script, while
# generating an output file containing OpenSCAD code to represent the
# structure described by the scripts.  See `pypeVue.dox.odt` for
# examples and details of scripts and how to write them.  Note, in
# default form the generated SCAD code depends upon access to file
# pypeVue.codebase.scad .

#  A layout script tells where to locate posts.  It has entries with
#  type of pattern (polygon, rectangular grid, triangular grid),
#  numbers of corners, rows, or columns, and radius or spacing.

#  A cylinders script tells what post-to-post cylinders to make.  It
#  has entries with optional <color>, <diam>, <post>, and <level>
#  elements.  Each semicolon terminator invokes cylinder production.

# Optional parameters for this program are described in a table in
# pypeVue.dox.odt, which tells how to specify parameter values endGap,
# postHi, pDiam, qDiam, SF, and numerous other params.  Parameter
# settings can appear on the command line as well as in =P lines of a
# scripts file.  A table in the documentation shows parameter names,
# parameter usage, and default values.

#  Note, an end gap is a small gap between a post and a cylinder end.
#  With endGap=3 default value, a gap of about 6 units is drawn
#  between the ends of cylinders meeting at the same point.  With
#  endGap=0, there'd be no gap.

# When modifying this code:
# (a) At outset (ie once only), at command prompt say:
#           SCF=pypeVue.scad;  PYF=${SCF%.scad}.py
#           STF=${SCF%.scad}.stl; exec-on-change $PYF "python3 $PYF" &
#           echo $PYF, $SCF, $STF;  openscad $SCF &
#     [To avoid a 'Text file busy' shell error message, instead of
#      just saying  ./$PYF  the commands above use python to run $PYF]

#     [Bash script exec-on-change runs a command upon changes.  See
#      https://github.com/ghjwp7/plastics/blob/master/exec-on-change .
#      Or just run the program manually as needed.  When you are
#      working on a particular script S, you can replace the e-o-c
#      shown above with:    exec-on-change S "./pypeVue f=S" &
#      which will run pypeVue0 on S whenever you change the file S.]

# (b) After program changes that you want to see the effect of, save
#     the file.  The exec-on-change script will be informed of the
#     file change and will run $PYF.  Then [if openscad's `Design ->
#     Automatic Reload and Preview` option is on] openscad will see
#     that $SCF changed*, and re-render its image.

from sys import argv, exit, exc_info, stderr
from datetime import datetime
from math import sqrt, pi, cos, sin, asin, atan2

def ssq(x,y,z):    return x*x + y*y + z*z
def sssq(x,y,z):   return sqrt(ssq(x,y,z))

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
    def inner(self, q):         # Inner product of two 3-vectors
        return (self.x*q.x + self.y*q.y + self.z*q.z)
    def diff(self, q):
        return (self.x-q.x, self.y-q.y, self.z-q.z)
    def str(self, places):
        x,y,z = (round(k,places) for k in (self.x, self.y ,self.z))
        return f'{x}, {y}, {z}'
    def __str__( self):  return self.str(3)
    def __repr__(self):  return self.str(8)
    def __lt__(a, b):     # To sort points in x,y,z order
        return (a.x < b.x) or (a.x == b.x and a.y < b.y) or (a.x == b.x and a.y == b.y and a.z <= b.z)

class Post:
    def __init__(self, foot, top=0, diam=0, hite=0, yAngle=0, zAngle=0, num=0, data=0):
        self.foot = foot
        self.top  = top
        self.diam = diam
        self.hite = hite
        self.yAngle = yAngle
        self.zAngle = zAngle
        self.num  = num
        self.data = data 
    def __str__( self):
        return f'Post {self.num} ({self.foot}) ({self.top}) {round(self.yAngle,1)} {round(self.zAngle,1)}  '
    def __repr__(self):  return self.__str__()

class Cylinder:
    def __init__(self, post1, post2, lev1, lev2, colo, thix, gap, data=0, num=0):
        self.put9 (post1, post2, lev1, lev2, colo, thix, gap, data, num)
        
    def get9(self):
        return self.post1, self.post2, self.lev1, self.lev2, self.colo, self.thix, self.gap, self.data, self.num
    
    def put9(self, post1, post2, lev1, lev2, colo, thix, gap, data, num):
        self.post1, self.post2 = post1, post2
        self.lev1,  self.lev2  = lev1,  lev2
        self.colo,  self.thix  = colo,  thix
        self.gap,   self.data =  gap,   data
        self.num = num
        
    def __str__( self):
        return f'Cylinder {self.num} ({self.post1},{self.post2}) {self.colo}{self.thix}{self.lev1}{self.lev2} {round(self.gap,2)}'
    def __repr__(self):  return self.__str__()
    
class Layout:
    def __init__(self, BP=Point(0,0,0), OP=Point(0,0,0),
                 posts=[], cyls=[],  edgeList={}):
        self.BP = BP  # Current basepoint value
        self.OP = OP  # Origin point of net
        self.posts = posts
        self.cyls  = cyls
        self.edgeList = edgeList
    def get4(self):
        return  self.BP, self.OP, self.posts, self.cyls
    def __str__( self):
        return f'Layout: BP({self.BP})  OP({self.OP});  {len(self.posts)} posts, {len(self.cyls)} cyls'
    
#---------------------------------------------------------
def arithmetic(line, xTrace):
    # Remove =A and any leading whitespace, to avoid indentation error
    code = line.lstrip()
    userLocals['gg']=globals() # Allow access to globals, albeit awkward
    if xTrace:   print (f'Code to exec:  {code}')
    try:
        if code:    # have we got any user-code?
            exec (code, userLocals)   # yes, try to execute it
    except SystemExit:
        exit(0)             # allow code to exit
    except Exception:           # catch normal exceptions
        er = exc_info()
        print(f'Got error:  {er[1]}   from  {code.strip()}')
#---------------------------------------------------------

def rotate2(a,b,theta):
    st = sin(theta)
    ct = cos(theta)
    return  a*ct-b*st, a*st+b*ct

def isTrue(x):
    '''Return false if x is None, or False, or an empty string, or a
    string beginning with f, F, N, or n.  Else, return True.    '''
    return not str(x)[:1] in 'fFNn'

# Compute coords of a letter-point on post p
def levelAt(lev, p):
    a = (ord(lev)-ord(levels[0]))/(len(levels)-1) # Get portion-of-post
    b = 1-a                                       # Get unused portion
    pf, pt = p.foot, p.top
    return Point(round(b*pf.x+a*pt.x, 2), round(b*pf.y+a*pt.y, 2), round(b*pf.z+a*pt.z, 2))

def thickLet(thix):
    if thix=='p':
        return SF*pDiam
    else: # diameters q, r, s, t... scale geometrically
        expo = max(0, ord(thix)-ord('q'))
        return round(SF * qDiam * pow(dRatio,expo), 2)

def addEdge(v,w, layout):
    if v in layout.edgeList:
        if w not in layout.edgeList[v]:
            layout.edgeList[v].append(w)
    else:
        layout.edgeList[v] = [w]

def addEdges(v,w, layout):
    addEdge(v,w,layout); addEdge(w,v,layout)
#---------------------------------------------------------
def generatePosts(code, numberTexts):
    '''Modify layout LO according to provided code and numbers'''
    B = LO.BP
    bx, by, bz = B.x, B.y, B.z
    nn = len(numberTexts)
    
    def getNums(j, k): # Get list of values from list of number strings
        nums = []
        try:
            if j > nn or nn > k :
                raise ValueError;
            for ns in numberTexts:
                nums.append(float(ns))
        except ValueError:
            print (f'Anomaly: code {code}, {numberTexts} has wrong count or format')
            return None
        return nums

    def postAt(x,y,z): LO.posts.append(Post(Point(x,y,z)))
    
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
                postAt(nums[0]+bx,  nums[1]+by, nums[2]+bz)
                nums = nums[3:]
            if len(nums)>0:
                print (f'Anomaly: code {code}, {numberTexts} has {nums} left over')
            return

    if code=='L':               # Create a line of posts
        nums = getNums(4,4)     # Need exactly 4 numbers
        if nums:
            n, dx, dy, dz = int(nums[0]), nums[1], nums[2], nums[3]
            x, y, z = bx, by, bz
            for k in range(n):
                postAt(x+dx, y+dy, z+dz)
            return

    if code=='P':               # Create a polygon of posts        
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:
            n, r, a0 = int(nums[0]), nums[1], nums[2]
            theta = 2*pi/n
            x, y = rotate2(r, 0, a0*pi/180) # a0 in degrees
            for post in range(n):
                postAt(bx+x, by+y, bz)
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
                    postAt(x,y,z)
                    x += dx
                y += dy
            return
    return                      # We might fail or fall thru

def postTop(p, OP):   # Given post location p, return loc. of post top
    x, y, z = p.foot.x, p.foot.y, p.foot.z
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
    return Point(tx,ty,tz), round(yAxisAngle,2), round(zAxisAngle,2)

#==================================
def scriptCyl(ss, preCyl):
    post1, post2, lev1, lev2, colo, thix, gap, nonPost, num = preCyl.get9()
    mode = 0                    # mode 0 = comments at start
    pc, code = '?', '?'
    for cc in ss:            
        if pc == '#':       # Insert a simple variable's value
            post1, post2 = post2, userLocals[cc]
            nonPost = False
        elif cc in colors: colo = cc
        elif cc in thixx: thix  = cc
        elif cc in levels:
            lev1, lev2 = lev2, cc
        elif cc in digits:
            if pc in digits:  post2 = post2 + cc
            else:             post1, post2 = post2, cc
            nonPost = False
        elif cc=='/':
            lev1, lev2 = lev2, lev1
        elif cc==';':
            p1, p2 = int(post1), int(post2)
            if nonPost:
                p1, p2 = p1+1, p2+1
                post1, post2 = str(p1), str(p2)
            num = len(LO.cyls)
            cyl = Cylinder(p1, p2, lev1, lev2, colo, thix, gap, 0, num)
            LO.cyls.append(cyl)
            addEdges(p1, p2, LO) # Add edges p1,p2 and p2,p1 to edges list
            nonPost = True
        pc = cc
    preCyl.put9(post1, post2, lev1, lev2, colo, thix, gap, nonPost, num)
    return preCyl
#==================================
def scriptPost(ss, prePost):
    pc, code, numbers = '?', '?', prePost.data
    for cc in ss:   # Process characters of script
        # Add character to number, or store a number, or what?
        if pc == '#':       # Set or use a simple variable
            if code=='?':
                userLocals[cc] = len(LO.posts)
            else:      # Substitute value into list of numbers
                numbers.append(userLocals[cc])
        elif cc in digits:
            num = num + cc if pc in digits else cc
        elif pc in digits:
            numbers.append(num) # Add number to list of numbers
        # Process a completed entry, or start a new entry?
        if cc==';':
            generatePosts(code, numbers)
            code = '?'
        elif cc in codes:
            pc, code, numbers = '?', cc, []
        pc = cc             # Prep to get next character
    prePost.data = numbers
    return
#==================================
def runScript(scripts):
    preCyl = Cylinder(0, 1, 'c','c', 'G', 'p', endGap, True, 0)
    prePost = Post(0, data=[])
    mode = 0                    # mode 0 = comments at start
    numbers = []
    for line in scripts:
        l1, l2, ss, ll = line[:1], line[:2], line[2:], line
        if   l2=='=C': mode = 'C'; ll=ss # Cylinders
        elif l2=='=L': mode = 'L'; ll=ss # Layout
        
        elif l2=='=P':          # Process Parameters line
            installParams((ss,paramTxt));
            continue
        elif l2=='=A':          # Process Arithmetic line
            arithmetic(ss, traceExec);
            continue
        elif l1=='=':           # Process comment line
            continue
        
        if   mode == 'L':       # Process Posts line
            scriptPost(ll, prePost)
        elif mode == 'C':       # Process Cylinders line
            scriptCyl (ll, preCyl)
#==================================
def writePosts(fout):
    LO.OP.scale(SF)      # Get ready to orient the posts: scale the OP
    # Scale the set of posts, and compute their tops and angles 
    pHi, pDi = SF*postHi, SF*postDiam
    for k, p in enumerate(LO.posts):
        p.num = k
        if isTrue(zSpread):
            zrat = 2/(1+p.foot.z/zSize) # assumes z centers at z==0
            p.foot.scalexy(zrat)
        p.foot.scale(SF)
        if isTrue(postList):
            print (f'p{k:<2}=Point( {p.foot})')
        p.diam, p.hite = pDi, pHi
        p.top, p.yAngle, p.zAngle = postTop(p, LO.OP)
    fout.write(f'//  onePost(num, diam, hite, yAngle, zAngle, foot x y z, top x y z);\n')
    for p in LO.posts:
        fout.write(f'    onePost({p.num}, {p.diam}, {p.hite}, {p.yAngle}, {p.zAngle}, {p.foot}, {p.top});\n')

#==================================
def writeLabels(fout):
    if not isTrue(postLabel):
        return
    fout.write(f'\n//  oneLabel(num, diam, hite, yAngle, zAngle, foot x y z, label x y z, tColor, tSize, txt);\n')
    for p in LO.posts:
        cName = colorSet['B']
        thik  = thickLet('t')
        lxyz  = levelAt('e', p)
        for cc in postLabel:
            if cc in colors: cName = colorSet[cc]
            if cc in thixx:  thik  = thickLet(cc)
            if cc in levels: lxyz  = levelAt(cc, p)
        txt = f'{str(p.num)}'
        fout.write(f'    oneLabel({p.num}, {p.diam}, {p.hite}, {p.yAngle}, {p.zAngle}, {p.foot}, {lxyz}, {cName}, {thik}, "{txt}");\n')

#=====================================
            
def writeCylinders(fout, clo, chi, listIt):
    fout.write('\n//  oneCyl (p1,2,cylDiam,cylLen,yAngle,zAngle,  c xyz,  e xyz, cColor)\n')
    posts = LO.posts
    nPosts = len(posts)
    for nCyl in range(clo, chi):
        cyl = LO.cyls[nCyl]     # Draw this cylinder
        post1, post2, lev1, lev2, colo, thix, gap, data, num = cyl.get9()
        gap = SF*gap            # gap needs scaling
        p1, p2 = min(post1,nPosts-1), min(post2,nPosts-1)
        try:
            pp, qq = posts[p1], posts[p2]
        except:
            print (f'Fatal Error with p1= {p1},   p2= {p2},  nPosts {nPosts}')
            exit(0)
        p = levelAt(lev1, pp)
        q = levelAt(lev2, qq)
        dx, dy, dz = q.diff(p)
        L = round(max(0.1, sssq(dx,  dy,  dz)), 2)
        cName = colorSet[colo]
        alpha = gap/L
        cc = Point(p.x+alpha*dx, p.y+alpha*dy, p.z+alpha*dz)
        cylNum= 1000*p1 + p2
        if isTrue(listIt):
            print (f'Make {cyl}  L {L:2.2f}  {cName}')
        yAngle = round((pi/2 - asin(dz/L)) * 180/pi, 2)
        zAngle = round( atan2(dy, dx)      * 180/pi, 2)
        diam = thickLet(thix)
        fout.write(f'    oneCyl({p1}, {p2}, {diam}, {round(L-2*gap,3)}, {yAngle}, {zAngle}, {cc}, {pp.foot}, {cName});\n')

#-------------------------------------------------------------
def autoAdder(fout):    # See if we need to auto-add cylinders
    cutoff = autoMax
    clo = len(LO.cyls) # Record how many cylinders are already processed
    nPosts = len(LO.posts)
    edgeList = LO.edgeList
    # in this version punt color, thix, levels ...
    colo, thix, lev1, lev2 = 'B', 'p', 'c','c'
            
    if cutoff > 0:     # See if any way for any more edges
        print (f'In auto-add, cutoff distance autoMax is {cutoff:7.3f}')
        cutoff2 = cutoff*cutoff
        for pn in range(nPosts):
            p = LO.posts[pn].foot
            for qn in range(1+pn, nPosts):
                q = LO.posts[qn].foot
                dx, dy, dz = p.diff(q)
                if abs(dx) > cutoff or abs(dy) > cutoff:
                    continue
                d2 = ssq(dx, dy, dz)
                if d2 > cutoff2: continue
                if pn not in edgeList or qn not in edgeList[pn]:
                    post1, post2 = str(pn), str(qn)          
                    cyl = Cylinder(pn,qn, lev1, lev2, colo, thix, endGap, 0,0)
                    LO.cyls.append(cyl)
                    addEdges(pn, qn, LO)
        writeCylinders(fout, clo, len(LO.cyls), autoList)
#-------------------------------------------------------------
def installParams(script):
    '''Given script lines that are Parameter-setting lines, this extracts
    variable names and values from it, like "var1=val1 var2=val2
    var3=val3 ...".  It converts the values to numeric forms, and
    stores them in globals() dict.  Note, white space is taboo within
    a var=val item.    '''
    flubs, glob = False,  globals()
    for parTxt in script:
        plist = parTxt.split()       # Split the list on white space
        for vev in plist:
            p = pq = vev.split('=')  # Split each equation on = sign
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
            if ok: glob[p] = v
            else:  flubs = True
        if flubs: print (f'Parameter-setting fail in "{parTxt}"')

def loadScriptFile(fiName):   # Read scripts from file
    with open(fiName) as fi:
        return fi.readlines()

if __name__ == '__main__':
    colors, levels = 'GYRBCMW',  'abcde'
    thixx,  digits = 'pqrstuvw', '01234356789'
    colorSet = {'G':'"Green"', 'Y':'"Yellow"', 'R':'"Red"', 'B':'"Blue"', 'C':'"Cyan"', 'M':'"Magenta"', 'W':'"White"'}   

    # Set initial values of main parameters
    pDiam,   qDiam,    dRatio   = 0.06, 0.02, sqrt(2)
    #endGap,  postHi,   postDiam = 0.03, 0.16, qDiam
    # *** debug issue: endGap not setting from params in file ***
    endGap,  postHi,   postDiam = 0.007, 0.16, qDiam
    f,       SF,    cylSegments = '', 100, 30
    paramTxt, postLabel= '', 'Bte' # Blue, size u, level e
    userCode = ''
    codeBase = f'pypeVue.codeBase.scad' # SCAD functions for posts, cyls, etc
    scadFile = f'pypeVue.scad'          # Name of scad output file
    postList = cylList = False # Control printing of post and cyl data
    autoMax, autoList  = 0, True
    zSpread, zSize, postAxial = False, 1, True
    userPar0 = userPar1 = userPar2 = '""'
    traceExec=False
    codes, digits = 'BCLOPRST', '01234356789+-.'
    script1 = '=P postDiam=.1 endGap=.05','=C Gpae 1,2;;;;1;Rea 1,2;;;;1;','=L C 0,0,0; P5,1,0;'
    for k in range(1,len(argv)):
        paramTxt = paramTxt + ' ' + argv[k]

    userLocals = {}               # Initialize empty user-space dict
    from os import path
    myname = path.splitext(path.basename(__file__))[0]
    exec(f'from {myname} import Point,Post,Layout', globals(), userLocals)
    LO = Layout()  # Start an empty layout.  LO is global.
    installParams((paramTxt,)) # To set f param from command line if specified
    scripts = script1 if f == '' else loadScriptFile(f)
    runScript(scripts)          # Create post locations layout
    date = datetime.today().strftime('%Y-%m-%d  %H:%M:%S')

    frontCode = f'''// File {scadFile}, generated  {date}
// by pypeVue from script "{f}"
$fn = {cylSegments};
userPar0 = {userPar0};
userPar1 = {userPar1};
userPar2 = {userPar2};
{'' if codeBase else '//'}use <{codeBase}>
{'' if userCode else '//'}include <{userCode}>
difference() {'{'}
  union() {'{'}
'''
    backCode = f'''\n
    addOn({SF}, {LO.BP}, {LO.OP}); // unscaled BP, OP
  {'}'}
  subOff ({SF}, {LO.BP}, {LO.OP}); // unscaled BP, OP
{'}'}'''

    with open(scadFile, 'w') as fout:
        fout.write(frontCode)
        writePosts    (fout)
        writeLabels   (fout)
        writeCylinders(fout, 0, len(LO.cyls), cylList)
        autoAdder     (fout)
        fout.write(backCode)
    print (f'For script "{f}", pypeVue wrote code to {scadFile} at {date}')
