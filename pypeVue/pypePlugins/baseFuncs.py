# jiw March 2020
'''Base functions for pypeVue, a program that generates OpenSCAD code
for tubes along selected edges between `posts` in a plane.  This
supports visualization of arrangements of edges in geodesic dome
structures.'''

from sys import argv, exit, exc_info, stderr
from datetime import datetime
from math import sqrt, pi, cos, sin, asin, atan2
from pypeVue import ssq, sssq, rotate2, isTrue 
from pypeVue import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList

#---------------------------------------------------------
def arithmetic(line, xTrace):
    # Remove =A and any leading whitespace, to avoid indentation error
    ref = FunctionList
    code = line.lstrip()
    if xTrace:   print (f'Code to exec:  {code}')
    try:
        if code:    # have we got any user-code?
            exec (code, ref.userLocals)   # yes, try to execute it
    except SystemExit:
        exit(0)             # allow code to exit
    except Exception:       # catch normal exceptions
        er = exc_info()
        print(f'Got error:  `{er[1]}`   from code:  `{code.strip()}`')
#---------------------------------------------------------
# Compute coords of a letter-point on post p
def levelAt(lev, p):
    ref = FunctionList
    a = (ord(lev)-ord(ref.levels[0]))/(len(ref.levels)-1) # Get portion-of-post
    b = 1-a                                       # Get unused portion
    pf, pt = p.foot, p.top
    return Point(round(b*pf.x+a*pt.x, 2), round(b*pf.y+a*pt.y, 2), round(b*pf.z+a*pt.z, 2))

def thickLet(thix):
    ref = FunctionList
    if type(thix)==float:
        return thix       # If thix is already real, return as is
    if thix=='p':
        return ref.SF*ref.pDiam
    else: # diameters q, r, s, t... scale geometrically
        expo = max(0, ord(thix)-ord('q'))
        return round(ref.SF * ref.qDiam * pow(ref.dRatio, expo), 2)

def addEdge(v,w, layout):
    if v in layout.edgeList:
        if w not in layout.edgeList[v]:
            layout.edgeList[v].append(w)
    else:
        layout.edgeList[v] = [w]

def addEdges(v,w, layout):
    ref = FunctionList
    ref.addEdge(v,w,layout); ref.addEdge(w,v,layout)
#---------------------------------------------------------
def generatePosts(code, numberTexts):
    '''Modify layout LO according to provided code and numbers'''
    ref = FunctionList
    B = ref.LO.BP
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

    def postAt(x,y,z): ref.LO.posts.append(Post(Point(x,y,z)))
    
    if code=='B':               # Set base point, BP
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:   ref.LO.BP = Point(*nums);  return

    if code=='C':               # Create a collection of posts
        nums = getNums(3,33333) # Need at least 3 numbers
        if nums:
            while len(nums) >= 3:
                postAt(nums[0]+bx,  nums[1]+by, nums[2]+bz)
                nums = nums[3:]
            if len(nums)>0:
                print (f'Anomaly: code {code}, {numberTexts} has {nums} left over')
            return

    if code=='D':      # Remove specified posts and references to them
        nums = getNums(1,33333)     # Accept any number of numbers
        if not nums: return
        lopo = ref.LO.posts;  nlop = len(lopo)
        nums = [int(x) for x in sorted(nums)]
        error = f'{nums[0]} < 0' if nums[0] < 0 else (f'{nums[-1]} >= {nlop}' if nums[-1] >= nlop else None)
        if error:
            print (f'Error: List of post numbers has a value {error} -- terminating.')
            exit(0)

        print (f'=  From {nlop} posts, deleting {len(nums)} of them: {nums}')
        # If we knew nums were distinct, we'd say polout = [0]*(nlop-len(nums))
        nin, nout, ninf, polout = 0, 0, nlop+13, [0]*nlop
        nums.append(ninf)
        transi = {} # Init post-number translation table
        # Copy posts, deleting some while making translation table
        for k in range(nlop):
            if k == nums[nin]:
                transi[k] = ninf;  nin += 1
            else:
                polout[nout] = lopo[k]
                transi[k] = nout;  nout += 1
        
        # We've deleted some posts; install into layout
        del lopo;  ref.LO.posts = polout[:nout]
        # Remove obsolete post numbers from cylinders list
        locy = ref.LO.cyls;   cylout = []
        for c in locy:
            if not(c.post1 in nums or c.post2 in nums):
                c.post1 = transi[c.post1]
                c.post2 = transi[c.post2]
                cylout.append(c)
        del locy;  ref.LO.cyls = cylout
        # If we wanted autoAdd to work ok after a D operation, at this
        # point we would clean up LO.edgeList.  But we don't care...
        return

    # This exclude-edge code is not so efficient but is good enough
    # for removing handfuls of edges, eg for doors or windows in
    # geodesics.
    #def edgecode(e,f): return max(e,f) + 262144*min(e,f)
    def edgecode(e,f): return max(e,f) + 1000*min(e,f)
    if code=='E':               # Exclude edges -- remove some cylinders
        nums = getNums(1,33333)     # Accept any number of numbers
        if not nums: return  
        pairs = [(int(x),int(y)) for x,y in zip(nums[::2],nums[1::2])]
        print (f'To remove: {pairs}')
        pairs = [edgecode(x,y) for x,y in pairs]
        drops = [];   locy = ref.LO.cyls
        for k in range(len(locy)):
            c = locy[k]
            if edgecode(c.post1, c.post2) in pairs:
                drops.append(k)
        if drops:
            drops.reverse()
            for k in drops: del locy[k]
        else:
            print (f'=  Error: None of edges {nums} found')
        return
        
    if code=='G':               # Create geodesic posts and cylinders
        from makeIcosaGeo import genIcosahedron
        nums = getNums(2,2)     # Need exactly 2 numbers
        if not nums: return
        Vfreq, gScale = int(round(nums[0])), nums[1]
        elo = Layout(posts=[], cyls=[],  edgeList={}) # Init an empty layout
        rlo = ref.LO
        # Rotation in following is not yet as advertised -- ie is normalizer not opt
        genIcosahedron(elo, Vfreq, rlo.clip1, rlo.clip2, rlo.rotavec.y, rlo.rotavec.z)
        epo = elo.posts;  eel = elo.edgeList;  nLoPo = len(rlo.posts)
        # Scale the generated posts by given scale; and copy to LO
        for p in epo:
            p.scale(gScale)
            rlo.posts.append(Post(p))
        # Generate sets of cylinders in various colors.
        colorTrans = {'Y':ref.geoColors[0], 'B':ref.geoColors[1], 'R':ref.geoColors[2], 'C':ref.geoColors[3] }
        for co in ('Y', 'B', 'R', 'C'):
            for j in sorted(eel.keys()):
                for k in sorted(eel[j]):
                    if j<k:   # Both of j,k and k,j are in the list
                        p, q = epo[j], epo[k]
                        oB = p.rank == q.rank
                        oY = p.nnbrs==5 or q.nnbrs==5
                        oC = p.dupl>1 and q.dupl>1 and not (oB or oY)
                        oR = not (oB or oY or oC)
                        if (co=='B' and oB and not oY) or (co=='Y' and oY) or (co=='R' and oR) or (co=='C' and oC):
                            cyl = Cylinder(j+nLoPo,k+nLoPo, 'c', 'c', colorTrans[co], 'p', ref.endGap, 0, 0)
                            rlo.cyls.append(cyl)
        return
        
    if code=='H':               # Create a clip box (particularly for geodesics)
        nums = getNums(6,6)     # Need exactly 6 numbers        
        if nums:
            ref.LO.clip1 = Point(*nums[:3]);
            ref.LO.clip2 = Point(*nums[3:]);
        return
        
    if code=='L':               # Create a line of posts
        nums = getNums(4,4)     # Need exactly 4 numbers
        if nums:
            n, dx, dy, dz = int(nums[0]), nums[1], nums[2], nums[3]
            x, y, z = bx, by, bz
            for k in range(n):
                postAt(x+dx, y+dy, z+dz)
            return

    if code=='O':               # Set origin point, OP
        nums = getNums(3,3)     # Need exactly 3 numbers
        if nums:   ref.LO.OP = Point(*nums);  return

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
    ref = FunctionList
    x, y, z = p.foot.x, p.foot.y, p.foot.z
    ox, oy, oz = (x, y, z-99) if ref.postAxial else (OP.x, OP.y, OP.z)
    u = ref.SF*ref.postHi       # Distance from p to post-top
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

#=========================================================
def scriptCyl(ss, preCyl):
    ref = FunctionList
    post1, post2, lev1, lev2, colo, thix, gap, nonPost, num = preCyl.get9()
    mode = 0                    # mode 0 = comments at start
    pc, code = '?', '?'
    for cc in ss:            
        if pc == '#':       # Insert a simple variable's value
            post1, post2 = post2, ref.userLocals[cc]
            nonPost = False
        elif cc in ref.colors: colo = cc
        elif cc in ref.thixx: thix  = cc
        elif cc in ref.levels:
            lev1, lev2 = lev2, cc
        elif cc in ref.digits:
            if pc in ref.digits:  post2 = post2 + cc
            else:             post1, post2 = post2, cc
            nonPost = False
        elif cc=='/':
            lev1, lev2 = lev2, lev1
        elif cc==';':
            p1, p2 = int(post1), int(post2)
            if nonPost:
                p1, p2 = p1+1, p2+1
                post1, post2 = str(p1), str(p2)
            num = len(ref.LO.cyls)
            cyl = Cylinder(p1, p2, lev1, lev2, colo, thix, gap, 0, num)
            ref.LO.cyls.append(cyl)
            addEdges(p1, p2, ref.LO) # Add edges p1,p2 and p2,p1 to edges list
            nonPost = True
        pc = cc
    preCyl.put9(post1, post2, lev1, lev2, colo, thix, gap, nonPost, num)
    return preCyl
#==================================
def scriptPost(ss, prePost):
    ref = FunctionList
    codes = 'BCDEGHILOPRST'
    pc, code, numbers = '?', '?', prePost.data
    for cc in ss:   # Process characters of script
        # Add character to number, or store a number, or what?
        if pc == '#':       # Set or use a simple variable
            if code=='?':
                ref.userLocals[cc] = len(ref.LO.posts)
            else:      # Substitute value into list of numbers
                numbers.append(ref.userLocals[cc])
        elif cc in ref.digits:
            num = num + cc if pc in ref.digits else cc
        elif pc in ref.digits:
            numbers.append(num) # Add number to list of numbers
        # Process a completed entry, or start a new entry?
        if cc==';':
            ref.generatePosts(code, numbers)
            code = '?'
        elif cc in codes:
            pc, code, numbers = '?', cc, []
        pc = cc             # Prep to get next character
    prePost.data = numbers
    return
#==================================
def runScript(scripts):
    ref = FunctionList
    preCyl = Cylinder(0, 1, 'c','c', 'G', 'p', ref.endGap, True, 0)
    prePost = Post(0, data=[])
    mode = 0                    # mode 0 = comments at start
    numbers = []
    for line in scripts:
        l1, l2, ss, ll = line[:1], line[:2], line[2:], line
        if   l2=='=C': mode = 'C'; ll=ss # Cylinders
        elif l2=='=L': mode = 'L'; ll=ss # Layout
        
        elif l2=='=P':          # Process Parameters line
            ref.installParams((ss, ref.paramTxt));
            continue
        elif l2=='=A':          # Process Arithmetic line
            ref.arithmetic(ss, ref.traceExec);
            continue
        elif l1=='=':           # Process comment line
            continue
        
        if   mode == 'L':       # Process Posts line
            ref.scriptPost(ll, prePost)
        elif mode == 'C':       # Process Cylinders line
            ref.scriptCyl (ll, preCyl)
#==================================
def writePosts(fout):
    ref = FunctionList
    try:
        ref.LO.OP.scale(ref.SF) # Get ready to orient the posts: scale the OP
    except:
        print ('In exception, dir(ref): ', [x for x in dir(PD) if not x.startswith('__')])
    # Scale the set of posts, and compute their tops and angles 
    pHi, pDi = ref.SF*ref.postHi, ref.SF*ref.postDiam
    for k, p in enumerate(ref.LO.posts):
        p.num = k
        if isTrue(ref.zSpread):
            zrat = 2/(1+p.foot.z/ref.zSize) # assumes z centers at z==0
            p.foot.scalexy(zrat)
        p.foot.scale(ref.SF)
        if isTrue(ref.postList):
            print (f'p{k:<2}=Point( {p.foot})')
        p.diam, p.hite = pDi, pHi
        p.top, p.yAngle, p.zAngle = postTop(p, ref.LO.OP)
    fout.write(f'//  onePost(num, diam, hite, yAngle, zAngle, foot x y z, top x y z);\n')
    for p in ref.LO.posts:
        fout.write(f'    onePost({p.num}, {p.diam}, {p.hite}, {p.yAngle}, {p.zAngle}, {p.foot}, {p.top});\n')

#==================================
def writeLabels(fout):
    ref = FunctionList
    if not isTrue(ref.postLabel):
        return
    fout.write(f'\n//  oneLabel(num, diam, hite, yAngle, zAngle, foot x y z, label x y z, tColor, tSize, txt);\n')
    for p in ref.LO.posts:
        cName = ref.colorSet['B']
        thik  = ref.thickLet('t')
        lxyz  = ref.levelAt('e', p)
        for cc in ref.postLabel:
            if cc in ref.colors: cName = ref.colorSet[cc]
            if cc in ref.thixx:  thik  = ref.thickLet(cc)
            if cc in ref.levels: lxyz  = ref.levelAt(cc, p)
        txt = f'{str(p.num)}'
        fout.write(f'    oneLabel({p.num}, {p.diam}, {p.hite}, {p.yAngle}, {p.zAngle}, {p.foot}, {lxyz}, {cName}, {thik}, "{txt}");\n')

#=====================================
            
def writeCylinders(fout, clo, chi, listIt):
    ref = FunctionList
    fout.write('\n//  oneCyl (p1,2,cylDiam,cylLen,yAngle,zAngle,  c xyz,  e xyz, cColor)\n')
    posts = ref.LO.posts
    nPosts = len(posts)
    for nCyl in range(clo, chi):
        cyl = ref.LO.cyls[nCyl]     # Draw this cylinder
        post1, post2, lev1, lev2, colo, thix, gap, data, num = cyl.get9()
        gap = ref.SF*gap            # gap needs scaling
        p1, p2 = min(post1,nPosts-1), min(post2,nPosts-1)
        try:
            pp, qq = posts[p1], posts[p2]
        except:
            print (f'Fatal Error with p1= {p1},   p2= {p2},  nPosts {nPosts}')
            exit(0)
        p = ref.levelAt(lev1, pp)
        q = ref.levelAt(lev2, qq)
        dx, dy, dz = q.diff(p)
        L = round(max(0.1, sssq(dx,  dy,  dz)), 2)
        cName = ref.colorSet[colo]
        alpha = gap/L
        cc = Point(p.x+alpha*dx, p.y+alpha*dy, p.z+alpha*dz)
        cylNum= 1000*p1 + p2
        if isTrue(listIt):
            print (f'Make {cyl}  L {L:2.2f}  {cName}')
        yAngle = round((pi/2 - asin(dz/L)) * 180/pi, 2)
        zAngle = round( atan2(dy, dx)      * 180/pi, 2)
        fout.write(f'    oneCyl({p1}, {p2}, {cyl.diam}, {round(L-2*gap,3)}, {yAngle}, {zAngle}, {cc}, {pp.foot}, {cName});\n')

#-------------------------------------------------------------
def autoAdder(fout):    # See if we need to auto-add cylinders
    ref = FunctionList
    rlo = ref.LO
    cutoff = ref.autoMax
    clo = len(rlo.cyls) # Record how many cylinders are already processed
    nPosts = len(rlo.posts)
    edgeList = rlo.edgeList
    # in this version punt color, thix, levels ...
    colo, thix, lev1, lev2 = 'B', 'p', 'c','c'
            
    if cutoff > 0:     # See if any way for any more edges
        print (f'In auto-add, cutoff distance autoMax is {cutoff:7.3f}')
        cutoff2 = cutoff*cutoff
        for pn in range(nPosts):
            p = rlo.posts[pn].foot
            for qn in range(1+pn, nPosts):
                q = rlo.posts[qn].foot
                dx, dy, dz = p.diff(q)
                if abs(dx) > cutoff or abs(dy) > cutoff:
                    continue
                d2 = ssq(dx, dy, dz)
                if d2 > cutoff2: continue
                if pn not in edgeList or qn not in edgeList[pn]:
                    post1, post2 = str(pn), str(qn)          
                    cyl = Cylinder(pn,qn, lev1, lev2, colo, thix, ref.endGap, 0,0)
                    rlo.cyls.append(cyl)
                    ref.addEdges(pn, qn, ref.LO)
        writeCylinders(fout, clo, len(rlo.cyls), ref.autoList)
#-------------------------------------------------------------
def installParams(script):
    '''Given script lines that are Parameter-setting lines, this extracts
    variable names and values from it, like "var1=val1 var2=val2
    var3=val3 ...".  It converts the values to numeric forms, and
    stores them in globals() dict.  Note, white space is taboo within
    a var=val item.    '''
    ref = FunctionList
    flubs = False
    for parTxt in script:
        plist = parTxt.split()       # Split the list on white space
        for vev in plist:
            p = pq = vev.split('=')  # Split each equation on = sign
            q, ok = '', False
            if len(pq)==2:
                p, q = pq
                #if p in glob.keys():
                if p in dir(ref):
                    t, v = type(getattr(ref, p)), q
                    try:
                        if   t==int:   v = int(q);     ok=True
                        elif t==float: v = float(q);   ok=True
                        elif t==bool:  v = isTrue(q);  ok=True
                        elif t==str:   v = q;          ok=True
                    except:  pass
            if ok: setattr(ref,p,v)
            else:  flubs = True
        if flubs: print (f'Parameter-setting fail in "{parTxt}"')
#-------------------------------------------------------------
def loadScriptFile(fiName):   # Read scripts from file
    with open(fiName) as fi:
        return fi.readlines()

def tell():
    return (addEdge, addEdges, arithmetic, autoAdder, generatePosts,
            installParams, levelAt, loadScriptFile, postTop,
            runScript, scriptCyl, scriptPost, thickLet,
            writeCylinders, writeLabels, writePosts)
