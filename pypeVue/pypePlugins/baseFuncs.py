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
        lopo = LO.posts;  nlop = len(lopo)
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
        del lopo;  LO.posts = polout[:nout]
        # Remove obsolete post numbers from cylinders list
        locy = LO.cyls;   cylout = []
        for c in locy:
            if not(c.post1 in nums or c.post2 in nums):
                c.post1 = transi[c.post1]
                c.post2 = transi[c.post2]
                cylout.append(c)
        del locy;  LO.cyls = cylout
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
        drops = [];   locy = LO.cyls
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
        # Rotation in following is not yet as advertised -- ie is normalizer not opt
        genIcosahedron(elo, Vfreq, LO.clip1, LO.clip2, LO.rotavec.y, LO.rotavec.z)
        epo = elo.posts;  eel = elo.edgeList;  nLoPo = len(LO.posts)
        # Scale the generated posts by given scale; and copy to LO
        for p in epo:
            p.scale(gScale)
            LO.posts.append(Post(p))
        # Generate sets of cylinders in various colors.
        colorTrans = {'Y':geoColors[0], 'B':geoColors[1], 'R':geoColors[2], 'C':geoColors[3] }
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
                            cyl = Cylinder(j+nLoPo,k+nLoPo, 'c', 'c', colorTrans[co], 'p', endGap, 0, 0)
                            LO.cyls.append(cyl)
        return
        
    if code=='H':               # Create a clip box (particularly for geodesics)
        nums = getNums(6,6)     # Need exactly 6 numbers        
        if nums:
            LO.clip1 = Point(*nums[:3]);
            LO.clip2 = Point(*nums[3:]);
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
        if nums:   LO.OP = Point(*nums);  return

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

#=========================================================
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
    codes = 'BCDEGHILOPRST'
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
        fout.write(f'    oneCyl({p1}, {p2}, {cyl.diam}, {round(L-2*gap,3)}, {yAngle}, {zAngle}, {cc}, {pp.foot}, {cName});\n')

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

def tell():
    return (addEdge, addEdges, arithmetic, autoAdder, generatePosts,
            installParams, levelAt, loadScriptFile, postTop,
            runScript, scriptCyl, scriptPost, thickLet,
            writeCylinders, writeLabels, writePosts)
