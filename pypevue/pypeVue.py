#!/usr/bin/env python3

# jiw 26 Dec 2018
'''A program that generates OpenSCAD code for tubes along selected
edges between `posts` in a plane.  This supports visualization of
arrangements of edges in geodesic dome structures.

development:
  ln -s <path to project>/pypevue ~/.local/lib/python3.6/site-packages/.
  ln -s <path to project>/pypevue/pypeVue.py ~/.local/bin/pypeVue
  create a separate working directory to run in and:
  copy an example file from ~/.local/lib/python3.6/site-packages/pypevue/
  pypeVue <example>

'''

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

# to debug add:
#import web_pdb; web_pdb.set_trace()

from sys import argv, exit, exc_info, stderr
import time, datetime
from math import sqrt, pi, cos, sin, asin, atan2
from pypevue.pypePlugins import FunctionList
from os import path

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
    def cross(self, q):         # cross product of the two 3-vectors
        return (self.y*q.z - self.z*q.y, self.z*q.x - self.x*q.z, self.x*q.y - self.y*q.x)
    def norm(self):
        mag = sssq(self.x, self.y, self.z)
        return (self.x / mag, self.y / mag, self.z / mag)
    def diff(self, q):
        return (self.x-q.x, self.y-q.y, self.z-q.z)
    def add(self, q):
        return (self.x+q.x, self.y+q.y, self.z+q.z)
    def mag(self):
        return sssq(self.x, self.y, self.z)
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
        diam = FunctionList.thickLet(thix)
        self.put9 (post1, post2, lev1, lev2, colo, diam, gap, data, num)
        
    def get9(self):
        return self.post1, self.post2, self.lev1, self.lev2, self.colo, self.diam, self.gap, self.data, self.num
    
    def put9(self, post1, post2, lev1, lev2, colo, diam, gap, data, num):
        self.post1, self.post2 = post1, post2
        self.lev1,  self.lev2  = lev1,  lev2
        self.colo,  self.diam  = colo,  diam
        self.gap,   self.data  = gap,   data
        self.num = num
        
    def __str__( self):
        return f'Cylinder {self.num} ({self.post1},{self.post2}) {self.colo}{self.diam:0.2f}{self.lev1}{self.lev2} {round(self.gap,2)}'
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
def rotate2(a,b,theta):
    st = sin(theta)
    ct = cos(theta)
    return  a*ct-b*st, a*st+b*ct
#---------------------------------------------------------
def isTrue(x):
    '''Return false if x is None, or False, or an empty string, or a
    string beginning with f, F, N, or n.  Else, return True.    '''
    return not str(x)[:1] in 'fFNn'
#---------------------------------------------------------
def setupData(c):
    ref = FunctionList
    c.levels, c.thixx,  c.digits = 'abcde', 'pqrstuvw', '01234356789+-.'
    c.colorSet = {'G':'"Green"', 'Y':'"Yellow"', 'R':'"Red"', 'B':'"Blue"', 'C':'"Cyan"', 'M':'"Magenta"', 'W':'"White"', 'P':'[.5,0,.5]', 'A':'"Coral"'}
    c.colors = c.colorSet.keys()

    # Set initial values of main parameters
    c.pDiam,   c.qDiam,    c.dRatio   = 0.06, 0.02, sqrt(2)
    c.endGap,  c.postHi,   c.postDiam = 0.03, 0.16, c.qDiam
    c.f,       c.SF,    c.cylSegments = '', 100, 30
    c.paramTxt, c.postLabel= '', 'Bte' # Blue, size u, level e
    c.userCode = ''
    c.codeBase = f'pypeVue.codeBase.scad' # SCAD functions for posts, cyls, etc
    c.scadFile = f'pypeVue.scad'          # Name of scad output file
    c.postList = c.cylList = False # Control printing of post and cyl data
    c.Plugins, c.autoMax, c.autoList  = '', 0, True
    c.zSpread, c.zSize,   c.postAxial = False, 1, True
    c.userPar0 = c.userPar1 = c.userPar2 = '""'
    c.traceExec=False
    c.geoColors = 'YBRC'          # Colors for pentagons, rings, rays, seams
    c.script1 = '=P postDiam=.1 endGap=.05','=C Gpae 1,2;;;;1;Rea 1,2;;;;1;','=L C 0,0,0; P5,1,0;'
    for k in range(1,len(argv)):
        c.paramTxt = c.paramTxt + ' ' + argv[k]
    c.userLocals = {}               # Initialize empty user-space dict
    import os.path
    myname = os.path.splitext(os.path.basename(__file__))[0]
    # Add classes Point, Post, Layout, FunctionList to userLocals, and ref.
    exec(f'from {myname} import Point,Post,Layout\nfrom pypePlugins import FunctionList\nref=FunctionList', c.userLocals)
#---------------------------------------------------------
def makePluginsList(ref):
    pll = ''
    for lin in list(ref.scripts) + [ref.paramTxt]: # For each line in script,
        for s in lin.split():            # split the line on white space.
            if s.startswith('Plugins='): # If it is a plugins param,
                pll = pll + ',' + s[8:]  # add its list to the plugins list.
    return pll
#---------------------------------------------------------

if __name__ == '__main__':
    t0 = time.time()
    FunctionList.registrar('')
    ref = FunctionList
    setupData(ref) 
    ref.installParams([ref.paramTxt]) # Should set f, script-name parameter
    if ref.f == '':
        ref.scripts = ref.script1
    else:
        with open(ref.f) as fi:
            ref.scripts = fi.readlines()
    # If our command line or script names any plugins, get them registered
    FunctionList.registrar(makePluginsList(ref))
    ref.setClipAndRota(ref)   # Create LO and its clip1, clip2, rotavec vals
    ref.runScript(ref.scripts)    # Run selected script
    ref.setCodeFrontAndBack(ref)  # Set up beginning and ending SCAD code
    with open(ref.scadFile, 'w') as fout:
        fout.write(ref.frontCode)
        ref.writePosts    (fout)
        ref.writeLabels   (fout)
        ref.writeCylinders(fout, 0, len(ref.LO.cyls), ref.cylList,
                           1 if ref.autoMax>0 else 3)
        ref.autoAdder     (fout)
        fout.write(ref.backCode)
    t1 = time.time()-t0
    print (f'For script "{ref.f}", pypeVue wrote code to {ref.scadFile} at {ref.date} in {t1:0.3f} seconds')

'''
failing examples
eg-fatpost-redstar, eg-arith-test
'''
