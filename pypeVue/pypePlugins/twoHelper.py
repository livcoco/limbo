# twoHelper.py, a python replacement for eg-two-helper.scad, with its
# OpenSCAD code for onePost and oneCyl examples.  In the plugins
# milieu, we override writeCylinders and writePosts instead of onePost
# and oneCyl.  - jiw 27 Mar 2020

from math import sqrt, pi, cos, sin, asin, atan2
from pypeVue import ssq, sssq, rotate2, isTrue 
from pypeVue import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList as ref
#---------------------------------------------------------
def writePosts(fout):
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
        p.top, p.yAngle, p.zAngle = ref.postTop(p, ref.LO.OP)

    fout.write('''
// Make a three-color post with specified colors
module flagPost (diam, high, co1, co2, co3) {
  eps=0.05;  h3=high/3;
  color(co1) cylinder(d=diam-eps, h=3*h3);
  color(co2) cylinder(d=diam,     h=2*h3);
  color(co3) cylinder(d=diam+eps, h=1*h3);
}

module canePost (diam, high, twirl, co1, co2) {
  eps=0.05;  sides=8*ceil(diam)+3;
  color(co1)
    linear_extrude(height=high, twist=twirl, $fn=sides)
      difference() {
        circle(d=diam, $fn=sides);
        square([diam/2,diam*2], center=true);
      }
  color(co2) cylinder(d=diam-eps, h=high+eps, $fn=sides);
}

// Decorate posts according to post number
module decoratedPost(num, diam, high) {
  if   (num<2) color("Red") {
                  cylinder(d1=0, d2=diam*2, h=high*1.1, $fn=5);
                  cylinder(d2=0, d1=diam*2, h=high*1.1, $fn=5);
               }
  else if (num<4) color("Magenta") cylinder(d=diam, h=high, $fn=4);
  else if (num<6) flagPost(diam, high, "Green", "Cyan", "Blue"); 
  else if (num<9) flagPost(diam, high, "Red", "White", "Green");
  else            canePost(diam, high, 300, "Red", "White");
}

module onePost (postNum, postDiam, postHi, yAngle, zAngle, px, py, pz) {
  translate (v=[px, py, pz]) {
    rotate(a=[0, yAngle, zAngle]) {
      decoratedPost(postNum, postDiam, postHi);
    }
  }
}

module makePosts() {
''')
    for p in ref.LO.posts:
        fout.write(f'''  onePost({p.num}, {p.diam}, {p.hite}, {p.yAngle:7.3f}, {p.zAngle:7.3f},   {p.foot});
''')
    fout.write('}\n')           # close the module

#---------------------------------------------------------
def writeCylinders(fout, clo, chi, listIt):
    posts = ref.LO.posts
    nPosts = len(posts)
    fout.write('''
function hash3(x) =
  let (a = (x+PI-3)*PI+userPar0, b=a*PI, c=b*PI)
  [a-floor(a), b-floor(b), c-floor(c)];

module roundEnds (cylDiam, cylLen, trans, rota) {
  translate(trans)
    rotate(rota) {
      cylinder(d=cylDiam, h=cylLen);
      sphere(d=cylDiam);
      translate ([0,0,cylLen]) sphere(d=cylDiam);
    }
}

// Add a ball at each end of pipe, and color some pipes per angles
module oneCyl (p1, p2, cylDiam, cylLen, yAngle, zAngle,
               cx, cy, cz, cColor="Blue") {
  if (p1==10)
    translate([cx, cy, cz])
      rotate([0, yAngle, zAngle])
        canePost (cylDiam, cylLen, 700, cColor, "Red");
  else
      color(c = hash3(yAngle-zAngle))
      roundEnds (cylDiam, cylLen, [cx, cy, cz], [0, yAngle, zAngle]);
}

module makeCylinders() {\n''')
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
        fout.write(f'    oneCyl({p1}, {p2}, {cyl.diam}, {round(L-2*gap,3)}, {yAngle}, {zAngle}, {cc}, {cName});\n')

    fout.write('}\n')           # close the module

#---------------------------------------------------------
def tell(): return (writePosts, writeCylinders)
#---------------------------------------------------------