# Helper code for use by eg-auto-test-2 example.  The function here
# overrides setCodeFrontAndBack from baseFuncs.py by adding a
# hemisphere centered at [0,0,0].  It is to prevent posts on the far
# side of a geodesic dome showing through when we view the near side.

import datetime
#-------------------------------------------------------------
def setCodeFrontAndBack(c):
    c.date = datetime.datetime.today().strftime('%Y-%m-%d  %H:%M:%S')
    c.frontCode = f'''// File {c.scadFile}, generated  {c.date}
// by pypeVue from script "{c.f}"
$fn = {c.cylSegments};
userPar0 = {c.userPar0};
userPar1 = {c.userPar1};
'''
    c.backCode = f'''
// Make an origin-centered hemisphere, of radius r=userPar0. 
// userPar0 can be varied in OpenSCAD's Customizer panel. 
// We make the hemisphere by removing a sphere's bottom half.
module makeHemi()
  intersection() {'{'}
    sphere(r = userPar0);
    translate([0,0,2*userPar0]) cube(size=4*userPar0, center=true);
  {'}'}

union() {'{'}
  makePosts();
  makeLabels();
  makeCylinders(); 
  makeHemi();
{'}'}
'''
#-------------------------------------------------------------
def tell():
    return (setCodeFrontAndBack,)
