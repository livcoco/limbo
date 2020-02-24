// eg-roundend.scad -- OpenSCAD code for oneCyl example.
// Modules here override some modules in pipeVue.codeBase.scad

// Add a ball at each end of cylinder

module oneCyl (p1=0, p2=0, cylDiam=0, cylLen=0, yAngle=0, zAngle=0,
               cx=0, cy=0, cz=0, ex=0, ey=0, ez=0, cColor="Blue") {
  translate (v=[cx, cy, cz]) {
    rotate(a=[0, yAngle, zAngle]) {
      color(c=cColor) {
        cylinder(d=cylDiam, h=cylLen);
        sphere(d=cylDiam);
        translate ([0,0,cylLen]) {sphere(d=cylDiam);}
      }
    }
  }
}
