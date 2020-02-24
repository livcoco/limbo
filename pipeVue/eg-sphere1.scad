// eg-sphere1.scad -- OpenSCAD code for addOn and subOff examples.
// Modules here can override some modules in pipeVue.codeBase.scad

// Add a sphere centered at OP, of radius 380
// (We want r=userPar0 but userPar0 apparently isn't in scope here -
// might need to have a function in frontCode that delivers them?)
  
module addOn (SF=0, Bx=0, By=0, Bz=0, Ox=0, Oy=0, Oz=0) {
  radi=380;
  translate([Ox*SF, Oy*SF, Oz*SF]) { sphere(r = radi);
  }
}

// Deduct stuff with negative z values
module subOff (SF=0, Bx=0, By=0, Bz=0, Ox=0, Oy=0, Oz=0) {
  translate([0,0,-2*SF]) { cube(size=4*SF , center=true);
  }
}
