// scad helper code for use by eg-auto-test-2's addOn and subOff
// examples.  Modules here override 2 modules in pipeVue.codeBase.scad

// Add a sphere centered at OP, of radius r=userPar0, which can be
// varied in OpenSCAD's Customizer panel  
module addOn (SF=0, Bx=0, By=0, Bz=0, Ox=0, Oy=0, Oz=0) {
  translate([Ox*SF, Oy*SF, Oz*SF]) { sphere(r = userPar0);
  }
}

// Deduct stuff with negative z values
module subOff (SF=0, Bx=0, By=0, Bz=0, Ox=0, Oy=0, Oz=0) {
  
  translate([0,0,-userPar1/2]) { cube(size=userPar1 , center=true);
  }
}
