// eg-roundend.scad -- OpenSCAD code for onePost and oneCyl examples.

// In this example, modules onePost and oneCyl override modules in
// pipeVue.codeBase.scad.  All the other modules and functions here
// are convenience functions.  jiw 24 Feb 2020

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

module onePost (postNum=0, postDiam=0, postHi=0, yAngle=0, zAngle=0,
                px=0, py=0, pz=0, tx=0, ty=0, tz=0) {
  translate (v=[px, py, pz]) {
    rotate(a=[0, yAngle, zAngle]) {
      decoratedPost(postNum, postDiam, postHi);
    }
  }
}

function pipix(x,k) = (x+PI-3)*pow(PI,k);
function fracn(x,k) = pipix(x,k) - floor(pipix(x,k));
function hash3(x) = [fracn(x,1), fracn(x,2), fracn(x,3)];

module roundEnds (cylDiam, cylLen, trans, rota) {
  translate(trans)
    rotate(rota) {
      cylinder(d=cylDiam, h=cylLen);
      sphere(d=cylDiam);
      translate ([0,0,cylLen]) sphere(d=cylDiam);
    }
}

// Add a ball at each end of pipe, and color some pipes per angles
module oneCyl (p1=0, p2=0, cylDiam=0, cylLen=0, yAngle=0, zAngle=0,
               cx=0, cy=0, cz=0, ex=0, ey=0, ez=0, cColor="Blue") {
  if (p1==10)
    translate([cx, cy, cz])
      rotate([0, yAngle, zAngle])
        canePost (cylDiam, cylLen, 700, cColor, "Red");
  else
    color(c = yAngle<92? cColor : hash3(zAngle))
      roundEnds (cylDiam, cylLen, [cx, cy, cz], [0, yAngle, zAngle]);
}
