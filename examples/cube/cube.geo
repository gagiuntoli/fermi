lc=100;
nn=10;

Point(1) = {0, 0, 0, 1};
Point(2) = {lc, 0, 0, 1};
Point(3) = {0, lc, 0, 1};
Point(4) = {lc, lc, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Transfinite Line {1} = nn Using Progression 1;
Transfinite Line {2} = nn Using Progression 1;
Transfinite Line {3} = nn Using Progression 1;
Transfinite Line {4} = nn Using Progression 1;

Transfinite Surface {6};
Recombine Surface{6};
out[] = Extrude{0,0,lc}{
   Surface {6};  
   Layers{nn-1}; 
   Recombine;
};

Physical Volume("MAT1") = {out[1]};

Physical Surface("SURF") = {6,23,27,28,19,15}; 

Color Goldenrod {Surface {6};}

Mesh.Light = 0;
General.SmallAxes = 0;

