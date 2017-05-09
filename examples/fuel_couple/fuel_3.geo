lc=100;
nn=10;
r =25;
Mesh.CharacteristicLengthMin = 5.00;

Point(1) = {0, 0, 0, 1};
Point(2) = {lc, 0, 0, 1};
Point(3) = {0, lc, 0, 1};
Point(4) = {lc, lc, 0, 1};

Point(5) = {0+lc/2, lc/2, 0, 1};
Point(6) = {0+lc/2+r, lc/2, 0, 1};
Point(7) = {0+lc/2, lc/2+r, 0, 1};
Point(8) = {0+lc/2-r, lc/2, 0, 1};
Point(9) = {0+lc/2, lc/2-r, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
Line Loop(9) = {7, 8, 5, 6};
Plane Surface(10) = {9};

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {9, 11};
Extrude {0, 0, 100} {
  Surface{12,10};
}
Physical Volume("FLUID") = {1};
Physical Volume("FUEL")  = {2};
Physical Surface("LAT_SURF") = {53, 41, 45, 49};
Physical Surface("BOT_SURF") = {12, 10};
Physical Surface("TOP_SURF") = {54, 76};
