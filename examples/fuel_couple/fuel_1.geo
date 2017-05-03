lc=100;
nn=10;
r =10;
Mesh.CharacteristicLengthMin = 4.00;

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

Translate {35, 0, 0} {
  Duplicata { Surface{10}; }
}
Translate {-35, 0, 0} {
  Duplicata { Surface{10}; }
}
Translate {-35, 35, 0} {
  Duplicata { Surface{10}; }
}
Translate {0, 35, 0} {
  Duplicata { Surface{10}; }
}
Translate {35, 35, 0} {
  Duplicata { Surface{10}; }
}
Translate {-35, -35, 0} {
  Duplicata { Surface{10}; }
}
Translate {0, -35, 0} {
  Duplicata { Surface{10}; }
}
Translate {35, -35, 0} {
  Duplicata { Surface{10}; }
}


Line Loop(51) = {4, 1, 2, 3};
Line Loop(52) = {22, 23, 24, 25};
Line Loop(53) = {19, 20, 17, 18};
Line Loop(54) = {27, 28, 29, 30};
Line Loop(55) = {45, 42, 43, 44};
Line Loop(56) = {47, 48, 49, 50};
Line Loop(57) = {12, 13, 14, 15};
Line Loop(58) = {32, 33, 34, 35};
Line Loop(59) = {39, 40, 37, 38};
Plane Surface(60) = {9, 51, 52, 53, 54, 55, 56, 57, 58, 59};
Extrude {0, 0, 100} {
  Surface{60, 21, 16, 36, 41, 46, 11, 31, 26, 10};
}

Physical Volume("FUEL_Z1") = {6, 2, 3, 4, 5, 7, 10, 9, 8};
Physical Volume("FLUID_Z1") = {1};
Physical Surface("BOT_SURF") = {262, 306, 284, 438, 416, 394, 372, 350, 460, 328};
Physical Surface("TOP_SURF") = {60, 21, 26, 31, 11, 16, 36, 41, 46, 10};
Physical Surface("LAT_SURF") = {121, 125, 133, 129};
