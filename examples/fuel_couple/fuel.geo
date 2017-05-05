lc=100;
nn=10;
r =5;
Mesh.CharacteristicLengthMin = 15.00;

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

Translate {25, 0, 0} {
  Duplicata { Surface{10}; }
}
Translate {-25, 0, 0} {
  Duplicata { Surface{10}; }
}
Translate {-25, 25, 0} {
  Duplicata { Surface{10}; }
}
Translate {0, 25, 0} {
  Duplicata { Surface{10}; }
}
Translate {25, 25, 0} {
  Duplicata { Surface{10}; }
}
Translate {-25, -25, 0} {
  Duplicata { Surface{10}; }
}
Translate {0, -25, 0} {
  Duplicata { Surface{10}; }
}
Translate {25, -25, 0} {
  Duplicata { Surface{10}; }
}

//Extrude {0, 0, 100} {
//  Surface{21, 26, 31, 11, 46, 41, 36, 16, 10};
//}
//Physical Volume("fuel_z1") = {7, 6, 5, 8, 9, 4, 1, 2, 3};
//Extrude {0, 0, 100} {
//  Surface{72, 226, 204, 182, 160, 138, 248, 94, 116};
//}
//Physical Volume("fuel_z2") = {10, 17, 18, 11, 16, 15, 12, 13, 14};
//Extrude {0, 0, 100} {
//  Surface{446, 424, 270, 292, 402, 380, 358, 336, 314};
//}
//Physical Volume("fuel_z3") = {21, 20, 19, 24, 25, 26, 23, 22, 27};
//Extrude {0, 0, 100} {
//  Surface{468, 490, 512, 534, 556, 578, 600, 622, 644};
//}
//Physical Volume("fuel_z4") = {30, 31, 36, 35, 32, 29, 28, 33, 34};
//Extrude {0, 0, 100} {
//  Surface{666, 776, 798, 688, 754, 820, 842, 732, 710};
//}
//Physical Volume("fuel_z5") = {37, 40, 45, 44, 43, 41, 42, 38, 39};
//Line Loop(1041) = {1, 2, 3, 4};
//Line Loop(1042) = {40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39, 40, 37, 38, 39};
//Plane Surface(1043) = {1041, 1042};
//
//fluid_1[]=Extrude {0, 0, 100} {
//  Surface{1043};
//};
//Physical Volume("fluid_z1") = {fluid_1[1]};
//
//Line Loop(1246) = {1048, 1045, 1046, 1047};
//Line Loop(1247) = {55, 52, 53, 54};
//Line Loop(1248) = {76, 77, 74, 75};
//Line Loop(1249) = {99, 96, 97, 98};
//Line Loop(1250) = {120, 121, 118, 119};
//Line Loop(1251) = {142, 143, 140, 141};
//Line Loop(1252) = {206, 207, 208, 209};
//Line Loop(1253) = {187, 184, 185, 186};
//Line Loop(1254) = {164, 165, 162, 163};
//Line Loop(1255) = {229, 230, 231, 228};
//Plane Surface(1256) = {1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255};
//fluid_2[]=Extrude {0, 0, 100} {
//  Surface{1245};
//};
//Physical Volume("fluid_z2") = {fluid_2[1]};
//
//Line Loop(1459) = {1258, 1259, 1260, 1261};
//Line Loop(1460) = {297, 294, 295, 296};
//Line Loop(1461) = {272, 273, 274, 275};
//Line Loop(1462) = {250, 251, 252, 253};
//Line Loop(1463) = {407, 404, 405, 406};
//Line Loop(1464) = {316, 317, 318, 319};
//Line Loop(1465) = {385, 382, 383, 384};
//Line Loop(1466) = {426, 427, 428, 429};
//Line Loop(1467) = {363, 360, 361, 362};
//Line Loop(1468) = {341, 338, 339, 340};
//Plane Surface(1469) = {1459, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468};
//fluid_3[]=Extrude {0, 0, 100} {
//  Surface{1458};
//};
//Physical Volume("fluid_z3") = {fluid_3[1]};
//Line Loop(1672) = {1472, 1473, 1474, 1471};
//Line Loop(1673) = {449, 450, 451, 448};
//Line Loop(1674) = {560, 561, 558, 559};
//Line Loop(1675) = {582, 583, 580, 581};
//Line Loop(1676) = {603, 604, 605, 602};
//Line Loop(1677) = {625, 626, 627, 624};
//Line Loop(1678) = {537, 538, 539, 536};
//Line Loop(1679) = {515, 516, 517, 514};
//Line Loop(1680) = {471, 472, 473, 470};
//Line Loop(1681) = {493, 494, 495, 492};
//Plane Surface(1682) = {1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 1680, 1681};
//fluid_4[]=Extrude {0, 0, 100} {
//  Surface{1671};
//};
//Physical Volume("fluid_z4") = {fluid_4[1]};
//Line Loop(1885) = {1685, 1686, 1687, 1684};
//Line Loop(1886) = {780, 781, 778, 779};
//Line Loop(1887) = {802, 803, 800, 801};
//Line Loop(1888) = {824, 825, 822, 823};
//Line Loop(1889) = {647, 648, 649, 646};
//Line Loop(1890) = {670, 671, 668, 669};
//Line Loop(1891) = {691, 692, 693, 690};
//Line Loop(1892) = {757, 758, 759, 756};
//Line Loop(1893) = {735, 736, 737, 734};
//Line Loop(1894) = {713, 714, 715, 712};
//Plane Surface(1895) = {1885, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894};
//fluid_5[]=Extrude {0, 0, 100} {
//  Surface{1884};
//};
//Physical Volume("fluid_z5") = {fluid_5[1]};
//Coherence;
//Coherence;
//Physical Surface("TOP_SURF") = {2097, 1018, 1040, 930, 864, 886, 908, 974, 996, 952};
//Physical Surface("BOTTOM_SURF") = {36, 41, 46, 10, 11, 1043, 16, 21, 26, 31};
//Physical Surface("LAT_SURF") = {1526, 1100, 1313, 1739, 1514, 1088, 1727, 1301, 1952, 1940, 1944, 1731, 1518, 1305, 1092, 1096, 1309, 1522, 1735, 1948};
//Coherence;
Line Loop(51) = {4, 1, 2, 3};
Line Loop(52) = {17, 18, 19, 20};
Line Loop(53) = {22, 23, 24, 25};
Line Loop(54) = {27, 28, 29, 30};
Line Loop(55) = {32, 33, 34, 35};
Line Loop(56) = {12, 13, 14, 15};
Line Loop(57) = {37, 38, 39, 40};
Line Loop(58) = {42, 43, 44, 45};
Line Loop(59) = {47, 48, 49, 50};
Plane Surface(60) = {9, 51, 52, 53, 54, 55, 56, 57, 58, 59};
Extrude {0, 0, 100} {
  Surface{60, 16, 21, 26, 31, 10, 11, 36, 41, 46};
}
Coherence;
Physical Volume("fuel_z1") = {6, 2, 3, 4, 5, 7, 10, 9, 8};
Physical Volume("fluid_z1") = {1};
Physical Surface("BOTTOM_SURF") = {11, 60, 46, 41, 36, 16, 10, 31, 26, 21};
Physical Surface("TOP_SURF") = {262, 306, 328, 350, 394, 460, 372, 438, 416, 284};
Physical Surface("LAT_SURF") = {121, 133, 129, 125};
