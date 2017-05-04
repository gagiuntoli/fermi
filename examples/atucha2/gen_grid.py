nx=62
ny=62
linex=1
lineloop=10000
surf=0
for i in range(ny):
    for j in range(nx):
	print "Point(", i*nx+j,")={", j-int(nx/2) ,"*dx/2,", i-int(ny/2) , "*dy/2," , "0, 1.0};"
        if (j!=0):
            print "Line(", linex,")={", i*nx+j-1 ,",", i*nx+j, "};"
            linex=linex+1
    if (i!=0):
        for j in range(nx):
            print "Line(", linex,")={", (i-1)*nx+j ,",", i*nx+j, "};"
            linex=linex+1
            if (j!=0 and i==1):
                print "Line Loop(", lineloop , ") = {", linex-1, "," , -(linex-nx-1) , "," , -linex+2 , "," , (linex-(2*nx)) , "};"
                print "Plane Surface(", surf, ") = {" , lineloop , "};"
                print "Color White {Surface {",surf,"};}"
                lineloop=lineloop+1
                surf=surf+1
            if (j!=0 and i!=1):
                print "Line Loop(", lineloop , ") = {", linex-1, "," , -(linex-nx-1) , "," , -linex+2 , "," , (linex-(3*nx)) , "};"
                print "Plane Surface(", surf, ") = {" , lineloop , "};"
                print "Color White {Surface {",surf,"};}"
                lineloop=lineloop+1
                surf=surf+1

linex=1
for i in range(ny):
    for j in range(nx):
        if (j!=0):
            print "Transfinite Line {", linex,"} = nx Using Progression 1 ;"
            linex=linex+1
    if (i!=0):
        for j in range(nx):
            print "Transfinite Line {", linex,"} = ny Using Progression 1 ;"
            linex=linex+1

for i in range(surf):
    print "Transfinite Surface { ", i,"};" 
    print "Recombine Surface {   ", i,"};" 
