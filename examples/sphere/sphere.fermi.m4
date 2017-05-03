$Mesh 
    mesh_file   sphere.msh
$EndMesh 

$Mode
  timedep QSTATIC 
  p0 2.5e6  
$EndMode

$Xs
  egn 1
  # cross sections
  #         F D      XA    nXF    eXF    CHI
  "FUEL"    0 1.72   XA_m4 nXF_m4 5.4e-6 1.0  
$EndXs

$Boundary
  "SURF" 1 0
$EndBoundary
