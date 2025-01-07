set title "3D Surface Plot"
set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Value"

set dgrid3d 30,30
set hidden3d
set pm3d
set palette rgb 33,13,10

splot "data.txt" using 1:2:3 with pm3d

pause -1 "Press Enter to exit"

#set terminal pngcairo
#set output "plot.png"
#replot