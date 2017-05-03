#!/bin/bash
# Script for determing the critical mass of a U235 sphere

rho=19.05    # [g/cm3]
MU=235.0      # [g/mol]
Na=6.022E23  # avogadro
xsA=1.5E-24  # microscopic fission cross section
nu=2.4       # neutrons emmited from fission
R=(2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0)
Es=0.2 

rm crit.dat
touch crit.dat
echo "# R    Mc    Keff" >> crit.dat

for i in ${R[@]}
do
  V=`echo "4 3 3.1415 $i" | awk '{printf "%e\n", $1/$2*$3*$4*$4*$4}'`
  M=`echo "$V $rho" | awk '{printf "%e\n", $1 * $2}'`
  N=`echo "$rho $Na $MU" | awk '{printf "%e\n", $1 * $2 / $3}'`
  XSA=`echo "$N $xsA" | awk '{printf "%e\n", $1 * $2}'`
  nXSF=`echo "$N $nu $xsA" | awk '{printf "%e\n", $1 * $2 * $3}'`
  echo "R="$i "M="$M "N="$N "XSA="$XSA "nXSF="$nXSF
  E=`echo "$Es $i" | awk '{printf "%e\n", $1 * $2}'`
  m4 -DR_m4=$i -DE_m4=$E sphere.geo.m4 > sphere.geo
  m4 -DXA_m4=$XSA -DnXF_m4=$nXSF sphere.fermi.m4 > sphere.fermi
  gmsh -3 sphere.geo > gmsh.out
  ../../fermi  sphere.fermi > sphere.out
  k=`awk '/Keff/ {print $2} ' sphere.out`
  echo "$i $M $k" >> crit.dat
done
pyxplot pyxplot.dat

