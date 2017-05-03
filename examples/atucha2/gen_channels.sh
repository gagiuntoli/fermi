#!/bin/bash


#directions=("nw" "ne" "se" "sw");
#
#for i in {1..451}; do
#    for j in ${directions[@]}; do
#       m4 -Dcanalm4=c$i-$j -Dsurfacem4=$i*4+$jj
#    done
#done

. channel_name.sh 
. surface_ID.sh

m4 -Dcanalm4=${channel_name[0]} -Dsurfacem4=${surface_ID[0]}  canal.geo.m4 > canal.geo
gmsh canal.geo -3 2&> null
echo `echo 0.0/${#channel_name[@]} | bc -l | awk '{printf "%.2lf", $0}'` % ${channel_name[0]} done
mv canal.msh old.msh

for i in $(seq 1 1 $((${#channel_name[@]}-1)));do  
   m4 -Dcanalm4=${channel_name[$i]} -Dsurfacem4=${surface_ID[$i]}  canal.geo.m4 > canal.geo
   gmsh canal.geo -3 2&> null
   echo `echo $i*100/${#channel_name[@]} | bc -l | awk '{printf "%.2lf", $0}'` % ${channel_name[$i]} done
   mv canal.msh new.msh
   gmshjoin old.msh new.msh
   mv new_result.msh old.msh
done

mv old.msh core.msh 