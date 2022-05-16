#!/bin/bash

# Get the distances between the lipid COM and the AH domain

echo 0 | gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o conf.gro -sep

# compute distances
for (( i=0; i<501; i++ ))
do
    gmx distance -s step7_umbrella_v1.tpr -f conf${i}.gro -n index.ndx -select 'com of group "MEMB" plus com of group "SOLU"' -oall dist${i}.xvg 
done

# compile summary
touch summary_distances.dat
for (( i=0; i<501; i++ ))
do
    d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
    echo "${i} ${d}" >> summary_distances.dat
    rm dist${i}.xvg
done

exit;
