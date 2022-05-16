#!/bin/bash

startstep=$1
endstep=$2

rm -f list_trr_oneline.txt
rm -f times_trr.txt
rm -f list_edr_oneline.txt
rm -f times_edr.txt

for ((i=$startstep;i<=$endstep;i++))
do
  echo -n " step7_${i}.trr" >> list_trr_oneline.txt;
  echo "c" >> times_trr.txt;
  echo -n " step7_${i}.edr" >> list_edr_oneline.txt;
  echo "c" >> times_edr.txt;
done

# Now concatenate the files automatically for the user
gmx trjcat -f `cat list_trr_oneline.txt` -o traj_continuous_v$startstep\_$endstep.xtc -settime < times_trr.txt
gmx eneconv -f `cat list_edr_oneline.txt` -o edr_continuous_v$startstep\_$endstep.edr -settime < times_edr.txt
