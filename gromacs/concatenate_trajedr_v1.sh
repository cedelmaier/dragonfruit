#!/bin/bash

rm -f list_trr_oneline.txt
rm -f times_trr.txt
rm -f list_edr_oneline.txt
rm -f times_edr.txt

for i in {901..1000}; do
  echo -n " step7_${i}.trr" >> list_trr_oneline.txt;
  echo "c" >> times_trr.txt;
  echo -n " step7_${i}.edr" >> list_edr_oneline.txt;
  echo "c" >> times_edr.txt;
done
