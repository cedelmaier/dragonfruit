#!/usr/bin/env bash

mpath=$1
vfinal=$2

scp -r popeye:$1/toppar .
scp -r popeye:$1/topol.top .
scp -r popeye:$1/index.ndx .
scp -r popeye:$1/step5* .
scp -r popeye:$1/step6*.tpr .
scp -r popeye:$1/step6*.gro .
scp -r popeye:$1/step6*.cpt .
scp -r popeye:$1/step7_production.mdp .
scp -r popeye:$1/reference.pdb .

scp -r popeye:$1/*v1_$vfinal* .
scp -r popeye:$1/step7_$vfinal.tpr .
scp -r popeye:$1/step7_$vfinal.cpt .
scp -r popeye:$1/step7_$vfinal.log .
scp -r popeye:$1/step7_$vfinal.gro .

scp -r popeye:$1/*.dat .
