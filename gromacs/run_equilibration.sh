#!/bin/bash

# Comments for running on the GPU cluster
#SBATCH --job-name=membrane_perf
#SBATCH --partition=volta-gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_access
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=15:00:00
#SBATCH --output=equil_run-%j.log

unset OMP_NUM_THREADS

module purge
module load gcc/9.1.0
module load cuda/9.2
source /nas/longleaf/apps/gromacs/2020.3/root-cuda9.2/bin/GMXRC.bash

# This is an attempt at getting the GROMACS equilibration of the membrane right

minit=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration
prod_prefix=step7_production
prod_step=step7

# Have to set the number of MPI threads to 1 so that it doesn't blow up
#ntmpi=1

lscpu

# See if we can run the precursor steps
gmx_gpu_mpi grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${minit}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx
gmx_gpu_mpi mdrun -v -deffnm ${mini_prefix}

# Equilibration only!
let cnt=1
let cntmax=6

while [ $cnt -le $cntmax ] ; do
  pcnt=$((cnt-1))
  istep=`printf ${equi_prefix} ${cnt}`
  pstep=`printf ${equi_prefix} ${pcnt}`

  if [ $cnt -eq 1 ]
  then
    pstep=${mini_prefix}
  fi
  gmx_gpu_mpi grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx
  gmx_gpu_mpi mdrun -v -deffnm ${istep} -ntmpi 3 -ntomp 8 -nb gpu -bonded gpu -pme gpu -npme 1
  ((cnt++))
done






