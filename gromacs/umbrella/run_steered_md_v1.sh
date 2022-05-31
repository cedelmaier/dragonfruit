#!/bin/bash

# Comments for running on the GPU cluster
#SBATCH --job-name=memb_smd
#SBATCH --partition=volta-gpu
#SBATCH --gres=gpu:3
#SBATCH --qos=gpu_access
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=smd_run-%j.log
#SBATCH --constraint=rhel8

unset OMP_NUM_THREADS

module purge
module load gcc/9.1.0
module load cuda/11.4
#source /nas/longleaf/apps/gromacs/2021.5/avx2_256-cuda11.4/bin/GMXRC.bash
source /nas/longleaf/apps/gromacs/2021.5/avx_512-cuda11.4/bin/GMXRC.bash

# This is an attempt at getting the GROMACS equilibration of the membrane right

minit=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration
prod_prefix=step7_production
prod_step=step7

# Get the number of cores and the linux version/name
lscpu
uname -a

# Run the steered MD simulations in preparation for umbrela sampling
gmx_gpu grompp -f step7_smd_cylinder.mdp -o step7_smd_cylinder.tpr -c step6.6_equilibration.gro -p topol.top -n index.ndx -maxwarn 1
gmx_gpu mdrun -v -deffnm step7_smd_cylinder -ntmpi 3 -ntomp 8 -nb gpu -bonded gpu -pme gpu -npme 1

