#!/usr/bin/env python3

import argparse
import os
import stat
import sys

import numpy as np

'''
Name: create_run_scripts_beta.py
Description: Create files for running membranes on longleaf
Input: To view options type create_run_scripts_beta.py
Output: Bash files for submission of scripts
'''

def parse_args():
    parser = argparse.ArgumentParser(prog='create_run_scripts.py')

    # Which files to create for submission
    parser.add_argument('--start', type=int, default=1,
            help='Starting number of files')
    parser.add_argument('--end', type=int, default=10,
            help='Ending number')
    parser.add_argument('--stride', type=int, default=10,
            help='Stride number per submission')

    # Threading and GPU settings for submission
    parser.add_argument('--ntmpi', type=int, default=4,
            help='Number of MPI tasks')
    parser.add_argument('--ntomp', type=int, default=20,
            help='Number of OMP threads per MPI task')
    parser.add_argument('--ngpu', type=int, default=4,
            help='Number of GPUs to use')


    # Input timing information for the number of ns/day possible
    parser.add_argument('--nsday', type=float, default=50.0,
            help='ns/day throughput')

    opts = parser.parse_args()

    return opts

###############################
if __name__ == "__main__":
    opts = parse_args()
    # Figure out the number of files
    sidx = np.int32(opts.start)
    eidx = np.int32(opts.end)
    stride = np.int32(opts.stride)
    ntmpi = np.int32(opts.ntmpi)
    ntomp = np.int32(opts.ntomp)
    ngpu = np.int32(opts.ngpu)

    # Calculate the total time needed
    nsday = np.float32(opts.nsday)
    total_time_hr = np.int32(np.ceil(np.float32(stride)*24.0/nsday*1.2))

    nfiles = np.int32(np.int32(eidx - sidx + 1) / stride)

    print(f"Creating run with the following parameters")
    print(f"  Start: {sidx}, End: {eidx}, Stride: {stride}")
    print(f"  ntmpi: {ntmpi}, ntomp: {ntomp}, ngpu: {ngpu}")
    print(f"  nsday: {nsday} --> time: {total_time_hr}")

    for ifile in range(nfiles):
        start_idx = sidx + ifile*stride
        end_idx = start_idx + stride - 1
        next_start_idx = sidx + (ifile+1)*stride
        next_end_idx = next_start_idx + stride - 1

        # Write to a bash file
        fname = 'run_prod_stride{}_v{}_{}.sh'.format(stride, start_idx, end_idx)
        fname_next = 'run_prod_stride{}_v{}_{}.sh'.format(stride, next_start_idx, next_end_idx)

        with open(fname, 'w') as rsh:
            rsh.write('''\
#!/bin/bash

# Comments for running on the GPU cluster
#SBATCH --job-name=membrane_prod_beta
#SBATCH --gres=gpu:{}
#SBATCH --ntasks={}
#SBATCH --cpus-per-task={}
#SBATCH --mem=8G
#SBATCH --time={}:00:00
#SBATCH --partition=beta-gpu
#SBATCH --output=log.%x.%j
#SBATCH --qos=gpu_access
#SBATCH --constraint=rhel8

unset OMP_NUM_THREADS

module purge
module load gcc/9.1.0
module load cuda/11.4
source /nas/longleaf/apps/gromacs/2021.5/avx2_256-cuda11.4/bin/GMXRC.bash

# This is an attempt at getting the GROMACS equilibration of the membrane right

minit=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration
prod_prefix=step7_production
prod_step=step7

lscpu

# Production, do another 10 ns
let cnt={}
let cntmax={}

while [ $cnt -le $cntmax ] ; do
  pcnt=$((cnt-1))
  istep=${{prod_step}}_${{cnt}}
  pstep=${{prod_step}}_${{pcnt}}

  if [ $cnt -eq 1 ]
  then
    pstep=`printf ${{equi_prefix}} 6`
    gmx_gpu grompp -f ${{prod_prefix}}.mdp -o ${{istep}}.tpr -c ${{pstep}}.gro -p topol.top -n index.ndx -maxwarn 1
  else
    gmx_gpu grompp -f ${{prod_prefix}}.mdp -o ${{istep}}.tpr -c ${{pstep}}.gro -t ${{pstep}}.cpt -p topol.top -n index.ndx -maxwarn 1
  fi
  gmx_gpu mdrun -v -deffnm ${{istep}} -ntmpi {} -ntomp {} -nb gpu -bonded gpu -pme gpu -npme 1
  ((cnt++))
done

# Set up the daisy-chain to next run
sbatch {}
'''.format(ngpu, ntmpi, ntomp, total_time_hr, start_idx, end_idx, ntmpi, ntomp, fname_next))

        st = os.stat(fname)
        os.chmod(fname, st.st_mode | stat.S_IEXEC)



