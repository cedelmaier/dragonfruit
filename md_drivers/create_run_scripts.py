#!/usr/bin/env python3

import argparse
import os
import stat
import sys

import numpy as np

from pathlib import Path

'''
Name: create_run_scripts.py
Description: Create files for running membranes on different systems
Input: To view options type create_run_scripts.py
Output: Bash files for submission of scripts

Example useage for rusty at FI:
    python3 <create_run_scripts.py> --engine gromacs --runtype cpuplumed --partition rome,ib --cluster rusty --start 1 --end 20 --stride 1 --nnodes 10 --ntmpi 120 --ntomp 1 --gfile step7_production.mdp --pfile plumed.dat

Example useage for popeye at FI:
    python3 <create_run_scripts.py> --engine gromacs --runtype cpuplumed --partition icelake,ib --cluster popeye --start 1 --end 20 --stride 1 --nnodes 10 --ntmpi 60 --ntomp 1 --gfile step7_production.mdp --pfile plumed.dat
'''

def parse_args():
    parser = argparse.ArgumentParser(prog='create_run_scripts.py')

    # What MD engine are we running, and what system are we using?
    parser.add_argument('--engine', type=str, required=True, choices=['gromacs', 'namd'],
            help='MD Engine')
    # What is the runtype we are using?
    parser.add_argument('--runtype', type=str, required=True, choices=['cpu', 'gpu', 'cpuplumed', 'gpuplumed'],
            help='Run type')
    # What is the partition that we are using?
    parser.add_argument('--partition', type=str, required=True, choices=['volta-gpu', 'beta-gpu', 'rome,ib', 'icelake,ib'],
            help='Partition to use (unique), constraint on FI computers')
    # What is the location we are running?
    parser.add_argument('--cluster', type=str, required=True, choices=['rusty', 'popeye', 'longleaf'],
            help='Cluster name')

    # How much runtime to create
    parser.add_argument('--start', type=int, default=1,
            help='Starting number of files')
    parser.add_argument('--end', type=int, default=10,
            help='Ending number')
    parser.add_argument('--stride', type=int, default=10,
            help='Stride number per submission')

    # Architecture information
    parser.add_argument('--nnodes', type=int, default=1,
            help='Number of Nodes')
    parser.add_argument('--ntmpi', type=int, default=1,
            help='Number of MPI tasks')
    parser.add_argument('--ntomp', type=int, default=1,
            help='Number of OMP threads per MPI task')
    parser.add_argument('--ngpu', type=int, default=1,
            help='Number of GPUs to use')

    # Control settings
    parser.add_argument('--gfile', type=str, default='step7_production.mdp',
            help='GROMACS .mdp file')
    parser.add_argument('--pfile', type=str, default='plumed.dat',
            help='PLUMED .dat file')

    # Input timing information for the number of ns/day possible
    parser.add_argument('--nsday', type=float, default=50.0,
            help='ns/day throughput')

    # Add verbosity control
    parser.add_argument('-v', '--verbose', action="store_true",
            help='Verbose output')

    opts = parser.parse_args()

    return opts

def configure_timing_information(stride, nsday):
    r"""Get timing information for simulation runs
    """
    return np.int32(np.ceil(np.float32(stride)*24.0/nsday*1.2))

def configure_cluster(mcluster, mruntype, mpartition, ntmpi, ntomp):
    r"""Configure cluster parameters for this run
    """
    gmx_exec    = ""    # GMX executable
    gmx_grompp  = ""    # GMX grompp executable
    gmx_src     = ""    # GMX source path
    gmx_options = ""    # GMX options (after commands, for like MD run specifics)
    module_list = []    # Module list for loading after purge

    # Check on the different combinations of things, like plumed, etc
    if mcluster == "longleaf":
        if mruntype == "gpu":
            gmx_exec = "gmx_gpu"
            gmx_grompp = "gmx_gpu"
            gmx_src  = "/nas/longleaf/apps/gromacs/2021.5/avx2_256-cuda11.4/bin/GMXRC"
            gmx_options = "-ntmpi {} -ntomp {} -nb gpu -bonded gpu -pme gpu -npme 1".format(ntmpi, ntomp)
            module_list.append('module load gcc/9.1.0')
            module_list.append('module load cuda/11.4')

        elif mruntype == "gpuplumed":
            gmx_exec = "gmx_gpu_plumed"
            gmx_grompp = "gmx_gpu_plumed"
            gmx_src  = "/nas/longleaf/apps/gromacs/2021.5/avx2_256-cuda11.4/bin/GMXRC"
            gmx_options = "-ntmpi {} -ntomp {}".format(ntmpi, ntomp) # This version doesn't like MPI with plumed
            module_list.append('module load gcc/9.1.0')
            module_list.append('module load cuda/11.4')

    elif mcluster == "rusty":
        if mruntype == "cpu":
            gmx_exec = "mpirun --map-by socket:pe=$OMP_NUM_THREADS -np {} gmx_mpi".format(ntmpi*nnodes)
            gmx_grompp = "mpirun -np 1 gmx_mpi"
            gmx_options = "-ntomp $OMP_NUM_THREADS"
            module_list.append('module purge')
            module_list.append('export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH')
            module_list.append('module load modules/2.1-alpha2')
            module_list.append('module load cuda')
            module_list.append('module load openmpi/cuda-4')
            module_list.append('module load plumed/mpi-2.8.1')
            module_list.append('module load fftw')
            module_list.append('module load gromacs/mpi-plumed-2022.3')
        elif mruntype == "cpuplumed":
            gmx_exec = "mpirun gmx_mpi"
            gmx_grompp = "mpirun -np 1 gmx_mpi"
            module_list.append('module load modules/2.0-20220630')
            module_list.append('module load openmpi/4.0.7')
            module_list.append('module load gromacs/mpi-plumed-2021.4')
            module_list.append('module load plumed/mpi-2.8.0')

    elif mcluster == 'popeye':
        if mruntype == "cpu":
            # Popeye command has the mpi number specified overall, so use a multiple
            gmx_exec = "mpirun --map-by socket:pe=$OMP_NUM_THREADS -np {} gmx_mpi".format(ntmpi*nnodes)
            gmx_grompp = "mpirun -np 1 gmx_mpi"
            gmx_options = "-ntomp $OMP_NUM_THREADS"
            module_list.append('module load modules/2.0-20220630')
            module_list.append('module load openmpi/4.0.7')
            module_list.append('module load gromacs/mpi-plumed-2021.4')
            module_list.append('module load plumed/mpi-2.8.0')
        elif mruntype == "cpuplumed":
            # Popeye command has the mpi number specified overall, so use a multiple
            gmx_exec = "mpirun --map-by socket:pe=$OMP_NUM_THREADS -np {} gmx_mpi".format(ntmpi*nnodes)
            gmx_grompp = "mpirun -np 1 gmx_mpi"
            gmx_options = "-ntomp $OMP_NUM_THREADS"
            module_list.append('module load modules/2.0-20220630')
            module_list.append('module load openmpi/4.0.7')
            module_list.append('module load gromacs/mpi-plumed-2021.4')
            module_list.append('module load plumed/mpi-2.8.0')

    return [gmx_exec, gmx_grompp, gmx_src, gmx_options, module_list]

# Create a header for UNC longleaf cluster
def create_header_longleaf(mpartition, ntmpi, ntomp, ngpu, total_time_hr, module_list, gmx_src):
    r"""Create a header for UNC longleaf cluster
    """
    # Slurm control parameters
    outstr = '''\
#!/bin/bash

# Comments for running on the UNC longleaf cluster
#SBATCH --job-name=grun
#SBATCH --partition={}
#SBATCH --ntasks={}
#SBATCH --cpus-per-task={}
#SBATCH --gres=gpu:{}
#SBATCH --qos=gpu_access
#SBATCH --mem=8G
#SBATCH --time={}:00:00
#SBATCH --constraint=rhel8

unset OMP_NUM_THREADS

module purge
'''.format(mpartition, ntmpi, ntomp, ngpu, total_time_hr)

    # Module information
    for module in module_list:
        outstr = outstr + "{}\n".format(module)

    # Get the source information for longleaf correct
    outstr = outstr + "source {}\n".format(gmx_src)

    return outstr

# Create a header for FI rusty cluster
def create_header_rusty(mpartition, nnodes, ntomp, module_list, gmx_src):
    r"""Create a header for FI rusty cluster
    """
    # Slurm control parameters
    outstr = '''\
#!/bin/bash

# Comments for running on the FI rusty cluster
#SBATCH --job-name=grun
#SBATCH --nodes={}
#SBATCH --cpus-per-task={}
#SBATCH --constraint={}
#SBATCH --partition=ccb

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
'''.format(nnodes, ntomp, mpartition)

    # Module information
    for module in module_list:
        outstr = outstr + "{}\n".format(module)

    return outstr

# Create a header for FI popeye cluster
def create_header_popeye(mpartition, nnodes, ntomp, module_list, gmx_src):
    r"""Create a header for FI popeye cluster
    """
    # Slurm control parameters
    outstr = '''\
#!/bin/bash

# Comments for running on the FI popeye cluster
#SBATCH --job-name=grun
#SBATCH --partition=ccb
#SBATCH --nodes={}
#SBATCH --cpus-per-task={}
#SBATCH --constraint={}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
'''.format(nnodes, ntomp, mpartition)

    # Module information
    for module in module_list:
        outstr = outstr + "{}\n".format(module)

    return outstr

# Write a block of plumed data files
# XXX: Note, this is specific to an analysis from CJE on 7/27/2022, probably need to modify for self
def write_plumed_files_block(start_idx, end_idx, pfile_strip):
    r"""Write plumed data files
    """
    for idx in range(start_idx, end_idx+1):
        fname_plumed = '{}_v{}.dat'.format(pfile_strip, idx)
        with open(fname_plumed, 'w') as stream:
            stream.write('''\
# vim:ft=plumed
MOLINFO STRUCTURE=reference.pdb

# Figure out the lipid and helix COM coordinates
lipid_com: COM ATOMS=313-99999
helix_com: COM ATOMS=1-312

# Get the Z distance
z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS

# Get the alpha value
alpha: ALPHARMSD RESIDUES=1-18

# Print to a file
PRINT ARG=* FILE=colvar_measure_{}.dat STRIDE=100
'''.format(idx))

    return

###############################
if __name__ == "__main__":
    opts = parse_args()
    if opts.verbose:
        print(f"    INFO: opts: {opts}")

    # Get the engine we are using (just gromacs for now)
    mdengine = opts.engine
    if mdengine == 'namd':
        print("Currently NAMD is just in testing mode, exiting!")
        sys.exit(1)

    # Which cluster are we configuring for?
    mcluster = opts.cluster

    # Run type information, for if we are doing plumed or not
    mruntype = opts.runtype

    # Partition information as well, to decide where we run
    mpartition = opts.partition

    # Convert opts to values
    sidx = np.int32(opts.start)
    eidx = np.int32(opts.end)
    stride = np.int32(opts.stride)
    nsday = np.float32(opts.nsday)

    nnodes = np.int32(opts.nnodes)
    ntmpi = np.int32(opts.ntmpi)
    ntomp = np.int32(opts.ntomp)
    ngpu = np.int32(opts.ngpu)

    gfile = opts.gfile
    gfile_strip = Path(gfile).stem
    gfile_base = gfile_strip.split('_')[0]

    pfile = opts.pfile
    pfile_strip = Path(pfile).stem

    # Configure the cluster information
    [gmx_exec, gmx_grompp, gmx_src, gmx_options, module_list] = configure_cluster(mcluster, mruntype, mpartition, ntmpi, ntomp)
    if opts.verbose:
        print(f"    INFO: gmx_exec: {gmx_exec}")
        print(f"    INFO: gmx_grompp: {gmx_grompp}")
        print(f"    INFO: gmx_src: {gmx_src}")
        print(f"    INFO: gmx_options: {gmx_options}")
        print(f"    INFO: module_list: {module_list}")
    # Also configure if we are doing plumed or not
    do_plumed = False
    plumed_str = ""
    if mruntype == "gpuplumed" or mruntype == "cpuplumed":
        do_plumed = True
        plumed_str = "-plumed ${{pfilenow}}".format()

    total_time_hr = configure_timing_information(stride, nsday)

    # Figure out how many files we are making
    nfiles = np.int32(np.int32(eidx - sidx + 1) / stride)
    nplumed = np.int32(np.int32(eidx -sidx + 1))

    # Print the information for the runs
    print(f"Creating run(s) with the following parameters")
    print(f"  Cluster:   {mcluster}")
    print(f"  Engine:    {mdengine}")
    print(f"  Run type:  {mruntype}")
    print(f"  Partition: {mpartition}")
    print(f"  ######## Run step information ########")
    print(f"  Start: {sidx}, End: {eidx}, Stride: {stride}")
    print(f"  ######## Architecture information ########")
    print(f"  N nodes:  {nnodes}")
    print(f"  N MPI:    {ntmpi}")
    print(f"  N OpenMP: {ntomp}")
    if mruntype == "gpu" or mruntype == "gpuplumed":
        print(f"  N GPU:    {ngpu}")

    ################
    # Actual creation of the bash script(s)
    ################
    rsh = ''
    if mcluster == "longleaf":
        rsh = create_header_longleaf(mpartition, ntmpi, ntomp, ngpu, total_time_hr, module_list, gmx_src)
    elif mcluster == "rusty":
        rsh = create_header_rusty(mpartition, nnodes, ntomp, module_list, gmx_src)
    elif mcluster == "popeye":
        rsh = create_header_popeye(mpartition, nnodes, ntomp, module_list, gmx_src)

    # The middle of the script isn't different
    rsh = rsh + '''\

minit=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration
prod_prefix={}
prod_step={}

'''.format(gfile_strip, gfile_base)

    # Figure out how many files we are going to generate
    for ifile in range(nfiles):
        start_idx       = sidx + ifile*stride
        end_idx         = start_idx + stride - 1
        next_start_idx  = sidx + (ifile+1)*stride
        next_end_idx    = next_start_idx + stride - 1

        # Write to a bash file
        fname = 'run_prod_stride{}_v{}_{}.sh'.format(stride, start_idx, end_idx)
        fname_next = 'run_prod_stride{}_v{}_{}.sh'.format(stride, next_start_idx, next_end_idx)

        # Set the start/end frames, switch to a loop interal string
        final_rsh = rsh + '''\
# Production, do another 10 ns
let cnt={}
let cntmax={}

'''.format(start_idx, end_idx)

        # If we need plumed, include it
        if do_plumed:
            final_rsh = final_rsh + '''\
# Do the plumed file as well
pfile={}

'''.format(pfile_strip)

        # Create the internal loop
        final_rsh = final_rsh + '''\
while [ $cnt -le $cntmax ] ; do
  pcnt=$((cnt-1))
  istep=${{prod_step}}_${{cnt}}
  pstep=${{prod_step}}_${{pcnt}}
'''.format()

        # If we need plumed, include it
        if do_plumed:
            final_rsh = final_rsh + '''\
  pfilenow=${{pfile}}_v${{cnt}}.dat
'''.format()

        final_rsh = final_rsh + '''\

  if [ $cnt -eq 1 ]
  then
    pstep=`printf ${{equi_prefix}} 6`
    {} grompp -f ${{prod_prefix}}.mdp -o ${{istep}}.tpr -c ${{pstep}}.gro -p topol.top -n index.ndx
  else
    {} grompp -f ${{prod_prefix}}.mdp -o ${{istep}}.tpr -c ${{pstep}}.gro -t ${{pstep}}.cpt -p topol.top -n index.ndx
  fi
  {} mdrun {} -deffnm ${{istep}} {}
  ((cnt++))
done

# Set up the daisy-chain to next run
sbatch {}
'''.format(gmx_grompp, gmx_grompp, gmx_exec, plumed_str, gmx_options, fname_next)

        if opts.verbose:
            print("-------- Next file --------")
            print(final_rsh)

        with open(fname, 'w') as stream:
            stream.write(final_rsh)
        st = os.stat(fname)
        os.chmod(fname, st.st_mode | stat.S_IEXEC)

        # Now, if we need to, write out the plumed files
        if do_plumed:
            write_plumed_files_block(start_idx, end_idx, pfile_strip)
