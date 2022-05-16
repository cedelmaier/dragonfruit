#!/usr/bin/env python3
__description__ = \
"""
Take a file containing the distance between the COM of the lipid
bilayer and the protein and generate the proper shell scripts
to process this for umbrella sampling
"""

import argparse
import os
import stat
import sys

import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(prog='setup_umbrella.py')

    # General options
    parser.add_argument('-dfile', '--distance', required = True, type = str,
            help = 'Summary file of distances')
    parser.add_argument('-i', '--interval', required = True, type = np.float32,
            help = 'Sampling interval (nm)')

    # Threading and GPU settings
    parser.add_argument('--ntmpi', type=int, default=4,
            help='Number of MPI tasks')
    parser.add_argument('--ntomp', type=int, default=20,
            help='Number of OMP threads per MPI task')
    parser.add_argument('--ngpu', type=int, default=4,
            help='Number of GPUs to use')

    # Input timing information
    parser.add_argument('--nsday', type=float, default=50.0,
            help='ns/day throughput')

    opts = parser.parse_args()
    return opts

class UmbrellaSampler(object):
    r""" Umbrella sampling class
    """
    def __init__(self, opts):
        self.opts = opts
        self.cwd = os.getcwd()

        self.ReadOpts()

        self.ProgOpts()

    def ReadOpts(self):
        # Read the distance file
        self.distance_table = self.readDistanceFile(self.opts.distance)
        self.sample_interval = self.opts.interval


    def ProgOpts(self):
        r""" Run the program options
        """
        self.sampled_indices = self.sampleDistances()

        # Now create the output files for volta-gpu
        self.generateOutputFiles()

    def readDistanceFile(self, distance_file):
        r""" Read the distance file that has been generated
        """
        f = open(distance_file, "r")
        lines = f.readlines()
        f.close()

        # Read the data from the bottom
        out_dict = {}
        for i in range(len(lines)-1,-1,-1):

            # Split on white space, grab frame/distance
            columns = lines[i].split()
            key = np.int32(columns[0])
            value = np.float32(columns[1])

            if key in out_dict:
                break
            else:
                out_dict[key] = value

        # Now reprocess
        keys = out_dict.keys()
        out = [(k,out_dict[k]) for k in sorted(out_dict)]

        return out

    def sampleDistances(self):
        r""" Sample the distances and give the frames that correspond to this
        """
        distances = [d[1] for d in self.distance_table]

        current_index = 0
        sampled_indices = [current_index]
        while current_index < len(distances):

            target_distance = distances[current_index] + self.sample_interval

            # Walk through and find where this occurs
            onward = [np.abs(target_distance-d) for d in distances[current_index:]]
            next_index = onward.index(min(onward)) + current_index

            if current_index == next_index:
                break
            else:
                sampled_indices.append(next_index)
                current_index = next_index

        return sampled_indices

    def generateOutputFiles(self):
        r""" Generate output files for ourselves for the umbrella sampling
        """
        # Get the information from opts for information numbers
        ntmpi = np.int32(self.opts.ntmpi)
        ntomp = np.int32(self.opts.ntomp)
        ngpu = np.int32(self.opts.ngpu)

        nsday = np.float32(opts.nsday)
        # Need a fake 10 runs (10ns) and then plus 2 things
        total_time_hr = np.int32(np.ceil(np.float32(12)*24.0/nsday*1.2))

        # Loop through all the possible files for this
        fname_list = []
        for i in range(len(self.sampled_indices)):
            frame = self.distance_table[self.sampled_indices[i]][0]
            dist = self.distance_table[self.sampled_indices[i]][1]

            if i >= len(self.sampled_indices)-1:
                next_frame = 'null'
            else:
                next_frame = self.distance_table[self.sampled_indices[i+1]][0]

            # Write to a bash script
            fname = 'run_umbrella_frame{}.sh'.format(frame)
            fname_next = 'run_umbrella_frame{}.sh'.format(next_frame)

            # Setup to write to a single batch file and submit all at once, as they are independent
            fname_list.append(fname)

            with open(fname, 'w') as rsh:
                rsh.write('''\
#!/bin/bash

# Comments for running on the GPU cluster
#SBATCH --job-name=membrane_umbrella
#SBATCH --partition=volta-gpu
#SBATCH --gres=gpu:{}
#SBATCH --qos=gpu_access
#SBATCH --ntasks={}
#SBATCH --cpus-per-task={}
#SBATCH --mem=8G
#SBATCH --time={}:00:00
#SBATCH --output=umbrella_run-%j.log
#SBATCH --constraint=rhel8

unset OMP_NUM_THREADS

module purge
module load gcc/9.1.0
module load cuda/11.4
source /nas/longleaf/apps/gromacs/2021.5/avx_512-cuda11.4/bin/GMXRC.bash

# Get the number of cores and the linux version/name
lscpu
uname -a

# Setup for a single umbrella run
# First, setup the NPT runs
gmx_gpu grompp -f npt_umbrella_v1.mdp -c conf{} -p topol.top -r conf{}.gro -n index.ndx -o npt{}.tpr -maxwarn 1
gmx_gpu mdrun -deffnm npt{} -ntmpi {} --ntomp {} -nb gpu -bonded gpu -pme gpu -npme 1

# Now run the actual simulations
gmx_gpu grompp -f md_umbrella_v1.mdp -c npt{}.gro -t npt{}.cpt -p topol.top -r npt{}.gro -n index.ndx -o umbrella{}.tpr
gmx_gpu mdrun -deffnm umbrella{} -ntmpi {} --ntomp {} -nb gpu --bonded gpu -pme gpu -npme 1
'''.format(ngpu, ntmpi, ntomp, total_time_hr, frame, frame, frame, frame, ntmpi, ntomp, frame, frame, frame, frame, frame, ntmpi, ntomp))

                st = os.stat(fname)
                os.chmod(fname, st.st_mode | stat.S_IEXEC)

        # Now, put everything into a single submission script
        with open('submit_all.sh', 'w') as rsh:
            rsh.write('#!/bin/bash\n\n')
            for iname in fname_list:
                rsh.write('sbatch {}\n'.format(iname))
            rsh.write('\n')

        st = os.stat('submit_all.sh')
        os.chmod('submit_all.sh', st.st_mode | stat.S_IEXEC)


###############
if __name__ == "__main__":
    opts = parse_args()
    x = UmbrellaSampler(opts)
