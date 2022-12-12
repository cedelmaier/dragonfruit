# XXX: Put a license here

""" Class to keep track of cluster topology """

import os
import random
import sys
import yaml

import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

class ClusterTopology(object):
    def __init__(self, yaml, opts):
        self.verbose = opts.verbose

        if self.verbose: print(f"ClusterTopology::init")

        # We aren't a 'seed' so have to setup the default yaml ourselves
        self.default_yaml = yaml
        # Read the data from said YAML file
        self.ReadData()

        if self.verbose: print(f"ClusterTopology::init return")

    def ReadData(self):
        r""" Read the YAML data for the cluster configuration
        """
        if self.verbose: print(f"ClusterTopology::ReadData")

        self.jobname    = self.default_yaml['jobname']
        self.cluster    = self.default_yaml['cluster']
        self.partition  = self.default_yaml['partition']
        self.nnodes     = np.int32(np.float64(self.default_yaml['n_nodes']))
        self.ntmpi      = np.int32(np.float64(self.default_yaml['n_mpi']))
        self.ntomp      = np.int32(np.float64(self.default_yaml['n_omp']))
        self.time       = self.default_yaml['total_time']

        if self.verbose: print(f"ClusterTopology::ReadData return")

    def PrintInformation(self):
        r""" Print out the cluster information
        """
        if self.verbose: print(f"ClusterTopology::PrintInformation")

        print(f"--------")
        print(f"Cluster topology")
        print(f"Job name                = {self.jobname}")
        print(f"Cluster                 = {self.cluster}")
        print(f"Partition               = {self.partition}")
        print(f"N Nodes                 = {self.nnodes}")
        print(f"N MPI threads           = {self.ntmpi}")
        print(f"N OpenMP threads/tMPI   = {self.ntomp}")
        print(f"Total time              = {self.time}")

        if self.verbose: print(f"ClusterTopology::PrintInformation return")

    def ConfigureCluster(self):
        r""" Configure the cluster paramters for runtime

        This keeps track of the type of executable, and some common options that we are going to be using. Also keeps
        track of the 'runtype' for the time being, might move this out later!
        """
        if self.verbose: print(f"ClusterTopology::ConfigureCluster")

        # This is probably not the best way to do this, but come up with something in the future for how
        # to control the large combination of paramters that we have seen for configuring the
        # cluster, etc

        self.sbatch_options = []

        if self.cluster == "longleaf":
            if self.partition == "snp":
                self.sbatch_options.append(f"#SBATCH --job-name={self.jobname}")
                self.sbatch_options.append(f"#SBATCH --partition={self.partition}")
                self.sbatch_options.append(f"#SBATCH --ntasks={self.ntmpi}")
                self.sbatch_options.append(f"#SBATCH --cpus-per-task={self.ntomp}")
                self.sbatch_options.append(f"#SBATCH --time={self.time}")
                self.sbatch_options.append(f"#SBATCH --qos=snp_access")
            else:
                print(f"ERROR: Partition {self.partition} not yet implemented for cluster: {self.cluster}, exiting")
                sys.exit(1)

        # Configure the cluster, header and all
        sh_str = (
            f'#!/bin/bash\n\n'
        )
        for line in self.sbatch_options:
            sh_str = sh_str + "{}\n".format(line)

        self.sbatch_header = sh_str

        if self.verbose: print(f"ClusterTopology::ConfigureCluster return")
