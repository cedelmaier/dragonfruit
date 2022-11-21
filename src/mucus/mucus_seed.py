# XXX: Put a license here

"""Class for a single mucus seed"""

import os
import random
import sys
import yaml

import gsd.hoomd
import numpy as np
import pandas as pd

# Analysis imports
import freud
from freud import box

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from seed_base import SeedBase
from common import radial_average, ragged_mean

class MucusSeed(SeedBase):
    def __init__(self, path, opts):
        self.verbose = opts.verbose
        print(f"MucusSeed::__init__")

        SeedBase.__init__(self, path, opts)

        self.ReadData()

        # Set the RNG state
        random.seed(self.nseed)

        print(f"MucusSeed::__init__ return")

    def ReadData(self):
        r""" Read the data from a YAML file for a single seed
        """
        print(f"MucusSeed::ReadData")
        # What engine are we going to be using?
        self.engine             = self.default_yaml['simulation']['engine']
        self.partition          = self.default_yaml['simulation']['partition']

        # Simulation parameters
        self.kT                 = np.float64(self.default_yaml['simulation']['kT'])
        self.t_damp             = np.float64(self.default_yaml['simulation']['t_damp'])
        self.bead_size          = np.float64(self.default_yaml['simulation']['bead_size'])
        self.deltatau           = np.float64(self.default_yaml['simulation']['deltatau'])
        self.nsteps_equilibrate = np.int32(np.float64(self.default_yaml['simulation']['nsteps_equilibrate']))
        self.nwrite_equilibrate = np.int32(np.float64(self.default_yaml['simulation']['nwrite_equilibrate']))
        self.nsteps             = np.int32(np.float64(self.default_yaml['simulation']['nsteps']))
        self.nwrite             = np.int32(np.float64(self.default_yaml['simulation']['nwrite']))
        self.lbox               = np.float64(self.default_yaml['simulation']['lbox'])
        self.nseed              = np.int32(np.float64(self.default_yaml['simulation']['seed']))

        self.nrowy              = np.int32(np.float64(self.default_yaml['mucus']['nrowy']))
        self.nrowz              = np.int32(np.float64(self.default_yaml['mucus']['nrowz']))
        self.lbond              = np.float64(self.default_yaml['mucus']['bond_length'])
        self.monomer_length     = np.int32(np.float64(self.default_yaml['mucus']['monomer_length']))
        self.n_term_length      = np.int32(np.float64(self.default_yaml['mucus']['n_term_length']))
        self.backbone_length    = np.int32(np.float64(self.default_yaml['mucus']['backbone_length']))
        self.n_cysteine         = np.int32(np.float64(self.default_yaml['mucus']['n_cysteine']))
        self.cysteine_locations = self.default_yaml['mucus']['cysteine_locations']
        self.mucin_charges      = self.default_yaml['mucus']['charges']
        self.dimers_per_poly    = np.int32(np.float64(self.default_yaml['mucus']['dimers_per_polymer']))

        self.nhist              = np.int32(np.float64(self.default_yaml['histones']['n']))
        self.histone_charges    = np.float64(self.default_yaml['histones']['charge'])

        if self.engine != "LAMMPS":
            print(f"ERROR: Only LAMMPS implementation currently supported for mucus, exiting!")
            sys.exit(1)

        print(f"MucusSeed::ReadData return")

    def PrintInformation(self, snap = None):
        r""" Print information about sim
        """
        print(f"--------")
        print(f"System information")
        print(f"kBT                     = {self.kT}")
        print(f"Temperature damping     = {self.t_damp}")
        print(f"Bead size               = {self.bead_size}")
        print(f"Delta tau               = {self.deltatau}")
        print(f"Nsteps equilibrate      = {self.nsteps_equilibrate}")
        print(f"Nwrite equilibrate      = {self.nwrite_equilibrate}")
        print(f"Nsteps                  = {self.nsteps}")
        print(f"Nwrite                  = {self.nwrite}")
        print(f"Simulation time (tau)   = {self.deltatau*self.nsteps}")
        print(f"Box size                = {self.lbox}")
        print(f"Seed                    = {self.nseed}")
        print(f"--------")
        print(f"Mucus polymer information")
        print(f"N polymers              = {self.n_mucins}")
        print(f"Monomer length          = {self.monomer_length}")
        print(f"Start N-terminus        = {self.start_n_terminus}")
        print(f"Start Backbone          = {self.start_backbone}")
        print(f"Start C-terminus        = {self.start_c_terminus}")
        print(f"N cysteines             = {self.n_cysteine}")
        print(f"  Cysteine locations      = {self.cysteine_locations}")
        print(f"Mucin charges           = {self.mucin_charges}")
        print(f"Nbeads per dimer (adj)  = {self.nper_dimer}")
        print(f"Nbeads per poly (adj)   = {self.nper_poly}")
        print(f"Total mucin beads       = {self.nbeads_mucin}")
        print(f"--------")
        print(f"Free histone (PCLS) information")
        print(f"N histones              = {self.nhist}")
        print(f"Histone charges         = {self.histone_charges}")
        print(f"--------")
        print(f"Total configured system")
        print(f"N types                 = {self.ntypes}")
        print(f"N beads                 = {self.natoms}")
        print(f"N bonds                 = {self.nbonds}")
        print(f"N angles                = {self.nangles}")
        print(f"Charges (per type)      = {self.charges}")

    def Configure(self, snap = None):
        r""" Configure self for running
        """
        print(f"MucusSeed::Configure")

        # Figrue out how many mucins to put in the system
        self.n_mucins           = self.nrowy * self.nrowz

        self.start_n_terminus   = 0
        self.start_backbone     = self.start_n_terminus + self.n_term_length
        self.start_c_terminus   = self.start_backbone + self.backbone_length
        self.c_term_length      = self.monomer_length - self.start_c_terminus
        self.nper_dimer         = 2*self.monomer_length - self.c_term_length

        self.nper_poly          = self.nper_dimer * self.dimers_per_poly
        self.nbeads_mucin       = self.nper_poly * self.n_mucins

        # Calculate the total number of atoms
        self.natoms             = self.nbeads_mucin + self.nhist
        self.nbonds             = self.n_mucins*(self.nper_poly - 1)
        self.nangles            = self.n_mucins*(self.nper_poly - 2)
        self.charges            = self.mucin_charges.copy()
        self.ntypes             = 3
        if self.nhist > 0:
            self.ntypes += 1
            self.charges.append(self.histone_charges)
        self.nbondstypes        = 1
        self.nanglestypes       = 1

        # Create the configuration file
        self.CreateLAMMPSConfig()

        # Create the equilibration runfile
        self.CreateEquilibrationScript()

        print(f"MucusSeed::Configure return")

    def CreateLAMMPSConfig(self):
        r""" Create a LAMMPS position file for this run
        """
        print(f"MucusSeed::CreateLAMMPSConfig")
        # Generate the LAMMPS configuration file
        lammps_str = (
            f'Autogenerated from dragonfruit\n'
            f'\n'
            f'{self.natoms} atoms\n'
            f'{self.ntypes} atom types\n'
            f'{self.nbonds} bonds\n'
            f'{self.nbondstypes} bond types\n'
            f'{self.nangles} angles\n'
            f'{self.nanglestypes} angle types\n'
            f'\n'
            f'  {-self.lbox/2.0} +{self.lbox/2.0} xlo xhi\n'
            f'  {-self.lbox/2.0} +{self.lbox/2.0} ylo yhi\n'
            f'  {-self.lbox/2.0} +{self.lbox/2.0} zlo zhi\n'
            f'\n'
            f'Masses\n\n'
        )
        # Loop to get masses
        for itype in range(self.ntypes):
            lammps_str += f'  {itype+1} 1.0\n'
        lammps_str += '\n'

        # Generate all of the positions, etc, for the beads (yuck)
        lammps_str += f'Atoms\n\n'
        dy1 = self.lbox/(2.0*self.nrowy)
        dz1 = self.lbox/(2.0*self.nrowz)
        icount = 0
        ichaincount = 0
        for iz in range(self.nrowz):
            z = dz1*(iz-self.nrowz/2.0 + 0.5)
            for jy in range(self.nrowy):
                y = dy1*(jy-self.nrowy/2.0 + 0.5)
                x = -90.0
                itype = 3
                icount += 1
                ichaincount += 1

                lammps_str += f'{icount} {ichaincount} {itype} {self.charges[itype-1]} {x} {y} {z} 0 0 0\n'

                for k in range(1, self.nper_poly):
                    icount += 1
                    x = x + self.lbond
                    itype = 1
                    if (k % self.nper_dimer < self.n_term_length): itype = 3
                    if (k % self.nper_dimer > (self.nper_dimer - self.n_term_length - 1)): itype = 3
                    if (k % self.nper_dimer < self.monomer_length) and (k % self.nper_poly >= self.monomer_length - self.c_term_length): itype = 3
                    for l in range(self.n_cysteine):
                        if (k % self.nper_dimer == self.cysteine_locations[l]): itype = 2
                        if ((self.nper_poly-k-1)%(self.nper_dimer) == self.cysteine_locations[l]): itype = 2

                    lammps_str += f'{icount} {ichaincount} {itype} {self.charges[itype-1]} {x} {y} {z} 0 0 0\n'

        # Now generate the histone positions
        beyondchaincount = ichaincount
        for ihist in range(self.nhist):
            icount += 1
            beyondchaincount += 1
            x = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            y = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            z = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            itype = 4
            lammps_str += f'{icount} {beyondchaincount} {itype} {self.charges[itype-1]} {x} {y} {z} 0 0 0\n'

        lammps_str += f'\n'

        # Bonds
        lammps_str += f'Bonds\n\n'
        icount1 = 1
        icount2 = 1
        for k in range(self.n_mucins):
            for n in range(self.nper_poly):
                if (n < self.nper_poly - 1):
                    lammps_str += f'{icount2} 1 {icount1} {icount1+1}\n'
                    icount2 += 1
                icount1 += 1

        lammps_str += f'\n'

        # Angles
        lammps_str += f'Angles\n\n'
        icount1 = 1
        icount2 = 1
        for k in range(self.n_mucins):
            for n in range(self.nper_poly):
                if (n < self.nper_poly - 2):
                    lammps_str += f'{icount2} 1 {icount1} {icount1+1} {icount1+2}\n'
                    icount2 += 1
                icount1 += 1

        lammps_str += f'\n'

        lammps_filename = f'nmucin{self.n_mucins}_{self.nrowy}x{self.nrowz}_mucus.lammps_config'
        with open(lammps_filename, 'w') as stream:
            stream.write(lammps_str)

        print(f"MucusSeed::CreateLAMMPSConfig return")

    def CreateEquilibrationScript(self):
        r""" Create an equilibration script
        """
        print(f"MucusSeed::CreateEquilibrationScript")

        sh_str = (
            f'#!/bin/bash\n\n'
            f'#SBATCH --job-name=lmp_mucus_equilibration\n'
        )

        # Get the other run options
        if self.partition == 'snp':
            sh_str += self.LongleafSNPHeader("04:00:00")

        print(sh_str)

        print(f"MucusSeed::CreateEquilibrationScript return")

    def LongleafSNPHeader(self, runtime):
        r""" Create an SNP header for ourselves on longleaf
        """
        sh_str =  (
            f'#SBATCH --partition=snp\n'
            f'#SBATCH --ntasks=128\n'
            f'#SBATCH --cpus-per-task=1\n'
            f'#SBATCH --time={runtime}\n'
            f'#SBATCH --qos=snp_access\n'
            f'\n'
        )
        return sh_str


