# XXX: Put a license here

"""Class for a single mucus seed"""

import os
import random
import stat
import sys
import yaml

#import gsd.hoomd
import numpy as np
import pandas as pd

# Analysis imports
#import freud
#from freud import box

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'cluster'))
from seed_base import SeedBase
from common import radial_average, ragged_mean
from cluster_topology import ClusterTopology

class MucusSeed(SeedBase):
    def __init__(self, path, opts):
        self.verbose = opts.verbose
        if self.verbose: print(f"MucusSeed::__init__")

        SeedBase.__init__(self, path, opts)

        # Initialize cluster topology to nothing
        self.cluster = None

        self.ReadData()

        # Set the RNG state
        random.seed(self.nseed)

        # Set up if we are initialized or not
        self.is_init = False
        self.getTypebyName = {}

        if self.verbose: print(f"MucusSeed::__init__ return")

    def ReadData(self):
        r""" Read the data from a YAML file for a single seed
        """
        if self.verbose: print(f"MucusSeed::ReadData")
        # What engine are we going to be using?
        self.engine             = self.default_yaml['simulation']['engine']

        # If we are using HOOMD, are we CPU or GPU computing?
        # XXX Go back and make the lammps stuff use the init_type as well!
        if self.engine == "HOOMD":
            self.compute_mode   = self.default_yaml['simulation']['compute_mode']
            self.init_type      = self.default_yaml['simulation']['init_type']
            self.mucusbend      = np.float64(self.default_yaml['interactions']['mucus_bend'])
            self.trajectory_file= self.default_yaml['simulation']['trajectory_file']
            self.init_filename  = self.default_yaml['simulation']['init_filename']

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
        if 'radius' in self.default_yaml['histones']:
            self.r_histone      = np.float64(self.default_yaml['histones']['radius'])
        else:
            self.r_histone      = 0.5

        if 'lennard_jones_ee' in self.default_yaml['interactions']:
            self.lennard_jones_ee   = np.float64(self.default_yaml['interactions']['lennard_jones_ee'])
        else:
            self.lennard_jones_ee   = 0.0
        self.lennard_jones      = np.float64(self.default_yaml['interactions']['lennard_jones'])
        self.bmh                = np.float64(self.default_yaml['interactions']['born_mayer_huggins'])

        if self.engine != "LAMMPS" and self.engine != "HOOMD":
            print(f"ERROR: Only LAMMPS and HOOMD implementation currently supported for mucus, exiting!")
            sys.exit(1)

        # Set up the cluster topology
        self.cluster = ClusterTopology(self.default_yaml['cluster'], self.opts)

        if self.verbose: print(f"MucusSeed::ReadData return")

    def PrintInformation(self, snap = None):
        r""" Print information about sim
        """
        print(f"--------")
        print(f"Engine information")
        print(f"Engine                  = {self.engine}")
        if self.engine == "HOOMD":
            print(f"Compute mode            = {self.compute_mode}")
            print(f"Init type               = {self.init_type}")
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
        print(f"N dimers per polymer    = {self.dimers_per_poly}")
        print(f"Monomer length          = {self.monomer_length}")
        print(f"Start N-terminus        = {self.start_n_terminus}")
        print(f"N-term length           = {self.n_term_length}")
        print(f"Start Backbone          = {self.start_backbone}")
        print(f"Backbone length         = {self.backbone_length}")
        print(f"Start C-terminus        = {self.start_c_terminus}")
        print(f"C-term length           = {self.c_term_length}")
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
        print(f"Histone radius          = {self.r_histone}")
        print(f"--------")
        print(f"Interactions")
        print(f"Electrostatics (brush)  = {self.lennard_jones_ee}")
        print(f"Electrostatics (LJ)     = {self.lennard_jones}")
        print(f"Hydrophobic (BMH)       = {self.bmh}")
        print(f"--------")
        print(f"Total configured system")
        print(f"N types                 = {self.ntypes}")
        print(f"N beads                 = {self.natoms}")
        print(f"N bonds                 = {self.nbonds}")
        print(f"N angles                = {self.nangles}")
        print(f"Charges (per type)      = {self.charges}")

        # Print out cluster topolgy information
        self.cluster.PrintInformation()

    def Configure(self, snap = None):
        r""" Configure self for running
        """
        if self.verbose: print(f"MucusSeed::Configure")

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

        # Master switch for running on LAMMPS vs HOOMD
        # XXX Really the correct way to do this is to use inheritance or something, but for now, just do via a switch
        if self.engine == "LAMMPS":
            # Create the configuration file
            self.CreateLAMMPSConfig()

            # Create the header for the cluster
            self.cluster.ConfigureCluster()

            # Create the equilibration runfile
            self.CreateEquilibrationScript()
            self.CreateEquilibrationControl()

            # Create the production runs
            self.CreateProductionScript()
            self.CreateProductionControl()

            if self.verbose: print(f"MucusSeed::Configure return (LAMMPS)")
            return

        elif self.engine == "HOOMD":
            # Do the HOOMD stuff!
            import gsd.hoomd
            import freud
            print(f"Configure detected HOOMD")

            # Create a default configuration that is the size of the box
            snap.configuration.box = [self.lbox, self.lbox, self.lbox, 0, 0, 0]
            if self.init_type == 'create_equilibration':
                self.InitMucus(snap)
            elif self.init_type == 'production':
                self.ReadMucus(snap)
            self.is_init = True

        if self.verbose: print(f"MucusSeed::Configure return")

    def CreateLAMMPSConfig(self):
        r""" Create a LAMMPS position file for this run
        """
        if self.verbose: print(f"MucusSeed::CreateLAMMPSConfig")
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
                    if ((k % self.nper_dimer) < self.n_term_length):
                        itype = 3
                    if ((k % self.nper_dimer) > (self.nper_dimer - self.n_term_length - 1)):
                        itype = 3
                    if (((k % self.nper_dimer) < self.monomer_length) and ((k % self.nper_dimer) >= (self.monomer_length - self.c_term_length))):
                        itype = 3
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
        self.lammps_filename = lammps_filename
        with open(lammps_filename, 'w') as stream:
            stream.write(lammps_str)

        if self.verbose: print(f"MucusSeed::CreateLAMMPSConfig return")

    def CreateEquilibrationScript(self):
        r""" Create an equilibration script
        """
        if self.verbose: print(f"MucusSeed::CreateEquilibrationScript")

        sh_str = self.cluster.sbatch_header

        # Add the module options for LAMMPS
        sh_str += (
            f'\n'
            f'module purge\n'
            f'module load lammps\n'
            f'\n'
            f'export OMP_NUM_THREADS={self.cluster.ntomp}\n'
            f'\n'
            f'CTRL_FILE={self.cluster.jobname+".Equilibrate.lammps_control"}\n'
            f'DATA_FILE={self.lammps_filename}\n'
            f'\n'
            f'SEED={self.nseed}\n'
            f'\n'
            f'ELECTROSTATIC=0.0\n'
            f'BMHA=0.0\n'
            f'BEND_K=0.0\n'
            f'\n'
            f'N_STEPS={self.nsteps_equilibrate}\n'
            f'OUTPUT_PREFIX="equilibrate"\n'
            f'\n'
        )

        # Add the actual command for lammps
        sh_str += (
            f'mpirun --map-by ppr:{self.cluster.ntmpi}:node:pe=$OMP_NUM_THREADS --bind-to hwthread --report-bindings \\\n'
            f'  lmp_mpi -sf omp -pk omp ${{OMP_NUM_THREADS}} -in ${{CTRL_FILE}} \\\n'
            f'  -var OUTPUT_PREFIX ${{OUTPUT_PREFIX}} \\\n'
            f'  -var DATA_FILE     ${{DATA_FILE}} \\\n'
            f'  -var TRAJ_EVERY    {self.nwrite_equilibrate} \\\n'
            f'  -var BMHA          ${{BMHA}} \\\n'
            f'  -var LJELECTROSTATIC  ${{ELECTROSTATIC}} \\\n'
            f'  -var BEND             ${{BEND_K}} \\\n'
            f'  -var DELTA_T          {self.deltatau} \\\n'
            f'  -var T_DAMP           {self.t_damp} \\\n'
            f'  -var SEED             ${{SEED}} \\\n'
            f'  -var N_STEPS          ${{N_STEPS}}\n'
            f'\n'
        )

        # Daisy chain the production run
        sh_str += (
            f'sbatch run_production.sh\n'
            f'\n'
        )

        with open('run_equilibrate.sh', 'w') as stream:
            stream.write(sh_str)
        st = os.stat('run_equilibrate.sh')
        os.chmod('run_equilibrate.sh', st.st_mode | stat.S_IEXEC)

        if self.verbose: print(f"MucusSeed::CreateEquilibrationScript return")

    def CreateProductionScript(self):
        r""" Create the production script for this run
        """
        if self.verbose: print(f"MucusSeed::CreateProductionScript")

        # A lot of this is duplicated from the equilibrate script
        sh_str = self.cluster.sbatch_header

        # Add the module options for LAMMPS
        sh_str += (
            f'\n'
            f'module purge\n'
            f'module load lammps\n'
            f'\n'
            f'export OMP_NUM_THREADS={self.cluster.ntomp}\n'
            f'\n'
            f'CTRL_FILE={self.cluster.jobname+".Production.lammps_control"}\n'
            f'DATA_FILE=equilibrate.final.data\n'
            f'\n'
            f'SEED={self.nseed}\n'
            f'\n'
            f'ELECTROSTATIC={self.lennard_jones}\n'
            f'BMHA={self.bmh}\n'
            f'BEND_K=0.0\n'
            f'\n'
            f'N_STEPS={self.nsteps}\n'
            f'OUTPUT_PREFIX="production"\n'
            f'\n'
        )

        # Add the actual command for lammps
        sh_str += (
            f'mpirun --map-by ppr:{self.cluster.ntmpi}:node:pe=$OMP_NUM_THREADS --bind-to hwthread --report-bindings \\\n'
            f'  lmp_mpi -sf omp -pk omp ${{OMP_NUM_THREADS}} -in ${{CTRL_FILE}} \\\n'
            f'  -var OUTPUT_PREFIX ${{OUTPUT_PREFIX}} \\\n'
            f'  -var DATA_FILE     ${{DATA_FILE}} \\\n'
            f'  -var TRAJ_EVERY    {self.nwrite} \\\n'
            f'  -var BMHA          ${{BMHA}} \\\n'
            f'  -var LJELECTROSTATIC  ${{ELECTROSTATIC}} \\\n'
            f'  -var BEND             ${{BEND_K}} \\\n'
            f'  -var DELTA_T          {self.deltatau} \\\n'
            f'  -var T_DAMP           {self.t_damp} \\\n'
            f'  -var SEED             ${{SEED}} \\\n'
            f'  -var N_STEPS          ${{N_STEPS}}\n'
            f'\n'
        )

        with open('run_production.sh', 'w') as stream:
            stream.write(sh_str)
        st = os.stat('run_production.sh')
        os.chmod('run_production.sh', st.st_mode | stat.S_IEXEC)

        if self.verbose: print(f"MucusSeed::CreateProductionScript return")

    def CreateEquilibrationControl(self):
        r""" Create the equilibration control
        """
        if self.verbose: print(f"MucusSeed::CreateEquilibrationControl")

        eq_str = r'''
log ${OUTPUT_PREFIX}.log

################
# Variables that need to be passed in via the command line
################
# OUTPUT_PREFIX     : Name of the simulation we are running
# DATA_FILE         : Data file (either original coordinates or restart)
# TRAJ_EVERY        : Trajectory output frequency
# DELTA_T           : Timestep
# T_DAMP            : Thermostat coupling strength
# N_STEPS           : Simulation length
# SEED              : Simulation RNG seed
#
# ELECTROSTATIC     : Electrostatic interactions epsilon value
# BMH               : Born-Mayer-Huggins prefactor A
# BEND              : Chain bending energy

################
# Setup variables for the whole system
################

units       lj
atom_style  full

boundary    p p p
read_data   ${DATA_FILE}

################
# Basic types and groups
################

variable muc_e  equal   1   # electrostatic regions
variable muc_c  equal   2   # cysteine knots
variable muc_h  equal   3   # hydrophobic termini
variable muc_p  equal   4   # positive junk

group muc_group type    ${muc_e} ${muc_c} ${muc_h} ${muc_p}
group muc_e_group type  ${muc_e}
group muc_p_group type  ${muc_p}

################
# Simulation output
################

thermo                          100
thermo_style                    custom step temp press etotal pe enthalpy epair pxy pxz pyz

dump                            muc_traj    muc_group   custom ${TRAJ_EVERY} ${OUTPUT_PREFIX}.muc.lammpstrj id type mol x y z
dump_modify                     muc_traj    sort id

################
# Born-Mayer-Huggins model of chained-mucus, equilibration run
################

variable SIGMA_16       equal   $(1.12246152962189)
variable Temperature    equal   $(1.0)

variable BMH_A      equal   $(v_BMHA)
variable BMH_RHO    equal   $(1.0 / 9.0)

variable LJ_E       equal   $(v_LJELECTROSTATIC)

variable FENE_K     equal   $(30.0)
variable FENE_R0    equal   $(1.5)
variable FENE_E     equal   $(1.0)
variable FENE_R     equal   $(1.0)

variable BEND_K     equal   $(v_BEND)

print ""
print "#"
print "# General"
print "# Sigma                      = $(1.0)"
print "# Sigma^{1/6}                = $(v_SIGMA_16)"
print "# Excluded Volume (WCA)"
print "# Sigma^{1/6}                = $(v_SIGMA_16)"
print "# Hydrophobic (BMH)"
print "# BMH A                      = $(v_BMH_A)"
print "# BMH rho                    = $(v_BMH_RHO)"
print "# Electrostatics"
print "# LJ epsilon                 = $(v_LJ_E)"
print "# Polymer properties"
print "# FENE K                     = $(v_FENE_K)"
print "# FENE R0                    = $(v_FENE_R0)"
print "# FENE E                     = $(v_FENE_E)"
print "# FENE R                     = $(v_FENE_R)"
print "# Bending K                  = $(v_BEND_K)"
print "#"
print ""

################
# Non-bonded interactions
################

pair_style soft 3.5
pair_coeff * * 10.0 ${SIGMA_16}
variable prefactor equal ramp(0,100)
fix var_soft_fix all adapt 1 pair soft a * * v_prefactor

################
# Bonded interactions (FENE)
################

bond_style fene
special_bonds fene
bond_coeff 1 ${FENE_K} ${FENE_R0} ${FENE_E} ${FENE_R}

################
# Angle interactions
################

angle_style harmonic
angle_coeff 1 ${BEND_K} 180.0

################
# Integration parameters
################

timestep        ${DELTA_T}
neighbor        1.0 bin
neigh_modify    every 1 delay 1 check yes

fix nve_fix     muc_group nve
fix nvt_fix     muc_group langevin ${Temperature} ${Temperature} ${T_DAMP} ${SEED}

################
# Final setup for run
################

comm_style tiled
comm_modify cutoff 2.5
fix rcb_balance all balance 10000 1.1 rcb

minimize 1.0e-4 1.0e-6 100 1000

info system communication

run ${N_STEPS}
write_data ${OUTPUT_PREFIX}.final.data nocoeff

################
# Post-processing
################
'''

        with open(self.cluster.jobname + ".Equilibrate.lammps_control", 'w') as stream:
            stream.write(eq_str)

        if self.verbose: print(f"MucusSeed::CreateEquilibrationControl return")

    def CreateProductionControl(self):
        r""" Write the Production control file
        """
        if self.verbose: print(f"MucusSeed::CreateProductionControl")

        eq_str = r'''
log ${OUTPUT_PREFIX}.log

################
# Variables that need to be passed in via the command line
################
# OUTPUT_PREFIX     : Name of the simulation we are running
# DATA_FILE         : Data file (either original coordinates or restart)
# TRAJ_EVERY        : Trajectory output frequency
# DELTA_T           : Timestep
# T_DAMP            : Thermostat coupling strength
# N_STEPS           : Simulation length
# SEED              : Simulation RNG seed
#
# ELECTROSTATIC     : Electrostatic interactions epsilon value
# BMH               : Born-Mayer-Huggins prefactor A
# BEND              : Chain bending energy

################
# Setup variables for the whole system
################

units       lj
atom_style  full

boundary    p p p
read_data   ${DATA_FILE}

################
# Basic types and groups
################

variable muc_e  equal   1   # electrostatic regions
variable muc_c  equal   2   # cysteine knots
variable muc_h  equal   3   # hydrophobic termini
variable muc_p  equal   4   # positive junk

group muc_group type    ${muc_e} ${muc_c} ${muc_h} ${muc_p}
group muc_e_group type ${muc_e}
group muc_p_group type ${muc_p}

################
# Simulation output
################

thermo                          1000
thermo_style                    custom step temp press etotal pe enthalpy epair pxy pxz pyz

dump                            muc_traj    muc_group   custom ${TRAJ_EVERY} ${OUTPUT_PREFIX}.muc.lammpstrj id type mol x y z
dump_modify                     muc_traj    sort id

################
# Born-Mayer-Huggins model of chained-mucus, equilibration run
################

variable SIGMA_16       equal   $(1.12246152962189)
variable Temperature    equal   $(1.0)

variable BMH_A      equal   $(v_BMHA)
variable BMH_RHO    equal   $(1.0 / 9.0)

variable LJ_E       equal   $(v_LJELECTROSTATIC)

variable FENE_K     equal   $(30.0)
variable FENE_R0    equal   $(1.5)
variable FENE_E     equal   $(1.0)
variable FENE_R     equal   $(1.0)

variable BEND_K     equal   $(v_BEND)

print ""
print "#"
print "# General"
print "# Sigma                      = $(1.0)"
print "# Sigma^{1/6}                = $(v_SIGMA_16)"
print "# Excluded Volume (WCA)"
print "# Sigma^{1/6}                = $(v_SIGMA_16)"
print "# Hydrophobic (BMH)"
print "# BMH A                      = $(v_BMH_A)"
print "# BMH rho                    = $(v_BMH_RHO)"
print "# Electrostatics"
print "# LJ epsilon                 = $(v_LJ_E)"
print "# Polymer properties"
print "# FENE K                     = $(v_FENE_K)"
print "# FENE R0                    = $(v_FENE_R0)"
print "# FENE E                     = $(v_FENE_E)"
print "# FENE R                     = $(v_FENE_R)"
print "# Bending K                  = $(v_BEND_K)"
print "#"
print ""

################
# Non-bonded interactions
################

# Do pair styles assuming everythingturned on for now
# Because of how LAMMPS does things, have to specify all atoms
# Order is e, c, h, p
if "${LJ_E} > 0.0" then &
    "pair_style hybrid lj/cut ${SIGMA_16} lj/cut 2.5 born 2.5" &
    "pair_modify shift yes" 

if "${LJ_E} > 0.0" then &
    "pair_coeff ${muc_e} ${muc_e} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_e} ${muc_c} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_e} ${muc_h} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_c} ${muc_c} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_c} ${muc_h} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_c} ${muc_p} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_h} ${muc_p} lj/cut 1 1.0 1.0 ${SIGMA_16}" &
    "pair_coeff ${muc_p} ${muc_p} lj/cut 1 1.0 1.0 ${SIGMA_16}" &

if "${LJ_E} > 0.0" then &
    "pair_coeff ${muc_e} ${muc_p} lj/cut 2 ${LJ_E} 1.0 2.5" &

if "${BMH_A} < 0.0" then &
    "pair_coeff ${muc_h} ${muc_h} born ${BMH_A} ${BMH_RHO} 1.0 1.0 1.0 2.5" &


################
# Bonded interactions (FENE)
################

bond_style fene
special_bonds fene
bond_coeff 1 ${FENE_K} ${FENE_R0} ${FENE_E} ${FENE_R}

################
# Angle interactions
################

angle_style harmonic
angle_coeff 1 ${BEND_K} 180.0

################
# Integration parameters
################

timestep        ${DELTA_T}
neighbor        1.0 bin
neigh_modify    every 1 delay 1 check yes

fix nve_fix     muc_group nve
fix nvt_fix     muc_group langevin ${Temperature} ${Temperature} ${T_DAMP} ${SEED}

################
# Pre-processing
################

compute msd_comp all msd
compute P all pressure thermo_temp
compute S all entropy/atom 0.25 2.5
compute STotal all reduce sum c_S
fix MSD all ave/time 100 1 100 c_msd_comp[4] file MSD.dat
fix Pvals all ave/time 100 1 100 c_P[4] c_P[5] c_P[6] file pressure-xy-xz-yz.dat
fix Svals all ave/time 100 1 100 c_STotal file entropy.dat
variable fe equal pe-(temp*c_STotal)

################
# Final setup for run
################

comm_style tiled
fix rcb_balance all balance 10000 1.1 rcb

#minimize 1.0e-4 1.0e-6 100 1000

info system communication

timer full

restart 50000000 restart
run ${N_STEPS}
write_data ${OUTPUT_PREFIX}.final.data nocoeff

################
# Post-processing
################

'''

        with open(self.cluster.jobname + ".Production.lammps_control", 'w') as stream:
            stream.write(eq_str)

        if self.verbose: print(f"MucusSeed::CreateProductionControl return")

    # Initialize a snapshot for the mucus simulations (based on membrane.py and the ah domain)
    def InitMucus(self, snap):
        if self.verbose: print(f"MucusSeed::InitMucus")

        print(f"  Creating mucus simulation from scratch!") 
        # Electrostatic ,cysteine, hydrophobic, positive(histonne)
        self.getTypebyName = {'muc_e': 0, 'muc_c': 1, 'muc_h': 2, 'muc_p': 3}
        self.getTypebyName['mucusbond'] = 0
        self.getTypebyName['mucusbend'] = 0

        # Start setting up the snapshot
        snap.particles.N        = self.natoms
        snap.particles.types    = ['muc_e', 'muc_c', 'muc_h', 'muc_p']
        snap.bonds.N            = self.nbonds
        snap.bonds.types        = ['mucusbond']
        snap.bonds.typeid[:]    = 0
        snap.angles.N           = self.nangles
        snap.angles.types       = ['mucusbend']
        snap.angles.typeid[:]   = 0

        # Create a freud box to do some fancy wrapping stuff
        import freud
        box = freud.Box.cube(self.lbox)

        # Now we can do the same thing as the LAMMPS config section
        # The -1 are for keeping track with the lammps code as well
        icount = -1
        ichaincount = -1
        dy1 = self.lbox/(2.0*self.nrowy)
        dz1 = self.lbox/(2.0*self.nrowz)
        for iz in range(self.nrowz):
            z = dz1*(iz-self.nrowz/2.0 + 0.5)
            for jy in range(self.nrowy):
                y = dy1*(jy-self.nrowy/2.0 + 0.5)
                x = random.uniform(0.0, self.lbox)
                itype = 3-1
                icount += 1
                ichaincount += 1

                # Every time we see a vector, wrap it using freud into a box
                r = [x, y, z]
                r_periodic = box.wrap(r)

                snap.particles.position[icount] = r_periodic
                snap.particles.typeid[icount] = itype
                snap.particles.mass[icount] = 1.0

                for k in range(1, self.nper_poly):
                    icount += 1
                    x = x + self.lbond
                    itype = 1-1
                    if ((k % self.nper_dimer) < self.n_term_length):
                        itype = 3 -1
                    if ((k % self.nper_dimer) > (self.nper_dimer - self.n_term_length - 1)):
                        itype = 3 -1
                    if (((k % self.nper_dimer) < self.monomer_length) and ((k % self.nper_dimer) >= (self.monomer_length - self.c_term_length))):
                        itype = 3 -1
                    for l in range(self.n_cysteine):
                        if (k % self.nper_dimer == self.cysteine_locations[l]): itype = 2-1
                        if ((self.nper_poly-k-1)%(self.nper_dimer) == self.cysteine_locations[l]): itype = 2-1

                    r = [x, y, z]
                    r_periodic = box.wrap(r)
                    snap.particles.position[icount] = r_periodic
                    snap.particles.typeid[icount] = itype
                    snap.particles.mass[icount] = 1.0

        beyondchaincount = ichaincount
        for ihist in range(self.nhist):
            icount += 1
            beyondchaincount += 1
            x = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            y = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            z = random.uniform(-self.lbox/2.0, self.lbox/2.0)
            itype = 4-1
            snap.particles.position[icount] = [x, y, z]
            snap.particles.typeid[icount] = itype
            snap.particles.mass[icount] = 1.0

        # Bond information
        icount1 = 1-1
        icount2 = 1-1
        for k in range(self.n_mucins):
            for n in range(self.nper_poly):
                if (n < self.nper_poly - 1):
                    snap.bonds.typeid[icount2] = self.getTypebyName['mucusbond']
                    snap.bonds.group[icount2] = [icount1, icount1+1]
                    icount2 += 1
                icount1 += 1

        # Angle (bend) information
        icount1 = 1-1
        icount2 = 1-1
        for k in range(self.n_mucins):
            for n in range(self.nper_poly):
                if (n < self.nper_poly - 2):
                    snap.angles.typeid[icount2] = self.getTypebyName['mucusbend']
                    snap.angles.group[icount2] = [icount1, icount1+1, icount1+2]
                    icount2 += 1
                icount1 += 1

        if self.verbose: print(f"MucusSeed::InitMucus return")

    def ReadMucus(self, snap):
        if self.verbose: print(f"MucusSeed::ReadMucus")

        print(f"  Reading mucus simulation from file {self.init_filename}!") 
        for itype,ntype in enumerate(snap.particles.types):
            self.getTypebyName[ntype] = itype

        if self.verbose: print(f"MucusSeed::ReadMucus return")
