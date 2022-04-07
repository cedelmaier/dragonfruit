# Here are some tips and tricks for GROMACS

# GROMACS running
Here are the tips and tricks for running GROMACS itself. These include the
ability to concatenate and process the trajectory files after a certain amount
of running gromacs.

## Running GROMACS
I have created run scripts for the membrane simulations that can be generated with
the create_run_scripts.py python script. This script is invoked in the following way,
for instance, if you have 100 steps that you want to simulate (100 ns, each 'step' is
1 ns worth of time). In this case, we start at step 1 after equilibration, then
run some number of scripts that generate 10ns worth of data (the stride value). The nsday
value is how efficient we think this is, to allocate the correct amount of time on the 
supercomputer.

    python3 create_run_scripts.py --start 1 --end 100 --stride 10 --ntmpi 3 --ntomp 8 --ngpu 1 --nsday 50

Then, start the first run script via something like this,

    sbatch run_prod_stride10_v1_10.sh

which will start the process of running all of the script you created. Eventually this will generate
100 output trajectory files from gromacs, each worth 1 ns of data. Then, you have to combine
them together in some meaningful way to see the full trajectory.

## Postprocessing GROMACS
We want to combine the trajectories together (usually) in a compressed data format. To do this, we
can use another of the packaged scripts here

    concatenate_trajedr_v1.sh 1 100

which will concatenate the first 1 to 100 files together into both a trajectory file (XTC) and an 
energy file (EDR). These files contain all of the information on the system, which includes the
water molecules. For many visualization purposes, this isn't ideal, as it includes lots of extra
information that isn't needed!

## CHARMM GUI membrane builder
Here is the recipe that I have for building the proper AH domain plus
bilayer online. You have to orient the proper AH dimensions for yourself
and make sure they're okay. I choose to align ILE4 and LEU11 along the 
z-axis.

    Align a vector (Two Atoms) along Z
    Choose: ILE4 and LEU11

Then you can choose to place the protein. If one is using the AH domain
that I built using Chimera, then you can just align the Y axis by 270 deg.

    Y rotate by 270 deg
    Z displace by <N> angstroms

Then build the system out of the relevant components.

To run equilibration, use the run_equilibration.sh script provided. This should
properly setup the environment and run equilibration steps 1-6

# VMD tips and tricks (groan)

## Simple display
If you want a really easy quick and dirty way to display the molecules in VMD, you can
use the following command to load a gromacs topology and trajectory into VMD from the
command line (OSX) like

    /Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/MacOS/startup.command step7_1.gro traj_continuous_v1_100.xtc

as this will load the topology and the frames from the continuous trajectory. Then, go into the graphical
representation of the molecules, and set one to protein (NewRibbons works well for me), then the following
two set to resname DOPC and resname PLPI, then set them to VDW radii or something similar. Set them to transparent, and voila,
you have something to show to people as a movie!

## Center a protein in the simulation box after simulation

