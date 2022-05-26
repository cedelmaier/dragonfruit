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

## Additional programs/feelings about GROMACS
Using just gromacs doesn't get you everything. For instance, you probably want to install the xmgrace
tool so that you can visualize the XVG files that it sometimes spits out (ugh). To do this, the one
packaged with brew seems okay.

    brew install grace

which you can then use for things like viewing these stupid XVG files

    xmgrace scount.xvg

To view the xpm files, you need to run them through something to convert to EPS format, like

    gmx xpm2ps -f dssp_test.xpm -o dssp_test.eps

which will give an eps file that can be converted into a plot for what kind of residues are turns/helices,
for instance.

## Umbrella sampling
**WARNING**: This is setup for pulling of COM between the AH domain and the membrane, without a reference
coordinate other than the COM of the lipid bilayer and the AH-domain. This may or may not be correct. So, 
make sure to think about what you are submitting!

Umbrella sampling is done similar to what is done in the gromacs tutorial 3 online. This is based off of
previous work that has also done umbrella sampling (see papers on Overleaf for more details). The basics of
this are to run a simulation and then pull the AH domain out of the membrane, measuring it at different intervals.
To do this, we need to run several commands in order, some of which do require GROMACS to run, and may not work
dependong on the version synchornizations between computers.

First, one must start from a simulation of the AH-domain and membrane. I've used the coiled AH inserted at 
+00 angstroms for this example. First, generate a configuration with the CHARMM-GUI, but include a *lot* of
water on both sides (NOTE: This would be better if it were one-sided), in order to give enough room to pull
the AH domain out of the membrane. Equilibration is done the same as the other types of simulations.

    run_equilibration.sh

This will generate the configuration to use going forward with the steered MD. This steered MD attempts to rip
the AH domain out of the membrane at some rate using a potential to restrict where the AH domain wants to be. The
equilibrated simulation can be fed into the steered MD via

    vim step7_umbrella_v1.mdp
    run_steered_md_v1.sh

to edit and then configure the exact pulling parameters for the simulation. Once this is done, we need to do 
some post-processing in order to get out the distances that we are going to use for our WHAM routine. This
can be done with the script get_distances.sh, which generates a distance file and configuration for EVERY
single frame in the previous trajectory.

    get_distances.sh

Once you have the summary distance file, you can generate the umbrella sampling MD scripts and submit them 
en mass to the cluster. This can be done as they are all independent simulations, and so don't need to run
one after the other. For example, for 0.1 nm spacing you can run

    vim npt_umbrella_v1.mdp
    vim md_umbrella_v1.mdp
    python3 setup_umbrella.py -dfile summary_distances.dat -i 0.1 --ntmpi 3 --ntomp 8 --ngpu 3 --nsday 50.0

this will also check the two setup files that you'll need for the umbrella sampling (always check these to make sure
that everything is okay!). Then you can just submit the final script and wait for the answer.

    ./submit_all.sh


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

# PyMOL tip and tricks (also groan)
Well, if you don't want to use VMD, this is another option. Not necessarily a better one, but it is different.

## Centering and ripping a protein out of a trajectory for visualization
PyMOL has some eccentricities, so let's deal with those first. First step is to recenter the protein in the box, and
then we can dump the information without the waters included, which significantly reduces the size. One probably wants 
to create a new index map that just has the lipid and protein members, which
you need to generate a new index file for (in this case Protein, PLPI, DOPC lipids)
    
    gmx make_ndx -f step7_1.gro -o extra_groups.ndx
    1 | 13 | 14
    13 | 14
    q

    gmx trjconv -s step7_1.gro -f traj_continuous_v1_200.xtc -o traj_nopbc.xtc -boxcenter rect -pbc mol -center -n extra_groups.ndx

    Select 1 for the centering, then 18 (or whatever) for what to actually write out

Either way, then you need the first frame of the trajectory as a topology file.

    gmx trjconv -s step7_1.gro -f traj_continuous_v1_200.xtc -o cfg_0.gro -boxcenter rect -pbc mol -center -dump 0 -n extra_groups.ndx 

    Select 1 for the centering, then 18 (or whatever) for what we actually write out

So really, just run these to get the trajectories and the initial setup

    gmx trjconv -s step7_1.gro -f traj_continuous_v1_200.xtc -o traj_nopbc.xtc -boxcenter rect -pbc mol -center -n extra_groups.ndx
    gmx trjconv -s step7_1.gro -f traj_continuous_v1_200.xtc -o xyz.pdb -boxcenter rect -pbc mol -center -n extra_groups.ndx -dump 0

Now you should have two files, and then we can rename them so that PyMOL doesn't throw a hissy fit when we load them (hopefully). Note, these
have to have the exact same name

    mv traj_nopbc.xtc ahpymol.xtc
    mv cfg_0.gro ahpymol.gro
    mv xyz.pdb ahpymol.pdb

20220520 Update: Just use the following. Center the box on the lipid layer, so that we can see the protein moving, and then run the
analysis. You still have to make the new groups as before (see above).

    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.pdb -n extra_groups.ndx -boxcenter rect -pbc mol -center -dump 0
    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.xtc -n extra_groups.ndx -boxcenter rect -pbc mol -center

## Some sort of workflow for GROMACS simulations
This is done using version 1 of the umbrella sampling that I did. This is to visualize the entire sequence in some sort of rational form.
First, preprocess this shit to make it easier in PyMOL.

    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.pdb -n extra_groups.ndx -boxcenter rect -pbc mol -center -dump 0
    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.xtc -n extra_groups.ndx -boxcenter rect -pbc mol -center

In pyMol, found these to be somewhat effective commands

    cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/coiled_umbrella/v1
    load step7_umbrella_v1_reduced.pdb
    load_traj step7_umbrella_v1_reduced.xtc
    smooth all, 30, 3
    select DOPC, resname DOPC
    select PLPI, resname PLPI
    select lipids, resname DOPC or resname PLPI
    select protein, chain A
    set sphere_transparency=0.9, lipids

    hide all
    show spheres, lipids
    show cartoon, protein
    select dopc_headgroups, resname DOPC and (name N or name C12 or name C13 or name C14 or name C15 or name C11 or name P or name O13 or name O12 or name O11)
    select dopc_notheadgroups, (resname DOPC) and (not name N and not name C12 and not name C13 and not name C14 and not name C15 and not name P and not name O13 and not name O12 and not name O11)
    

# Analysis

## CD measurements
We can make predicted CD measurements of our structures with the online tool pdbmd2cd. This allows us to take a bunch of pdb files, submit
then to the server, and then have an analysis of the helix nature of the simulated samples. One can generate a list of the PDB files
via the following commands from GROMACS.

    gmx trjconv -s step7_1.tpr -f traj_continuous_v1_600.xtc -o xyz.pdb -boxcenter rect -pbc mol -center -dt 10000 -sep
    tar -czvf xyz.tar.gz xyz*.pdb

This allows multiple structure files to be submitted at once to the online tool. Notice that there is a separation and a dt in ters of (ps)
that is used to control how many structures we submit to the server. To dump the last frame, you can also do something like the following.

    gmx trjconv -s step7_1.tpr -f traj_continuous_v1_600.xtc -o frame_6000.pdb -boxcenter rect -pbc mol -center -dump 600000

## DSSP
This is a way to assign helix-like states to the different things in PDB files and structures. It's also super, super finnicky, so this
is gonna suck. In order to install a version that works on OSX, I found that you can install it via anaconda in it's own environment

    conda create --name dssp_osx
    conda activate dssp_osx
    conda install -c salilab dssp

This will create the damn thing, as the online versions (hssp, xssp, whatever else) don't seem to work properly. You can run this in standalone mode
on a single PDB file with somethig like

    mkdssp -i final_frame.pdb

And this will generate the records. To run through gromacs, you need to run something like

    gmx do_dssp -s step7_1.tpr -f traj_continuous_v1_1000.xtc -n index.ndx -dt 100000 -ssdump -o unfolded_zdepth30_s1_dssp.xpm -a -aa -ta

One can then turn this into an image (or examine the different stuff in the files), you need to run something like

    gmx xpm2ps -f dssp_test.xpm -o dssp_test.eps