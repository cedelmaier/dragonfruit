# Here are some tips and tricks for GROMACS

# Running GROMACS
Here are the tips and tricks for running GROMACS itself. These include the
ability to concatenate and process the trajectory files after a certain amount
of running gromacs.

## UNC: Longleaf: Loading GROMACS
Longleaf has problems with the modules correctly working with gromacs. To get around this, load
the modules needed, and then directly source the correct setup script in the gromacs directory
that you want to use.

For CUDA-based gromacs, use the following commands.
	
	module load gcc/9.1.0
	module load cuda/11.4
	source /nas/longleaf/apps/gromacs/2021.5/avx_512-cuda11.4/bin/GMXRC
	
For CUDA-based-PLUMED gromacs, use the following commands.

	module load gcc/9.1.0
	module load cuda/9.2
	source /nas/longleaf/apps/gromacs/???


## Flatiron Institute: Rome: Loading GROMACS
For the Flatiron computers, to load gromacs with the plumed compilation, this seems
to work the best.

	module load modules/2.0-20220630
	module load openmpi/4.0.7
	module load gromacs/skylake-mpi-plumed-2021.4
	module load plumed/mpi-2.8.0
	
Make sure that we have a conda environment loaded with the proper MDAnalysis/MDTraj dependencies. At the moment, this doesn't work. Figure out why another time.
	
At the moment, have a special create_runscripts_gromacsplumed.py that generates the plumed input data file for measuring the CVs of the system. This requires a separate plumed file for each run, so keep track of this here! Try to merge the create_run_scripts.py files in the future to work for UNC and FI, gromacs, gromacs+plumed, and whatever else I can come up with.


## Running GROMACS
I have created run scripts for the membrane simulations that can be generated with
the create\_run\_scripts.py python script. This script is invoked in the following way,
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

Note, some of these do not care about the ns/day (like at Flatiron), because they have a very generous allocation. Still, it should be used if you have some idea about how long things run.

## Preprocessing GROMACS
In order to use the structure files from Ronit's lab, or really just as a sanity check, you need to make sure that you're properly patching the end termini of the protein under study. For instance, from Ronit's lab, we want to add NH2 onto both the Nterminal and Cterminal of the AH domain. To do this, several steps must occur.

### Freeman lab instructions

First, remove the hydrogens in the configuration file, as they will get in the way.

	vmd
	mol load pdb your.pdb
	set sel [atomselect top "not hydrogen"]
	$sel writepdb monomer_preprocessed.pdb
	
Then, go and edit this new file (whatever you named it) to remove the trailing nitrogen from the PDB file, like LEU18 has an N at N01 near the end. The peptide is now ready for the terminal patching as done by GROMACS. Temrinal patching is done in GROMACS.

	gmx pdb2gmx -f monomer_preprocessed.pdb -o monomer_processed.gro -water tip3p -ignh -ter
	
Where you interactively choose the correct options. This will genereate a gromacs file for use. To use the same in the CHARMM gui, you need to do use the preprocessed monomer with the correct termianl patching, as done in the CHARMM gui. This corresponds to options NNEU and CT2 when importing your protein. In the end, you should wind up with PDB files

	monomer_preprocessed.pdb
	monomer_charmmgui.pdb
	
and gromacs file

	monomer_processed.gro
	
Really, overall, the whole point is to properly get the terminal patching done with either GROMACS or the CHARMM-gui, and this can then be used as inputs to the run system. Hopefully, overall, you use the CHARMM-gui to get things configured, and then can just use my canned run scripts to submit everything to whatever supercomputer you are using.

## Postprocessing GROMACS
We want to combine the trajectories together (usually) in a compressed data format. To do this, we
can use another of the packaged scripts here

    concatenate_trajedr_v1.sh 1 100

which will concatenate the first 1 to 100 files together into both a trajectory file (XTC) and an 
energy file (EDR). These files contain all of the information on the system, which includes the
water molecules. For many visualization purposes, this isn't ideal, as it includes lots of extra
information that isn't needed!

## GROMACS performance
Here is an easy awk script ot get the performance of a bunch of gromacs simulations in the same directory.

	grep -ir '^Performance' step7_*.log | awk '{sum+=$2; n++}; END {printf "%.3f\n", sum/n}'

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

Then build the system out of the relevant components. To generate an asymmetric (or at least big enough) water box, use something like 50 angstroms of space above and below the bilayer.

	50.0 angstroms water

To run equilibration, use the run_equilibration.sh script provided. This should
properly setup the environment and run equilibration steps 1-6,

A dual AH domain simulation is prepared by aligning the ILE4 residues so that the proteins are oriented the proper way, and displacing the second AH-domain by some number of angstroms in pymol using the translate command. There is probably a better way to do this, by creating the peptides in PyMOL to start, copying them, and then generating the water/etc around the peptides in the membrane.

## CHARMM GUI preparing asymmetric water box
In the case of PMF calculations we need to prepare an asymmetric water box, in order to have the
correct amount of space in the direction we are pulling the protein. In this case, we can
start with the same process as the membrane builder above. However, at the step where the amount
of water on top/below the bilayer is chosen, choose something large. The bilayer will still be
centered, but we will fix that in a moment.

Try this type of thing first. This should use the native GROMACS magic to translate the
entire system, and then create the input files.

    gmx editconf -f step5_input.gro -o step5_pbcshift.gro -translate 0 0 -3
    gmx trjconv -f step5_pbcshift.gro -o step5_pbcshift_res.gro -pbc res -s step5_input.gro

However, it has issues with splitting up residues. Maybe use pymol instead?

## Additional programs/feelings about GROMACS
Using just gromacs doesn't get you everything. For instance, you probably want to install the xmgrace
tool so that you can visualize the XVG files that it sometimes spits out (ugh). To do this, the one
packaged with brew seems okay.

    brew install grace

which you can then use for things like viewing these stupid XVG files

    xmgrace scount.xvg

To view the xpm files, you need to run them through something to convert to EPS format, like

    gmx xpm2ps -f dssp_test.xpm -o dssp_test.eps

which will give an eps file that can be converted into a plot for what kind of residues are turns/helices, for instance.

## Umbrella sampling via GROMACS
**WARNING**: This is setup for pulling of COM between the AH domain and the membrane, without a reference coordinate other than the COM of the lipid bilayer and the AH-domain. This may or may not be correct. So, make sure to think about what you are submitting!

Umbrella sampling is done similar to what is done in the gromacs tutorial 3 online. This is based off of previous work that has also done umbrella sampling (see papers on Overleaf for more details). The basics of this are to run a simulation and then pull the AH domain out of the membrane, measuring it at different intervals. To do this, we need to run several commands in order, some of which do require GROMACS to run, and may not work dependong on the version synchornizations between computers.

First, one must start from a simulation of the AH-domain and membrane. I've used the coiled AH inserted at +00 angstroms for this example. First, generate a configuration with the CHARMM-GUI, but include a *lot* of water on both sides (NOTE: This would be better if it were one-sided), in order to give enough room to pull the AH domain out of the membrane. Equilibration is done the same as the other types of simulations.

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

Then, after they have all run, you can combine them together into data files for the gmx WHAM routine to 
combine them together. Note, this requires some maniuplation of the files to get them in the correct
order, but that's fine.

    ls -1 umbrella*.tpr >> tpr_files.dat
    ls -1 umbrella*_pullf* >> pull_files.dat
    gmx wham -it tpr_files.dat -if pull_files.dat -o -hist -unit kCal

This generates some files that you can then look into for the binding energy!

## Umbrella Sampling via PLUMED
Plumed and GROMACS work slightly differently. You can look at the online tutorial for plumed (any of the master-classes should work) to understand how this is working. Currently I'm running a hybrid design. First, the peptide is steered either outside to inside (or the reverse) using PLUMED, and then the actual sampling is done with GROMACS directly. This is because the GROMACS system seems to be slightly easier to interpret the 1D umbrella sampling and can be done direclty with gromacs, rather than using PLUMED to generate the final free energy surface files.

## Metadynamics
We can also use metadynamics to figure out the free energy of the peptide. In this case, I started from an 'equilibrated' state where the peptide is at the surface of the membrane. Then, this can be run through PLUMED in a metadynamics simulation.

Currently, my plumed file looks like this for a neutral-cap monomer.

	# vim:ft=plumed
	MOLINFO STRUCTURE=reference.pdb

	# Figure out the lipid and helix COM coordinates
	lipid_com: COM ATOMS=313-48624
	helix_com: COM ATOMS=1-312

	# Get the Z distance
	z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS

	# Get the alpha value
	alpha: ALPHARMSD RESIDUES=1-18
	
	# Get some end to end alpha helix distance
	alpha_dist: DISTANCE ATOMS=4,308

	# Set up walls to prevent the helix from escaping
	uwall: UPPER_WALLS ARG=z_dist.z AT=5.0  KAPPA=200.0 EXP=2 EPS=1 OFFSET=0
	lwall: LOWER_WALLS ARG=z_dist.z AT=-5.0 KAPPA=200.0 EXP=2 EPS=1 OFFSET=0

	# Set up the simple 1D metadynamics
	metad: METAD ARG=z_dist.z ...
    	# Deposit a gaussian every 500 time steps, with an initial height
    	PACE=500 HEIGHT=0.5 BIASFACTOR=15
    	# Try a gaussian width of 0.2 nm
    	SIGMA=0.2
    	# Gaussian written to file and stored on grid
    	FILE=HILLS GRID_MIN=-5.0 GRID_MAX=5.0
	...

	# Print to a file
	PRINT ARG=* FILE=colvar_metad.dat STRIDE=500
	
### Getting the masses and charges out of GROMACS
We usually need the masses and charges out of a gromacs simulation for plumed to use the XTC files. To do this, use the following file and commands.

`plumed_dumpmasscharge.dat`:

	DUMPMASSCHARGE FILE=mcfile.dat

Then you can run the command.

	mpirun -np 1 gmx_mpi mdrun -s step6.0_minimization.tpr -nsteps 1 -plumed plumed_dumpmasscharge.dat

and you will have a file `mcfile.dat`. Note this has to be done with a version of gromacs that is compiled with PLUMED (yuck).

### Harvesting plumed info from previous simulations
After you get out the masses/charges from a previously run simulation, you can then runt he plumed driver on it. If you use the following file, then you can harvest things like the Z distance and the alpha variable.

`plumed_metad_measure.dat`:

	# vim:ft=plumed
	MOLINFO STRUCTURE=reference.pdb

	# Figure out the lipid and helix COM coordinates
	lipid_com: COM ATOMS=313-48624
	helix_com: COM ATOMS=1-312

	# Get the Z distance
	z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS

	# Get the alpha value
	alpha: ALPHARMSD RESIDUES=1-18

	# Get some end to end alpha helix distance
	alpha_dist: DISTANCE ATOMS=4,308

	# Print to a file
	PRINT ARG=* FILE=colvar_metad.dat STRIDE=1

You can then run the plumed driver with the file to harvest the data.

	plumed driver --mf_xtc traj_continuous_v1_100.xtc --mc mcfile.dat --plumed plumed_metad_measure.dat --kt 2.494339
	
### Reweighting the original simulation
We can also re-bias the original simulation to see if the free energies are the same (converged) or not (hint: they're not, as we don't have nearly the statistics we need!). First, you are going to need the results of the metadynamics simulation in a HILLS file (see above). Then, you can use this to "unbias" the original simulation. First, you need a new plumed file that doesn't deposit new gaussian kernels and can then be used to get out the free energy contributions. NOTE: This is done on the 'original' simulation, not the metadynamics version.

`plumed_metad_reweight.dat`:

	# vim:ft=plumed
	MOLINFO STRUCTURE=reference.pdb

	# Figure out the lipid and helix COM coordinates
	lipid_com: COM ATOMS=313-48624
	helix_com: COM ATOMS=1-312

	# Get the Z distance
	z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS

	# Get the alpha value
	alpha: ALPHARMSD RESIDUES=1-18

	# Get some end to end alpha helix distance
	alpha_dist: DISTANCE ATOMS=4,308

	# Set up walls to prevent the helix from escaping
	#uwall: UPPER_WALLS ARG=z_dist.z AT=5.0  KAPPA=200.0 EXP=2 EPS=1 OFFSET=0
	#lwall: LOWER_WALLS ARG=z_dist.z AT=-5.0 KAPPA=200.0 EXP=2 EPS=1 OFFSET=0

	# Set up the simple 1D metadynamics, increase time and don't deposit hills
	metad: METAD ARG=z_dist.z ...
    	# Deposit a gaussian (NEVER) with zero height
    	PACE=1000000000 HEIGHT=0.0 BIASFACTOR=15
    	# Try a gaussian width of 0.2 nm
    	SIGMA=0.2
    	# Gaussian written to file and stored on grid
    	FILE=/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/enhanced_sampling/rfmonomer_aglipid_11x11_zdepth00_50mMKCl_metad_v2/HILLS GRID_MIN=-5.0 GRID_MAX=5.0
    	# Make sure we restart this
    	RESTART=YES
	...

	# Use the metadynamics bias as argument
	as: REWEIGHT_BIAS ARG=metad.bias

	# Calculate histograms of dz and alpha every 50 steps using the weights from the bias potential
	hhdz: HISTOGRAM ARG=z_dist.z STRIDE=50 GRID_MIN=-5.0 GRID_MAX=5.0 GRID_BIN=50 BANDWIDTH=0.05 	LOGWEIGHTS=as
	hhalpha: HISTOGRAM ARG=alpha STRIDE=50 GRID_MIN=0.0 GRID_MAX=18.0 GRID_BIN=50 BANDWIDTH=0.05 	LOGWEIGHTS=as
	#Convert histograms h(s) to free energies
	ffdz: CONVERT_TO_FES GRID=hhdz
	ffalpha: CONVERT_TO_FES GRID=hhalpha
	#Print to the file
	DUMPGRID GRID=ffdz FILE=ffdz.dat
	DUMPGRID GRID=ffalpha FILE=ffalpha.dat

	# Print to a file
	PRINT ARG=z_dist.z,alpha,metad.bias FILE=colvar_reweight.dat STRIDE=1

One can then see how hilariously off the free energy estimates are by using the plumed driver.

	plumed driver --mf_xtc traj_continuous_v1_100.xtc --mc mcfile.dat --plumed plumed_metad_reweight.dat --kt 2.494339


### Selecting a subset of atoms for the CVs
Sometimes we might want to select some subset of atoms for our CVs, as they are quite large and computationally intensive. Here are some tricks, then.

Phosphorus only from lipids:

	grep -Ew "P" reference.pdb.bak | awk '{printf "%s,",$2}'
    
Carbon backbone only from helix:

	grep -Ew "CA" reference.pdb.bak | awk '{printf "%s,",$2}'
	
# Metadynamics Version 2
This is an updated version of the metadynamics information that I've found. So this is a work in progress.

## FES via sum\_hills
We can use the native plumed sum_hills utility to generate a free energy landscape that should be the corrected version (well-tempered or standard metadynamics).

Use the following command, for instance, to generate a 2D FES of the landscape from the hills file.

	mpirun -np 1 plumed sum_hills --bin 255,255 --min 0,0 --max 5.0,13.1 --hills HILLS_SINGLEWALKER --mintozero --outfile fes_sumhills_original.dat
	
One can also generate the negative bias alone.

	mpirun -np 1 plumed sum_hills --bin 255,255 --min 0,0 --max 5.0,13.1 --hills HILLS_SINGLEWALKER --mintozero --negbias --outfile fes_sumhills_negbias.dat
	
## Unbias a multiwalker simulation
If we want access CVs that are not part of the original biasing routine with metadynamics, we need to come up with a way to take the biased trajectories and unbias them. This is complicated further when we use multiple walker metadynamics, as each walker's trajectory is biased by a potential coming from all of the walkers. Here are the steps that I was able to put together in order to unbias a trajectory, giving a result that is similar between getting the FES from the HILLS file(s), and one that goes through the unbiasing routine.

The following is done assuming there are 4 walkers (or really however many) that each read from a common HILLS file defining the metadynamics potential as it is written.

### Step 1: Reduce single walker trajectory file contents and concatenate trajectories
First one needs to combine all of the XTC (or TRR) trajectory files from GROMACS into a single entity. It is very advisable to reduce the number of atoms that are recorded at this stage. If you're using a `density_groups.ndx` file (similar to this document's version), then the internals of the following script will center the box on the lipid group and write out the lipid+helix group. This has to be done **independently** for each walker (very important!).

`reduce_concatenate_metadblocks.sh`

    #!/bin/bash
    
    startstep=$1
    endstep=$2

    # Remove previous control groups (for density)
    rm control_groups.txt

    echo "18" >> control_groups.txt
    echo "19" >> control_groups.txt

    # First reduce the content of all of the files
    for ((i=$startstep;i<=$endstep;i++))
    do
      mpirun -np 1 gmx_mpi trjconv -s metad_2d_block${i}.tpr -f metad_2d_block${i}.xtc -o metad_2d_block${i}reduced.pdb -n density_groups.ndx -boxcenter rect -pbc mol -center -dump 0 < control_groups.txt
      mpirun -np 1 gmx_mpi trjconv -s metad_2d_block${i}.tpr -f metad_2d_block${i}.xtc -o metad_2d_block${i}reduced.xtc -n density_groups.ndx -boxcenter rect -pbc mol -center < control_groups.txt
    done

    # Now concatenate into a singular file
    rm -f list_trr_oneline.txt
    rm -f times_trr.txt
    rm -f list_edr_oneline.txt
    rm -f times_edr.txt

    for ((i=$startstep;i<=$endstep;i++))
    do
      echo -n " metad_2d_block${i}reduced.xtc" >> list_trr_oneline.txt;
      echo "c" >> times_trr.txt;
    done

    # Now concatenate the files automatically for the user
    mpirun -np 1 gmx_mpi trjcat -f `cat list_trr_oneline.txt` -o traj_continuous_v$startstep\_$endstep\_reduced.xtc -settime < times_trr.txt

This will result in a file (for blocks 1 to 25) of `traj_continuous_v1_25_reduced.xtc`, as well as a number of PDB files that are useful for visualizing the results.

### Step 2: Combine trajectories from all walkers together
Now we need to create a single trajectory from the different walkers to properly read into the plumed driver.

    mpirun -np 1 gmx_mpi trjcat -cat -f walker?/traj_continuous_v1_25_reduced.xtc -o traj_multi_reduced.xtc
    
### Step 3: Create a consistent reference and mass/charge files for PLUMED
Since we have stripped out atoms/whatever from our trajectory, we need to make sure that both the reference and mass/charge files are consistent. The definition of the reference file is probably fine, assuming that one has stripped out the water and salt ions from the original PDF file, and have a reference file `reference.pdb`. Briefly, here is what you do to generate this file (and keep a useful backup for things like finding atoms in all cases later).

`reference.pdb`

    mpirun -np 1 gmx_mpi editconf -f metad_2d_block1.tpr -o reference.pdb
    cp reference.pdb reference.pdb.bak

Then you have to edit the file and remove the salt/water molecules. To get the masses and charges out of the original file for PLUMED, you can use the following commands.

`plumed_dumpmasscharge.dat`:

	DUMPMASSCHARGE FILE=mcfile.dat

Then you can run the command.

	mpirun -np 1 gmx_mpi mdrun -s metad_2d_block1.tpr -nsteps 1 -plumed plumed_dumpmasscharge.dat

and you will have a file `mcfile.dat`. Note this has to be done with a version of gromacs that is compiled with PLUMED. After this is done, you need to go back through and remove the extra molecules/atoms that were stripped out for the updated reference file as well.

Really, in the end, you just need to have consistent `reference.pdb` and `mcfile.dat` files.

### Step 4: Determine what collective variables you want to analyze
I recommend including the bias variables in this, as then you have a way to test the FES against the version one can get from the HILLS file. Here is an example file that generates FES grids based on the CVs of dz, alpha, and a new 'tilt' angle.  As a note, this does the reweighting via both the exponential and Tiwary weights.

`plumed_multiwalker_reweight.dat`

    # vim:ft=plumed
    MOLINFO STRUCTURE=reference.pdb

    # Get detailed timing
    DEBUG DETAILED_TIMERS

    # Figure out the lipid and helix COM coordinates
    # lipid: center of mass of phosphorus only
    lipid_com: COM ATOMS=332,470,608,746,884,1022,1160,1298,1436,1574,1712,1850,1988,2126,2264,2402,2540,2678,2816,2954,3092,3230,3368,3506,3644,3782,3920,4058,4196,4334,4472,4610,4748,4886,5024,5162,5300,5438,5576,5714,5852,5990,6128,6266,6404,6542,6680,6818,6956,7094,7232,7370,7508,7646,7784,7922,8060,8198,8336,8474,8612,8750,8888,9026,9164,9302,9440,9578,9716,9854,9992,10130,10268,10406,10544,10682,10820,10958,11096,11234,11372,11510,11648,11786,11924,12062,12200,12338,12476,12614,12752,12890,13028,13166,13304,13442,13580,13718,13856,13994,14132,14270,14408,14546,14684,14822,14960,15098,15236,15374,15512,15650,15788,15926,16064,16202,16340,16478,16616,16754,16892,17030,17168,17306,17444,17582,17720,17858,17996,18134,18272,18410,18548,18686,18824,18962,19100,19238,19376,19514,19652,19790,19928,20066,20204,20342,20480,20618,20756,20894,21032,21170,21308,21446,21584,21722,21860,21998,22136,22274,22412,22550,22688,22826,22964,23102,23240,23378,23516,23654,23792,23930,24068,24206,24344,24482,24620,24758,24896,25034,25172,25310,25448,25586,25724,25862,26000,26138,26276,26414,26552,26690,26828,26966,27104,27242,27380,27518,27656,27794,27932,28070,28208,28346,28484,28622,28760,28898,29036,29174,29312,29450,29588,29726,29864,30002,30140,30278,30416,30554,30692,30830,30968,31106,31244,31382,31520,31658,31796,31934,32072,32210,32348,32486,32624,32762,32900,33038,33176,33314,33452,33590,33728,33866,34004,34142,34280,34418,34556,34694,34832,34970,35108,35246,35384,35522,35660,35798,35936,36074,36212,36350,36488,36626,36767,36902,37037,37172,37307,37442,37577,37712,37847,37982,38117,38252,38387,38522,38657,38792,38927,39062,39197,39332,39467,39602,39737,39872,40007,40142,40277,40412,40547,40682,40817,40952,41087,41222,41357,41492,41627,41762,41897,42032,42167,42302,42437,42572,42707,42842,42977,43112,43247,43382,43517,43652,43787,43922,44057,44192,44327,44462,44597,44732,44867,45002,45137,45272,45407,45542,45677,45812,45947,46082,46217,46352,46487,46622,46757,46892,47027,47162,47297,47432,47567,47702,47837,47972,48107,48242,48377,48512
    helix_com: COM ATOMS=4,23,38,53,72,89,96,118,134,156,178,197,212,227,244,260,282,293

    # Get the Z distance
    z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS NOPBC
    dz: COMBINE ARG=z_dist.z PERIODIC=NO

    # Get the alpha value
    alpha: ALPHARMSD RESIDUES=1-18

    # Construct a tilt via proxy between CA of ILE4 and LEU11
    rvec: DISTANCE ATOMS=53,178 COMPONENTS NOPBC
    tilt: CUSTOM ...
        ARG=rvec.x,rvec.y,rvec.z
        VAR=rx,ry,rz
        FUNC=acos((rz)/sqrt(rx*rx+ry*ry+rz*rz))
        PERIODIC=NO
    ...

    # Set up parallel replica metadynamics
    metad: METAD ... 
        ARG=dz,alpha
        SIGMA=0.05,0.1
        HEIGHT=0.0
        BIASFACTOR=15
        PACE=10000000
        GRID_MIN=-7.0,0.0
        GRID_MAX=7.0,18.0
        FILE=HILLS_MULTIWALKER.block25
        RESTART=YES
        WALKERS_MPI
    ...

    # Try the different reweighting schemes (including none). In theory, the metadynamcis reweight maps onto the Tiwary reweights
    as: REWEIGHT_BIAS   ARG=metad.bias
    tw: REWEIGHT_METAD  ARG=metad.bias

    # Biased histograms just to see what is going on
    hh_z_alpha_biased:      HISTOGRAM ARG=dz,alpha      STRIDE=1 GRID_MIN=0.0,0.0  GRID_MAX=5.0,13.1 GRID_BIN=255,255 BANDWIDTH=0.05,0.1
    hh_z_tilt_biased:       HISTOGRAM ARG=dz,tilt       STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=5.0,pi   GRID_BIN=255,255 BANDWIDTH=0.05,0.1
    hh_alpha_tilt_biased:   HISTOGRAM ARG=alpha,tilt    STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=13.1,pi  GRID_BIN=255,255 BANDWIDTH=0.05,0.1

    # Unbiasing via reweight_bias (maybe the exponential bias talked about in the reweighting paper?)
    hh_z_alpha_as:          HISTOGRAM ARG=dz,alpha      STRIDE=1 GRID_MIN=0.0,0.0  GRID_MAX=5.0,13.1 GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=as
    hh_z_tilt_as:           HISTOGRAM ARG=dz,tilt       STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=5.0,pi   GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=as
    hh_alpha_tilt_as:       HISTOGRAM ARG=alpha,tilt    STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=13.1,pi  GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=as

    # Unbiasing via reweight_bias (maybe the exponential bias talked about in the reweighting paper?)
    hh_z_alpha_tw:          HISTOGRAM ARG=dz,alpha      STRIDE=1 GRID_MIN=0.0,0.0  GRID_MAX=5.0,13.1 GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=tw
    hh_z_tilt_tw:           HISTOGRAM ARG=dz,tilt       STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=5.0,pi   GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=tw
    hh_alpha_tilt_tw:       HISTOGRAM ARG=alpha,tilt    STRIDE=1 GRID_MIN=0.0,-pi  GRID_MAX=13.1,pi  GRID_BIN=255,255 BANDWIDTH=0.05,0.1 LOGWEIGHTS=tw

    # Convert to FES
    ff_z_alpha_biased:      CONVERT_TO_FES GRID=hh_z_alpha_biased
    ff_z_tilt_biased:       CONVERT_TO_FES GRID=hh_z_tilt_biased
    ff_alpha_tilt_biased:   CONVERT_TO_FES GRID=hh_alpha_tilt_biased

    ff_z_alpha_as:          CONVERT_TO_FES GRID=hh_z_alpha_as
    ff_z_tilt_as:           CONVERT_TO_FES GRID=hh_z_tilt_as
    ff_alpha_tilt_as:       CONVERT_TO_FES GRID=hh_alpha_tilt_as

    ff_z_alpha_tw:          CONVERT_TO_FES GRID=hh_z_alpha_tw
    ff_z_tilt_tw:           CONVERT_TO_FES GRID=hh_z_tilt_tw
    ff_alpha_tilt_tw:       CONVERT_TO_FES GRID=hh_alpha_tilt_tw

    # Dump all of the grids
    DUMPGRID GRID=ff_z_alpha_biased     FILE=ff_z_alpha_biased.dat
    DUMPGRID GRID=ff_z_tilt_biased      FILE=ff_z_tilt_biased.dat
    DUMPGRID GRID=ff_alpha_tilt_biased  FILE=ff_alpha_tilt_biased.dat

    DUMPGRID GRID=ff_z_alpha_as         FILE=ff_z_alpha_as.dat
    DUMPGRID GRID=ff_z_tilt_as          FILE=ff_z_tilt_as.dat
    DUMPGRID GRID=ff_alpha_tilt_as      FILE=ff_alpha_tilt_as.dat

    DUMPGRID GRID=ff_z_alpha_tw         FILE=ff_z_alpha_tw.dat
    DUMPGRID GRID=ff_z_tilt_tw          FILE=ff_z_tilt_tw.dat
    DUMPGRID GRID=ff_alpha_tilt_tw      FILE=ff_alpha_tilt_tw.dat

With this file, you can then run the PLUMED driver on the combined trajectory file (reduced) and get several free energy grid outputs.

### Step 5: Unbias the simulation
This command will run the plumed driver on the combined multiwalker trajectories to unbias the simulation with the collective variables defined in the `plumed_multiwalker_reweight.dat` file.

    mpirun -np 1 plumed driver --ixtc traj_multi_reduced.xtc --plumed plumed_multiwalker_reweight.dat --kt 2.494339 --mc mcfile.dat
    
### Step 6: Santy check on free energies
This command will generate several `ff*.dat` files. For our specific example, the `ff_z_alpha_biased.dat` file contains the *biased* version of the free energy landscape for our collective variables. Metadynamics should be flattening this out so that all states are equally probable, and so this is a good sanity check, along with comparing the unbiased landscape will will have in a moment. One can look at this in a jupyter notebook using the following (for example). **NOTE**: The following was taken directly from a jupyter notebook, and so you need to make sure you break this up to actually use it effectively, as a warning.

    import subprocess
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import plumed
    import os
    
    data = plumed.read_as_pandas("unbias_attempt1/ff_z_alpha_biased.dat")
    npoints = 256
    dz    = np.array(data["dz"]).reshape(npoints,npoints)
    alpha = np.array(data["alpha"]).reshape(npoints,npoints)
    fes   = np.array(data["ff_z_alpha_biased"]).reshape(npoints,npoints)

    fes_adj = fes - np.min(fes)

    print(f"Biased max after subtraction: {np.nanmax(fes_adj[fes_adj != np.inf])}")
    
    # Plot the free energy version
    from matplotlib import ticker, cm
    plt.contour(dz, alpha, fes_adj, levels=range(0,60,5), linewidths=0.5, colors='k')
    cntr = plt.contourf(dz, alpha, fes_adj, levels=range(0,60), cmap=cm.rainbow)
    plt.colorbar(cntr, label="FES [kJ/mol]")
    # If I'm interpreting this correctly, this one should eventually be flat!
    
    # Normal reweighing (exponential?)
    data_as = plumed.read_as_pandas("unbias_attempt1/ff_z_alpha_as.dat")
    npoints = 256
    dz_as    = np.array(data_as["dz"]).reshape(npoints,npoints)
    alpha_as = np.array(data_as["alpha"]).reshape(npoints,npoints)
    fes_as   = np.array(data_as["ff_z_alpha_as"]).reshape(npoints,npoints)

    fes_adj_as = fes_as - np.min(fes_as)

    print(f"AS max after subtraction: {np.nanmax(fes_adj_as[fes_adj_as != np.inf])}")
    
    plt.contour(dz_as, alpha_as, fes_adj_as, levels=range(0,150,5), linewidths=0.5, colors='k')
    cntr = plt.contourf(dz_as, alpha_as, fes_adj_as, levels=range(0,150), cmap=cm.rainbow)
    plt.colorbar(cntr, label="FES [kJ/mol]")
    
This should give you two plots, one for the biased energy landscape, and one for the unbiased energy landscape.  You can then use the same machinery to look at all of the free energy grids.


# VMD tips and tricks (groan)

## Simple display
If you want a really easy quick and dirty way to display the molecules in VMD, you can
use the following command to load a gromacs topology and trajectory into VMD from the
command line (OSX) like

	alias vmd=/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/MacOS/startup.command
	vmd step7_1.gro traj_continuous_v1_100.xtc

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

Just use the following. Center the box on the lipid layer, so that we can see the protein moving, and then run the
analysis. You still have to make the new groups as before (see above).

    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.pdb -n extra_groups.ndx -boxcenter rect -pbc mol -center -dump 0
    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.xtc -n extra_groups.ndx -boxcenter rect -pbc mol -center

## Some sort of workflow for GROMACS simulations
This is done using version 1 of the umbrella sampling that I did. This is to visualize the entire sequence in some sort of rational form. First, preprocess this shit to make it easier in PyMOL. I like to use group 19 (DOPC + PLPI) for centering, and then write out all of the lipid and helix coordinates (usually group protein + DOPC + PLPI). This centers the box around the bilayer, but keeps the protein coordinates as well in the end.

    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.pdb -n extra_groups.ndx -boxcenter rect -pbc mol -center -dump 0
    gmx trjconv -s step7_umbrella_v1.tpr -f step7_umbrella_v1.xtc -o step7_umbrella_v1_reduced.xtc -n extra_groups.ndx -boxcenter rect -pbc mol -center

Then, you can run the pymol script that we have here, and it will generate the correct viewing of the AH-domain and membrane simulations! You should set it to load into the correct directory, and load the correct trajectory.

    pymol setup_graphics_pymol.pml

If you generate a movie, you can combine it togther with the following command.

    ffmpeg -y -f image2 -i frame_%4d.png -vcodec libx264 -profile baseline -pix_fmt yuv420p -r 15 -qscale -0.8 test.mov

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

## Density
We also probably want the density of various systems. GROMACS provides a way to extract this information for us. We can create density groups for various portions of the lipids and the amphipathic helix, in theory. To create the density groups, we need to go into GROMACS and generate the information for each lipid in the system.

	gmx make_ndx -f step7_20.tpr -o density_groups.ndx
	
The, create a density group for the following. The phosphate groups are the PO4- selection, the carbonyls are the C double bonded to O selection(s), and the terminal methyl groups are just the end of the two chains (neglect the hydrogens, just take the carbons as the position). You should also always make a combined lipid bilayer set of information. Following are the examples for a DOPC/PLPI bilayer with a helix, otherwise, the numbers will be off.

### General densities
	r DOPC PLPI
	name 18 lipids
	
	1 | 18
	name 19 helix_lipids

### DOPC

	13 & a P | a O11 | a O12 | a O13 | a O14
	name 20 dopc_phosphates
	
	13 & a C21 | a O22 | a C31 | a O32
	name 21 dopc_carbonyls
	
	13 & a C218 | a C318
	name 22 dopc_terminalmethyls
	
### PLPI

	14 & a P | a O11 | a O12 | a O13 | a O14
	name 23 plpi_phosphates
	
	14 & a C21 | a O22 | a C31 | a O32
	name 24 plpi_carbonyls
	
	14 & a C218 | a C316
	name 25 plpi_terminalmethyls
	
### Combined densities

	20 | 23
	name 26 all_phosphates
	
	21 | 24
	name 27 all_carbonyls
	
	22 | 25
	name 28 all_terminalmethyls
	
### Cheat file
This file allows us to cheat the make\_ndx interactive portion, as it sucks to do. Put this into create\_density\_map.txt.

	r DOPC PLPI
	name 18 lipids
	1 | 18
	name 19 helix_lipids
	13 & a P | a O11 | a O12 | a O13 | a O14
	name 20 dopc_phosphates
	13 & a C21 | a O22 | a C31 | a O32
	name 21 dopc_carbonyls
	13 & a C218 | a C318
	name 22 dopc_terminalmethyls
	14 & a P | a O11 | a O12 | a O13 | a O14
	name 23 plpi_phosphates
	14 & a C21 | a O22 | a C31 | a O32
	name 24 plpi_carbonyls
	14 & a C218 | a C316
	name 25 plpi_terminalmethyls
	20 | 23
	name 26 all_phosphates
	21 | 24
	name 27 all_carbonyls
	22 | 25
	name 28 all_terminalmethyls
	q
	
Then you can run the command.

	gmx make_ndx -f step7_20.tpr -o density_groups.ndx < create_density_map.txt

### Cheat file 2
This file is for the more complicated Ronit Freeman group monomer.

	r CHL1 POPC DOPE DOPS SAPI24 SAPI25
	name 22 lipids
	1 | 22
	name 23 helix_lipids
	14 & a P | a O11 | a O12 | a O13 | a O14
	name 24 popc_phosphates
	14 & a C21 | a O22 | a C31 | a O32
	name 25 popc_carbonyls
	14 & a C218 | a C316
	name 26 popc_terminalmethyls
	15 & a P | a O11 | a O12 | a O13 | a O14
	name 27 dope_phosphates
	15 & a C21 | a O22 | a C31 | a O32
	name 28 dope_carbonyls
	15 & a C218 | a C318
	name 29 dope_terminalmethyls
	16 & a P | a O11 | a O12 | a O13 | a O14
	name 30 dops_phosphates
	16 & a C21 | a O22 | a C31 | a O32
	name 31 dops_carbonyls
	16 & a C218 | a C318
	name 32 dops_terminalmethyls
	17 & a P | a O11 | a O12 | a O13 | a O14
	name 33 sapi25_phosphates
	17 & a C21 | a O22 | a C31 | a O32
	name 34 sapi25_carbonyls
	17 & a C220 | a C318
	name 35 sapi25_terminalmethyls
	18 & a P | a O11 | a O12 | a O13 | a O14
	name 36 sapi24_phosphates
	18 & a C21 | a O22 | a C31 | a O32
	name 37 sapi24_carbonyls
	18 & a C220 | a C318
	name 38 sapi24_terminalmethyls
	24 | 27 | 30 | 33 | 36
	name 39 all_phosphates
	25 | 28 | 31 | 34 | 37
	name 40 all_carbonyls
	26 | 29 | 32 | 35 | 38
	name 41 all_terminalmethyls
	q

	
### Density calculation itself
You need to figure out which group you are going to center around. For instance, in the above, I would center around group 18 (the lipids). Then, you can write out the density map for everybody else. In this case, I wrote out the density for 9 groups, which includes the ions. I also only do the density after 100,000 ps (halfway through our 200 ns simulations), as they should have equilibrated better by then. I also divide the box into 100 slices, rather than 50, to give finer details.

	gmx density -s step7_20.tpr -f traj_continuous_v1_20.xtc -n density_groups.ndx -o density_map.xvg -d Z -center -ng 6 -b 100000 -sl 100
	
### Charge density
Charge densities can be extracted in a similar manner to this, by using the gmx potential function. For example, to calculate the potential of the entire system, one can use a setup similar to the gmx density command. Note, this should be done with the complete trajectory, as we need information about the ions and water to make sure we get everything correct.

	gmx potential -s step7_20.tpr -f traj_continuous_v1_20.xtc -n density_groups.ndx -o potential_map.xvg -d Z -center -b 100000 -sl 100 -ng 5
	
# Lipid preparation

## Ronit Freeman bilayer
Based off of Beber 2019 (or earlier work, see Zotero). Complicated bilayer, with the following conditions (ratios are preserved, not percentages).

	Cholesterol 30
	POPC 113 (stand-in for Egg PC)
	DOPE 20
	DOPS 20
	SAPI24 8
	SAPI25 8

Everything here is done at 50 mM NaCl unless stated explicitly otherwise.
	
	
	
	




