# Masterclass 3

## Preparing files for viewing
Once you have a simulation trajectory, you need to prepare the files to view. There are several options for how to do this, and so this is just one possible workflow you can use. First, you need to concatenate your trajectories together, strip out the unnecessary atoms for visualization, and then finally visualize the trajectory.

### Concatenating the files together
I made a script called `concatenate_trajedr_v1.sh` found in the repository. To use this, make sure that you have the proper modules loaded on your machine, and then you can run the file on the server to merge the trajectories together. For longleaf, you should just be able to load the modules relevant for our gromacs and then go from there.

	module purge
	module load gcc/9.1.0
	source /nas/longleaf/apps/gromacs/2021.5/avx2_256/bin/GMXRC
	
This will load gromacs (not the GPU version!) for you to use. Then, in the directory you have your output files in, you can simply run the following command.

	<path to concatenate>/concatenate_trajedr_v1.sh 1 100

Again, this just stiches all of the files together for yourself into a `.xtc` file (GROMACS trajectory format) for you to use. We gave it the command to do 100 files, so this should match what we've been doing, but it will take some number of sequential files and stitch them together based on what numbers you give it. Now, we just need to get a reduced set of the atoms in the system to make visualization easier.

### Reducing the number of atoms in the trajectory
Usually we don't want all of the atoms in the trajectory for visualization, and so just a subset of them will be good. To do this, on the server, we need to strip out the extra atoms. First, we need to make a new index file so that the next steps are easier. An Index file just sorts the atoms into different groups for ease. Create the file `create_density_map.txt` and place the following into it.

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

We can then run the following command to generate a new index file.

	gmx make_ndx -f step7_1.tpr -o density_groups.ndx < create_density_map.txt
	
Now you can actually created reduced trajectory files. Issue the following commands.

	gmx trjconv -s step7_1.tpr -f traj_continuous_v1_100.xtc -o traj_continuous_v1_100reduced.pdb -n density_groups.ndx -boxcenter rect -pbc mol -center -dump 0
	gmx trjconv -s step7_1.tpr -f traj_continuous_v1_100.xtc -o traj_continuous_v1_100reduced.xtc -n density_groups.ndx -boxcenter rect -pbc mol -center
	
For each of these, choose the lipids as the centering option (group 18) and write out the combined helix/lipid group (group 19) (this is if the density creation command worked). You now have these two files, `traj_continuous_v1_100reduced.pdb` and `traj_continuous_v1_100reduced.xtc` for use. Download these somewhere convenient on your computer.

### Visualizing the trajectory
There is another file in this directory of the repository `setup_graphics_pymol.pml`. This file sets up some nice graphics specifically for peptide+membrane simulations. Also, if you didn't have anything to follow along with in the last section, there are 2 trajectory files in this directory that you can use. Simply run PyMOL, and then execute the script when you're in the same directory.

## Analyzing the trajectory (still under construction)
This is going to be interesting, as I'm not exactly sure how we will do this on the windows machines. What I normally do is bring back the trajectories (the full trajectories, not the reduced ones) for analysis. A HUGE part of dragonfruit is devoted to the analysis of the trajectories, but I have never used it on a windows machine. It might run on the cluster in the correct conditions.

The basics of the analysis is in the file `src/GromacsAnalysis.py`. This is the master file that will take a single random number 'seed' and proceed through my *current* full analysis. The big problem is that it requires a BUNCH of extra python installs to work correctly.

### Installing analysis prerequisites
If you have anaconda installed, you can create a new environment to run this analysis in. I always start with a fresh installation.

	conda create -c conda-forge --name mda python=3.10.4
	
Then add the modules that we will need (there are a lot).

	conda install -c conda-forge gsd freud fresnel jupyter notebook pytest pybind11 matplotlib pytables cereal pandas scipy pyyaml numpy mdanalysis mdtraj plumed
	
These take a while to install. Then, you also need to separately (via pip) install membrane-curvature since it isn't in conda-forge at the moment.

	pip install membrane-curvature
	
Now you can run the analysis (in theory).

### Running the analysis
The analysis is, again, in `src/GromacsAnalysis.py`. We need a `config.yaml` file to run the analysis that tells it what structure, trajectory, etc files to use, so create that now with the following.

	structure: step7_1.tpr
	trajectory: traj_continuous_v1_100.xtc
	gromacs: step7_1.gro

Then you can run the analysis. To see what the analysis options are, you can run the following command.

	python3 <pathtorepo>/src/GromacsAnalysis.py -sd --yaml config.yaml -A -G -W -F

This runs the full analysis with the 'Force' option (-F). There is a trick inside of the analysis that if the analysis has already been run, it won't rerun the whole thing. The 'Force' option forces this to happen, so that if you change something in the analysis, you can then rerun everything. The graph option (-G) will force graphing of the resultant variables.






