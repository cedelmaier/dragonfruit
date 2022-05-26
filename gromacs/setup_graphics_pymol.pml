# Load the PDB file and the trajectory
# Run via the command run setup_graphics_pymol.pml

#cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/coiled_umbrella/v1
#load "step7_umbrella_v1_reduced.pdb"
#load_traj "step7_umbrella_v1_reduced.xtc"

cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/coiled/zdepth_00angstroms/s1
load "traj_continuous_v1_1000_reduced.pdb"
load_traj "traj_continuous_v1_1000_reduced.xtc", start=9500, stop=10000

#cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/coiled/zdepth_30angstroms/s1
#load "traj_continuous_v1_1000_reduced.pdb"
#load_traj "traj_continuous_v1_1000_reduced.xtc", start=0, stop=500

#cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/unfolded/zdepth_00angstroms/s1
#load "traj_continuous_v1_1000_reduced.pdb"
#load_traj "traj_continuous_v1_1000_reduced.xtc", start=9500, stop=10000

#cd /Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/unfolded/zdepth_10angstroms/s1
#load "traj_continuous_v1_1000_reduced.pdb"
#load_traj "traj_continuous_v1_1000_reduced.xtc", start=1000, stop=1500

# Smooth the trajectory and set up options
smooth all, 30, 3
select DOPC, resname DOPC
select PLPI, resname PLPI
select lipids, resname DOPC or resname PLPI
select helix, chain A

# Set some sort of transparency
set sphere_transparency=0.7, lipids

# Set up selections
select dopc_headgroups, (resname DOPC) and (name N or name C11 or name C12 or name C13 or name C14 or name C15 or name P or name O11 or name O12 or name O13 or name O14) and (not name H*)
select dopc_notheadgroups, (resname DOPC) and (not name N and not name C11 and not name C12 and not name C13 and not name C14 and not name C15 and not name P and not name O11 and not name O12 and not name O13 and not name O14 and not name H*)

select plpi_headgroups, (resname PLPI) and (name P or name C11 or name C12 or name C13 or name C14 or name C15 or name C16 or name O1 or name O2 or name O3 or name O4 or name O5 or name O11 or name O12 or name O13 or name O14) and (not name H*)
select plpi_notheadgroups, (resname PLPI) and (not name P and not name C11 and not name C12 and not name C13 and not name C14 and not name C15 and not name C16 and not name O1 and not name O2 and not name O3 and not name O4 and not name O5 and not name O11 and not name O12 and not name O13 and not name O14 and not name H*)

select hydrophobes, (resn ala+gly+val+ile+leu+phe+met)

# Setup colors
color blue, dopc_headgroups
color blue, plpi_headgroups
color cyan, dopc_notheadgroups
color cyan, plpi_notheadgroups

#color orange, hydrophobes
color magenta, hydrophobes

# Display options
hide all
show spheres, dopc_headgroups
show spheres, plpi_headgroups
show spheres, dopc_notheadgroups
show spheres, plpi_notheadgroups

show cartoon, helix
show sticks, (hydrophobes and (!name c+n+o))

# Final unselect?
disable hydrophobes

# Change view options
rotate x, angle=270, state=0

# Set up movie code
#intra_fit name ca, 1
set ray_trace_frames=1
set cache_frames=0
