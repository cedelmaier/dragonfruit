# Load the PDB file and the trajectory

#load "traj_continuous_v1_100_reduced.pdb"
#load_traj "traj_continuous_v1_100_reduced.xtc", start=0, stop=1000
load "traj_continuous_v1_200_reduced.pdb"
load_traj "traj_continuous_v1_200_reduced10.xtc"

########
# Actual program
########

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
select dopc_phosphates, (resname DOPC) and (name P or name O11 or name O12 or name O13 or name O14) and (not name H*)
select dopc_notheadgroups, (resname DOPC) and (not name N and not name C11 and not name C12 and not name C13 and not name C14 and not name C15 and not name P and not name O11 and not name O12 and not name O13 and not name O14 and not name H*)

select plpi_headgroups, (resname PLPI) and (name P or name C11 or name C12 or name C13 or name C14 or name C15 or name C16 or name O1 or name O2 or name O3 or name O4 or name O5 or name O6 or name O11 or name O12 or name O13 or name O14) and (not name H*)
select plpi_phosphates, (resname PLPI) and (name P or name O11 or name O12 or name O13 or name O14) and (not name H*)
select plpi_notheadgroups, (resname PLPI) and (not name P and not name C11 and not name C12 and not name C13 and not name C14 and not name C15 and not name C16 and not name O1 and not name O2 and not name O3 and not name O4 and not name O5 and not name O6 and not name O11 and not name O12 and not name O13 and not name O14 and not name H*)

#select hydrophobes, (resn ala+gly+val+ile+leu+phe+met)
select hydrophobes, (resn ala+ile+leu+met+val+phe+trp+tyr+gly)
select cterm, (resi 18)

# Setup colors
color blue, dopc_headgroups
color blue, plpi_headgroups
color cyan, dopc_notheadgroups
color cyan, plpi_notheadgroups
color orange, dopc_phosphates
color orange, plpi_phosphates

#color orange, hydrophobes
color magenta, hydrophobes
color red, cterm

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
disable cterm

# Change view options
rotate x, angle=270, state=0

#########
## Movie generation
#########
## Set up movie code
##intra_fit name ca, 1
#set ray_trace_frames=1
#set cache_frames=0
#
## This is what actually djmps a bunch of PNG files for a movie
#mpng frame_

