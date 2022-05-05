# Run this script
# source ~/Projects/Biophysics/septin_project/updated_performance/bi_75dopc_25plpi_ahv3/setup_graphics.tcl
# source /Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/gromacs/setup_graphics.tcl

# Before doing anything else, unwrap the PBCs that we might be dealing with for the lipids and the protein
#pbc unwrap -all -sel "protein or resname DOPC or resname PLPI"
#pbc wrap -compound fragment -center com -centersel "resname DOPC or resname PLPI" -sel "resname DOPC or resname PLPI or protein" -all

# Go back to 0
animate goto 0

# First, setup the different representations that we want
# Get the protein setup in the proper circumstances
mol modselect 0 0 protein
mol modstyle 0 0 NewRibbons 0.6 12.0 3.0 0
mol modcolor 0 0 ResType
mol smoothrep 0 0 3

# Create new representations for the other lipids and stuff
mol addrep 0
mol modselect 1 0 user 1.0
mol modcolor 1 0 ColorID 0
mol modstyle 1 0 VDW 1.0 12.0
mol modmaterial 1 0 Transparent
mol smoothrep 0 1 3

mol addrep 0
mol modselect 2 0 user 2.0
mol modcolor 2 0 ColorID 0
mol modstyle 2 0 VDW 1.0 12.0
mol modmaterial 2 0 Transparent
mol smoothrep 0 2 3

mol addrep 0
mol modselect 3 0 user 1.1
mol modcolor 3 0 ColorID 10
mol modstyle 3 0 VDW 1.0 12.0
mol modmaterial 3 0 Transparent
mol smoothrep 0 3 3

mol addrep 0
mol modselect 4 0 user 2.1
mol modcolor 4 0 ColorID 10
mol modstyle 4 0 VDW 1.0 12.0
mol modmaterial 4 0 Transparent
mol smoothrep 0 4 3

# This sets the headgroups to the proper thing
set dopc_headgroups [atomselect top "resname DOPC and name N C12 C13 C14 C15 C11 P O13 O12 O11 "]
set dopc_notheadgroups [atomselect top "(resname DOPC) and (not name N C12 C13 C14 C15 C11 P O13 O12 O11 \"H.*\") "]
$dopc_headgroups set user 1.0
$dopc_notheadgroups set user 1.1

set plpi_headgroups [atomselect top "resname PLPI and name C12 O2 C13 O3 C14 O4 C15 O5 C16 O6 P O11 O12 O13 O14 "]
set plpi_notheadgroups [atomselect top "(resname PLPI) and (not name C12 O2 C13 O3 C14 O4 C15 O5 C16 O6 P O11 O12 O13 O14 \"H.*\") "]
$plpi_headgroups set user 2.0
$plpi_notheadgroups set user 2.1

# Choose your rotation
scale by 1.2
scale by 1.2
scale by 1.2
scale by 1.2

# Rotate into new X position
#rotate x by 270
