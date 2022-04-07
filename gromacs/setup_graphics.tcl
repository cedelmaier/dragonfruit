# Run this script
# source ~/Projects/Biophysics/septin_project/updated_performance/bi_75dopc_25plpi_ahv3/setup_graphics.tcl

# Before doing anything else, unwrap the PBCs that we might be dealing with for the lipids and the protein
#pbc unwrap -all -sel "protein or resname DOPC or resname PLPI"
#pbc wrap -compound fragment -center com -centersel "resname DOPC or resname PLPI" -sel "resname DOPC or resname PLPI or protein" -all

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
rotate x by 270
