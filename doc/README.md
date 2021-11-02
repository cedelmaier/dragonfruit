# Septin project

## Running verion of different systems

### Bilayer alone
This is just a bilayer on it's own, 75:25 DOPC:PLPI, 150 mM NaCl, 300K.
    
    bi_75dopc_25plpi_base
    93311 atoms

### Bilayer with AH
Created a bilayer with AH from frame 51 of the melting process. All other variables
remain the same.

    bi_75dopc_25plpi_ahv1
    136801 atoms

### Bilayer with close AH
Created a bilayer with AH from frame 50 of the metling process. Placed the AH domain very
close in Z to one of the bilayers to promote interactions.

    bi_75dopc_25plpi_ahv2
    113223 atoms

## Recipes
General recipes for taking notes on the septin project.

### Interactively running on nodes on longleaf
I put these into my .bashrc so that I could use them, but you can just run the commands.

    alias sinteractive="srun -t 8:00:00 -p interact -N 1 --mem=6G --x11=first --pty /bin/bash"
    alias sinteractivecompile="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --pty /bin/bash"
    alias sinteractivegpu="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --partition=gpu --gres=gpu:1 --qos=gpu_access --pty /bin/bash"
    alias sinteractivevolta="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --partition=volta-gpu --gres=gpu:1 --qos=gpu_access --pty /bin/bash"

### Pulling down information to analyze locally with rsync

    rsync -avh -e ssh /source/path/ host:/destination/path

### Installing HOOMD-blue on longleaf cluster (UNC)
This is the first set of information for getting HOOMD and switching to the most recent (v3.0.0-beta.10) tag. This will
also load the correct modules on longleaf. Warning, this is not complete. You should also be following the
instructions for hoomd installation online.

    https://github.com/glotzerlab/hoomd-blue
 
    git clone --recursive https://github.com/glotzerlab/hoomd-blue
    cd hoomd-blue
    git fetch --all --tags
    git checkout tags/v3.0.0-beta.10 -b v3.0.0-mahsa
    
    module load git
    module load cmake
    module load python/3.8.8
    module load cuda/11.2
    module load gcc/9.1.0

You also need the custom code that CJE has for the soft sine-based lipid potential.

    https://github.com/cedelmaier/dragonfruit/

Create a python virtual environment and source it.

    python3 -m venv ~/virtual_envs/hoomd300beta9 --system-site-packages
    source ~/virtual_envs/hoomd300beta9/bin/activate
    python3 install-prereq-headers.py
   
Compile with a virtual envinronment python, as well as CUDA. Note: this will take quite a while, and I did it on the
native GPU node. Be careful with GNinja, as it can cause memory errors.

    CC=gcc CXX=g++ cmake -B build/hoomd hoomd-blue
    ccmake .
    cmake --build ./
    cmake --install ./

### Gromacs tips and tricks
To concatenate a bunch of different files into gromacs, for instance, to create a continuous trajectory, run the following
commands. The first command will create the list of files to run on for the latter command.

    ./concatenate_traj_v1.sh
    gmx trjcat -f `cat list_trr_oneline.txt` -o traj_continuous.xtc -settime < times_trr.txt

To do the same with energy files, run the following commands

    ./concatenate_edr_v1.sh
    gmx eneconv -f `cat list_edr_oneline.txt` -o edr_continuous_v1.edr -settime < times_edr.txt

For our simulations, GROMACS seems to run the best with commands similar to the following (performance testing is an
ongoing process).

    gmx_gpu mdrun -v -deffnm ${istep} -ntomp 8 -ntmpi 1 -nb gpu -bonded gpu -pme gpu

### PackMem
This is an attempt to get PackMem to work with GROMACS 2021.1. There is a mismatch between the named atoms in the old PDB file format, and
the new ones. The current workaround is to do the analysis on longleaf, and then transfer the files for analysis on a local
machine.

To dump the GROMACS files in a way that is consistent with the online tutorial, use this command. This is based on the example of
the lipid bilayer generated for a DOPC and PLPI set of layers with GROMACS 2021.2, combined on a mac (see below).

    gmx trjconv -f traj_continuous.xtc -s step7_1.tpr -o pdb/bilayertest.pdb -sep -n index.ndx -dt 100

To use PackMem, you have to first make sure to have the lipid properly encoded into the text files associated with PackMem. This means
that you might have to convert the lipid and force-field files into the VDW radius. I created a helper script to do this. The only 
addition this script needs is to go in by hand and label the amphipathic parts of the lipid chain, as that is not done automatically.
Additionally, you need to be able to go in and change the order of the lipid atom names for the PDB, as that is what causes the 
mismatch issues between GROMACS 5 and later versions. These files that need changing are the following.

    param_Charmm.txt
    vdw_radii_Charmm.txt

You can then run the PackMem script on the files that have been generated, via the following command.

    bash ScriptPackMem.sh

Note that this will take a while to run (about 2 minutes per PDB file on my laptop). To use R in the Terminal on a mac, 
run the following command, which will compute the statistics of the different defects per bilayer.

    export PATH="/Library/Frameworks/R.framework/Resources:$PATH"
    R --vanilla < Script_fit_and_plot.R

### Visualizing lipids from GROMACS with VMD

First, you can access VMD on the command line via this command.

    /Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/MacOS/startup.command

The lipids are usually selected by the following commands. The second command undoes the
boundary conditions so that you don't have to carry all the waters together to unwrap
the system.

    resname DOPC and resname PLPI
    pbc unwrap -all -sel "protein or resname DOPC or resname PLPI"
    pbc wrap -compound fragment -center com -centersel "resname DOPC or resname PLPI" -sel "resname DOPC or resname PLPI or protein" -all

It gets quite a bit more complicated from here. To visualize the separate headgroup for DOPC and PLPI, we
must go through several different stages of parsing with VMD. Here are some of the more salient points. This gets
the headgroups (or not headgroups), and then gives them a custom color so that we can visualize them.

    source ~/Projects/Biophysics/septin_project/updated_performance/bi_75dopc_25plpi_ahv3/setup_graphics.tcl

### Creating a PDB file from the last frame in GROMACS

For the AH simulations we need the final frame from the NPT AH runs in GROMACS, in PDB form. Read in the above XTC file (for example,
npt_production_isotropic.xtc) into VMD. Then, correct the periodic boundary conditions. Then, finally, use VMD to dump the molecule
in PDB format using Save Coordinates.

One can figure out the periodic boundary conditions and center the protein via using the following.

    pbc unwrap -all
    set sel [atomselect top "protein"]
    measure center $sel
    measure minmax $sel

### Density analysis in GROMACS

To do the density analysis, or just for good bookeeping, it is ideal to keep track of which groups are which. I have named these
as such (have Brandy double-check them). Note, the selection groups need to be changed if a different order of the groups
is used.

To create the density map itself, run this command.

    gmx make_ndx -f step7_1.tpr -o density_groups.ndx

Then, for each of the following, you can define groups to allow you to create density profiles.

    DOPC_Headgroups
    2 & a C12 | a C13 | a C14 | a C15 | a C11 | a P | a O13 | a O14 | a O12 | a O11
    (resname DOPC)and(name "N" "C12" "C13" "C14" "C15" "C11" "P" "O13" "O12" "O11")
    set dopc_headgroups [atomselect top "(resname DOPC)and(name 'N' 'C12' 'C13' 'C14' 'C15' 'C11' 'P' 'O13' 'O12' 'O11')"]

    DOPC_Notheads
    2 & ! a C12 & ! a C13 & ! a C14 & ! a C15 & ! a C11 & ! a P & ! a O13 & ! a O14 & !a O12 & !a O11
    (resname DOPC)and(not name "N" "C12" "C13" "C14" "C15" "C11" "P" "O13" "O12" "O11")and(not name "H.*")

    PLPI_Headgroups
    3 & a C12 | a O2 | a C13 | a O3 | a C14 | a O4 | a C15 | a O5 | a C16 | a O6 | a P | a O11 | a O12 | a O13 | a O14 

    PLPI_Notheads
    3 & ! a C12 & ! a O2 & ! a C13 & ! a O3 & ! a C14 & ! a O4 & ! a C15 & ! a O5 & ! a C16 & ! a O6 & ! a P & ! a O11 & ! a O12 & ! a O13 & ! a O14 

The density command itself is run by things like.

    gmx density -s step7_1.tpr -f traj_continuous_v1_200.xtc -n density_groups.ndx -o dens_headgroups.xvg -d Z -center -ng 3













