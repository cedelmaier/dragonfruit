# Dragonfruit (and old_chi_pet) tutorials

This is a generic tutorial for using the old_chi_pet and dragonfruit packages to
run a set of simulations in HOOMD, launch them on longleaf, and then harvest and
analyze the results.

# Getting old_chi_pet
First, you need to get the correct branch of old_chi_pet from cedelmaier on github.

    git clone git@github.com:cedelmaier/old_chi_pet.git
    git checkout hoomd

# Brief tutorial
This portion uses old_chi_pet. When you have multiple simulations that you might want to run,
you can use the example.yaml file found in this directory to run a membrane simulation using the
NPT ensemble on longleaf.

old_chi_pet uses YAML files to control the simulation space that is created. For this example, we will
be using a membrane simulation that requires HOOMD to be installed with the GrimeLipid pair potential. This
example runs the simulations at 2 different box sizes, each with 2 seeds, in order to understand the variability
of the outputs. Note that this isn't necessary if you know your system is ergodic (double check on this with others
later). You can see the seed and parameter behavior controlled in the YAML file by the commands

    seed: ChiSeed(bounds = [0,2])
    lbox: ChiParam(format_str = "lbox{}", exec_str = "[100.0, 103.0]") # Size of final box, simulation units (sigma)

There is also an associated args.yaml for all simulations. This is going to be modified in the future, as right
now the submission system cannot run different submission states. Currently, you need to modify the args.yaml
file to point to where you have LipidVolta.py and SeptinAnalysis.py installed.

You create a simulation run by invoking the command

    Chi.py -C *.yaml -a args.yaml

in the directory where you have your original YAML file. This will create a substructure based on the
parameters that you are choosing in your YAML file. You then launch the simulations with the following command:

    Chi.py -L simulations/*

and then follow the on-screen prompts. The first thing to do is use the command 'longleaf' when asked what partition
to use.

    List space separated states you wish to run (leave blank for all): <ENTER>
    Input supercomputer name (default disBatch): longleaf
    <HIT ENTER FOR THE REST>

This should submit all of the jobs to the longleaf cluster on the volta-gpu partition, run the simulations, and then
run the analysis (in whatever current form it is in) on the cluster.

# Directly running simulations
If you want to directly run simulations, you need to be setup in a HOOMD environment locally. For instance, for
an anaconda installation, you would need to get into the environment that you initially installed HOOMD into, such
as

    conda activate hoomd_300beta9

Then, you can directly import HOOMD into a python script for running. If you change the old_chi_pet version of the YAML
control file for a specific value (specific values), you can then run the command:

    python3 <path_to_dragonfruit>/src/LipidVolta.py --yaml example.membrane.npt.yaml




