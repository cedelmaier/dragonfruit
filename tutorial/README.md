# Dragonfruit (and old_chi_pet) tutorials

This is a generic tutorial for using the old_chi_pet and dragonfruit packages to
run a set of simulations in HOOMD, launch them on longleaf, and then harvest and
analyze the results.

# Brief tutorial
This portion uses old_chi_pet. When you have multiple simulations that you might want to run,
you can use the example.yaml file found in this directory to run a membrane simulation using the
NPT ensemble on longleaf.

old_chi_pet uses YAML files to control the simulation space that is created. For this example, we will
be using a membrane simulation that requires HOOMD to be installed with the GrimeLipid pair potential. This
example runs the simulations at 2 different box sizes, each with 2 seeds, in order to understand the variability
of the outputs. Note that this isn't necessary if you know your system is ergodic (double check on this with others
later). You can see the seed and parameter behavior controlled in the YAML file by the commands

  seed: ChiSeed(bounds = [0,1])
  lbox: ChiParam(format_str = "lbox{}", exec_str = "[100.0, 103.0]") # Size of final box, simulation units (sigma)

There is also an associated args.yaml for all simulations.

You create a simulation run by invoking the command

    Chi.py -C *.yaml -a args.yaml

in the directory where you have your original YAML file. This will create a substructure based on the
parameters that you are choosing in your YAML file.
