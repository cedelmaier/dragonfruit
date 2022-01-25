# Dragonfruit
Collection of software and packages for simulating membranes with HOOMD

# Structure
This is the simulation structure setup
```bash
-- Run
   |-- Simulation (parameter value 1)
       |-- Seed 1
       |-- Seed 2
   `-- Simulation (parameter value 2)
       |-- Seed 1
       |-- Seed 2
```
and the corresponding directory structure (for 2 choices of box_size and 2 seeds each).

```bash
-- simulation/
   |-- box_size100.0/
       |-- s0/
       |-- s1/
   `-- box_size100.0/
       |-- s0/
       |-- s1/
```


These are easiest to describe in a bottom-up approach (for now).

## Simulation Seed
Each simulation seed is a specific choice of initial random number state for the simulation. This is done
for non-ergodic systems, where the time-average and ensemble-average are not guaranteed to be the same.
Multiple random number seeds are created in order to understand the ensemble behavior of the simulation under
question. So that is what is done here. For example, if you had a parameter box_size, and two values of box_size,
each one might have two associated seeds for each box_size, leading to 4 "real" simulations in total.

## Simulation Simulation
The name is confusing, so maybe change at some point. A simulation is the collection of all seeds run at a specific
point in the parameter space.

## Simulation Run
The collection of all Simulation Simulations for a set of parameter choices (sweeps) in parameter space.

# Setting up dragonfruit and old_chi_pet
You need to do some special setup after cloning the repositories for dragonfruit and old_chi_pet. Usually the
best way to get access to the programs of interest is to symlink them into a bin/ directory somewhere that your
path can access. I have mind setup as

    ~/bin/cedelmaier

and then symlink in the relevant files via

    ln -s <install location>/dragonfruit/src/LipidVolta.py LipidVolta
    ln -s <install location>/old_chi_pet/Chi.py Chi

# Related packages
You need to also get old_chi_pet (under development) to correctly use the
distributed launch system. Eventually the two packages will be combined.

    git clone git@github.com:cedelmaier/old_chi_pet.git
    git checkout hoomd
