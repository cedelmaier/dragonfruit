# dragonfruit
Collection of software and packages for simulating membranes with HOOMD

# Basic structure of code
```
-- Run
   |-- Simulation
       |-- Seed 1
       |-- Seed 2
   `-- Simulation
       |-- Seed 1
       |-- Seed 2
```

# Related packages
You need to also get old_chi_pet (under development) to correctly use the
distributed launch system. Eventually the two packages will be combined.

  git clone git@github.com:cedelmaier/old_chi_pet.git
  git checkout hoomd
