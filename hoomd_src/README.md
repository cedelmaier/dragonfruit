# source code
This is a repository for source code that might go into HOOMD, kept separate, until I can get that
branch/fork working correctly.

# Compiling HOOMD
Here are the commands (somewhat up to date) for compiling HOOMD in different environments.

## OSX
Here, we can setup an anaconda environment to take care of everything, which makes stuff a lot easier. For
example, for build 3.2.0, we can use the following to generate a new conda environment.

    conda create -c conda-forge --name hoomd_320 python=3.9.6

and then load the environment with the proper stuff

    conda install -c conda-forge gsd freud fresnel jupyter notebook pytest pybind11 matplotlib pytables cereal pandas scipy pyyaml numpy

then things get more complicated. Because of the ever-evolving nature of HOOMD, including the special potentials requires
going in and setting them up by hand, as the way in which files are defined keeps changing. The EvalulatorPairGrimeLipid.h
file should be roughly in the right shape, and this can be used as a prototype for everything else.

## longleaf (UNC)
I like to login to a compile node (you can use one without graphics for compilation) for these steps, although
you can also login to a GPU node. I put these commands in my .bashrc so that I have easy access to different
options.

    alias sinteractive="srun -t 8:00:00 -p interact -N 1 --mem=6G --x11=first --pty /bin/bash"
    alias sinteractivecompile="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --pty /bin/bash"
    alias sinteractivegpu="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --partition=gpu --gres=gpu:1 --qos=gpu_access --pty /bin/bash"
    alias sinteractivevolta="srun --ntasks=1 --cpus-per-task=8 --mem=8G --time=8:00:00 --partition=volta-gpu --gres=gpu:1 --qos=gpu_access --pty /bin/bash"

Load the correct modules to compile HOOMD.

    module load git
    module load cmake
    module load python/3.9.6
    module load gcc/11.2.0
    module load cuda/11.8

Create a python virtual environment to load everything into with 'install' command later.
NOTE: You can get away with doing this various ways. Just be sure you know what exactly
you are doing! I put these in a special place

    python3 -m venv virtual_envs/hoomd381
    source ~/virtual_envs/hoomd381/bin/activate

Get the HOOMD source code and put it somewhere convenient (also install the prerequisites for 
compilation).

    git clone --recursive https://github.com/glotzerlab/hoomd-blue
    cd hoomd-blue
    git fetch --all --tags
    git checkout tags/v3.8.1
    python3 ./install-prereq-headers.py

Any other prerequisites needed?

    python3 -m pip install gsd
    python3 -m pip install freud-analysis

Copy in the files that are necessary from this repository for the additional features in HOOMD.
These features include the 'soft' excluded volume potential as well as a prototype Born-Mayer-Huggins
potential for another project. These are easier to compile directly into the HOOMD potentials, rather
than use an external module.

I do my compilation in the scratch space.

    cd /work/<youpathhere>
    CC=gcc CXX=g++ cmake -B build/hoomd -S ~/hoomd-blue/
    cd build/hoomd

This requires specifying the C++ compiler for some reason, hence the CC and CXX commands for the
cmake section. Now edit the cmake configuration to enable GPU computing.

    ccmake .
    <EDIT THE ENABLE_GPU FEATURE, THEN PRESS C>

Now you should be ready for the compilation of CUDA itself for GPU processing. I do this with only 2 threads
beacuse with CUDA, sometimes you can run out of memory on the compile nodes and get weird, non-specific
errors that are hard to track down. NOTE: This will take a long time to complete.

    cmake --build ./ -j2
    cmake --install ./
    
After this, it is a good idea to run the tests recommended by the HOOMD website in order to make sure that everything is okay. Both the C++ and python tests are useful. See the main website for details.

# Data I/O
Here is an example of how to transfer entire directories between longleaf and a local machine.

    scp -r edelmaie@longleaf.unc.edu:/pine/scr/e/d/edelmaie/supra_cg/dragonfruit/data/20220204 .
    scp -r <uname>@longleaf.unc.edu:<pine location>/<directory> .

In order to get stuff off of longleaf that doesn't contain the .gsd files (as they are huge), one
can use the following command to download everything else.

    rsync -avh -e ssh --exclude='*.gsd' edelmaie@longleaf.unc.edu:/pine/scr/e/d/edelmaie/dragonfruit/20220227 .

# Flatiron
Modules you might want to compile with both GPU and MPI support!

As of 20230209:

	module load modules/2.1
	module load cuda/11.8
	module load gcc/11.3.0
	module load python/3.10.8
	module load openmpi/cuda-4.0.7
	module load cmake/3.25.1




    
    
    


