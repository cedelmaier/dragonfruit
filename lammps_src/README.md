# Compiling LAMMPS
Here are the commands to compile LAMMPS with python support. Note, this may not be complete, so
take care to make sure things like the proper dependencies are installed.

    git clone -b develop https://github.com/lammps/lammps.git mylammps

Set up cmake to work correctly, and then enable a bunch of options. From what I remember, it is like this
although then more options need to be enabled. You also need the pair_grime.h and pair_grime.cpp code
found in this directory in the proper location.

    cmake -C ../cmake/presets/most.cmake -DBUILD_MPI=yes -DBUILD_OPENMP=yes -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_C_COMPILER=gcc-11 -DBUILD_SHARED_LIBS=yes

and then use ccmake . to adjust the values to get the python package enabled. Build with

    cmake --build .
    cmake --install .

into the proper environment. Then you also need to make sure that the shared library is found on OSX

    export DYLD_LIBRARY_PATH=$HOME/.local/lib:$DYLD_LIBRARY_PATH

# Compiling LAMMPS on longleaf
Can't use MPI, sad face. Use a build environment that is set up on the pine system. First, however, make a virtual environment.

    python3 -m venv ~/virtual_envs/lammps_20220309
    source ~/virtual_envs/lammps_20220309/bin/activate

for example. Then, make sure to set the CMAKE_INSTALL_PREFIX for the installation, as this will otherwise screw up where it sends the files for the pythong installation.

    cmake -C ~/mylammps/cmake/presets/most.cmake -DBUILD_OPENMP=yes -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_Fortran_COMPILER=gfortran -DBUILD_SHARED_LIBS=yes ~/mylammps/cmake

Then make sure to run ccmake . and check everything. Enable the following, and also make sure we have a virtual environment set up.

    Python
    LAMMPS_Exceptions

Then you can build the setup

    cmake --build . -j8
    cmake --install .

You also need to add a line to the activate script for the python stuff. This can easily be accomplished by

    echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> ~/virtual_envs/lammps_20220309/bin/activate
