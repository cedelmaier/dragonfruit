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
