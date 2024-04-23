#!/bin/bash

# from {$maindir}/sundials/builddir, run
# cp ../../summa/build/cmake_external/build_cmakeSundials.bash build_cmake
# run script from the builddir directory with ./build_cmake
# run `make`, then `make install`
# Note, using -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF, if want to run examples should change

cmake ../../sundials-software/ -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DCMAKE_Fortran_COMPILER=gfortran -DENABLE_MAGMA=ON -DMAGMA_DIR=  -DENABLE_CUDA=ON -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples
