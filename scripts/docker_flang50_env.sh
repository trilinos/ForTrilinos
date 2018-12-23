# Load Spack env
SPACK_ROOT=$1
source $SPACK_ROOT/share/spack/setup-env.sh

# Order matters! For some reason, llvm also provides
# "flang" as a symlink to clang
spack load llvm
spack load flang
spack load pgmath

# Use flang and clang instead of gfortran and gcc
export CC=mpicc
export CXX=mpicxx
export CC=clang
export CXX=clang++
export FC=flang
export OMPI_CC=clang
export OMPI_CXX=clang++
export OMPI_FC=flang

