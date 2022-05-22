#!/bin/sh
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_PROC_BIND=true

mkdir openmp_build; cd openmp_build
cmake -DKokkos_ENABLE_OPENMP=ON \
      ..
make -j 8

./ShallowWater
./ShallowWater_par1
./ShallowWater_par2
./ShallowWater_par3
./ShallowWater_par4
