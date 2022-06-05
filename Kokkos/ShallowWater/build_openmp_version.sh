#!/bin/sh
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_PROC_BIND=true

mkdir openmp_build; cd openmp_build
cmake -DKokkos_ENABLE_OPENMP=ON \
      ..
make -j 8

for exec in ShallowWater ShallowWater_par?
do
   echo ""
   echo "Running $exec version"
   ./$exec
   echo "==========================="
done
