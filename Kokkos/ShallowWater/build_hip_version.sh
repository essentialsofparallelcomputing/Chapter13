#!/bin/sh
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_PROC_BIND=true

mkdir hip_build; cd hip_build
cmake -DKokkos_ENABLE_HIP=ON -DCMAKE_CXX_COMPILER=hipcc ..
make -j 8

for exec in ShallowWater ShallowWater_par?
do
   echo ""
   echo "Running $exec version"
   ./$exec
   echo "==========================="
done
