#!/bin/sh
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_PROC_BIND=true
mkdir cuda_build; cd cuda_build
cmake -DKokkos_ENABLE_OPENMP=ON \
      -DKokkos_ARCH_VOLTA70=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
      ..
make -j 8

./ShallowWater
./ShallowWater_par1
./ShallowWater_par2
./ShallowWater_par3
./ShallowWater_par4
