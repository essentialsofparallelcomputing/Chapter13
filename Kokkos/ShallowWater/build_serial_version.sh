#!/bin/sh
mkdir serial_build; cd serial_build
cmake \
      ..
make -j 8

for exec in ShallowWater ShallowWater_par?
do
   echo ""
   echo "Running $exec version"
   ./$exec
   echo "==========================="
done
