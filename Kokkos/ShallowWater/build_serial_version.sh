#!/bin/sh
mkdir serial_build; cd serial_build
cmake \
      ..
make -j 8

./ShallowWater
./ShallowWater_par1
./ShallowWater_par2
./ShallowWater_par3
./ShallowWater_par4
