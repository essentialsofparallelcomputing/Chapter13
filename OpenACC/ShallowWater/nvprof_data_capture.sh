#!/bin/sh

set -vex

module load pgi/19.7 cuda/9.1 cmake/3.12.4

rm -rf prof
mkdir prof
for prog in ShallowWater.c ShallowWater_par1.c ShallowWater_par2.c ShallowWater_par3.c ShallowWater_par4.c
do
   sed -e '1,$s/ntimes = 2000/ntimes = 4/' ${prog} > prof/${prog}
done
cd prof
cp ../ShallowWater.h ../CMakeLists.txt ../memory.c ../timer.h ../timer.c .
cmake .
make

./ShallowWater |& tee ShallowWater.out
for prog in ShallowWater_par1 ShallowWater_par2 ShallowWater_par3 ShallowWater_par4
do
#  nvprof -o ${prog}.summary ./${prog} |& tee ${prog}_summary.out
   nvprof --export-profile ${prog}_timeline.prof ./${prog} |& tee ${prog}_timeline.out
done

tar -cvf ShallowProf.tgz *.out *.prof # *.summary
#rm -f ShallowWater_par1_timeline.prof
#tar -cvf ShallowProfSmall.tgz *.out *.prof # *.summary

mv ShallowProf.tgz ShallowProf_v100_9.1.tgz
mv ShallowProfSmall.tgz ShallowProfSmall_v100_9.1.tgz
