# Chapter 13 GPU profiling and tools
This is from Chapter 13 of Parallel and High Performance Computing, Robey and Zamora,
Manning Publications, available at http://manning.com

The book may be obtained at
   http://www.manning.com/?a_aid=ParallelComputingRobey

Copyright 2019-2020 Robert Robey, Yuliana Zamora, and Manning Publications
Emails: brobey@earthlink.net, yzamora215@gmail.com

See License.txt for licensing information.

# Shallow Water example
   Directory OpenACC/ShallowWater
   Build with cmake
      mkdir build && cd build
      cmake -DENABLE_GRAPHICS=1 ..
      make
   Run the code with
      ./ShallowWater

# Profile the Shallow Water example (Book: listing 13.1 - 13.3)
   Directory OpenACC/ShallowWater_profiled
   Optimizations are implemented in ShallowWater_par1.c to _par4.c
   Build with cmake
      mkdir build && cd build
      cmake -DENABLE_GRAPHICS=1 ..
      make
   Run the code with
      ./ShallowWater
   Post-processing profile output
      gprof -l -pg ./ShallowWater

# Virtual Box
  
