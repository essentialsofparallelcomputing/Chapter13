#include <stdio.h>
#include "timer.h"
#include <Kokkos_Core.hpp>

int main(int argc, char *argv[]){
   struct timespec tstart;
   // initializing data and arrays
   int nsize=80000000; int ntimes=16;
   double scalar = 3.0, time_sum = 0.0;

   Kokkos::initialize(argc, argv); {

   // initializing data and arrays
   Kokkos::View<double *> a( "a", nsize);
   Kokkos::View<double *> b( "b", nsize);
   Kokkos::View<double *> c( "c", nsize);

   Kokkos::parallel_for(nsize, KOKKOS_LAMBDA (int i) {
      a[i] = 1.0;
      b[i] = 2.0;
   });

   for (int k=0; k<ntimes; k++){
      cpu_timer_start(&tstart);
      // stream triad loop
      Kokkos::parallel_for(nsize, KOKKOS_LAMBDA (int i) {
         c[i] = a[i] + scalar*b[i];
      });
      time_sum += cpu_timer_stop(tstart);
      // to keep the compiler from optimizing out the loop
      c[1]=c[2];
   }

   printf("Average runtime is %lf msecs\n", time_sum/ntimes);

   } Kokkos::finalize();
}
