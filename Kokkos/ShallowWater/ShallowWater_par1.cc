#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "timer.h"

#include <Kokkos_Core.hpp>

/*********************************************************************************************
 * WAVE -- 2D Shallow Water Equation Model
 *                            Bob Robey, Los Alamos National Laboratory
 * ******************************************************************************************/

//define macro for squaring a number
#define SQ(x) ((x)*(x))

int main(int argc, char *argv[])
{
  int      nx = 500, ny = 200;
  int      ntimes = 2000, nburst = 100;
  double   deltaX=1.0, deltaY=1.0;         //size of cell
  double   g = 9.80;                           // gravitational constant
  double   sigma = 0.95;
  double   time=0.0;                           //computer simulation time
  double   totaltime;   //variables to calculate time taken for the program to run
  struct   timespec starttime;
  double   TotalMass, origTM;    //variables for checking conservation of mass

  Kokkos::initialize(argc, argv);
  
  /* allocate the memory dynamically for the matrix */
  // state variables
  Kokkos::View<double **> H("H", ny+2, nx+2);
  Kokkos::View<double **> U("U", ny+2, nx+2);
  Kokkos::View<double **> V("V", ny+2, nx+2);

  Kokkos::View<double **> Hnew("Hnew", ny+2, nx+2);
  Kokkos::View<double **> Unew("Unew", ny+2, nx+2);
  Kokkos::View<double **> Vnew("Vnew", ny+2, nx+2);

  // half-step arrays
  Kokkos::View<double **> Hx("Hx", ny, nx+1);
  Kokkos::View<double **> Ux("Ux", ny, nx+1);
  Kokkos::View<double **> Vx("Vx", ny, nx+1);

  Kokkos::View<double **> Hy("Hy", ny+1, nx);
  Kokkos::View<double **> Uy("Uy", ny+1, nx);
  Kokkos::View<double **> Vy("Vy", ny+1, nx);

  /*initialize matrix*/
  
  for(int j=0;j<=ny+1;j++){
    for(int i=0;i<=nx+1;i++){
      H(j,i)=2.0;
      U(j,i)=0.0;
      V(j,i)=0.0;
    }
  }
  for(int j=0;j<=ny+1;j++){
    for(int i=0;i<=(nx+1)/2;i++){
      H(j,i)=10.0 - ((10.0 - 2.0)/ (double)((nx+1)/2))*(double)(i);
    }
  }

  origTM=0.0;
  for(int j=1;j<=ny;j++){
    for(int i=1;i<=nx;i++){
        origTM+=H(j,i);
    }
  }
  
  //print iteration info
  double deltaT = 1.0e30;
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      double wavespeed = sqrt(g*H(j,i));
      double xspeed = (fabs(U(j,i))+wavespeed)/deltaX;
      double yspeed = (fabs(V(j,i))+wavespeed)/deltaY;
      double my_deltaT = sigma/(xspeed+yspeed);
      if (my_deltaT < deltaT) deltaT = my_deltaT;
    }
  }
  printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", 0, time, deltaT, origTM);

  cpu_timer_start(&starttime);

  /* run the simulation for given number of iterations */
  for (int n = 0; n < ntimes; ) {

    for (int ib=0; ib<nburst; ib++){

      //set boundary conditons
      for(int j=1;j<=ny;j++){
        H(j,0)=H(j,1);
        U(j,0)=-U(j,1);
        V(j,0)=V(j,1);
        H(j,nx+1)=H(j,nx);
        U(j,nx+1)=-U(j,nx);
        V(j,nx+1)=V(j,nx);
      }
      for(int i=0;i<=nx+1;i++){
        H(0,i)=H(1,i);
        U(0,i)=U(1,i);
        V(0,i)=-V(1,i);
        H(ny+1,i)=H(ny,i);
        U(ny+1,i)=U(ny,i);
        V(ny+1,i)=-V(ny,i);
      }
    
      //set timestep
      deltaT = 1.0e30;
      for (int j = 1; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
          double wavespeed = sqrt(g*H(j,i));
          double xspeed = (fabs(U(j,i))+wavespeed)/deltaX;
          double yspeed = (fabs(V(j,i))+wavespeed)/deltaY;
          double my_deltaT = sigma/(xspeed+yspeed);
          if (my_deltaT < deltaT) deltaT = my_deltaT;
        }
      }

      //first pass
      //x direction
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i<=nx;i++) {
          //density calculation
          Hx(j,i)=0.5*(H(j+1,i+1)+H(j+1,i  )) - deltaT/(2.0*deltaX)*
                              (U(j+1,i+1)-U(j+1,i  ));
          //momentum x calculation
          Ux(j,i)=0.5*(U(j+1,i+1)+U(j+1,i  )) - deltaT/(2.0*deltaX)*
                             ((SQ(U(j+1,i+1))/H(j+1,i+1) + 0.5*g*SQ(H(j+1,i+1))) -
                              (SQ(U(j+1,i  ))/H(j+1,i  ) + 0.5*g*SQ(H(j+1,i  ))));
          //momentum y calculation
          Vx(j,i)=0.5*(V(j+1,i+1)+V(j+1,i  )) - deltaT/(2.0*deltaX)*
                             ((U(j+1,i+1)*V(j+1,i+1)/H(j+1,i+1)) -
                              (U(j+1,i  )*V(j+1,i  )/H(j+1,i  )));
        }
      }
    
      //y direction
      for(int j = 0; j<=ny; j++){
        for(int i=0; i<nx; i++){
          //density calculation
          Hy(j,i)=0.5*(H(j+1,i+1)+H(j  ,i+1)) - deltaT/(2.0*deltaY)*
                              (V(j+1,i+1)-V(j  ,i+1));
          //momentum x calculation
          Uy(j,i)=0.5*(U(j+1,i+1)+U(j  ,i+1)) - deltaT/(2.0*deltaY)*
                             ((V(j+1,i+1)*U(j+1,i+1)/H(j+1,i+1)) -
                              (V(j  ,i+1)*U(j  ,i+1)/H(j  ,i+1)));
          //momentum y calculation
          Vy(j,i)=0.5*(V(j+1,i+1)+V(j  ,i+1)) - deltaT/(2.0*deltaY)*
                             ((SQ(V(j+1,i+1))/H(j+1,i+1) + 0.5*g*SQ(H(j+1,i+1))) -
                              (SQ(V(j  ,i+1))/H(j  ,i+1) + 0.5*g*SQ(H(j  ,i+1))));
        }
      }

      //second pass
      for (int j = 1; j <=ny; j++) {
        for (int i = 1; i<=nx;i++) {
          //density calculation
          Hnew(j,i) = H(j,i) - (deltaT/deltaX)*(Ux(j-1,i  )-Ux(j-1,i-1))
                                     - (deltaT/deltaY)*(Vy(j  ,i-1)-Vy(j-1,i-1));
          //momentum x calculation
          Unew(j,i) = U(j,i) - (deltaT/deltaX)*
                                  ((SQ(Ux(j-1,i  ))/Hx(j-1,i  ) +0.5*g*SQ(Hx(j-1,i  ))) -
                                   (SQ(Ux(j-1,i-1))/Hx(j-1,i-1) +0.5*g*SQ(Hx(j-1,i-1))))
                               - (deltaT/deltaY)*
                                  ((Vy(j  ,i-1)*Uy(j  ,i-1)/Hy(j  ,i-1)) -
                                   (Vy(j-1,i-1)*Uy(j-1,i-1)/Hy(j-1,i-1)));
          //momentum y calculation
          Vnew(j,i) = V(j,i) - (deltaT/deltaX)*
                                  ((Ux(j-1,i  )*Vx(j-1,i  )/Hx(j-1,i  )) -
                                   (Ux(j-1,i-1)*Vx(j-1,i-1)/Hx(j-1,i-1)))
                               - (deltaT/deltaY)*
                                  ((SQ(Vy(j  ,i-1))/Hy(j  ,i-1) +0.5*g*SQ(Hy(j  ,i-1))) -
                                   (SQ(Vy(j-1,i-1))/Hy(j-1,i-1) +0.5*g*SQ(Hy(j-1,i-1))));
        }
      }
    
      // Need to replace swap with copy
      for(int j=1;j<=ny;j++){
        for(int i=1;i<=nx;i++){
           H(j,i) = Hnew(j,i);
           U(j,i) = Unew(j,i);
           V(j,i) = Vnew(j,i);
        }   
      }   

    } // burst loop
      
    TotalMass=0.0;
    for(int j=1;j<=ny;j++){
      for(int i=1;i<=nx;i++){
          TotalMass+=H(j,i);
      }
    }

    if(((fabs(TotalMass-origTM)>1.0E-6)||isnan(TotalMass))){
       printf("Conservation of mass\nMass difference:%e\n", TotalMass-origTM);
       printf("Problem occured on iteration %5.5d at time %f.\n", n, time);
       //exit(0);
    }
    time+=deltaT;
    n+=nburst;

    //print iteration info
    printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", n, time, deltaT, TotalMass);

  }  // End of iteration loop
  
  /* Compute the average time taken/processor */
  totaltime = cpu_timer_stop(starttime);
  
  /* Print the total time taken */
  printf(" Flow finished in %lf seconds\n",totaltime);//*/

  Kokkos::finalize();

  exit(0);
}
