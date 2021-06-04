#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"


void v_cycle(N,b,x0,L,nu1,nu2)
{

  int i,j;
  // Parameter fuer Jacobi
  double w = 0.6;
  double tol = 0.0;

  // Array mit N und den Gittergroessen
  double *N_vec = malloc(L*sizeof(double));
  double *N_pot = malloc(L*sizeof(double));

  N_vec[0] = N;
  N_pot[0] = N*N;
  for (i=1; i<L; i++)
  {
    N_vec[i]=(N_vec[i-1]-1)/2;
    N_pot[i] = N_vec[i]*N_vec[i];
    if ((N_vec[i]%2)==0)
    {
      printf("N darf auf keinem Level gerade sein! Dem Programm fehlen die Fallunterscheidungen.\n");
    }
  }

  // Arrays fuer die Gitter
  // TODO b nur ein Array
  double *ArrOfx = (double**) malloc (L*sizeof(double*));
  double *ArrOfb = (double**) malloc (L*sizeof(double*));
  double *ArrOfr = (double**) malloc (L*sizeof(double*));

  for (int i = 0; i < L; i++)
  {
      ArrOfx[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfb[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfr[i] = (double*) malloc (N_pot[i]*sizeof(double));
  }

  // Belege erstes x und b
  duplicate(ArrOfx[0],x0,N_pot[0]);
  duplicate(ArrOfb[0],b,N_pot[0]);


  for (i=0;i<L-1;i++)
  {
    // Vorglaettung mit Jacobi
    jac(ArrOfx[i], ArrOfb[i], ArrOfr[i], w, nu1, tol, N_vec[i]);

    // Residuum berechnen (Jetzt schon in Jacobi)
    //residual(ArrOfr[i], ArrOfx[i], ArrOfb[i], N);

    //Restriktion
    restr(ArrOfr[i],ArrOfr[i+1], N_vec[i]);

    duplicate(ArrOfb[i+1], ArrOfr[i+1], N_pot[i]);

    zeroes(ArrOfx[i+1], N_pot[i+1]);

  }


  // Loesen mit LU
  lu_decomp(...)


  for (i=2;i<L+1;i++)
  {
    // Prolongieren und korrigieren
    prol(ArrOfx[L-i+1], ArrOfx[L-i], N_vec[L-i]);

    // Nachglaettung
    jac(ArrOfx[L-i], ArrOfb[L-i], ArrOfr[L-i], w, nu2, tol, N_vec[L-i]);
  }

  // Schreibe das Ergebnis in x0
  duplicate(x0,ArrOfx[0],N_pot[0]);

}
