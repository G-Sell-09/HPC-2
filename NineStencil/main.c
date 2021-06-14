#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "alu.h"
#include "stencils.h"
#include "functions.h"
#include "jacobi.h"
#include "gaussseidel.h"
#include "sor.h"
#include "smoothing.h"
#include "prol.h"
#include "restr.h"
#include "ludecomp.h"
#include "v_cycle.h"
#include "w_cycle.h"

int main (int argc , char **argv)
{

  if( argc < 3 )
  {
    printf("Not enough input parameters.\n");
    return 0;
  }

  int N = atoi(argv[1]);
  int L = atoi(argv[2]); // Give 0 if L shall be calculated
  char *GL = argv[3]; // Options: jac gau sor

  //int tolF = atoi(argv[4]); //


  int maxIter = 1000;
  int step = 0;
  int i,j, nu1, nu2;

  // Smoothing parameters
  nu1 = 20;
  nu2 = 20;
  int mu = 0;

  // tolerance for ending multi-grid
  double tol = 1.e-8;

  // if L is given as 0, calc fitting L for given grid size
  if (L == 0)
  {
    calculateL(N,L);
  }
  // step size
  double h  = 1./(double)(N+1);

  double nrm_x, r_norm, R, S, w;
  double *e_iter = malloc(maxIter*sizeof(double));

  // Function handle for smoothing
  myFuncDef smooth = functionHandle(GL, &w);

  // Allocate N_vec, N_pot

  int *N_vec = malloc(L*sizeof(int));
  int *N_pot = malloc(L*sizeof(int));

  // Calc Gridsize for every level
  int checkN = initN(N_vec, N_pot, L, N);
  if (checkN != 0)
    return 1;

  // Allocate x, b and r

  double **ArrOfx = (double**) malloc (L*sizeof(double*));
  double **ArrOfb = (double**) malloc (L*sizeof(double*));
  double **ArrOfr = (double**) malloc (L*sizeof(double*));

  for (int i = 0; i < L; i++)
  {
      ArrOfx[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfb[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfr[i] = (double*) malloc (N_pot[i]*sizeof(double));
  }

  // Allocate A, L and U

  double **A = (double**) malloc (N_pot[L-1]*sizeof(double*));
  double **lower = (double**) malloc (N_pot[L-1]*sizeof(double*));
  double **upper = (double**) malloc (N_pot[L-1]*sizeof(double*));

  for (int i = 0; i < N_pot[L-1]; i++)
  {
      A[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
      lower[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
      upper[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
  }

  double StartALU = getTimeStamp();
  // Initiate L and U for roughest grid size
  initLU9(A, lower, upper, N_vec[L-1], N_pot[L-1]);
  //initLU9(A, lower, upper, N_vec[L-1], N_pot[L-1]);
  double EndeALU = getTimeStamp();

  // Initiate inital values and right side
  initXB(ArrOfx[0], ArrOfb[0], N, h);


  double Q = getTimeStamp();
  for (i=1;i<=maxIter;i++)
  {

    step += 1;

    //v_cycle(N, ArrOfx, ArrOfb, ArrOfr, L, nu1, nu2, w, step, N_vec, N_pot, lower, upper, smooth);
    w_cycle(N, ArrOfx, ArrOfb, ArrOfr, L, mu, nu1, nu2, w, step, N_vec, N_pot, lower, upper, smooth);

    // Calc residual norm
    e_iter[i] = norm2(ArrOfr[0],N*N);
    S = getTimeStamp();


    printf("\n\n Fehler %e, Time: %f \n", e_iter[i], S-Q);
    printf("-----------------------------------------------------\n");

    // terminator
    if (e_iter[i]<tol)
      break;

  }
  printf("Time for ALU: %f\n", EndeALU-StartALU);

  // for (int i =0;i<N*N;i++)
  //   printf("x[%d]: %f\n", i, ArrOfx[0][i]);

  return 0;
}
