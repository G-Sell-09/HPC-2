#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "alu.h"
#include "stencils.h"
#include "functions.h"
#include "jacobi.h"
#include "sor.h"
#include "prol.h"
#include "restr.h"
#include "ludecomp.h"
#include "v_cycle.h"
#include "mu_cycle.h"
#include "functionHandles.h"
/*!
 * \brief Recursive multigrid method with pre- and post-smoothing.
 *
 *
 *
 * @param[in] N        grid size options: has to be an uneven number for all levels but the last
 * @param[in] L        multigrid level options: 0 gives an optimized value for level, 1 gives a direct solution with LU decomposition and inputs >1 set a static level
 * @param[in] GL       smoothing operator options: jac, gau or sor
 * @param[in] stencil  stencils options: 5 or 9
 * @param[in] cyc      cycle options: v or w
 * @param[in] mu       parameter for mu_cycle options: 1 for v_cycle, 2 for w_cycle
 *
 * @param[out] {0,1}   returns 0 if the method finished without errors, 1 otherwise
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
int main (int argc , char **argv)
{

  if( argc < 6 )
  {
    printf("Not enough input parameters.\n");
    return 1;
  }

  int N = atoi(argv[1]);         // gridsize
  int L = atoi(argv[2]);         // # of Levels of the multigrid method
  char *GL = argv[3];            // smoothing operator
  int stencil = atoi(argv[4]);   // stencil
  char *cyc = argv[5];           // cycle procedure
  int mu;                        // optional mu for mu_cycle

  // if L is given as 0, calc fitting L for given grid size
  if (L == 0)
  {
    L = calculateL(N,L);
  }

  if(inputs(L, N, GL, stencil, cyc, argc, argv, &mu))
  {
    return 1;
  }

  int i,j, nu1, nu2;
  double nrm_x, r_norm, R, Q, S, T, w;

  /// number of pre and post smoothings
  nu1 = 30;
  nu2 = 30;

  // tolerances and maximum iterations
  double tol = 1.e-8; // tolerance for the residual norm
  double tolSmoo = 0.0; // tolerance for smoothing
  int maxIter = 1000; // Maximum number of iterations

  // step size
  double h  = 1./(double)(N+1);

  // array with residual errors
  double *x_sol = malloc(N*N*sizeof(double));
  double *e_iter = malloc(maxIter*sizeof(double));
  double e_analytic;


  // function handle for smoothingreturn
  myInitLU initLU   = initLUHandle(stencil);
  myMult mfMult     = multHandle(stencil);
  myInitXB initXB   = initXBHandle(stencil);
  mySmoother smooth = smoothHandle(GL, &w, stencil);
  myCycle cycle     = cycleHandle(cyc, &mu);

  // allocate N_vec, N_pot
  int *N_vec = malloc(L*sizeof(int));
  int *N_pot = malloc(L*sizeof(int));

  // calc. gridsize for every level
  if (initN(N_vec, N_pot, L, N) != 0)
    return 1;

  // allocate x, b and r
  double **ArrOfx = (double**) malloc (L*sizeof(double*));
  double **ArrOfb = (double**) malloc (L*sizeof(double*));
  double **ArrOfr = (double**) malloc (L*sizeof(double*));

  for (int i = 0; i < L; i++)
  {
      ArrOfx[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfb[i] = (double*) malloc (N_pot[i]*sizeof(double));
      ArrOfr[i] = (double*) malloc (N_pot[i]*sizeof(double));
  }

  // allocate A, L and U
  double **A = (double**) malloc (N_pot[L-1]*sizeof(double*));
  double **lower = (double**) malloc (N_pot[L-1]*sizeof(double*));
  double **upper = (double**) malloc (N_pot[L-1]*sizeof(double*));

  for (int i = 0; i < N_pot[L-1]; i++)
  {
      A[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
      lower[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
      upper[i] = (double*) malloc (N_pot[L-1]*sizeof(double));
  }


  R = getTimeStamp();
  // initiate L and U for roughest grid size
  initLU(A, lower, upper, N_vec[L-1], N_pot[L-1]);

  Q = getTimeStamp();


  // initiate inital values and right side
  initXB(ArrOfx[0], x_sol, ArrOfb[0], N, h);

  S = getTimeStamp();


  for (i=1;i<=maxIter;i++)
  {

    cycle(ArrOfx, ArrOfb, ArrOfr, L, L, nu1, nu2, mu, w, tolSmoo, N_vec, N_pot, lower, upper, smooth, mfMult);

    // calculate residual norm
    e_iter[i] = norm2(ArrOfr[0],N*N);
    // tolerance for the residual norm
    if (e_iter[i]<tol)
      break;

  }

  // error to the analytical solution
  vec2(x_sol,ArrOfx[0],-1,N*N);
  // ||x-u|| / ||x||
  e_analytic = norm2(x_sol,N*N) / norm2(ArrOfx[0],N*N);

  T = getTimeStamp();

  printf("-----------------------------------------------\n");
  printf("N: %d     L: %d    S: %s    sten: %d   cyc: %s\n", N,L,GL,stencil,cyc);
  printf("-----------------------------------------------\n");
  printf("Iterations : %d\ne_residual      : %e\ne_analytic      : %e\n", i, e_iter[i], e_analytic);
  printf("-----------------------------------------------\n");
  printf("Time ALU   : %f Time multigrid : %f\nTime total : %f\n", Q-R, T-S, T-R);
  printf("-----------------------------------------------\n");
  printf("-----------------------------------------------\n");

  return 0;
}
