#define _GNU_SOURCE
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
#include "w_cycle.h"
#include "functionHandles.h"

void ILU0_Iter( int N, double* LL, double* UU);
void mfMult_alt(int N, double* r, double* y);
double dot_alt(int N, double* x, double* y);
void substitionen(int N, double* x, double* z, double* r, double* LL, double* UU);
void pcg(int N, double* b, double* x, double* LL, double* UU, double eps, double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int L_temp, double nu1, double nu2, double w, double tolSmoo, int *N_vec,int *N_pot, \
  double **lower, double **upper, mySmoother smooth, myMult mfMult, myCycle cycle, int verfahren);

int main (int argc , char **argv)
{

  // Anfang unsere main
  int N = atoi(argv[1]);         // gridsize
  int L = atoi(argv[4]);         // # of Levels of the multigrid method
  char *GL = argv[5];            // smoothing operator
  int stencil = atoi(argv[6]);   // stencil
  char *cyc = argv[7];           // cycle procedure
  int verfahren = atoi(argv[8]); // 0 for ILU0 1 for multigrid

  if(inputs(GL, stencil, cyc))
  {
    return 1;
  }


  int i,j, nu1, nu2;
  double nrm_x, r_norm, R, Q, S, T, w;

  /// number of pre and post smoothings
  nu1 = 50;
  nu2 = 50;

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
  myInitLU initLU = initLUHandle(stencil);
  myMult mfMult = multHandle(stencil);
  myInitXB initXB = initXBHandle(stencil);
  mySmoother smooth = smoothHandle(GL, &w, stencil);
  myCycle cycle = cycleHandle(cyc);

  // if L is given as 0, calc fitting L for given grid size
  if (L == 0)
    L = calculateL(N,L);

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



  // Ende unsere main
	//int N = atoi(argv[1]);

	//double h  = 1./(double)(N+1);
	double hh = h*h;

	double *x = malloc(N*N*sizeof(double));
	double *b = malloc(N*N*sizeof(double));

	double *LL = malloc(2*N*N*sizeof(double));
	double *UU = malloc(3*N*N*sizeof(double));

	// f(x,y) = 2 * pi^2 * sin(pi * x) sin(pi * y)
	for ( int i = 0; i < N; i++  )
		for ( int j = 0; j < N ; j++ )
		{
			x[i*N + j] = 1; // start vector
			b[i*N + j] = hh*2*M_PI*M_PI * sin( M_PI * (j+1)*h ) * sin( M_PI * (i+1)*h  ); // b_{iN+j} = hh*f_ji = hh*f(jh,ih)
		}

	double eps = 1;
	// for (int i = 0; i < atoi(argv[2]); i++ )
	// 	eps = eps/10.0;
  eps = tol;
	ILU0_Iter(N, LL, UU);
	double time = getTimeStamp();

	pcg(N, b, x, LL, UU, eps, ArrOfx, ArrOfb, ArrOfr, L, L, nu1, nu2, w, tolSmoo, N_vec, N_pot, lower, upper, smooth, mfMult, cycle, verfahren);
	time = getTimeStamp() - time;
	//printf("Runtime: %f\n", time );

	if( atoi(argv[3]) )
	{
		printf("approx.  &&  exact solution\n");
		for ( int i = 0; i < N; i++  )
			for ( int j = 0; j < N ; j++  )
			{
				printf("%lf     ", x[i*N +j] );
				printf("%lf\n", sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h  ) ); //u_ij = u(ih, jh)
			}
	}

	free (x);
	free (b);

	return 0;

}

void ILU0_Iter( int N, double* LL, double* UU){

#pragma omp parallel
{
    // set startvalue L+U = A
#pragma omp for schedule(static)
    for (int i = 0; i < N*N; i++)
    {
        LL[      + i] = -1;
        LL[  N*N + i] = -1;
        UU[        i] =  4;
        UU[  N*N + i] = -1;
        UU[2*N*N + i] = -1;
    }
#pragma omp for nowait schedule(static)
    for (int i = 0; i < N; i++)
    {
        LL[i] = 0;
        UU[2*N*N + N*N-i-1] = 0;
    }

#pragma omp for nowait schedule(static)
    for (int i = 0; i < N*N; i+=N)
    {
        LL[N*N + i    ] = 0;
        UU[N*N + i+N-1] = 0;
    }

}


	int ij;
	for (int k = 0; k < 1000; k++)
    {
#pragma omp parallel
{
        // bottom left corner, no update needed

		// lower row
#pragma omp for nowait schedule(static)
		for ( int i = 1; i < N; i++)
		{
			UU[        i]  = 4 - LL[N*N + i]*UU[N*N + i-1];
			LL[  N*N + i] = -1 / UU[i-1];
		}

		// left column
#pragma omp for nowait schedule(static)
		for (int j = 1; j < N; j++)
		{
			UU[      j*N] =  4 - LL[j*N] * UU[2*N*N + (j-1)*N];
			LL[      j*N] = -1 / UU[(j-1)*N];
		}

		// core
#pragma for collapse(2)
		for ( int i = 1; i < N; i++)
		{
			for ( int j = 1; j < N; j++)
			{
				ij = i*N + j;
				// LL(i,1) = (-1)/UU(i-N,1)
				LL[        ij] = -1 / UU[ij-N];
				// LL(i,2) = (-1)/UU(i-1,1)
				LL[  N*N + ij] = -1 / UU[ij-1];
				// UU(i,1) = 4-LL(i,1)*UU(i-N,3)-LL(i,2)*UU(i-1,2)
				UU[        ij] =  4 - LL[ij] * UU[2*N*N + ij-N] - LL[N*N + ij] * UU[N*N + ij-1];
				// UU(i,2) = -1 no update needed
				// UU(i,3) = -1 no update needed
			}
		}

}

    }

}

void mfMult_alt(int N, double* r, double* y)
{
	// A * r = y
#pragma omp parallel
{
	// core
#pragma omp for nowait schedule(static)
	for ( int i = 1; i < N-1; i++ )
		for ( int j = 1; j < N-1; j++ )
			y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

	// margins
#pragma omp for nowait schedule(static)
	for ( int i = 1; i < N-1; i++ )
	{
		y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // lower row
		y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // upper row
		y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // left column
		y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // right column
	}

	// corners
#pragma omp single nowait
{
	y[0]		 = 4*r[0]                 - r[N]             - r[1];  		        // bottom left
	y[N-1]		 = 4*r[N-1]               - r[2*N-1]         - r[N - 2];		    // bottom right
	y[(N-1)*N] 	 = 4*r[(N-1)*N]           - r[(N-1)*N + 1]   - r[(N-2)*N];	    	// top left
	y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];	// top right
}

}

}

double dot_alt(int N, double* x, double* y)
{
	double result = 0;
#pragma omp parallel for schedule(static) reduction(+:result)
	for ( int i = 0; i < N*N; i++ )
		result += x[i]*y[i];
	return result;
}

void substitionen(int N, double* x, double* z, double* r, double* LL, double* UU)
{

    // LU z = r

    // L x = r
        x[0] = r[0];
    for (int i = 1; i < N; i++)
        x[i] = r[i] - LL[N*N + i] * x[i-1];
    for (int i = N; i < N*N; i++)
        x[i] = r[i] - LL[N*N + i] * x[i-1] - LL[i] * x[ i-N ];

    // U z = x
	z[N*N-1] = x[N*N-1] / UU[N*N-1];
    for (int i = N*N-2; i > N*(N-1)-1; i--)
        z[i] = ( x[i] - UU[N*N + i] * z[i+1] ) / UU[i];
    for (int i = N*(N-1)-1; i >= 0; i--)
        z[i] = ( x[i] - UU[N*N + i] * z[i+1] - UU[2*N*N + i]* z[i+N] ) / UU[i];

}

void pcg(int N, double* b, double* x, double* LL, double* UU, double eps, double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int L_temp, double nu1, double nu2, double w, double tolSmoo, int *N_vec,int *N_pot, \
  double **lower, double **upper, mySmoother smooth, myMult mfMult, myCycle cycle, int verfahren)
{

	int i;
	double alpha, beta, r_dot, z_dot;
	double *r = malloc(N*N*sizeof(double));
  double *z = malloc(N*N*sizeof(double));
	double *p = malloc(N*N*sizeof(double));
	double *m = malloc(N*N*sizeof(double));
	double* t = malloc(N*N*sizeof(double));

	mfMult_alt(N, x, r); // A*x0 = r0

  double e_mg;

duplicate(ArrOfb[0],b,N*N);

#pragma omp parallel for schedule(static)
	for ( i = 0; i < N*N ; i++ )
		r[i] = r[i] - b[i];	// r0 = r0 - b
    if(verfahren==0)
      substitionen(N, t, z, r, LL, UU); // z_0 = (LU)^(-1) * r_o
    else{
      duplicate(ArrOfx[0],z,N*N);
      duplicate(ArrOfr[0],r,N*N);
      v_cycle(ArrOfx, ArrOfb, ArrOfr, L, L, nu1, nu2, w, tolSmoo, N_vec, N_pot, lower, upper, smooth, mfMult);
      duplicate(z, ArrOfx[0],N*N);
      duplicate(r, ArrOfr[0],N*N);
    }
#pragma omp parallel for schedule(static)
	for ( i = 0; i < N*N ; i++)
		p[i] = -z[i];		// po = -r0

	double z0_norm = sqrt( dot_alt(N, z, z) ); // for breaking condition

	double rz_dot = dot_alt(N, r, z);
	int count = 0;
    for(int k = 0; k < 10000; k++)
	{

		// alpha_k = < r_k, r_k > / < p_k, m_k >
		mfMult_alt(N, p, m);
		r_dot = rz_dot;
		alpha = r_dot / dot_alt(N, p, m);

		// x_{k+1} = x_k + alpha_k * p_k
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++ )
			x[i] += alpha * p[i];

		// r_{k+1} = r_k + alpha_k * m_k
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++ )
			r[i] += alpha * m[i];

      if (verfahren==0)
      {
        // z_{k+1} = (LU)^(-1) * r_k <=> LU z = r
        substitionen(N, t, z, r, LL, UU);
      }
      else
      {
        duplicate(ArrOfx[0],z,N*N);
        duplicate(ArrOfr[0],r,N*N);
        cycle(ArrOfx, ArrOfb, ArrOfr, L, L, nu1, nu2, w, tolSmoo, N_vec, N_pot, lower, upper, smooth, mfMult);
        duplicate(z, ArrOfx[0],N*N);
        duplicate(r, ArrOfr[0],N*N);

      }
    if (sqrt( dot_alt(N, r, r) ) < eps ){
  			printf("Err: %.2e\n", sqrt( dot_alt(N, r, r)));
  			break;
		}

    rz_dot = dot_alt(N, r, z);
		// beta_{k+1} = < r_{k+1}, r_{k+1} > / < r_k, r_k  >
		beta = rz_dot / r_dot;

		// p_{k+1} = beta_{k+1} * p_k - r_{k=1}
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++)
			p[i] =  beta * p[i] - z[i];

		count ++;
	}

	free(r);
  free(z);
	free(p);
	free(m);
	free(t);
	printf("Iterations: %i\n", count);

}
