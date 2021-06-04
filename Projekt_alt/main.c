#define _GNU_SOURCE
#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>

/// Berechnet die Loesungsfunktion in  einem Gitterpunkt
double u(double x, double y)
{
  x = sin(M_PI*x) * sin(M_PI*y);
  return x;
}

void duplicate(double* x, double * y, int N)
{
  for (int i = 0; i < N*N; i++  )
  {
    y[i] = x[i];
  }
}




int main (int argc , char **argv)
{

  if( argc < 3 )
  {
    printf("Not enough input parameters.\n");
    return 0;
  }

  int N = atoi(argv[1]);
  int L = atoi(argv[2]);

  int maxIter = 500;
  int step = 0;
  int i,j;

  double tol = 1.e-8;

  if (L == 0)
  {
    //Berechne L anhand Größe der Matrix
  }
  double h  = 1./(double)(N+1);

  double *x = malloc(N*N*sizeof(double));
	double *b = malloc(N*N*sizeof(double));
  double *x_sol = malloc(N*N*sizeof(double));
  double *e_iter = malloc(maxIter*sizeof(double));
  double *x0 = malloc(N*N*sizeof(double));

  for (i = 0; i < N; i++  )
  {
		for (j = 0; j < N ; j++ )
		{
			x[i*N + j] = 0; // start vector
			b[i*N + j] = 1; // right side
      x_sol[i*N + j] = u((i+1)*h, (j+1)*h);
    }
  }

  for (i=1;i<maxIter;i++){
    e_iter[i] = 0;
  }


  for (i=1;i<maxIter;i++)
  {

    step += 1;

    duplicate(x, x0);

    v_cycle(N,b,x,L,nu1,nu2);
    
    if (N<1000){
      e_iter[i] = norm2(vec2(x_sol,x,1,N)) / norm2(x);
    }
    else{
      e_iter[i] = norm2(vec2(x,x0,1,N)) / norm2(x);
    }

    if (e_iter[i]<tol){
      break
    }
  }







}
