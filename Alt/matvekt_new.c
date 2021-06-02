#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <likwid.h>



double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}


int main(int argc, char const *argv[])
{
  if( argc < 2 )
	{
		printf("Nicht genuegend Eingabeparameter.\n");
		return 0;
	}

  int N = atoi(argv[1]);
  int count = N*N;
  //double *a = malloc((N)*sizeof(double));
  //double *x = malloc((N)*sizeof(double));
  //const double **b = malloc((N)*(N)*sizeof(double));
  double **M;
  double *x;
  double *b;

  M  = malloc(N*sizeof(double));
  x =  malloc(N*sizeof(double));
  b =  malloc(N*sizeof(double));

  for (int i = 0; i < N; i++)
  {
    M[i] = (double*) malloc (N*sizeof(double));
  }


  LIKWID_MARKER_INIT;


  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      M[i][j] = (double)i / (double)j;
    }
    x[i] = (double)i;
  }

  double A,B;
  A = getTimeStamp();
  printf("Hallo \n");

  #pragma omp parallel
  {
  LIKWID_MARKER_START("MatVec");
  #pragma omp for
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      b[i] += M[i][j] * x[j];
    }
  }
  LIKWID_MARKER_STOP("MatVec");
  }

  B = getTimeStamp();


  printf("\n Rechenzeit Var 1: %lf \n", B-A );

  LIKWID_MARKER_CLOSE;

  return 0;
}
