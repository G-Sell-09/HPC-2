#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>


double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}

//Variante1
void var1(int N)
{
  double *a, *x;
  double **b;

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      printf("Fehler1\n");
      printf("a-Wert %f\n", &a);
      a[i]=a[i]+b[i][j]*x[j];
      printf("Fehler2\n");
    }
  }
}

//Variante2

// void var2(int N)
// {
//   double *a, *b_val, *b_col, *b_row_ptr, *x;
//   for(int i=0;i<N;i++)
//   {
//     for(int j=b_row_ptr[i]; j< b_row_ptr[i+1]; j++)
//     {
//       a[i] = a[i] + b_val[j] * x[b_col[j]];
//     }
//   }
// }

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
  //double **b = malloc((N)*(N)*sizeof(double));
  double *array_1 = malloc((N)*sizeof(double));
  double wert;

  for (int i = 0; i < count; i++)
  {
    wert = rand()%count;
    array_1[i] = wert;
    printf("Wert: %f\n", array_1[i]);
  }
  printf("Hallo \n");
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      double b[N][N];
      b[i][j] = array_1[i+j];
      printf("b-Wert %f\n", b[i][j]);
    }
  }
  for (int i = 0; i < N; i++)
  {
    double a[N];
    double x[N];
    a[i] = array_1[i];
    x[i] = array_1[i]/2*3*array_1[i+1];
    printf("a-Wert %f\n", a[i]);
    printf("x-Wert %f\n", x[i]);
  }
  double A,B,C,D;
  A = getTimeStamp();
  var1(N);
  B = getTimeStamp();
  //var2(N);
  D = getTimeStamp();

  printf("\n Rechenzeit Var 1: %lf \n", B-A );
  printf("\n Rechenzeit Var 2:%lf \n", C-B );

  return 0;
}
