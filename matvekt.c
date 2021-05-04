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
void var1(int a, int b, int x, int N) {
  double ∗a, ∗x;
  double ∗∗b;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      a[i]=a[i]+b[i][j]∗x[j];
    }
  }
}

//Variante2
void var2(int a, int b, int x, int N) {
  double ∗a, ∗b_val, ∗b_col, ∗b_row_ptr, ∗x;
  for(int i=0;i<N;i++){
    for(int j = b_row_ptr[i],j < b_row_ptr[i+1], j++){
      a[i] = a[i] + b_val[j] ∗ x[b_col[j]];
    }
  }
}

int main(int argc, char const *argv[]) {
  if( argc < 2 )
	{
		printf("Nicht genuegend Eingabeparameter.\n");
		return 0;
	}

  int N = atoi(argv[1]);
  int count = N*N;

  for (i = 0; i < count; i++);
  {
    double wert = rand()%count;
    double array_1 = wert;
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int a[i][j] = wert[i*j+j];
    }
  }
  for (int i = 0; i < N; i++) {
    int b[i] = wert[i];
    int x[i] = wert[i]/2*3*wert[i+1];
}
  double A,B,C,D;
  A = getTimeStamp();
  var1(a,b,N);
  B = getTimeStamp();
  C = getTimeStamp();
  var2(a,b,N);
  D = getTimeStamp();
  return 0;
}
