#include <math.h>
#include <sched.h>
#include <time.h>
#include <string.h>

typedef void (*myFuncDef)(double*, double*, double* , double, double, double, int);


void dummy(double *x, double *b, double *r, double w, double nu, double tol, int N)
{
  residual(r, x, b, N);
}

myFuncDef functionHandle(char* GL, double *w) {

  myFuncDef smooth;
  if (strcmp(GL, "jac") == 0)
	{
		smooth = &jac;
    *w = 0.6;
	}
	else if (strcmp(GL, "gau") == 0)
	{
		smooth = &gau;
    *w = 1.0;
	}
	else if (strcmp(GL, "sor") == 0)
	{
		smooth = &sor;
    *w = 1.95;
	}
  else
  {
    smooth = &dummy;
    printf("Kein Glaetter ausgewaehlt\n");
  }
	return smooth;
}
