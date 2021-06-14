#include <math.h>
#include <sched.h>
#include <time.h>
#include <string.h>

double min(double x, double y)
{
  if (x < y)
	 return x;
	else
	 return y;
}

/// Berechnet Vektor plus SKalar mal Vektor
void vec1(double * restrict a,
              double * restrict b,
              double * restrict c,
              double ska,
              int N
              )
{
	int i;
	#pragma omp parallel for shared(ska,a,b,N) private(i)
  for (int i=0;i<N;i++)
  {
    a[i] = b[i] + ska * c[i];
  }
}

/// Berechnet Skalar mal Vektor minus Vektor
void vec2(double * restrict a,
              double * restrict b,
              double ska,
              int N
              )
{
	int i;
	#pragma omp parallel for shared(ska,a,b,N) private(i)
  for (int i=0;i<N;i++)
  {
    a[i] += ska * b[i];
  }
}

/// Berechnet Skalar mal Vektor minus Vektor
void vec3(double * restrict a,
          double ska,
          int N
          )
{
	int i;
	#pragma omp parallel for shared(ska,a,b,N) private(i)
  for (int i=0;i<N;i++)
  {
    a[i] = ska * a[i];
  }
}

/// Berechnet das Skalarprodukt von zwei Vektoren
double dot(double * restrict x,
        double * restrict y,
        int N)
{
  double vproduct = 0.0; ///< Skalarprodukt

  /// Parallele Berechnung des Skalarprodukts und Ausgabe der jeweiligen Thread-ID
	int i;
  #pragma omp parallel for schedule(static) reduction(+:vproduct)
  for(int i=0;i<N;i++)
  {
    vproduct += x[i] * y[i];
  }
  return vproduct;
}

/// Berechnet die 2-Norm einen Vektors
double norm2(double * restrict x, int N)
{
  double wurzel = sqrt(dot(x,x,N));
  return wurzel;
}

/// Berechnet das geg. f in einem Gitterpunkt
double f(double x, double y)
{
  x = 2 * pow(M_PI, 2) * sin(M_PI*x) * sin(M_PI*y);
  return x;
}

/// Berechnet die Loesungsfunktion in  einem Gitterpunkt
double u(double x, double y)
{
  x = sin(M_PI*x) * sin(M_PI*y);
  return x;
}

double initial_val(double x, double y, int k, int l, double h)
{
	x = sin(k*M_PI*x*h) + sin(l*M_PI*y*h);
}

/// Gibt aktuellen Zeitstempel zurueck
double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}

// ToDo: Residuum für mfMult5 oder mfMult9 je nach Eingabewert?
void residual(double* res, double* x, double *b, int N)
{
  mfMult9(x, res, N);

  vec1(res, b, res, -1, N*N);
}

void duplicate(double *a, double *b, int N)
{
  for (int i=0;i<N;i++)
  {
    a[i] = b[i];
  }
}

void zeroes(double *a, int N)
{
  for (int i=0;i<N;i++)
  {
    a[i]=0.0;
  }
}

void calculateL(int N, int L)
{
	int N_help=N;
	for (int i = 0; i < 100; i++)
	{
		if (N_help > 50)
		{
			L += 1;
			N_help = (N-1)/2;
		}
		else
		break;
	}
}

int initN(int *N_vec,int *N_pot, int L, int N)
{
	N_vec[0] = N;
  N_pot[0] = N*N;
  for (int i=1; i<L; i++)
  {
    N_vec[i]=(N_vec[i-1]-1)/2;
    N_pot[i] = N_vec[i]*N_vec[i];
    if ((N_vec[i]%2)==0 && i!=L-1)
    {
      printf("N darf auf keinem Level gerade sein! Dem Programm fehlen die Fallunterscheidungen.\n");
    	return 1;
    }
  }
	return 0;
}

void initXB(double *x, double *b, int N, double h){
	// int k = 1;
	// int l = 1;
  for (int i = 0; i < N; i++  )
  {
    for (int j = 0; j < N ; j++ )
    {
      x[i*N + j] = 1.0; //initial_val(i+1, j+1, k, l, h); // start vector
      b[i*N + j] = h*h*f((i+1)*h, (j+1)*h); // right side
      //x_sol[i*N + j] = u((i+1)*h, (j+1)*h);
    }
  }
}
