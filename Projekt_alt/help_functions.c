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
              double ska,
              int N
              )
{
	int i;
	#pragma omp parallel for shared(ska,a,b,N) private(i)
  for (int i=0;i<N*N;i++)
  {
    a[i] += ska * b[i];
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
  for (int i=0;i<N*N;i++)
  {
    a[i] = ska * a[i] - b[i];
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
  for(int i=0;i<N*N;i++)
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

/// Gibt aktuellen Zeitstempel zurueck
double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}
