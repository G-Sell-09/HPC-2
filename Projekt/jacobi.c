#include "functions.h"
void jac(double *x, double *b, double *r, double w, double nu, double tol, int N)
{

  double h = 1/(N+1);

  //TODO ohne res_k?
  //res_k = (double*) malloc (N*N*sizeof(double));
  //ax = (double*) malloc (N*N*sizeof(double));

  // Erstes Residuum
  residual(r, x, b, N);

  for (int i = 0;i<nu;i++)
  {
    // Naechstes x_k
    vec2(x, r, w*pow(h,2)/2, N*N);

    // Naechstes Residuum
    residual(r, x, b, N)

    if (norm2(r) / norm2(b) < tol)
      break;

  }
  // Uebernehme letztes Residuum als aktuellles in v_cyclus
  //duplicate(r, res_k, N*N);
}
