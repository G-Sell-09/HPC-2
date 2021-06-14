void jac(double *x, double *b, double *r, double w, double nu, double tol, int N)
{
  // first residual
  residual(r, x, b, N);

  for (int i = 0;i<nu;i++)
  {

    // next x_k
    vec2(x, r, w/20, N*N);

    // next residual
    residual(r, x, b, N);

    if (norm2(r, N*N) / norm2(b,N*N) < tol)
      break;
  }
}
