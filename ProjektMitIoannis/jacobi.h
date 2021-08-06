/*!
 * \brief jacobi smoothing operator using the 5-point stencil
 *
 * @param[in] x       vector with current values
 * @param[in] b       right hand side
 * @param[in] r       residual
 * @param[in] w       relaxation factor
 * @param[in] nu      number of iterations
 * @param[in] tol     tolerance
 * @param[in] N       grid size
 * @param[in] mfMult  typdef myMult with correct stencil
 *
 * @param[out] x      updated values
 * @param[out] r      residual
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
 */
void jac5(double *x, double *b, double *r, double w, double nu, double tol, int N, myMult mfMult)
{
  // first redidual
  residual(r, x, b, N, mfMult);

  for (int i = 0;i<nu;i++)
  {
    // next x_k
    vec2(x, r, w/4, N*N);
    // next residiual
    residual(r, x, b, N, mfMult);

    if (norm2(r, N*N) / norm2(b,N*N) < tol)
      break;
  }
}

/*!
 * \brief jacobi smoothing operator using the 9-point stencil
 *
 * @param[in] x       vector
 * @param[in] b       right hand side
 * @param[in] r       residual
 * @param[in] w       relaxation factor
 * @param[in] nu      number of iterations
 * @param[in] tol     tolerance
 * @param[in] N       grid size
 * @param[in] mfMult  typdef myMult for correct stencil
 *
 * @param[out] x      updated values
 * @param[out] r      residual
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
 */
void jac9(double *x, double *b, double *r, double w, double nu, double tol, int N, myMult mfMult)
{
  // first residual
  residual(r, x, b, N, mfMult);

  for (int i = 0;i<nu;i++)
  {
    // next x_k
    vec2(x, r, w/20, N*N);

    // next residual
    residual(r, x, b, N, mfMult);

    if (norm2(r, N*N) / norm2(b,N*N) < tol)
      break;
  }
}
