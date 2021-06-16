/*!
 * \brief Performs a recursive multigrid method with pre- and post-smoothing using the w_cycle procedure.
 *
 * @param[in] ArrOfx    matrix with enough space to store x on all levels
 * @param[in] ArrOfb    matrix with enough space to store right hand sides on all levels
 * @param[in] ArrOfr    matrix with enough space to store residuals on all levels
 * @param[in] L         # of levels for multigrid method
 * @param[in] L_temp    current level for rekursive call
 * @param[in] nu1       number of iteration in the pre-smoothing operator
 * @param[in] nu2       number of iteration in the post-smoothing operator
 * @param[in] mu        # of iterations in each smoothing operator
 * @param[in] w         relaxation factor for smoothing operator
 * @param[in] tolSmoo   tolerance for the smoothing operator
 * @param[in] N_vec     array with values of N on each level
 * @param[in] N_pot     array with values of N squared on each level
 * @param[in] A         matrices L and U for the LU decomposition
 * @param[in] smooth    typedef that points to the chosen smoothing operator
 * @param[in] mfMult    typedef that points to the chosen matrixfree mat-vec multiplication
 *
 * @param[out] ArrOfx   matrix with enough space to store x on all levels
 * @param[out] ArrOfr   matrix with enough space to store residuals on all levels
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 */
void mu_cycle(double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int L_temp, double nu1, double nu2, int mu, double w, double tolSmoo, int *N_vec,int *N_pot, \
  double **A, mySmoother smooth, myMult mfMult)
{


  // solve the system Ax = r on the roughest grid with the prior calculated LU decomposition
  // after solving prolongate form the roughest grid to the next thinner grid
  if (L_temp == 1){
    lu_Decomp(ArrOfx[L-1],ArrOfb[L-1],N_pot[L-1], A);
    prol(ArrOfx[L-1], ArrOfx[L-2], N_vec[L-1]);
    // post-smoothing operator
    smooth(ArrOfx[L-2], ArrOfb[L-2], ArrOfr[L-2], w, nu2, tolSmoo, N_vec[L-2], mfMult);
  }

  else
  {
    for (int i=1;i<=mu;i++)
    {

      // pre-smoothing operator
      smooth(ArrOfx[L-L_temp], ArrOfb[L-L_temp], ArrOfr[L-L_temp], w, nu1, tolSmoo, N_vec[L-L_temp], mfMult);

      // restriction to a rougher grid
      restr(ArrOfr[L-L_temp], ArrOfr[L-L_temp+1], N_vec[L-L_temp]);

      duplicate(ArrOfb[L-L_temp+1], ArrOfr[L-L_temp+1], N_pot[L-L_temp+1]);
      zeroes(ArrOfx[L-L_temp+1], N_pot[L-L_temp+1]);

      // recursive call to w_cycle with one lower level
      mu_cycle(ArrOfx, ArrOfb, ArrOfr, L, L_temp-1, nu1, nu2, mu, w, tolSmoo, N_vec, N_pot, A, smooth, mfMult);
    }
    // prolongate to a thinner grid if not already on the thinnest
    if (L != L_temp){
      prol(ArrOfx[L-L_temp], ArrOfx[L-L_temp-1], N_vec[L-L_temp]);

      // post-smoothing operator
      smooth(ArrOfx[L-L_temp-1], ArrOfb[L-L_temp-1], ArrOfr[L-L_temp-1], w, nu2, tolSmoo, N_vec[L-L_temp-1], mfMult);
    }

  }
}
