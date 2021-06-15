/*!
 * \brief Performs a recursive multigrid method with pre- and post-smoothing using the v_cycle procedure.
 *
 * @param[in] N         grid size
 * @param[in] ArrOfx    matrix with enough space to store x on all levels
 * @param[in] ArrOfb    matrix with enough space to store right hand sides on all levels
 * @param[in] ArrOfr    matrix with enough space to store residuals on all levels
 * @param[in] L         # of levels for multigrid method
 * @param[in] L_temp    current level for rekursive call
 * @param[in] nu1       number of iteration in the pre-smoothing operator
 * @param[in] nu2       number of iteration in the post-smoothing operator
 * @param[in] w         relaxation factor for smoothing operator
 * @param[in] tolSmoo   tolerance for the smoothing operator
 * @param[in] N_vec     array with values of N on each level
 * @param[in] N_pot     array with values of N squared on each level
 * @param[in] lower     matrix L for the LU decomposition
 * @param[in] upper     matrix U for the LU decomposition
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
void v_cycle(double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int L_temp, double nu1, double nu2, int mu, double w, double tolSmoo, int *N_vec,int *N_pot, \
  double **lower, double **upper, mySmoother smooth, myMult mfMult)
{
  // solve the system Ax = r on the roughest grid with the prior calculated LU decomposition
  if (L_temp == 1){
    lu_Decomp(ArrOfx[L-1],ArrOfb[L-1],N_pot[L-1], lower, upper);
  }
  else
  {

    // pre-smoothing operator
    smooth(ArrOfx[L-L_temp], ArrOfb[L-L_temp], ArrOfr[L-L_temp], w, nu1, tolSmoo, N_vec[L-L_temp], mfMult);

    // restriction to a rougher grid
    restr(ArrOfr[L-L_temp], ArrOfr[L-L_temp+1], N_vec[L-L_temp]);

    duplicate(ArrOfb[L-L_temp+1], ArrOfr[L-L_temp+1], N_pot[L-L_temp+1]);
    zeroes(ArrOfx[L-L_temp+1], N_pot[L-L_temp+1]);

    // recursive call to the next lower level
    v_cycle(ArrOfx, ArrOfb, ArrOfr, L, L_temp-1, nu1, nu2, mu, w, tolSmoo, N_vec, N_pot, lower, upper, smooth, mfMult);

    // prolongate and correct
    prol(ArrOfx[L-L_temp+1], ArrOfx[L-L_temp], N_vec[L-L_temp+1]);

    // post-smoothing operator
    smooth(ArrOfx[L-L_temp], ArrOfb[L-L_temp], ArrOfr[L-L_temp], w, nu2, tolSmoo, N_vec[L-L_temp], mfMult);

  }
}
