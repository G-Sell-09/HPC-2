
/*!
 * \brief Performs a iterative multigrid method with pre and after smoothing.
 *
 * @param[in] N        dimension
 * @param[in] ArrOfx   matrix with enough space to store x in all levels
 * @param[in] ArrOfb   matrix with enough space to store right hand sides in all levels
 * @param[in] ArrOfr   matrix with enough space to store residuals in all levels
 * @param[in] L        # of levels for multigrid
 * @param[in] L_temp   current Level for rekursion
 * @param[in] nu1      number of iteration in the pre-smoother
 * @param[in] nu2      number of iteration in the post-smoother
 * @param[in] w        relaxation factor for smoother
 * @param[in] tolSmoo  tolerance for the smoother
 * @param[in] N_vec    array with values of N for each level
 * @param[in] N_vec    array with values of N squared for each level
 * @param[in] lower    matrix L for the LU decomposition
 * @param[in] upper    matrix U for the LU decomposition
 * @param[in] smooth   typedef that points to the chosen smoothing operator
 * @param[in] mfMult   typedef that points to the chosen matrixfree mat-vec multiplication
 *
 * @param[in] ArrOfx   matrix with enough space to store x in all levels
 * @param[in] ArrOfr   matrix with enough space to store residuals in all levels
 */
void v_cycle_iter(int N,double **ArrOfx, double **ArrOfb, double **ArrOfr, int L,int mu, double nu1, double nu2, double w, double tolSmoo, int *N_vec,int *N_pot, \
  double **lower, double **upper, mySmoother smooth, myMult mfMult)
{

  int i,j;

  for (i=0;i<L-1;i++)
  {
    // Pre-Smoothing with Jacobi, Gauss or SOR
    smooth(ArrOfx[i], ArrOfb[i], ArrOfr[i], w, nu1, tolSmoo, N_vec[i], mfMult);

    //restriction
    restr(ArrOfr[i], ArrOfr[i+1], N_vec[i]);

    duplicate(ArrOfb[i+1], ArrOfr[i+1], N_pot[i+1]);

    zeroes(ArrOfx[i+1], N_pot[i+1]);

  }

  // lu decomposition for the roughest grid
  lu_Decomp(ArrOfx[L-1],ArrOfb[L-1],N_pot[L-1], lower, upper);

  for (i=2;i<L+1;i++)
  {
    // prolong and correct
    prol(ArrOfx[L-i+1], ArrOfx[L-i], N_vec[L-i+1]);

    // Post-Smoothing with Jacobi, Gauss or SOR
    smooth(ArrOfx[L-i], ArrOfb[L-i], ArrOfr[L-i], w, nu2, tolSmoo, N_vec[L-i], mfMult);
  }

  //printf("Times in step %d: for1: %f, LU: %f, for2: %f\n", step, T-S, U-T, V-U);
}
