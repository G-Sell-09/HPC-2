/*!
 * \brief Calculates the 5-point poisson matrix for the roughest grid and its LU-decomposition
 *
 * @param[in] A       poisson matrix
 * @param[in] N_vec   size of the roughest grid
 * @param[in] N_pot   size of the roughest poisson matrix
 *
 * @param[out] A      L + U for LU decomposition
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void initLU5(double **A, int N_vec, int N_pot){

  // Insert the right values into the correspondix diagonales of the matrix A
  int i,j,k;
  #pragma omp parallel
  {
    #pragma omp for nowait
    for (int i = 0; i < N_pot; i++){
      for (int j = 0; j < N_pot; j++){
        if (i == j){
          A[i][j] = 4;
        }

        if (i==j+1 || i==(j+N_vec)){
          A[i][j] = -1;
        }
        if (i==j-1 || i==(j-N_vec)){
          A[i][j] = -1;
        }
      }
    }
  }

  #pragma omp parallel
  {
    // Setting some values on the diagonales to 0 to get the correct Poisson-matrix
    #pragma omp for nowait
    for (int k = 0; k < N_vec-1; k++){
      // below the diagonal
      A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // above the diagonal
      A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
    }
  }

  // Calcualte the L and U matrix of the Poisson Matrix A
  for(int k=0;k<N_pot-1;k++)
  {
    for(int i=k+1;i<N_pot;i++)
    {
      A[i][k] = A[i][k]/A[k][k];

      for(int j=k+1;j<N_pot;j++)
      {
        A[i][j] -= A[i][k]*A[k][j];
      }
    }
  }

}



//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------end of 5-point stencil--------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

/*!
 * \brief Calculates the 9-point poisson matrix for the roughest grid and its LU-decomposition
 *
 * @param[in] A       poisson matrix
 * @param[in] N_vec   size of the roughest grid
 * @param[in] N_pot   size of the roughest poisson matrix
 *
 * @param[out] A      L + U for LU decomposition
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void initLU9(double **A, int N_vec, int N_pot){

  // Insert the right values into the correspondix diagonales of the matrix A
  int i,j,k;
  #pragma omp parallel
  {
    #pragma omp for nowait
    for (int i = 0; i < N_pot; i++){
      for (int j = 0; j < N_pot; j++){
        if (i == j){
          A[i][j] = 20;
        }

        if (i==j+1 || i==(j+N_vec)){
          A[i][j] = -4;
        }
        if (i==(j+N_vec-1) || i==(j+N_vec+1)){
          A[i][j] = -1;
        }
        if (i==j-1 || i==(j-N_vec)){
          A[i][j] = -4;
        }
  			if (i==(j-N_vec-1) || i==(j-N_vec+1)){
  				A[i][j] = -1;
  			}
      }
    }
  }

  #pragma omp parallel
  {
    // Setting some values on the diagonales to 0 to get the correct Poisson-matrix
    #pragma omp for nowait
    for (int k = 0; k < N_vec-1; k++){
      // below the diagonal
      A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // above the diagonal
      A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
    }

    // Upper and lower blocks
    #pragma omp for nowait
		for (int k = 0; k < N_vec; k++){
			// in the lower block above the diagonal
			A[N_vec-1+k*N_vec][k*N_vec] = 0;
      // in the upper block above the diagonal
			A[k*N_vec][N_vec-1+k*N_vec] = 0;
		}
    #pragma omp for nowait
		for (int k = 0; k < N_vec-2; k++)
		{
			// in the lower block below the diagonal
			A[2*N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // in the upper block below the diagonal
			A[N_vec-1+k*N_vec][2*N_vec+k*N_vec] = 0;
		}
  }

  // Calcualte the L and U matrix of the Poisson Matrix A
  for(int k=0;k<N_pot-1;k++)
  {
    for(int i=k+1;i<N_pot;i++)
    {
      A[i][k] = A[i][k]/A[k][k];

      for(int j=k+1;j<N_pot;j++)
      {
        A[i][j] = A[i][j]-A[i][k]*A[k][j];
      }
    }
  }
}
