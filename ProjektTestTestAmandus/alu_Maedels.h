/*!
 * \brief Calculates the 5-point poisson matrix for the roughest grid and its LU-decomposition
 *
 * @param[in] A       matrix A
 * @param[in] lower   matrix L
 * @param[in] upper   matrix U
 * @param[in] N_vec   size of the grid in one dimension
 * @param[in] N_pot   size of the matrices in one dimension
 *
 * @param[out] A      poisson matrix
 * @param[out] lower  lower triangle matrix from the LU-decomposition
 * @param[out] upper  upper triangle matrix from the LU-decomposition
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void initLU5(double **A, double **lower,double **upper,int N_vec,int N_pot){

  printf("1");

  // Insert the right values into the correspondix diagonales of the matrix A
  int i,j,k;
  #pragma omp parallel
  {
    #pragma omp for nowait
    for(int i=0;i<N_pot-N_vec;i++)
    {
            A[i][i]    = 4.0;
            A[i][i+N_vec] = -1.0;
            A[i+N_vec][i] = -1.0;
    }

    #pragma omp for nowait
    for(int i=N_pot-N_vec;i<N_pot;i++)
    {
            A[i][i] = 4.0;
    }

    #pragma omp for collapse(2) nowait
    for(int k=0;k<N_vec;k++)
    {
            for(int i=0;i<N_vec-1;i++)
            {
                    A[k*N_vec+i][k*N_vec+i+1] = -1.0;
                    A[k*N_vec+i+1][k*N_vec+i] = -1.0;
            }
    }
  }

  #pragma omp parallel
  {
    // Setting some values on the diagonales to 0 to get the correct Poisson-matrix
    for (int k = 0; k < N_vec-1; k++){
      // below the diagonal
      A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // above the diagonal
      A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
    }
  }

  printf("2");

  #pragma omp parallel
  {
    // Create L and U of the Poisson Matrix for the smallest grid
    #pragma omp for nowait
    for(int j=0; j<N_pot; j++)
    {
      for(int i=0; i<N_pot; i++)
      {
        if(i<=j)
        {
          upper[i][j]=A[i][j];

          for(int k=0; k<=i-1; k++)
              upper[i][j]-=lower[i][k]*upper[k][j];
          if(i==j)
              lower[i][j]=1;
          else
              lower[i][j]=0;
        }
        else
        {
          lower[i][j]=A[i][j];
          for(int k=0; k<=j-1; k++)
              lower[i][j]-=lower[i][k]*upper[k][j];
          lower[i][j]/=upper[j][j];
          upper[i][j]=0;
        }
      }
    }
  }
  printf("3");
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------end of 5-point stencil--------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

/*!
 * \brief Calculates the 9-point poisson matrix for the roughest grid and its LU-decomposition
 *
 * @param[in] A       matrix A
 * @param[in] lower   matrix L
 * @param[in] upper   matrix U
 * @param[in] N_vec   size of the grid in one dimension
 * @param[in] N_pot   size of the matrices in one dimension
 *
 * @param[out] A      poisson matrix
 * @param[out] lower  lower triangle matrix from the LU-decomposition
 * @param[out] upper  upper triangle matrix from the LU-decomposition
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void initLU9(double **A, double **lower,double **upper,int N_vec,int N_pot){

  // Insert the right values into the correspondix diagonales of the matrix A
  int i,j,k;
  #pragma omp parallel
  {
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
    for (int k = 0; k < N_vec-1; k++){
      // below the diagonal
      A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // above the diagonal
      A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
    }

    // Upper and lower blocks
		for (int k = 0; k < N_vec; k++){
			// in the lower block above the diagonal
			A[N_vec-1+k*N_vec][k*N_vec] = 0;
      // in the upper block above the diagonal
			A[k*N_vec][N_vec-1+k*N_vec] = 0;
		}
		for (int k = 0; k < N_vec-2; k++)
		{
			// in the lower block below the diagonal
			A[2*N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
      // in the upper block below the diagonal
			A[N_vec-1+k*N_vec][2*N_vec+k*N_vec] = 0;
		}
  }

  // Calcualte the L and U matrix of the Poisson Matrix A
  for(int j=0; j<N_pot; j++)
  {
    for(int i=0; i<N_pot; i++)
    {
      if(i<=j)
      {
        upper[i][j]=A[i][j];

        for(int k=0; k<=i-1; k++)
            upper[i][j]-=lower[i][k]*upper[k][j];
        if(i==j)
            lower[i][j]=1;
        else
            lower[i][j]=0;
      }
      else
      {
        lower[i][j]=A[i][j];
        for(int k=0; k<=j-1; k++)
            lower[i][j]-=lower[i][k]*upper[k][j];
        lower[i][j]/=upper[j][j];
        upper[i][j]=0;
      }
    }
  }
}
