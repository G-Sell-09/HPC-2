void initLU5(double **A, double **lower,double **upper,int N_vec,int N_pot){

  // ToDo: Modulo Abfrage einbauen falls noch notwendig
  // Create the Poisson Matrix for the smallest grid
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

  for (int k = 0; k < N_vec-1; k++){
        A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
        A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
  }

  // Create L and U of the Poisson Matrix for the smallest grid
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

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------end of 5-point stencil--------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void initLU9(double **A, double **lower,double **upper,int N_vec,int N_pot){

  // ToDo: Modulo Abfrage einbauen falls noch notwendig
  // Create the Poisson Matrix for the smallest grid
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

  for (int k = 0; k < N_vec-1; k++){
        // unterhalb der Diagonalen
        A[N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;

        // oberhalb der Diagonalen
        A[N_vec-1+k*N_vec][N_vec+k*N_vec] = 0;
  }

  // Unteren Block die Nullen setzen
		for (int k = 0; k < N_vec; k++)
		{
					//printf("UNTEN OBERHALB x: %d y: %d\n", N_vec-1+k*N_vec, k*N_vec);
					// unterhalb der Diagonalen im unteren Block die obere Diagonale
					A[N_vec-1+k*N_vec][k*N_vec] = 0;
		}

		for (int k = 0; k < N_vec-2; k++)
		{
					//printf("UNTEN UNTERHALB x: %d y: %d\n", 2*N_vec+k*N_vec, N_vec-1+k*N_vec);
					// unterhalb der Diagonalen im unteren Block die untere Diagonale
					A[2*N_vec+k*N_vec][N_vec-1+k*N_vec] = 0;
		}

		// Oberen Block die Nullen setzen
		for (int k = 0; k < N_vec; k++)
		{
					//printf("OBEN UNTERHALB x: %d y: %d\n", k*N_vec, N_vec-1+k*N_vec);
					// unterhalb der Diagonalen im unteren Block die obere Diagonale
					A[k*N_vec][N_vec-1+k*N_vec] = 0;
		}

		for (int k = 0; k < N_vec-2; k++)
		{
					//printf("OBEN OBERHALB x: %d y: %d\n", N_vec-1+k*N_vec, 2*N_vec+k*N_vec);
					// unterhalb der Diagonalen im unteren Block die untere Diagonale
					A[N_vec-1+k*N_vec][2*N_vec+k*N_vec] = 0;
		}

  // Create L and U of the Poisson Matrix for the smallest grid
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
