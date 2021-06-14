
void lu_Decomp(double *X, double *B, int N, double **L, double **U)
{

		double *Y = malloc(N*sizeof(double));

    // forward propagation
    for(int i=0; i<N; i++)
    {
        Y[i]=B[i];
        for(int j=0; j<i; j++)
        {
            Y[i]-=L[i][j]*Y[j];
        }
    }

		// backwards propagation
    for(int i=N-1; i>=0; i--)
    {
        X[i]= Y[i];
        for(int j=i+1; j<N; j++)
        {
            X[i]-=U[i][j]*X[j];
        }
        X[i]/=U[i][i];
    }

		free(Y);
}
