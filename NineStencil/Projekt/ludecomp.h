/*
 * \brief LU forward backward propagation
 *
 * @param[in] x 			vector
 * @param[in] b 			right hand side
 * @param[in] N 			size of the roughest grid
 * @param[in] A				L + U for LU decomposition
 *
 * @param[out] x 			values after forward backward propagation
 * @param[out] b 			right hand side
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void lu_Decomp(double *x, double *b, int N, double **A)
{

		double *y = malloc(N*sizeof(double));

    // forward propagation
    for(int i=0; i<N; i++)
    {
        y[i]=b[i];
        for(int j=0; j<i; j++)
        {
            y[i]-=A[i][j]*y[j];
        }
    }

		// backwards propagation
    for(int i=N-1; i>=0; i--)
    {
        x[i]= y[i];
        for(int j=i+1; j<N; j++)
        {
            x[i]-=A[i][j]*x[j];
        }
        x[i]/=A[i][i];
    }
		free(y);
}
