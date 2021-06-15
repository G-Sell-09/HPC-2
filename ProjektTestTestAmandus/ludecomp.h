/*
 * \brief LU forward backward propogation
 *
 * @param[in] x 			vector x
 * @param[in] b 			right hand side
 * @param[in] N 			size of the matrices in one dimension
 * @param[in] lower 	matrix L
 * @param[in] upper 	matrix U
 *
 * @param[out] x 			vector x
 * @param[out] b 			right hand side
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void lu_Decomp(double *x, double *b, int N, double **lower, double **upper)
{

		double *y = malloc(N*sizeof(double));

    // forward propagation
    for(int i=0; i<N; i++)
    {
        y[i]=b[i];
        for(int j=0; j<i; j++)
        {
            y[i]-=lower[i][j]*y[j];
        }
    }

		// backwards propagation
    for(int i=N-1; i>=0; i--)
    {
        x[i]= y[i];
        for(int j=i+1; j<N; j++)
        {
            x[i]-=upper[i][j]*x[j];
        }
        x[i]/=upper[i][i];
    }
		free(y);
}
