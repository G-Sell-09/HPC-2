/*!
 * \brief Matrixfree matrix-vector multiplication with the Poisson Matrix of the 5 point stencil.
 *
 * @param[in] x   	vector x
 * @param[in] ax  	matrix-vector product
 * @param[in] N  		dimension N
 *
 * @param[out] ax		matrix-vector product
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
void mfMult5(double* x, double* ax, int N)
{

	int i,j;
	// inner grid
	for ( i = 1; i < N-1; i++ )
		for ( j = 1; j < N-1; j++ )
			ax[i*N + j] = 4*x[i*N + j] - x[i*N + j-1] - x[i*N + j+1] - x[(i-1)*N + j] - x[(i+1)*N + j];

	// edges
	for ( int i = 1; i < N-1; i++ )
	{
		ax[i]           = 4*x[i]           - x[i-1]           - x[i+1]           - x[N + i];           // bottom edge
		ax[(N-1)*N + i] = 4*x[(N-1)*N + i] - x[(N-1)*N + i-1] - x[(N-1)*N + i+1] - x[(N-2)*N + i];     // top edge
		ax[i*N]         = 4*x[i*N]         - x[i*N + 1]       - x[(i-1)*N]       - x[(i+1)*N];         // left edge
		ax[i*N + N-1]   = 4*x[i*N + N-1]   - x[i*N + N-2]     - x[(i-1)*N + N-1] - x[(i+1)*N + N-1];   // right edge
	}

	{
	// corners
	ax[0]		          = 4*x[0]             - x[N]             - x[1];  		        // bottom left
	ax[N-1]		        = 4*x[N-1]           - x[2*N-1]         - x[N - 2];		      // bottom right
	ax[(N-1)*N] 	    = 4*x[(N-1)*N]       - x[(N-1)*N + 1]   - x[(N-2)*N];		    // top left
	ax[(N-1)*N + N-1] = 4*x[(N-1)*N + N-1] - x[(N-1)*N + N-2] - x[(N-2)*N + N-1];	// top right
	}
}
//
// void mfMult5(double* x, double* ax, int N)
// {
//
// 	int i,j;
// 	#pragma omp parllel
// 	{
// 		// inner grid
// 		#pragma omp for private(i,j) nowait
// 		for ( i = 1; i < N-1; i++ )
// 			for ( j = 1; j < N-1; j++ )
// 				ax[i*N + j] = 4*x[i*N + j] - x[i*N + j-1] - x[i*N + j+1] - x[(i-1)*N + j] - x[(i+1)*N + j];
//
// 		// edges
// 		#pragma omp for private(i) nowait
// 		for ( int i = 1; i < N-1; i++ )
// 		{
// 			ax[i]           = 4*x[i]           - x[i-1]           - x[i+1]           - x[N + i];           // bottom edge
// 			ax[(N-1)*N + i] = 4*x[(N-1)*N + i] - x[(N-1)*N + i-1] - x[(N-1)*N + i+1] - x[(N-2)*N + i];     // top edge
// 			ax[i*N]         = 4*x[i*N]         - x[i*N + 1]       - x[(i-1)*N]       - x[(i+1)*N];         // left edge
// 			ax[i*N + N-1]   = 4*x[i*N + N-1]   - x[i*N + N-2]     - x[(i-1)*N + N-1] - x[(i+1)*N + N-1];   // right edge
// 		}
//
// 		#pragma omp single
// 		{
// 		// corners
// 		ax[0]		          = 4*x[0]             - x[N]             - x[1];  		        // bottom left
// 		ax[N-1]		        = 4*x[N-1]           - x[2*N-1]         - x[N - 2];		      // bottom right
// 		ax[(N-1)*N] 	    = 4*x[(N-1)*N]       - x[(N-1)*N + 1]   - x[(N-2)*N];		    // top left
// 		ax[(N-1)*N + N-1] = 4*x[(N-1)*N + N-1] - x[(N-1)*N + N-2] - x[(N-2)*N + N-1];	// top right
// 		}
// 	}
// }

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------end of 5-point mfMult--------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

/*!
 * \brief Matrixfree matrix-vector multiplication with the Poisson Matrix of the 9 point stencil.
 *
 * @param[in] x      Vector x
 * @param[in] ax     Matrix-vector product
 * @param[in] N  		 Dimension N
 *
 * @param[out]  ax   Matrix-vector product
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
void mfMult9(double* x, double* ax, int N)
{

	int i,j;
	#pragma omp parllel
	{
		// inner grid
		#pragma omp for private(i,j) nowait
		for ( i = 1; i < N-1; i++ )
			for ( j = 1; j < N-1; j++ )
				ax[i*N+j] = 20*x[i*N+j] - 4*x[i*N+j-1] - 4*x[i*N+j+1] - 4*x[(i-1)*N+j] - 4*x[(i+1)*N+j] - x[(i-1)*N+j-1] - x[(i-1)*N+j+1] - x[(i+1)*N+j+1]- x[(i+1)*N+j-1];

		// edges
		#pragma omp for private(i) nowait
		for ( int i = 1; i < N-1; i++ )
		{
			ax[i]         = 20*x[i]         - 4*x[i-1]         - 4*x[i+1]           - 4*x[i+N]        - x[i+N-1]       - x[i+N+1]   ;    // untere Kante
			ax[(N-1)*N+i] = 20*x[(N-1)*N+i] - 4*x[(N-1)*N+i-1] - 4*x[(N-1)*N+i+1]   - 4*x[(N-2)*N+i] - x[(N-2)*N+i-1] - x[(N-2)*N+i+1]; // obere Kante
			ax[i*N]       = 20*x[i*N]       - 4*x[i*N + 1]     - 4*x[(i-1)*N]       - 4*x[(i+1)*N]    - x[(i-1)*N+1]   - x[(i+1)*N+1]  ; // links Kante
			ax[i*N+N-1]   = 20*x[i*N+N-1]   - 4*x[i*N+N-2]     - 4*x[(i-1)*N+N-1] - 4*x[(i+1)*N+N-1] - x[(i-1)*N+N-2] - x[(i+1)*N+N-2];  // rechts Kante
		}

		#pragma omp single
		{
		// corners
		ax[0]             = 20*x[0]           - 4*x[1]           - 4*x[N]           - x[N+1];  		// unten links
		ax[N-1]           = 20*x[N-1]         - 4*x[N-2]         - 4*x[2*N-1]       - x[2*N-2];		// unten rechts
		ax[(N-1)*N]       = 20*x[(N-1)*N]     - 4*x[(N-1)*N+1]   - 4*x[(N-2)*N]     - x[(N-2)*N+1] ;		// oben links
		ax[(N-1)*N + N-1] = 20*x[(N-1)*N+N-1] - 4*x[(N-1)*N+N-2] - 4*x[(N-2)*N+N-1] - x[(N-2)*N+N-2];	// oben rechts
		}
	}
}
