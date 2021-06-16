/*
 * \brief SOR smoothing operator using the 5 point stencil.
 *
 * @param[in] x       initial values
 * @param[in] b       right hand side
 * @param[in] r       residual
 * @param[in] w       relaxation factor
 * @param[in] nu      number of iterations
 * @param[in] tol     tolerance
 * @param[in] N       grid size
 * @param[in] mfMult  typdef myMult for correct stencil
 *
 * @param[out] x      smoothed values
 * @param[out] b      right hand side
 * @param[out] r      residual, calculated with the new values
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void sor5(double *x, double *b, double *r, double w, double nu, double tol, int N, myMult mfMult)
{
  double *ax = malloc(N*N*sizeof(double));
  int i,j,k;

  for (k = 0;k<nu;k++)
  {
      // bottom left corner
  		ax[0] =  (1-w)*x[0] + (w/4.0) *(b[0] + x[N] + x[1]);

      // bottom edge
  		for (i = 1; i < N-1; i++ )
  		{
  			ax[i] = (1-w)*x[i] + (w/4.0) *(b[i] + ax[i-1] + x[i+1] + x[i+N]);
  		}

      // bottom right corner
  		ax[N-1] = (1-w)*x[N-1] + (w/4.0) *(b[N-1] + ax[N-2] + x[2*N-1]);

      //left edge
  		for (i = 1;i<N-1;i++)
  		{
  			ax[i*N] = (1-w)*x[i*N] + (w/4.0) *(b[i*N] + ax[(i-1)*N] + x[i*N+1] + x[(i+1)*N]);
  		}

  		// inner grid
  		for ( i = 1; i < N-1; i++ )
  			for ( j = 1; j < N-1; j++ )
  				ax[i*N+j] = (1-w)*x[i*N+j] + (w/4.0) *(b[i*N+j] + ax[(i-1)*N+j] + ax[i*N+j-1] + x[i*N+j+1] + x[(i+1)*N+j]);

  		// right edge
  		for (i = 1; i < N-1; i++ )
  		{
  			ax[i*N+N-1] = (1-w)*x[i*N+N-1] + (w/4.0) *(b[i*N+N-1] + ax[(i-1)*N + N-1] + ax[i*N+N-2] + x[(i+1)*N + N-1]);
  		}

      // upper left corner
  		ax[N*(N-1)] = (1-w)*x[N*(N-1)] + (w/4.0) *(b[N*(N-1)] + ax[N*(N-2)] + x[N*(N-1)+1]);

  		// upper edge
  		for (i = 1;i<N-1;i++)
  		{
  			ax[N*(N-1)+i] = (1-w)*x[N*(N-1)+i] + (w/4.0) *(b[N*(N-1)+i] + ax[N*(N-1)+i-1] + ax[N*(N-2)+i] + x[N*(N-1)+i+1]);
  		}

      // upper right corner
  		ax[N*N-1] = (1-w)*x[N*N-1] + (w/4.0) *(b[N*N-1] + ax[N*N-2] + ax[N*(N-1)-1]);

      duplicate(x,ax,N*N);

      residual(r, x, b, N, mfMult);

      if (norm2(r, N*N) / norm2(b,N*N) < tol)
        break;
    }
  free(ax);
  }

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  //----------------------end of 5-point stencil--------------------------------------------
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  /*
   * \brief SOR smoothing operator using the 9 point stencil.
   *
   * @param[in] x       initial values
   * @param[in] b       right hand side
   * @param[in] r       residual
   * @param[in] w       relaxation factor
   * @param[in] nu      number of iterations
   * @param[in] tol     tolerance
   * @param[in] N       grid size
   * @param[in] mfMult  typdef myMult for correct stencil
   *
   * @param[out] x      smoothed values
   * @param[out] b      right hand side
   * @param[out] r      residual, calculated with the new values
   *
   * \author Robin Sell
   * \author Neil Vetter
   *
   * \version 1.0
   * \copyright HPC-2 Team RSNV
   *
  */
  void sor9(double *x, double *b, double *r, double w, double nu, double tol, int N, myMult mfMult)
  {
    double *ax = malloc(N*N*sizeof(double));

    int i,j,k;

    for (k=0;k<nu;k++)
    {
      // bottom left corner
      ax[0] = (1-w)*x[0] + (w/20.0) *(b[0] + 4*x[1] + 4*x[N] + x[N+1]);

      // bottom edge
      for ( int i = 1; i < N-1; i++ )
        ax[i] = (1-w)*x[i] + (w/20.0) *(b[i] + 4*ax[i-1] + 4*x[i+1] + 4*x[i+N] + x[i+N-1] + x[i+N+1]);

      // bottom right corner
      ax[N-1] =  (1-w)*x[N-1] + (w/20.0) *(b[N-1] + 4*ax[N-2] + 4*x[2*N-1] + x[2*N-2]);

      // inner grid including the left and right edges
  		for ( i = 1; i < N-1; i++ ){
        // left edge
        ax[i*N] = (1-w)*x[i*N] + (w/20.0) *(b[i*N] + 4*x[i*N+1] + 4*ax[(i-1)*N] + 4*x[(i+1)*N] + ax[(i-1)*N+1] + x[(i+1)*N+1]);

        // inner grid
  			for ( j = 1; j < N-1; j++ ){
  				ax[i*N+j] = (1-w)*x[i*N+j] + (w/20.0) * (b[i*N+j] + 4*ax[i*N+j-1] + 4*x[i*N+j+1] + 4*ax[(i-1)*N+j] + 4*x[(i+1)*N+j] + ax[(i-1)*N+j-1] + ax[(i-1)*N+j+1] + x[(i+1)*N+j+1] + x[(i+1)*N+j-1]);
        }

        // right edge
        ax[i*N+N-1] = (1-w)*x[i*N+N-1] + (w/20.0) * (b[i*N+N-1] + 4*ax[i*N+N-2] + 4*ax[(i-1)*N+N-1] + 4*x[(i+1)*N+N-1] + ax[(i-1)*N+N-2] + x[(i+1)*N+N-2]);
      }

      // upper left corner
      ax[(N-1)*N] = (1-w)*x[(N-1)*N] + (w/20.0) *(b[(N-1)*N] + 4*x[(N-1)*N+1] + 4*ax[(N-2)*N] + ax[(N-2)*N+1]);

      // upper edge
  		for ( int i = 1; i < N-1; i++ )
  			ax[(N-1)*N+i] = (1-w)*x[(N-1)*N+i] + (w/20.0) *(b[(N-1)*N+i] + 4*ax[(N-1)*N+i-1] + 4*x[(N-1)*N+i+1] + 4*ax[(N-2)*N+i] + ax[(N-2)*N+i-1] + ax[(N-2)*N+i+1]);

      // upper right corner
      ax[(N-1)*N+N-1] = (1-w)*x[(N-1)*N+N-1] + (w/20.0) *(b[(N-1)*N+N-1] + 4*ax[(N-1)*N+N-2] + 4*ax[(N-2)*N+N-1] + ax[(N-2)*N+N-2]);

      duplicate(x,ax,N*N);

      residual(r, x, b, N, mfMult);

      if (norm2(r, N*N) / norm2(b,N*N) < tol)
        break;
    }
    free(ax);
  }
