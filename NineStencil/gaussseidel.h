void gau(double *x, double *b, double *r, double w, double nu, double tol, int N)
{
  double *ax = malloc(N*N*sizeof(double));

  int i,j,k;

  for (k=0;k<nu;k++)
  {
    ax[0] = (1-w)*x[0] + (w/20.0) *(b[0] + 4*x[1] + 4*x[N] + x[N+1]);  		// unten links

    // untere Kante
    for ( int i = 1; i < N-1; i++ )
      ax[i] = (1-w)*x[i] + (w/20.0) *(b[i] + 4*ax[i-1] + 4*x[i+1] + 4*x[i+N] + x[i+N-1] + x[i+N+1]);

    ax[N-1] =  (1-w)*x[N-1] + (w/20.0) *(b[N-1] + 4*ax[N-2] + 4*x[2*N-1] + x[2*N-2]);		// unten rechts

    // Mittlere Punkte und Kanten links u. rechts
		for ( i = 1; i < N-1; i++ ){
      ax[i*N] = (1-w)*x[i*N] + (w/20.0) *(b[i*N] + 4*x[i*N+1] + 4*ax[(i-1)*N] + 4*x[(i+1)*N] + ax[(i-1)*N+1] + x[(i+1)*N+1]); // links Kante

      // Mittlere Punkte
			for ( j = 1; j < N-1; j++ ){
				ax[i*N+j] = (1-w)*x[i*N+j] + (w/20.0) * (b[i*N+j] + 4*ax[i*N+j-1] + 4*x[i*N+j+1] + 4*ax[(i-1)*N+j] + 4*x[(i+1)*N+j] + ax[(i-1)*N+j-1] + ax[(i-1)*N+j+1] + x[(i+1)*N+j+1] + x[(i+1)*N+j-1]);
      }

      ax[i*N+N-1] = (1-w)*x[i*N+N-1] + (w/20.0) * (b[i*N+N-1] + 4*ax[i*N+N-2] + 4*ax[(i-1)*N+N-1] + 4*x[(i+1)*N+N-1] + ax[(i-1)*N+N-2] + x[(i+1)*N+N-2]);  // rechts Kante
    }

    ax[(N-1)*N] = (1-w)*x[(N-1)*N] + (w/20.0) *(b[(N-1)*N] + 4*x[(N-1)*N+1] + 4*ax[(N-2)*N] + ax[(N-2)*N+1]);		// oben links


		for ( int i = 1; i < N-1; i++ )
			ax[(N-1)*N+i] = (1-w)*x[(N-1)*N+i] + (w/20.0) *(b[(N-1)*N+i] + 4*ax[(N-1)*N+i-1] + 4*x[(N-1)*N+i+1] + 4*ax[(N-2)*N+i] + ax[(N-2)*N+i-1] + ax[(N-2)*N+i+1]); // obere Kante

    ax[(N-1)*N+N-1] = (1-w)*x[(N-1)*N+N-1] + (w/20.0) *(b[(N-1)*N+N-1] + 4*ax[(N-1)*N+N-2] + 4*ax[(N-2)*N+N-1] + ax[(N-2)*N+N-2]);	// oben rechts

    duplicate(x,ax,N*N);

    residual(r, x, b, N);

    if (norm2(r, N*N) / norm2(b,N*N) < tol)
      break;
  }
  free(ax);
}
