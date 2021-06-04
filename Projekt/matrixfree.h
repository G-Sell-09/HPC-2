/// Wendet den 5-Punkte Differenzenstern für Mat-Vek Multiplikation an
void mfMult(double* r, double* y, int N)
{
	// innere Punkte
	int i,j;
	#pragma omp parllel
	{

		#pragma omp for private(i,j) nowait
		for ( i = 1; i < N-1; i++ )
			for ( j = 1; j < N-1; j++ )
				y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

		// Kanten
		#pragma omp for private(i) nowait
		for ( int i = 1; i < N-1; i++ )
		{
			y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // untere Kante
			y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // obere Kante
			y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // links Kante
			y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // rrechts Kante
		}

		#pragma omp single
		{
		// Ecken
		y[0]		 = 4*r[0]             - r[N]             - r[1];  		// unten links
		y[N-1]		 = 4*r[N-1]           - r[2*N-1]         - r[N - 2];		// unten rechts
		y[(N-1)*N] 	 = 4*r[(N-1)*N]       - r[(N-1)*N + 1]   - r[(N-2)*N];		// oben links
		y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];	// oben rechts
		}
	}
}
