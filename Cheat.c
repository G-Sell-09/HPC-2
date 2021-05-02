#define _GNU_SOURCE
#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>


void mfMult(int N, double* r, double* y)
{
	// core
	for ( int i = 1; i < N-1; i++ )
		for ( int j = 1; j < N-1; j++ )
			y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

	// margins
	for ( int i = 1; i < N-1; i++ )
	{
		y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // lower row
		y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // upper row
		y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // left column
		y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // right column
	}

	// corners
	y[0]		 = 4*r[0]             - r[N]             - r[1];  		// bottom left
	y[N-1]		 = 4*r[N-1]           - r[2*N-1]         - r[N - 2];		// bottom right
	y[(N-1)*N] 	 = 4*r[(N-1)*N]       - r[(N-1)*N + 1]   - r[(N-2)*N];		// top left
	y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];	// top right

}




double dot(int N, double* x, double* y)
{
	double result = 0;
	for ( int i = 0; i < N*N; i++ )
		result += x[i]*y[i];
	return result;
}



void cg(int N, double* b, double* x, double eps)
{

	int i;
	double alpha, beta, r_dot;
	double *r = malloc(N*N*sizeof(double));
	double *p = malloc(N*N*sizeof(double));
	double *m = malloc(N*N*sizeof(double));

	mfMult(N, x, r); // A*x0 = r0

	for ( i = 0; i < N*N ; i++ )
		r[i] = r[i] - b[i];	// r0 = r0 - b

	for ( i = 0; i < N*N ; i++)
		p[i] = -r[i];		// po = -r0

	double r0_norm = sqrt( dot(N, r, r) ); // for breaking condition

	while(1)
	{

		// alpha_k = < r_k, r_k > / < p_k, m_k >
		mfMult(N, p, m);
		r_dot = dot(N, r, r);
		alpha = r_dot / dot(N, p, m);

		// x_{k+1} = x_k + alpha_k * p_k
		for ( i = 0; i < N*N; i++ )
			x[i] += alpha * p[i];

		// r_{k+1} = r_k + alpha_k * m_k
		for ( i = 0; i < N*N; i++ )
			r[i] += alpha * m[i];

		// breaking condition
		if ( sqrt( dot(N, r, r) )/r0_norm < eps )
			break;

		// beta_{k+1} = < r_{k+1}, r_{k+1} > / < r_k, r_k  >
		beta = dot(N, r, r) / r_dot;

		// p_{k+1} = beta_{k+1} * p_k - r_{k=1}
		for ( i = 0; i < N*N; i++)
			p[i] =  beta * p[i] - r[i];

	}

	free(r);
	free(p);
	free(m);

}

int main (int argc , char **argv)
{

	if( argc < 3 )
	{
		printf("Not enough input parameters.\n");
		return 0;
	}

	int N = atoi(argv[1]);

	double h  = 1./(double)(N+1);
	double hh = h*h;

	double *x = malloc(N*N*sizeof(double));
	double *b = malloc(N*N*sizeof(double));

	// f(x,y) = 2 * pi^2 * sin(pi * x) sin(pi * y)
	for ( int i = 0; i < N; i++  )
		for ( int j = 0; j < N ; j++ )
		{
			x[i*N + j] = 1; // start vector
			b[i*N + j] = hh*2*M_PI*M_PI * sin( M_PI * (j+1)*h ) * sin( M_PI * (i+1)*h  ); // b_{iN+j} = hh*f_ji = hh*f(jh,ih)
		}

	double eps = 1;
	for (int i = 0; i < atoi(argv[2]); i++ )
		eps = eps/10.0;


	double time = omp_get_wtime();

	cg(N, b, x, eps);

	time = omp_get_wtime() - time;
	printf("Runtime: %f\n", time );

	if( atoi(argv[3]) )
	{
		printf("approx.  &&  exact solution\n");
		for ( int i = 0; i < N; i++  )
			for ( int j = 0; j < N ; j++  )
			{
				printf("%lf     ", x[i*N +j] );
				printf("%lf\n", sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h  ) ); //u_ij = u(ih, jh)
			}
	}

	free (x);
	free (b);

	return 0;

}
