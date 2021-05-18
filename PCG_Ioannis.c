#define _GNU_SOURCE
#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>


void ILU0_Iter( int N, double* LL, double* UU);
void pcg(int N, double* b, double* x, double* LL, double* UU, double eps);
void mfMult(int N, const double* restrict r, double* restrict y);
double dot(int N, double* x, double* y);
void substitution(int N, double* z, double* r, double* LL, double* UU);

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

	double *LL = malloc(2*N*N*sizeof(double));
	double *UU = malloc(3*N*N*sizeof(double));

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

	ILU0_Iter(N, LL, UU);
	pcg(N, b, x, LL, UU, eps);

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

void ILU0_Iter( int N, double* LL, double* UU){

#pragma omp parallel
{
    // set startvalue L+U = A
#pragma omp for nowait schedule(static)
    for (int i = 0; i < N*N; i++)
    {
        LL[      + i] = -1;
        LL[  N*N + i] = -1;
        UU[        i] =  4;
        UU[  N*N + i] = -1;
        UU[2*N*N + i] = -1;
    }
#pragma omp for nowait schedule(static)
    for (int i = 0; i < N; i++)
    {
        LL[i] = 0;
        UU[2*N*N + N*N-i-1] = 0;
    }

#pragma omp for nowait schedule(static)
    for (int i = 0; i < N*N; i+=N)
    {
        LL[N*N + i    ] = 0;
        UU[N*N + i+N-1] = 0;
    }

}

    //Iter
    // wann soll ich aufhoeren?
    // abbruchsbedingung
    for (int k = 0; k < 5000; k++)
    {
        // bottom left corner
        // no update needed
#pragma omp parallel
		{
			// lower row + bottom right corner
#pragma omp for nowait schedule(static)
			for ( int i = 1; i < N; i++){
				LL[  N*N + i] = -1 / UU[i-1];
			}

			// core
#pragma omp for nowait schedule(static)
			for ( int i = N; i < N*N; i++)
			{
				// LL(i,1) = (-1)/UU(i-N,1)
				LL[        i] = -1 / UU[i-N];
				// LL(i,2) = (-1)/UU(i-1,1)
				LL[  N*N + i] = -1 / UU[i-1];
				// UU(i,1) = 4-LL(i,1)*UU(i-N,3)-LL(i,2)*UU(i-1,2)
				UU[        i] =  4 - LL[i] * UU[2*N*N + i-N] - LL[N*N + i] * UU[N*N + i-1];
				// UU(i,2) = -1 no update needed
				// UU(i,3) = -1 no update needed
			}
		}



    }

}

void mfMult(int N, const double* restrict r, double* restrict y)
{
	// A * r = y
#pragma omp parallel
{
	// core
#pragma omp for nowait schedule(static)
	for ( int i = 1; i < N-1; i++ )
		for ( int j = 1; j < N-1; j++ )
			y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

	// margins
#pragma omp for nowait schedule(static)
	for ( int i = 1; i < N-1; i++ )
	{
		y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // lower row
		y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // upper row
		y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // left column
		y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // right column
	}

	// corners
#pragma omp single nowait
{
	y[0]		 = 4*r[0]                 - r[N]             - r[1];  		        // bottom left
	y[N-1]		 = 4*r[N-1]               - r[2*N-1]         - r[N - 2];		    // bottom right
	y[(N-1)*N] 	 = 4*r[(N-1)*N]           - r[(N-1)*N + 1]   - r[(N-2)*N];	    	// top left
	y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];	// top right
}

}

}

double dot(int N, double* x, double* y)
{
	double result = 0;
#pragma omp parallel for schedule(static) reduction(+:result)
	for ( int i = 0; i < N*N; i++ )
		result += x[i]*y[i];
	return result;
}

void pcg(int N, double* b, double* x, double* LL, double* UU, double eps)
{

	int i;
	double alpha, beta, r_dot, z_dot;
	double *r = malloc(N*N*sizeof(double));
    double *z = malloc(N*N*sizeof(double));
	double *p = malloc(N*N*sizeof(double));
	double *m = malloc(N*N*sizeof(double));

	mfMult(N, x, r); // A*x0 = r0

#pragma omp parallel for schedule(static)
	for ( i = 0; i < N*N ; i++ )
		r[i] = r[i] - b[i];	// r0 = r0 - b

#pragma omp parallel for schedule(static)
	for ( i = 0; i < N*N ; i++)
		p[i] = -r[i];		// po = -r0

    substitution(N, z, r, LL, UU); // z_0 = (LU)^(-1) * r_o

	double z0_norm = sqrt( dot(N, z, z) ); // for breaking condition

	double rz_dot = dot(N, r, z);

    while(1)
	{

		// alpha_k = < r_k, r_k > / < p_k, m_k >
		mfMult(N, p, m);
		r_dot = rz_dot;
		alpha = r_dot / dot(N, p, m);

		// x_{k+1} = x_k + alpha_k * p_k
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++ )
			x[i] += alpha * p[i];

		// r_{k+1} = r_k + alpha_k * m_k
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++ )
			r[i] += alpha * m[i];

        // z_{k+1} = (LU)^(-1) * r_k <=> LU z = r
        substitution(N, z, r, LL, UU);

		// breaking condition
		if ( sqrt( dot(N, z, z) )/z0_norm < eps )
			break;

        rz_dot = dot(N, r, z);
		// beta_{k+1} = < r_{k+1}, r_{k+1} > / < r_k, r_k  >
		beta = rz_dot / r_dot;

		// p_{k+1} = beta_{k+1} * p_k - r_{k=1}
#pragma omp parallel for schedule(static)
		for ( i = 0; i < N*N; i++)
			p[i] =  beta * p[i] - z[i];

	}

	free(r);
    free(z);
	free(p);
	free(m);

}

void substitution(int N, double* z, double* r, double* LL, double* UU)
{

    // LU z = r
    double* x = malloc(N*N*sizeof(double));
    // L x = r
        x[0] = r[0];
    for (int i = 1; i < N; i++)
        x[i] = r[i] - LL[N*N + i] * x[i-1];
    for (int i = N; i < N*N; i++)
        x[i] = r[i] - LL[N*N + i] * x[i-1] - LL[i] * x[ i-N ];

    // U z = x
    z[N*N-1] = x[N*N-1] / UU[N*N-1];
    for (int i = N*N-2; i > N*(N-1)-1; i--)
        z[i] = ( x[i] - UU[N*N + i] * z[i+1] ) / UU[i];
    for (int i = N*(N-1)-1; i >= 0; i--)
        z[i] = ( x[i] - UU[N*N + i] * z[i+1] - UU[2*N*N + i]* z[i+N] ) / UU[i];

    free(x);

}
