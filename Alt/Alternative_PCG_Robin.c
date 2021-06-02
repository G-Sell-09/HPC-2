#define _GNU_SOURCE
#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>


void matrixfree_A(int N, const double* restrict r, double* restrict y)
{
    int i, j;
    #pragma omp parallel
    {
    	// Inside
        #pragma omp for private(i, j) nowait
    	for (i=1; i<N-1; i++)
    		for ( j=1; j<N-1; j++ )
    			y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

    	// Edges
        #pragma omp for private(i) nowait
    	for (int i=1; i<N-1; i++)
        {
    		y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // lower edge
    		y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // upper edge
    		y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // left edge
    		y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // right edge
    	}

        // Corners
        #pragma omp single
        {
        	y[0]		 = 4*r[0]             - r[N]             - r[1];                   // bottom left
        	y[N-1]		 = 4*r[N-1]           - r[2*N-1]         - r[N - 2];               // bottom right
        	y[(N-1)*N] 	 = 4*r[(N-1)*N]       - r[(N-1)*N + 1]   - r[(N-2)*N];		       // top left
        	y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];   // top right
        }
    }
}


double dot(int N, const double* restrict x, const double* restrict y)
{
	double result = 0;
    #pragma omp parallel for reduction(+: result)
	for ( int i=0; i<N*N; i++ )
		result += x[i]*y[i];
	return result;
}


int cg(int N, const double* restrict b, double* restrict x, double rel_tol,
    int maxIter)
{
    int i, cg_iter;
    double alpha, beta, r_zero_norm, rdot_old, rdot_new, pAp;
    double *r = malloc(N*N*sizeof(double));
    double *p = malloc(N*N*sizeof(double));
    double *Ap = malloc(N*N*sizeof(double));

    // Init
    cg_iter = 0;
    matrixfree_A(N, x, r);  // A * x_0
    #pragma omp parallel for
    for(i=0; i<N*N; i++)
    {
        r[i] = b[i] - r[i];  // r_0 = b - Ax_0
        p[i] = r[i];  // p_0 = r_0
    }
    rdot_old = dot(N, r, r);
    r_zero_norm = sqrt(rdot_old);

    printf("    r zero norm: %.2e\n", r_zero_norm);

    while(cg_iter < maxIter)
    {
        cg_iter +=1;
        // A*p_k
        matrixfree_A(N, p, Ap);

        // dot(p_k, A * p_k)
        pAp = 0;
        pAp = dot(N, p, Ap);

        // alpha = dot(r_k, r_k)/dot(p_k, A * p_k)
        alpha = rdot_old / pAp;

        // x_k+1 = x_k + alpha *p_k
        // r_k+1 = r_k - alpha_k * A * p_k
        #pragma omp parallel for
        for (i=0; i<N*N; i++)
        {
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * Ap[i];
        }

        // dot(r_k+1, r_k+1)
        rdot_new = dot(N, r, r);

        // Abbruchskriterium
        if ((sqrt(rdot_new) / r_zero_norm) < rel_tol)
        {
            printf("    Reached relative tolerance after %i iterations with %.2e\n", cg_iter, sqrt(rdot_new));
            break;
        }
        if (cg_iter == maxIter)
        {
            printf("    Reached max Iterations %i\n with %.2e", cg_iter, sqrt(rdot_new));
            break;
        }

        // beta = dot(r_k+1, r_k+1)/dot(r_k, r_k)
        beta = rdot_new / rdot_old;

        // p_k+1 = r_k+1 + beta * p_k
        #pragma omp parallel for
        for (i=0; i<N*N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }

        rdot_old = rdot_new;
    }

    free(r);
    free(p);
    free(Ap);
}




int iterative_ILU_0(double* LL, double*UU, int N)
{
  // Moeglichen Startwert aus Alg. 1.8 setzen
  for (int i = N; i < N*(N+1); i++)
  {
      LL[i]         = -1;
      LL[N*N+i]     = -1;
      UU[i-N]       =  4;
      UU[N*N+i-N]   = -1;
      UU[2*N*N+i-N] = -1;
  }

  // unter LL und oberhalb von UU Nullen-Layer
  for (int i = 0; i < N; i++)
  {
      LL[i]         = 0;
      UU[3*N*N-1-i] = 0;
  }

  for (int i = 0; i < 2*N*N; i+=N)
  {
      LL[N + i    ]   = 0;
      UU[N*N + i+N-1] = 0;
  }



  for ( int i = 0; i < N*N; i++)
  {
    // Da ein Layer unterhalb von LL den Index um +N erhöhen
    // LL(i,1) = (-1)/UU(i-N,1)
    LL[N+i] = -1/UU[i-N];
    // LL(i,2) = (-1)/UU(i-1,1)
    LL[N+N*N+i] = -1/UU[i-1];

    // Da ein Layer oberhalb von UU muss nichts beim Index angepasst werden
    // UU(i,1) = 4-LL(i,1)*UU(i-N,3)-LL(i,2)*UU(i-1,2)
    UU[N+i] =  4 - LL[i] * UU[2*N*N + i-N] - LL[N*N + i] * UU[N*N + i-1];
    // UU(i,2) und UU(i,3) bleiben unverändert (s. VL)
}








int iterative_ILU_0(double* LL, double*UU, int N)
{
  // Moeglichen Startwert aus Alg. 1.8 setzen
  for (int i = N+1; i < (N+1)*(N+1); i++)
  {
      LL[i]         = -1;
      LL[N*N+i]     = -1;
      UU[i-N]       =  4;
      UU[N*N+i-N]   = -1;
      UU[2*N*N+i-N] = -1;
  }

  // unter LL und oberhalb von UU Nullen-Layer
  for (int i = 0; i < N+1; i++)
  {
      LL[i]         = 0;
      UU[3*N*N-1-i] = 0;
  }

  for (int i = 0; i < (N+1)*(N+1); i+=N+1)
  {
      LL[i]   = 0;
      UU[N*N + i+N-1] = 0;
  }

  for ( int i = 0; i < N*N; i++)
  {
    // Da ein Layer unterhalb von LL den Index um +N erhöhen
    // LL(i,1) = (-1)/UU(i-N,1)
    LL[N+i] = -1/UU[i-N];
    // LL(i,2) = (-1)/UU(i-1,1)
    LL[N+N*N+i] = -1/UU[i-1];

    // Da ein Layer oberhalb von UU muss nichts beim Index angepasst werden
    // UU(i,1) = 4-LL(i,1)*UU(i-N,3)-LL(i,2)*UU(i-1,2)
    UU[N+i] =  4 - LL[i] * UU[2*N*N + i-N] - LL[N*N + i] * UU[N*N + i-1];
    // UU(i,2) und UU(i,3) bleiben unverändert (s. VL)
}














int main (int argc , char **argv)
{
    // relative tolerance
    double rel_tol= 1e-12;
    // grid size
    int N = 100;
    double h  = 1./(double)(N+1);
    double hh = h*h;

    // error
    double err;

    // solution
    double *x = malloc(N*N*sizeof(double));
    // right side
    double *b = malloc(N*N*sizeof(double));

    // LL und UU initialisieren
    double *LL = malloc(N+2*N*N*sizeof()double));
    double *UU = malloc(N+3*N*N*sizeof()double));

    // for timing
    double time;

    // define right side
    // f(x,y) = 2 * pi^2 * sin(pi * x) sin(pi * y)
	for (int i=0; i<N; i++ )
    {
        for (int j=0; j<N ; j++)
		{
			b[i*N + j] = hh * 2*M_PI*M_PI * sin(M_PI*(j+1)*h) * sin(M_PI*(i+1)*h);
		}
    }

    int max_threads = omp_get_max_threads();
    printf(" -- Maximum num. threads: %i\n", max_threads);
    printf(" -- Start \n");
    printf("\n");
    for (int threads=1; threads<max_threads+1; threads++)
    {
        // initial guess
        for (int i=0; i<N*N; i++)
                x[i] = 1;

        omp_set_num_threads(threads);
        time = omp_get_wtime();
    	cg(N, b, x, rel_tol, 5000);
    	time = omp_get_wtime() - time;

    	printf("    Runtime Parallel with %i threads: %f\n", threads, time);

        // cumulative error
        err = 0;
        for ( int i=0; i<N; i++  ){
            for ( int j=0; j<N ; j++ ){
                err += x[i*N + j] - sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h);
                // printf("x = %f,     true_x = %f\n", x[i*N + j], sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h));
            }
        }
        printf("    cum. error: %f\n", err);

        printf("\n");
    }
    printf(" -- End \n");
}
