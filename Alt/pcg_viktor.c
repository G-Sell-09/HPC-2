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
	for (int i=0; i<N*N; i++)
		result += x[i]*y[i];
	return result;
}

void forward_sub(int N, double** LL, double** UU, double* z, double* r)
{
    int k;
	z[0] = r[0];
	for (int j = 1; j <= N; j++)
	{
        //  1. Zeile
		if (j==1)
		{
			for( int i=2;i <=N;i++)
			{
				k = i-1;
                z[k]=r[k]-LL[k][1]*z[k-1];
			}
		}else{
			for( int i=1;i<=N;i++)
			{
				k=(j-1)*N+(i-1);
				z[k]=r[k]-LL[k][0]*z[k-N]-LL[k][1]*z[k-1];
			}
		}
	}
}


void backward_sub(int N, double** LL, double** UU, double* z, double* y)
{
    	int k;
    	z[(N-1)*N+(N-1)] = y[(N-1)*N+(N-1)] / UU[k][0];
    	for (int j = N; j >=1; j--)
    	{
            // N-te Zeile
    		if (j==N)
    		{
    			for(int i=N-1; i>=1;i--)
    			{
    				k = N*(N-1)+(i-1);
                    z[k] = (y[k] - UU[k][1]*z[k+1]) / UU[k][0];
    			}
    		}else{
    			for(int i=N; i>=1;i--)
    			{
    				k = (j-1) * N+(i-1);
                    z[k] = (y[k] - UU[k][1] *z[k+1] - UU[k][2]*z[k+N]) / UU[k][0];
    			}
    		}
    	}
}


void forward_sub_safe(int N, double** LL, double** UU, double* y, double* r)
{
	for (int k=0; k<N*N; k++){
		y[k] = r[k];
	}
}


void backward_sub_safe(int N, double** LL, double** UU, double* z, double* y)
{
	for (int k=N*N-1; k>=0; k--)
    {
	       z[k] = y[k];
    }
}


int pcg(int N, const double* restrict b, double* restrict x, double** LL,
    double** UU, double rel_tol, int maxIter)
{
    int i, cg_iter;
    double alpha, beta, z_zero_norm, zdot, pAp, rdotz_old, rdotz_new;
    double *r = malloc(N*N*sizeof(double));  // Residuum
    double *z = malloc(N*N*sizeof(double));  // M^(-1) * r
    double *y = malloc(N*N*sizeof(double));  // Zwischenprodukt forward_sub
    double *p = malloc(N*N*sizeof(double));  // Suchrichtung
    double *Ap = malloc(N*N*sizeof(double));

    // Init
    cg_iter = 0;
    matrixfree_A(N, x, r);  // A * x_0
    #pragma omp parallel for
    for(i=0; i<N*N; i++)
    {
        r[i] = b[i] - r[i];  // r_0 = b - Ax_0
    }

    // z_0 = M r_o
    forward_sub(N, LL, UU, y, r);
    backward_sub(N, LL, UU, z, y);

    // p_0 = h_0
    #pragma omp parallel for
    for(i=0; i<N*N; i++)
    {
        p[i] = z[i];   // p_0 = z_0
    }

    zdot = dot(N, z, z);
    rdotz_old = dot(N, r, z);
    z_zero_norm = sqrt(zdot);

    printf("    r zero norm: %.2e\n", z_zero_norm);

    while(cg_iter < maxIter)
    {
        cg_iter +=1;
        // A*p_k
        matrixfree_A(N, p, Ap);
        // dot(p_k, A * p_k)
        pAp = dot(N, p, Ap);
        // alpha = dot(r_k, z_k)/dot(p_k, A * p_k)
        alpha = rdotz_old / pAp;

        #pragma omp parallel for
        for (i=0; i<N*N; i++)
        {
            // x_k+1 = x_k + alpha *p_k
            x[i] = x[i] + alpha * p[i];
            // r_k+1 = r_k - alpha_k * A * p_k
            r[i] = r[i] - alpha * Ap[i];
        }

        // h_k+1 = M^(-1) r_k+1
        forward_sub(N, LL, UU, y, r);
        backward_sub(N, LL, UU, z, y);
        // dot(z_k+1, z_k+1)
        zdot = dot(N, z, z);
        // dot(r_k+1, z_k+1)
        rdotz_new = dot(N, r, z);
        // Abbruchskriterium
        if ((sqrt(zdot) / z_zero_norm) < rel_tol)
        {
            printf("    Reached relative tolerance after %i iterations with %.2e\n", cg_iter, sqrt(zdot));
            break;
        }
        if (cg_iter == maxIter)
        {
            printf("    Reached max Iterations %i with %.2e\n", cg_iter, sqrt(zdot));
            break;
        }
        // beta = dot(r_k+1, z_k+1)/dot(r_k, z_k)
        beta = rdotz_new / rdotz_old;
        // p_k+1 = - z_k+1 + beta * p_k
        #pragma omp parallel for
        for (i=0; i<N*N; i++)
        {
            p[i] = z[i] + beta * p[i];
        }
        rdotz_old = rdotz_new;
    }

    free(r);
    free(p);
    free(Ap);
}


int iterative_ILU_0 (double **LL, double **UU, int N, int iter)
{
	int i, j, k;

	#pragma parallel for collapse(2)
	for (j=0; j<N; j++) //
	{
		for (i=0; i<N; i++) //
		{
			k = j*N + i; // k-ter Gitterpunkt
            UU[k][0] = 4; // Gitterpunkt selbst
            if (j==0) {  // untere Zeile
                LL[k][0] = 0; // unterer Nachbar
            } else {
                LL[k][0] = -1; // unterer Nachbar
            }

            if (j==N-1){  // obere Zeile
                UU[k][2] = 0; // oberer Nachbar
            } else {
                UU[k][2] = -1; // oberer Nachbar
            }

            if (i==0){  // linke Spalte
                LL[k][1] = 0; // linker Nachbar
            } else {
                LL[k][1] = -1; // linker Nachbar
            }

            if (i==N-1){  // rechte Spalte
                UU[k][1] = 0; // rechter Nachbar
            } else {
                UU[k][1] = -1; // rechter Nachbar
            }
		}
	}


    /*
    for (int k=0; k<N*N; k++)
    {
        printf("%f     ", LL[k][0]);
        printf("%f     ", LL[k][1]);
        printf("\n");
    }
    for (int k=0; k<N*N; k++)
    {
        printf("%f     ", UU[k][0]);
        printf("%f     ", UU[k][1]);
        printf("%f     ", UU[k][2]);
        printf("\n");
    }
    */

	for (int l=0; l<iter; l++)
	{
        #pragma parallel for collapse(2)
		for (j=0; j<N; j++)
		{
			for (i=0; i<N; i++)
			{
				k = j*N+i; // k-ter Gitterpunkt
                UU[k][0] = 4; // Gitterpunkt selbst

                if (j>0){ // Nicht untere Zeile
                    UU[k][0] -= LL[k][0] * UU[k-N][2]; // Gitterpunkt selbst
                    LL[k][0] = -1 / UU[k-N][0]; // unterer Nachbar
                }
                if (j<N-1){ // Nicht obere Zeile
                    UU[k][2] = -1; // oberer Nachbar
                }
                if (i>0){ // Nicht linke Spalte
                    UU[k][0] -=  LL[k][1]*UU[k-1][1]; // Gitterpunkt selbst
    				LL[k][1] = -1 / UU[k-1][0]; // linker Nachbar
                }
                if (i<N-1){ // Nicht rechte Spalte
                    UU[k][1] = -1; // rechter Nachbar
                }
			}
		}
	}
	return 0;
}


int main (int argc , char **argv)
{
    // relative tolerance
    double rel_tol= 1e-12;
    // grid size
    int N = 200;
    double h  = 1./(double)(N+1);
    double hh = h*h;

    // error
    double err;

    // solution
    double *x = malloc(N*N*sizeof(double));
    // right side
    double *b = malloc(N*N*sizeof(double));

    double **LL = (double**) malloc (N*N*sizeof(double*));
    double **UU = (double**) malloc (N*N*sizeof(double*));
    for (int i = 0; i < N*N; i++)
    {
        LL[i] = (double*) malloc (2*sizeof(double));
        UU[i] = (double*) malloc (3*sizeof(double));
    }

    // for timing
    double time;

    #pragma omp parallel
    {
        // initialize initial guess
        #pragma omp for
        for (int i=0; i<N*N; i++)
                x[i] = 1;

        // define right side
        // f(x,y) = 2 * pi^2 * sin(pi * x) sin(pi * y)
        #pragma omp for collapse(2)
        for (int i=0; i<N; i++ )
        {
            for (int j=0; j<N ; j++)
    		{
    			b[i*N + j] = hh * 2*M_PI*M_PI * sin(M_PI*(j+1)*h) * sin(M_PI*(i+1)*h);
    		}
        }

        #pragma omp for
        for (int k=0; k<N*N; k++)
        {
            LL[k][0] = 0;
            LL[k][1] = 0;

            UU[k][0] = 1;
            UU[k][1] = 0;
            UU[k][2] = 0;
        }
    }

    printf(" -- Calling iterative ILU(0)\n");
    iterative_ILU_0(LL, UU, N, 15);

    int max_threads = omp_get_max_threads();
    printf(" -- Maximum num. threads: %i\n", max_threads);
    printf(" -- Start \n");
    printf("\n");
    for (int threads=1; threads<4+1; threads++)
    {
        omp_set_num_threads(threads);
        time = omp_get_wtime();
    	pcg(N, b, x, LL, UU, rel_tol, 1000);
    	time = omp_get_wtime() - time;

    	printf("    Runtime Parallel with %i threads: %f\n", threads, time);

        err = 0;
        #pragma omp parallel
        {
            // cumulative error
            #pragma omp for collapse(2) reduction(+: err) nowait
            for ( int i=0; i<N; i++  ){
                for ( int j=0; j<N ; j++ ){
                    err += x[i*N + j] - sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h);
                    // printf("x = %f,     true_x = %f\n", x[i*N + j], sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h));
                }
            }

            // re-initialize initial guess
            #pragma omp for nowait
            for (int i=0; i<N*N; i++)
                    x[i] = 1;
        }
        printf("    cum. error: %f\n", err);

        printf("\n");
    }
    printf(" -- End \n");
}
