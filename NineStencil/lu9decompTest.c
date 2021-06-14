#define _GNU_SOURCE
#include <stdio.h>
#include <omp.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>

void lu_Decomp(double *X, double *B, int N)
{
		int N_help = sqrt(N);
    int i,j,k;
    int n = N;
		// double *L = malloc(N*N*sizeof(double));
		// double *u = malloc(N*N*sizeof(double));
		double *Y = malloc(N*sizeof(double));

    double **A = (double**) malloc (N*sizeof(double*));
    double **L = (double**) malloc (N*sizeof(double*));
    double **U = (double**) malloc (N*sizeof(double*));

    for (int i = 0; i < N; i++)
    {
        A[i] = (double*) malloc (N*sizeof(double));
        L[i] = (double*) malloc (N*sizeof(double));
        U[i] = (double*) malloc (N*sizeof(double));
    }

		//double A[9][9] = {0};
		//double L[9][9] = {0};
		//double U[9][9] = {0};

    //double Y[9] = {0};




		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				if (i == j){
					A[i][j] = 20;
				}

				if (i==j+1 || i==(j+N_help)){
					A[i][j] = -4;
				}
				if (i==(j+N_help-1) || i==(j+N_help+1)){
					A[i][j] = -1;
				}
				if (i==j-1 || i==(j-N_help)){
					A[i][j] = -4;
				}
				if (i==(j-N_help-1) || i==(j-N_help+1)){
					A[i][j] = -1;
				}
			}
		}

		for (int k = 0; k < N_help-1; k++)
		{
					// unterhalb der Diagonalen
					A[N_help+k*N_help][N_help-1+k*N_help] = 0;
					// oberhalb der Diagonalen
					A[N_help-1+k*N_help][N_help+k*N_help] = 0;
		}

		// Unteren Block die Nullen setzen
		for (int k = 0; k < N_help; k++)
		{
					printf("UNTEN OBERHALB x: %d y: %d\n", N_help-1+k*N_help, k*N_help);
					// unterhalb der Diagonalen im unteren Block die obere Diagonale
					A[N_help-1+k*N_help][k*N_help] = 0;
		}

		for (int k = 0; k < N_help-2; k++)
		{
					printf("UNTEN UNTERHALB x: %d y: %d\n", 2*N_help+k*N_help, N_help-1+k*N_help);
					// unterhalb der Diagonalen im unteren Block die untere Diagonale
					A[2*N_help+k*N_help][N_help-1+k*N_help] = 0;
		}

		// Oberen Block die Nullen setzen
		for (int k = 0; k < N_help; k++)
		{
					printf("OBEN UNTERHALB x: %d y: %d\n", k*N_help, N_help-1+k*N_help);
					// unterhalb der Diagonalen im unteren Block die obere Diagonale
					A[k*N_help][N_help-1+k*N_help] = 0;
		}

		for (int k = 0; k < N_help-2; k++)
		{
					printf("OBEN OBERHALB x: %d y: %d\n", N_help-1+k*N_help, 2*N_help+k*N_help);
					// unterhalb der Diagonalen im unteren Block die untere Diagonale
					A[N_help-1+k*N_help][2*N_help+k*N_help] = 0;
		}

		for (int i=0;i<N;i++)
	  {
			for (int j=0;j<N;j++){
				if (j == 0){
					printf("\n");
				}
				printf("A: %f ", A[i][j]);
			}
	  }	printf("\n"); printf("\n");	printf("\n");





    for(j=0; j<n; j++)
        {
            for(i=0; i<n; i++)
            {
                if(i<=j)
                {
                    U[i][j]=A[i][j];
                    //for(k=0; k<i-1; k++)
                    for(k=0; k<=i-1; k++)
                        U[i][j]-=L[i][k]*U[k][j];
                    if(i==j)
                        L[i][j]=1;
                    else
                        L[i][j]=0;
                }
                else
                {
                    L[i][j]=A[i][j];
                    for(k=0; k<=j-1; k++)
                        L[i][j]-=L[i][k]*U[k][j];
                    L[i][j]/=U[j][j];
                    U[i][j]=0;
                }
            }
        }

        // printf("[L]: \n");
        // for(i=0; i<n; i++)
        // {
        //     for(j=0; j<n; j++)
        //         printf("%9.3f",L[i][j]);
        //     printf("\n");
        // }
        // printf("\n\n[U]: \n");
        // for(i=0; i<n; i++)
        // {
        //     for(j=0; j<n; j++)
        //         printf("%9.3f",U[i][j]);
        //     printf("\n");
        // }

        // Ly=b
        for(i=0; i<n; i++)
        {
            Y[i]=B[i];
            for(j=0; j<i; j++)
            {
                Y[i]-=L[i][j]*Y[j];
            }
        }

        // printf("\n\n[Y]: \n");
        // for(i=0; i<n; i++)
        // {
        //     printf("%9.3f",Y[i]);
        // }

        // Ux = y
        for(i=n-1; i>=0; i--)
        {
            X[i]= Y[i];
            for(j=i+1; j<n; j++)
            {
                X[i]-=U[i][j]*X[j];
            }
            X[i]/=U[i][i];
        }
        // printf("\n\n[X]: \n");
        // for(i=0; i<n; i++)
        // {
        //     printf("%9.3f",X[i]);
        // }

    }



    void main (int argc , char **argv)
    {
      int N = atoi(argv[1]); // für 4x4

    	double *x = malloc (N*sizeof(double));
    	double *b = malloc (N*sizeof(double));

    	// x und b Vektor der Länge 4
    	for (int i = 0; i < N; i++){
    		x[i]=0;
    		b[i]=1;
    	}

      lu_Decomp(x,b,N);
    }
