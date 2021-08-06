#include <stdio.h>
#include <stdlib.h>




// Anzahl Threads pro Blockdimension
const int blocksize = 32;




/////////////////////////////////////////////////////////////////
//                       Hilfsfunktionen                       //
/////////////////////////////////////////////////////////////////


/* Funktion zur Ausgabe eines Vektors. */
void print_vec (double *v, int dim, char *name)
{
        printf("\n%s:\n", name);
        for (int i = 0; i < dim; i++)
        {
                printf("% lf\n" , v[i]);
        }
        printf("\n");
}


/* Funktion zum Ueberfuehren eines Vektors auf den inneren Gitterpunkten
 auf das Gitter mit Ghost Layer. */
void ghostify (double *v, double *v_ghost, int N)
{
    int k; // Gitterpunkt des Gitteers mit Ghost Layer
    int ind = 0; // Aktueller Eintrag des Ausgangsvektors

    // Durchgehen des Gitters mit Ghost Layer
    for (int j=1; j <= N+2; j++) // Gitterzeile
    {
        for (int i=1; i <= N+2; i++) // Gitterspalte
        {
            k = (j-1)*(N+2)+(i-1); // Gitterpunkt

            if ( (k < N+2) || (k%(N+2) == 0) || (k%(N+2) == N+1) || (k >= (N+1)*(N+2)) ) // Gitterpunkt faellt in Ghost Layer
            {
                v_ghost[k] = 0;
            }
            else
            {
                v_ghost[k] = v[ind]; // Gitterpunkt ist innerer Gitterpunkt
                ind += 1;
            }
        }
    }
}


/* Funktion zum Ueberfuehren eines Vektors auf dem Gitter mit Ghost Layer
 auf die inneren Gitterpunkte. */
void deghostify (double *v_ghost, double *v, int N)
{
    int k; // Gitterpunkt des Gitters mit Ghost Layer
    int ind = 0; // Aktueller Eintrag des Zielvektors

    // Durchgehen des Gitters mit Ghost Layer nur auf den inneren Punkten
    for (int j=2; j <= N+1; j++) // Gitterzeile
    {
        for (int i=2; i <= N+1; i++) // Gitterspalte
        {
            k = (j-1)*(N+2) + (i-1); // Gitterpunkt

            v[ind] = v_ghost[k];
            ind += 1;
        }
    }
}


/* Funktion zur Bestimmung des Skalarprodukts zweier Vektoren fester Laenge. */
double dot (double *v1, double *v2, int dim)
{
    double res = 0; // Ergebnis

        for(int i = 0; i < dim; i++)
        {
                res += v1[i]*v2[i];
        }

        return res;
}


/* Funktion zur Berechnung der euklidischen Norm eines Vektors. */
double norm_cpu (double *v, int dim)
{
        double res = sqrt(dot(v,v,dim));
        return res;
}


/* Funktion zur Berechnung der Wurzel von einem Vektor mit einem Eintrag. */
void sq (double *v)
{
        v[0]=sqrt(v[0]);
}

/* Funktion zur Bestimmung des absoluten Fehlers der Naeherung eines Loesungsvektors als Norm der Differenz
 zwischen der Naeherung und der Loesung. */
double get_abs_err (double *u, double *u_exakt, int dim)
{
        // Differenzvektor
        double *diff = (double *)malloc(dim*sizeof(double));

        for (int i = 0; i < dim; i++)
        {
                diff[i] = u[i] - u_exakt[i];
        }

        // Absoluter Fehler als Norm des Differenzvektors
        double abs_err = norm_cpu(diff,dim);

        free(diff);

        return abs_err;
}

/////////////////////////////////////////////////////////////////
//            Funktionen zum Befuellen von Vektoren            //
/////////////////////////////////////////////////////////////////


/* Funktion zum Befuellen eines Vektors mit der exakten Loesung des Poisson-Problems an
 den inneren Gitterpunkten nach zeilenweiser Nummerierung. */
void get_solution (double *u_exakt, int N)
{
        // Schrittweite
        double h = 1/(double)(N+1);

        // Durchgehen des Gitters
        for (int j=1; j<=N; j++) // Gitterzeile
        {
                for (int i=1; i<=N; i++) // Gitterspalte
                {
                        u_exakt[(i-1)+(j-1)*N]=sin(M_PI*i*h)*sin(M_PI*j*h);
                }
        }
}


/* Funktion zum Befuellen eines Vektors mit der Auswertung der gegebenen rechten Seite des
 Poissons-Problems an den inneren Gitterpunkten nach zeilenweiser Nummerierung. */
void get_RHS (double *f, int N)
{
        // Schrittweite
        double h = 1/(double)(N+1);

        // Durchgehen des Gitters
        for (int j=1; j<=N; j++) // Gitterzeile
        {
                for (int i=1; i<=N; i++) // Gitterspalte
                {
                        f[(i-1)+(j-1)*N]=h*h*2*M_PI*M_PI*sin(M_PI*i*h)*sin(M_PI*j*h); // h^2 auf rechter Seite
                }
        }
}


/* Funktion zum Befuellen eines Vektors mit zufaelligen Eintraegen in einem spezifizierten Wertebereich. */
void get_random_vec (double *v, int dim, double min, double max)
{
        // Setzen eines Random-Seeds
        srand(time(0));

        // Fuellen des Vektors mit zufaelligen Eintraegen
        for (int i = 0; i < dim; i++)
        {
                v[i] = min + ((double)rand()/(double)RAND_MAX) * (max - min);
        }
}


/* Funktion zum Befuellen eines Vektors als Startvektor fuer das Jacobi-Verfahren. */
void get_u0 (double *u0, int N, char fill_type)
{
        // Zufaelliger Startvektor
        if (fill_type == 'r')
        {
                get_random_vec(u0,N*N,0,1);
        }
        // Testvektor nur aus Einsen
        else if (fill_type == 't')
        {
                for (int i = 0; i < N*N; i++)
                {
                        u0[i] = 1;
                }
        }
}


/////////////////////////////////////////////////////////////////
//                      Kernel-Funktionen                      //
/////////////////////////////////////////////////////////////////


/* Kernel-Funktion zur eintragsweisen einmaligen Jacobi-Iteration.
 Jeder Thread berechnet einen Eintrag der neuen Iterierten. */
__global__ void jacobi (double *u_ghost, double *f, int N)
{

    // Statische Allokation des Shared Memory; Daten von einem Teilgitter sollen geladen werden
    __shared__ double u_ghost_loc[blocksize*blocksize];
    __shared__ double f_loc[(blocksize-2)*(blocksize-2)];


    // Bestimmung der Indizes des Vektoreintrags, der zu aufrufendem Thread korrespondiert
    int k_ghost = blockIdx.y*(blocksize-2)*(N+2) + threadIdx.y*(N+2) + blockIdx.x*(blocksize-2) + threadIdx.x ; // Globaler Index auf Gitter mit Ghost Layer
    int k_ghost_loc = threadIdx.y*blocksize + threadIdx.x; // Lokaler Index auf Gitter mit Ghost Layer

    int k = (blockIdx.y*(blocksize-2) + threadIdx.y-1)*N + blockIdx.x*(blocksize-2) + threadIdx.x-1; // Globaler Index auf inneren Gitterpunkten
    int k_loc = (threadIdx.y-1)*(blocksize-2) + threadIdx.x-1; // Lokaler Index auf inneren Gitterpunkten

    // Achtung: Die Betrachtung von k und k_loc ist nur fuer Threads sinnvoll, die tatsaechlich zu inneren Gitterpunkten korrespondieren!


    // Lade Daten von Teilgitter in Shared Memory; jeder Thread uebernimmt einen Eintrag
    u_ghost_loc[k_ghost_loc] = u_ghost[k_ghost]; // Eintraege der vorherigen Iterierten
    if ( (k_ghost_loc >= blocksize) && (k_ghost_loc%blocksize != 0) && (k_ghost_loc%blocksize != blocksize-1) && (k_ghost_loc < (blocksize-1)*blocksize) )
    {
        f_loc[k_loc] = f[k]; // Eintraege der rechten Seite werden nur fuer innere Gitterpunkte geladen
    }

    __syncthreads();


    // Kernel: Eintragsweise Jacobi-Iteration; nur Threads zu inneren Gitterpunkten des Teilgitters rechnen; Ergebnis wird direkt in Global Memory geschrieben
    if ( (k_ghost_loc >= blocksize) && (k_ghost_loc%blocksize != 0) && (k_ghost_loc%blocksize != blocksize-1) && (k_ghost_loc < (blocksize-1)*blocksize) )
    {
        u_ghost[k_ghost] = (f_loc[k_loc] + u_ghost_loc[k_ghost_loc-blocksize] + u_ghost_loc[k_ghost_loc-1] + u_ghost_loc[k_ghost_loc+1] + u_ghost_loc[k_ghost_loc+blocksize])/4;
    }

}


/* Kernel-Funktion zur eintragsweisen Bestimmung des Residuums einer Iterierten des Jacobi-Verfahrens.
 Jeder Thread berechnet einen Eintrag des Residuums. */
__global__ void residual (double *res, double *u_ghost, double *f, int N)
{

    // Statische Allokation des Shared Memory; Daten von einem Teilgitter sollen geladen werden
    __shared__ double u_ghost_loc[blocksize*blocksize];
    __shared__ double f_loc[(blocksize-2)*(blocksize-2)];


    // Bestimmung der Indizes des Vektoreintrags, der zu aufrufendem Thread korrespondiert
    int k_ghost = blockIdx.y*(blocksize-2)*(N+2) + threadIdx.y*(N+2) + blockIdx.x*(blocksize-2) + threadIdx.x ; // Globaler Index auf Gitter mit Ghost Layer
    int k_ghost_loc = threadIdx.y*blocksize + threadIdx.x; // Lokaler Index auf Gitter mit Ghost Layer

    int k = (blockIdx.y*(blocksize-2) + threadIdx.y-1)*N + blockIdx.x*(blocksize-2) + threadIdx.x-1; // Globaler Index auf inneren Gitterpunkten
    int k_loc = (threadIdx.y-1)*(blocksize-2) + threadIdx.x-1; // Lokaler Index auf inneren Gitterpunkten

    // Achtung: Die Betrachtung von k und k_loc ist nur fuer Threads sinnvoll, die tatsaechlich zu inneren Gitterpunkten korrespondieren!


    // Lade Daten von Teilgitter in Shared Memory; jeder Thread uebernimmt einen Eintrag
    u_ghost_loc[k_ghost_loc] = u_ghost[k_ghost]; // Eintraege der vorherigen Iterierten
    if ( (k_ghost_loc >= blocksize) && (k_ghost_loc%blocksize != 0) && (k_ghost_loc%blocksize != blocksize-1) && (k_ghost_loc < (blocksize-1)*blocksize) )
    {
        f_loc[k_loc] = f[k]; // Eintraege der rechten Seite werden nur fuer innere Gitterpunkte geladen
    }

    __syncthreads();


    // Kernel: Eintragsweise Berechnung des Residuums; nur Threads zu inneren Gitterpunkten des Teilgitters rechnen; Ergebnis wird direkt in Global Memory geschrieben
    if ( (k_ghost_loc >= blocksize) && (k_ghost_loc%blocksize != 0) && (k_ghost_loc%blocksize != blocksize-1) && (k_ghost_loc < (blocksize-1)*blocksize) )
    {
        res[k] = (- u_ghost_loc[k_ghost_loc-blocksize] - u_ghost_loc[k_ghost_loc-1] + 4*u_ghost_loc[k_ghost_loc] - u_ghost_loc[k_ghost_loc+1] -u_ghost_loc[k_ghost_loc+blocksize]) - f_loc[k_loc];
    }

}



/* Kernel-Funktion zur punktweisen multiplikation zweier Vektoren */
__global__ void pmult (double *v1, double *v2)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    //Ergebnis liegt nun in v1
    v1[i]*=v2[i];
}



/* Kernel-Funktion zur Reduktion bzw. Fan-in */
__global__ void reduction (double *iv, double *ov)
{
    // Statische Allokation des Shared Memory; Groesse des Blocks
    __shared__ double sv[blocksize-2];

    //Hole Daten in shared memory
    int tid = threadIdx.x;
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    sv[tid]=iv[i];
    __syncthreads();

    //Reduktion in shared memory auf Block durchfuehren
    for ( int s=1; s<blockDim.x; s*=2)
    {
                // Nur gerade Threads rechnen
                if (tid%(2*s)==0)
                {
                        //Wenn Zugriff innerhalb des Vektors erfolgt
                        if (tid+s<blocksize-2)
                        {
                                sv[tid]+=sv[tid+s];
                        }
                }
                __syncthreads();
        }

        //Ergebnis herausschreiben
        if(tid==0)
        {
                ov[blockIdx.x]=sv[0];
        }

}
/* Optimierte Kernel-Funktion zur Reduktion bzw. Fan-in */
__global__ void reduction2 (double *iv, double *ov)
{
    // Statische Allokation des Shared Memory; Groesse des Blocks
    __shared__ double sv[blocksize];

    //Hole Daten in shared memory
    int tid = threadIdx.x;

    // Direkte Addition zweier Bloecke
    int i = blockIdx.x*(blockDim.x*2)+threadIdx.x;
    sv[tid]=iv[i]+iv[i+blockDim.x];
    __syncthreads();

    //Reduktion in shared memory auf Block durchfuehren
    for ( int s=blockDim.x/2; s>0; s>>=1)
    {
                if (tid<s)
                {
                        //Wenn Zugriff innerhalb des Vektors erfolgt
                        if (tid+s<blocksize)
                        {
                                sv[tid]+=sv[tid+s];
                        }
                }
                __syncthreads();
        }

        //Ergebnis herausschreiben
        if(tid==0)
        {
                ov[blockIdx.x]=sv[0];
        }
}

/* Kernel-Funktion zur Bestimmung der Norm eines Vektors auf den inneren Gitterpunkten.
__global__ void norm (double *v, int N)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    //Ergebnis liegt nun in v
    v[i]*=v[i];
    __syncthreads();
    int temp = v[0];
    v[0]=0;
    v[0]+=v[i];
    v[0]+=temp;

}*/



////////////////////////////////////////////////////////////////
//                        Main-Funktion                        //
/////////////////////////////////////////////////////////////////


int main()
{

    printf("\n*******************************************\n");
    printf("*      Jacobi-Verfahren auf der GPU       *\n");
    printf("*******************************************\n\n");



    // Deklaration benoetigter Variablen

    int N; // Problemgroesse
    int max_iter; // Maximale Iterationszahl
    double rel_tol; // Relative Toleranz
    char fill_type; // Art der Befuellung des Startvektors

    double *u; // Iterierte des Jacobi-Verfahrens
    double *u_ghost, *dev_u_ghost; // Iterierte des Jacobi-Verfahrens auf Gitter mit Ghost Layer
    double *f, *dev_f; // Rechte Seite
    double *sol; // Exakte Loesung

    double *dev_res; // Residuum => Nur auf GPU benoetigt!
    double *dev_res_init; // Initiales Residuum => Nur auf GPU benoetigt!

    double *norm_res, *dev_norm_res; // Residuumsnorm
    double *norm_res_init, *dev_norm_res_init; // Residuumsnorm des initalen Residuums

    int iter; // Iterationszahl

    // Zeitnahme
    clock_t start;
    clock_t stop;
    clock_t difference;
    int msec;



    // Nutzereingaben

    printf("Problemgroesse (Vielfaches von %d):  ",blocksize-2);
    scanf("%d",&N);

    // Achtung: N%(blocksize-2) == 0 erforderlich!

    printf("Maximale Iterationszahl: ");
    scanf("%d",&max_iter);

    printf("Relative Toleranz: ");
    scanf("%lf",&rel_tol);

    printf("Startvektor (r,t): ");
    scanf("%c",&fill_type);
    fill_type = getchar();

    printf("\n");



    // Beginn der Zeitnahme
    start = clock();



    // Speicherallokation und Befuellung der Vektoren auf der CPU

    printf("Erstelle Daten auf CPU ...\n");

    // Iterierte
    u = (double *)malloc(N*N*sizeof(double));
    get_u0(u,N,fill_type);

    // Iterierte auf Gitter mit Ghost Layer
    u_ghost = (double *)malloc((N+2)*(N+2)*sizeof(double));
    ghostify(u,u_ghost,N);

    // Rechte Seite
    f = (double *)malloc(N*N*sizeof(double));
    get_RHS(f,N);

    // Exakte Loesung
    sol = (double *)malloc(N*N*sizeof(double));
    get_solution(sol,N);

    // Residuumsnorm
    norm_res = (double *)malloc(sizeof(double));

    // Residuumsnorm des initialen Residuums
    norm_res_init = (double *)malloc(sizeof(double));




    // Speicherallokation im globalen Speicher der GPU

    printf("Allokiere Speicher auf GPU ...\n");

    // Iterierte auf Gitter mit Ghost Layer
    cudaMalloc((void**)&dev_u_ghost,(N+2)*(N+2)*sizeof(double));

    // Rechte Seite
    cudaMalloc((void**)&dev_f,N*N*sizeof(double));

    // Residuum
    cudaMalloc((void**)&dev_res,N*N*sizeof(double));

    // Initiales Residuum
    cudaMalloc((void**)&dev_res_init,N*N*sizeof(double));

    // Residuumsnorm
    cudaMalloc((void**)&dev_norm_res,sizeof(double));

    // Residuumsnorm des initialen Residuums
    cudaMalloc((void**)&dev_norm_res_init,sizeof(double));



    // Kopieren benoetigter Daten in den globalen Speicher der GPU

    printf("Kopiere Daten auf GPU ...\n");

    // Iterierte auf Gitter mit Ghost Layer
    cudaMemcpy(dev_u_ghost,u_ghost,(N+2)*(N+2)*sizeof(double),cudaMemcpyHostToDevice);

    // Rechte Seite
    cudaMemcpy(dev_f,f,N*N*sizeof(double),cudaMemcpyHostToDevice);



    // Baue Gitter von Blocks fuer die GPU-Threads; 2D-Gitter zur Erfassung eines Teilgitters der Diskretisierung

    printf("Baue Gitter fuer die GPU ...\n");

    int nblocks = N/(blocksize-2); // Anzahl Blocks pro Gitterdimension
    dim3 gridDim (nblocks,nblocks);
    dim3 blockDim (blocksize,blocksize);

    // Ggf. neues Gitter zur Normbestimmung ...
    //Bestimme Anzahl rekursiver Aufrufe fuer den Redukton-Kernel
        int l=1;
        int NN=N;
        while(NN>blocksize-2)
        {
                    l+=1;
                    NN=NN/(blocksize-2);
                    // Fuer Reduction2:
                    //NN=NN/(2*blocksize);
            }

            //Anzahl Bloecke je Level
            int *nb=(int*)malloc(l*sizeof(double));
            nb[0]=N/(blocksize-2);
            for (int i=0; i<l; i++)
            {
                    nb[i]=nb[i-1]/(blocksize-2);
                    // Fuer Reduction2:
                    //nb[i]=nb[i-1]/(2*blocksize);
            }

            //Zwei Arrays fuer Ergebnisse
            double *dev_r1;
            cudaMalloc((void**)&dev_r1,nb[0]*sizeof(double));
            double *dev_r2;
            cudaMalloc((void**)&dev_r2,nb[0]*sizeof(double));


        // Ausfuehrung des Jacobi-Verfahrens auf der GPU

        printf("Jacobi-Iterationen ...\n");

        // Bestimmung der Norm des initialen Residuums fuer die relative Abbruchbedingung
        residual<<<gridDim,blockDim>>>(dev_res_init,dev_u_ghost,dev_f,N);
        // norm <<<...,...>>>(dev_res_init,N);

        // Berechnen der Norm

        pmult<<<nb[0],blocksize-2>>>(dev_res_init,dev_res_init);
        // Fuer Reduction2:
        // pmult<<<nb[0]*2,blocksize>>>(dev_res_init,dev_res_init);

        //Reduktion 0, dev_res_init wir blockweise in dev_r1 reduziert
        reduction<<<nb[0],blocksize-2>>>(dev_res_init,dev_r1);
        // Fuer Reduction2:
        // reduction2<<<nb[0],blocksize>>>(dev_res_init,dev_r1);

        //Restliche Reduktionen, immer weniger Bloecke
        //Abwechselnd ist r1 Input und r2 Output und umgekehrt
        for (int i=1; i<l; i++)
        {
                    if (i%2==1)
                    {
                            reduction<<<nb[i],blocksize-2>>>(dev_r1,dev_r2);
                    }
                    else
                    {
                            reduction<<<nb[i],blocksize-2>>>(dev_r2,dev_r1);
                    }
            }

            // Fuer Reduction2:
            //Restliche Reduktionen, immer weniger Bloecke
        //Abwechselnd ist r1 Input und r2 Output und umgekehrt
        /*for (int i=1; i<l; i++)
        {
                    if (nb[i]>0)
                    {
                            if (i%2==1)
                            {
                                    reduction2<<<nb[0],blocksize>>>(dev_r1,dev_r2);
                            }
                            else
                            {
                                    reduction2<<<nb[i],blocksize>>>(dev_r2,dev_r1);
                            }
                    }
                    else
                    {
                            if (i%2==1)
                            {
                                    reduction2<<<nb[i-1]/2,blocksize>>>(dev_r1,dev_r2);
                            }
                            else
                            {
                                    reduction2<<<nb[i-1]/2,blocksize>>>(dev_r2,dev_r1);
                            }
                    }
            }*/


        // Ergebnis von GPU zu Host kopieren
        if (l%2==1)
            {
                    cudaMemcpy(norm_res_init,dev_r1,sizeof(double),cudaMemcpyDeviceToHost);
            }
            else
            {
                    cudaMemcpy(norm_res_init,dev_r2,sizeof(double),cudaMemcpyDeviceToHost);
            }
            printf("Residuumsnorm: %e\n", *norm_res_init);
            //Wurzelziehen vom Eintrag des Vektors, der auf die CPU kopiert wurde
            //sq(norm_res_init);

        // cudaMemcpy(norm_res_init,dev_norm_res_init,sizeof(double),cudaMemcpyDeviceToHost);

        // Jacobi-Schleife
        iter = 0;
        while (true)
        {
                // Ausfuehrung der Jacobi-Iteration
                jacobi<<<gridDim,blockDim>>>(dev_u_ghost,dev_f,N);
                iter += 1;

                // Bestimmung der Residuumsnorm
                residual<<<gridDim,blockDim>>>(dev_res,dev_u_ghost,dev_f,N);
                // norm <<<...,...>>>(dev_res,N);

                pmult<<<nb[0],blocksize-2>>>(dev_res,dev_res);
                // Fuer Reduction2:
                        // pmult<<<nb[0]*2,blocksize>>>(dev_res_init,dev_res_init);


                        //Reduktion 0, dev_res_init wir blockweise in dev_r1 reduziert
                        reduction<<<nb[0],blocksize-2>>>(dev_res,dev_r1);
                        // Fuer Reduction2:
                        // reduction2<<<nb[0],blocksize>>>(dev_res,dev_r1);


                        //Restliche Reduktionen, immer weniger Bloecke
                        //Abwechselnd ist r1 Input und r2 Output und umgekehrt
                        for (int i=1; i<l; i++)
                        {
                                if (i%2==1)
                                {
                                        reduction<<<nb[i],blocksize-2>>>(dev_r1,dev_r2);
                                }
                                else
                                {
                                        reduction<<<nb[i],blocksize-2>>>(dev_r2,dev_r1);
                                }
                        }

                        // Fuer Reduction2:
                        //Restliche Reduktionen, immer weniger Bloecke
                        //Abwechselnd ist r1 Input und r2 Output und umgekehrt
                        /*for (int i=1; i<l; i++)
                        {
                                if (nb[i]>0)
                                {
                                        if (i%2==1)
                                        {
                                                reduction2<<<nb[i],blocksize>>>(dev_r1,dev_r2);
                                        }
                                        else
                                        {
                                                reduction2<<<nb[i],blocksize>>>(dev_r2,dev_r1);
                                        }
                                }
                                else
                                {
                                        if (i%2==1)
                                        {
                                                reduction2<<<nb[i-1]/2,blocksize>>>(dev_r1,dev_r2);
                                        }
                                        else
                                        {
                                                reduction2<<<nb[i-1]/2,blocksize>>>(dev_r2,dev_r1);
                                        }
                                }
                        }*/


                        // Ergebnis von GPU zu Host kopieren
                        if (l%2==1)
                        {
                                cudaMemcpy(norm_res,dev_r1,sizeof(double),cudaMemcpyDeviceToHost);
                        }
                        else
                        {
                                cudaMemcpy(norm_res,dev_r2,sizeof(double),cudaMemcpyDeviceToHost);
                        }
                        //Wurzelziehen vom Eintrag des Vektors, der auf die CPU kopiert wurde
                        //sq(norm_res);
                // cudaMemcpy(norm_res_init,dev_norm_res_init,sizeof(double),cudaMemcpyDeviceToHost);

                // Abbruchbedingung
                if (iter >= max_iter || sqrt((*norm_res))/sqrt((*norm_res_init)) < rel_tol)
                {
                    break;
                }
            }



            // Kopieren des Ergebnisses von globalem Speicher der GPU zurueck auf die CPU

            printf("Kopiere Ergebnis auf CPU ...\n");

            // Iterierte auf Gitter mit Ghost Layer
            cudaMemcpy(u_ghost,dev_u_ghost,(N+2)*(N+2)*sizeof(double),cudaMemcpyDeviceToHost);
            deghostify(u_ghost,u,N);

            printf("\n  => Berechnungen abgeschlossen!\n\n");



            // Ende der Zeitnahme
            stop = clock();
            difference = stop-start;
            msec = difference*1000/CLOCKS_PER_SEC;



            // Ausgabe der Resultate

            printf("Absoluter Fehler: %e\n", get_abs_err(u,sol,N*N));
            printf("Residuumsnorm: %e\n", *norm_res);
            printf("Iterationen: %d\n", iter);
            printf("Benoetigte Zeit: %d.%ds\n\n", msec/1000,msec%1000);
            // Speicherfreigabe auf der CPU
                free(u);
                free(u_ghost);
                free(f);
                free(sol);
                free(norm_res);
                free(norm_res_init);
                free(nb);



                // Speicherfreigabe auf der GPU
                cudaFree(dev_u_ghost);
                cudaFree(dev_f);
                cudaFree(dev_res);
                cudaFree(dev_res_init);
                cudaFree(norm_res);
                cudaFree(norm_res_init);
                cudaFree(dev_r1);
                cudaFree(dev_r2);



                return 0;
            }
