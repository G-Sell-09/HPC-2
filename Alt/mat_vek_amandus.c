//Programm mat_vek auf amandus
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <likwid.h>

/* Funktion zur Berechnung der Matrix-Vektor-Multiplikation Ax einer nxn-Matrix A mit einem nx1-Vektor x */
double * mat_vek_mult (double **A, double *x, int n)
{
    double *res = (double *) malloc(n*sizeof(double));

    for (int i = 0; i < n; i++)
    {
        res[i] = 0.0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i] = res[i] + A[i][j]*x[j];
        }
    }

    return res;
}

double *getRandomVec (int n, double min, double max)
{
	// Zufallsvektor
	double *randomVec = (double *)malloc(n*sizeof(double));

	// Setzen eines Random-Seeds
	srand(time(0));

	//Fuellen des Vektors mit Zufallszahlen
	for (int i = 0; i < n; i++)
	{
		randomVec[i] = min + ((double)rand()/RAND_MAX) * (max - min); 
	}
	
	return randomVec;
}


int main (int argc, char **argv)
{
	LIKWID_MARKER_INIT;
	/* Deklaration der benoetigten Variablen*/
    char eingabe[6];
    char eingabeVektor[3];
    char ausgabe[3];
    int n;
    double *x;
    double *res;
    double **A;
	
	/* Abfangen der Faelle ohne oder mit zu vielen Eingabeparametern */
    if (argc != 1){
        printf("Es wird kein Eingabeparameter verlangt!\n");
        return 1;
    }
	
	/* Abfrage von n */
    printf("Wie groß soll die Dimension sein?\n");
    scanf("%d",&n);
    printf("\n");
	
	x = (double*) malloc (n*sizeof(double));

    A = (double**) malloc (n*sizeof(double*));
    
    res = (double *) malloc(n*sizeof(double));
    

    //Speicher fuer Matrix A allokieren
    for (int i = 0; i < n; i++)
    {
      A[i] = (double*) malloc (n*sizeof(double));
    }


	//Matrix und Vektor initialisieren.

    /* Abfrage Vektor*/
    printf("Soll der Vektor zufaellig sein? (yes, no)\n");
    scanf("%s",eingabeVektor);
    printf("\n");
    
    if (strcmp(eingabeVektor, "yes")==0)
    {		
		//Random Vektor
		double min;
		double max;
		printf("Wertebereich angeben\n");
		printf("Minimum:");
		scanf("%lf",&min);
		printf("Maximum:");
		scanf("%lf",&max);
		printf("\n");
		
		LIKWID_MARKER_START("Erstellung Random Vec");
		x=getRandomVec(n,(double)min,(double)max);
		LIKWID_MARKER_STOP("Erstellung Random Vec");
	}
	else if (strcmp(eingabeVektor, "no")==0)
	{
		//Vektor x mit 1en initialisieren
		for (int i=0; i < n;i++)
		{
			x[i]=1;
		}
	}else{
		printf("Andere Eingabeform waehlen.");
	}
	
	/* Abfrage Matrix Format*/
    printf("Wie soll die Matrix aufgebaut sein? (sparse, dense, random)\n");
    scanf("%s",eingabe);
    printf("\n");
    
	
	if (strcmp(eingabe, "sparse")==0)
    {
		//Sparsematrix erstellen
		LIKWID_MARKER_START("Sparse-Matrix aufstellen");
        for (int i= 0; i<n; i++)
		{
			for (int j= 0; j<n; j++)
			{
				if (j==i)
				{
					A[i][j]=4.0;
				}
				else if (j==i-1)
				{
					A[i][j]=-1.0;
				}else if (j-1==i)
				{
					A[i][j]=-1.0;
				}else{
					A[i][j]=0.0;
				}
	
			}
		}
		LIKWID_MARKER_STOP("Sparse-Matrix aufstellen");
    }
    else if (strcmp(eingabe, "dense")==0)
    {
		//Densematrix erstellen
		LIKWID_MARKER_START("Dense-Matrix aufstellen");
        for (int i= 0; i<n; i++)
		{
			for (int j= 0; j<n; j++)
			{
				if (j==i-1 || j-1==i)
				{
					A[i][j]=0.0;
				}else if (j==i)
				{
					A[i][j]=(double)n*(double)n;
				}else{
					A[i][j]=-1.0*(double)n;
				}
	
			}
		}
		LIKWID_MARKER_STOP("Dense-Matrix aufstellen");      
    }
    else if (strcmp(eingabe, "random")==0)
	{
		// Setzen eines Random-Seeds
		srand(time(0));
		
		double min;
		double max;
		printf("Wertebereich angeben\n");
		printf("Minimum:");
		scanf("%lf",&min);
		printf("Maximum:");
		scanf("%lf",&max);
		printf("\n");
		
		//Random-Matrix erstellen
		LIKWID_MARKER_START("Random-Matrix aufstellen");
        for (int i= 0; i<n; i++)
		{
			for (int j= 0; j<n; j++)
			{
				A[i][j]= min + ((double)rand()/RAND_MAX) * (max - min);	
			}
		}
		LIKWID_MARKER_STOP("Dense-Matrix aufstellen");  
	}else{
		printf("Andere Matrixeigenschaft waehlen.\n");
		return 1;
	}
	  
	// Berechnung
	LIKWID_MARKER_START("Berechnung Mat-Vek-Mult");
	res=mat_vek_mult(A,x,n);
	LIKWID_MARKER_STOP("Berechnung Mat-Vek-Mult");  
	
	
	/* Abfrage, ob Ausgabe*/
    printf("Sollen A, x und b ausgegeben werden? (yes, no)\n");
    scanf("%s",ausgabe);
    printf("\n");
    
    if (strcmp(ausgabe, "yes")==0)
    {	 
		// Ausgabe Matrix
		printf("A = \n");
		for (int i= 0; i<n; i++)
		{
			for (int j= 0; j<n; j++)
			{
				printf("%lf  ", A[i][j]);
			}
			printf("\n");
		}
		printf("\n\n");
	
		//Ausgabe Vektor x
		printf("x = ");
		for (int i = 0; i < n; i++)
		{
			printf("%lf ", x[i]);
		}
		printf("\n\n");
	
		// Ausgabe Loesungsvektor
		printf("Ergebnis Matrix-Vektor-Multiplikation b:\n");
		printf("b = ");
		for (int i = 0; i < n; i++)
		{
			printf(" %lf ", res[i]);
		}
		printf("\n\n");
	 }else{
		 printf("\n");
	 }
	
	LIKWID_MARKER_CLOSE;
	
    return 0;
}
