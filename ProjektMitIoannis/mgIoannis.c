/*!
 * @file mg.c
 * @brief Mehrgitterverfahren zur Loesung eines Randwertproblems
 * \author Ioannis Lilikakis
 * \author Eric Windler
 * @date 2021-06-15
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <sched.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
*	@brief Gauss-Verfahren bzw. LU-Zerlegung fuer Poissonmatrix.
*
*	Für den 5- und 9-Punkte-Stern wird zunächst die entsprechende Matrix aufgestellt und anschließend mit dem Gauss-kij Algorithmus L und U bestimmt.
*	Dabei wird ausgenutzt, dass L+U eine Bandbmatrix ist, wenn die zugrundeliegende Matrix ebenfalls eine Bandmatrix ist.
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] LU Matrix, die die Nichtnull Diagonalen, Spaltenweise, von L und U speichert
*/
void gauss( const int N, const int S, double* restrict LU );

/**
*	@brief Vorwaerts- und Rueckwaertssubstitution.
*
*	Mithilfe der in der Methode gauss berechneten LU-Matrix, wird über L z = b und U x = z die exakte Lösung ermittelt.
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] LU Matrix, in der die Nichtnull Diagonalen, Spaltenweise, von L und U gespeichert hat
*	@param[in] x Speicherort fuer Loesung
*	@param[in] b Rechte Seite des Gleichungssystems
*/
void substitutionen( const int N, const int S, const double* restrict LU, double* restrict x, const double* restrict b );


/**
*	@brief Matrixfreie Multiplikation mit der Poisson-Matrix des 5- oder 9-Punkte-Stern
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] r Zu multiplizierender Vektor
*	@param[in] y Speicherort fuer Loesung
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*/
void mfMult( const int N, const double* restrict r, double* restrict y, const int S );

/**
*	@brief Normberechnung
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] x Zu Normierender Vektor
*	@param[out] norm
*/
double norm( const int N, const double* restrict x );

/**
*	@brief Glätter auf Basis des Jacobi-Verfahren
*
*	Formel: x_{k+1} = x_{k} + w*D^(-1)*( b - Ax_{k} ) mit w = 0.6 \n
*	Hier wird die regelmaessige Struktur der Poissonmatrix ausgenutzt, sodass mithilfe der Methode mfMult matrixfrei gearbeitet werden kann.
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] x Speicherort fuer Loesung
*	@param[in] b Rechte Seite des Gleichungssystems
*	@param[in] maxIter Maximale Iterationsanzahl. Es existieren keine weiteren Abbruchskriterien
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*/
void jacobi( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S );

/**
*	@brief Glaetter auf Basis des Gauss-Seidel-Verfahrens
*
*	Formel: x_{k+1} = Bx_{k} + C mit B = (D − L)^(-1)U und C = ( D − L)^(−1)b \n
*	Hier wird die regelmaessige Struktur der Poissonmatrix ausgenutzt, sodass matrixfrei gerechnet werden kann.
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] x Speicherort fuer Loesung
*	@param[in] b Rechte Seite des Gleichungssystems
*	@param[in] maxIter Maximale Iterationsanzahl. Es existieren keine weiteren Abbruchskriterien
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*/
void gaussSeidel( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S );

/**
*	@brief Glaetter auf Basis des SOR-Verfahrens
*
*	Formel: x_{k+1} = Bx_{k} + C mit B = ( 1/w D − L)^(-1)( (1-w)/w D + U ) und C = ( 1/w D − L)^(−1)b mit w = 1.9 \n
*	Hier wird die regelmaessige Struktur der Poissonmatrix ausgenutzt, sodass matrixfrei gerechnet werden kann.
*
*	@param[in] N Anzahl Gitterpunkte entlang einer Kante
*	@param[in] x Speicherort fuer Loesung
*	@param[in] b Rechte Seite des Gleichungssystems
*	@param[in] maxIter Maximale Iterationsanzahl. Es existieren keine weiteren Abbruchskriterien
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*/
void SOR( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S );

/**
*	@brief linerae Prolongationsmethode
*
*	@param[in] l Level des zu prolongierenden Vektors x
*	@param[in] x Vektor aus groeberem Gitter
*	@param[in] p Zielvektor in feinerem Gitter
*	@param[in] sizeofL Speicherort der verschiedenen Gittergroessen entsprechend der Level
*/
void prol( const int l, const double* restrict x, double* restrict p, const int* restrict sizeofL );

/**
*	@brief lineare full-weight Restriktionsmethode
*
*	@param[in] l Level des zu restringierenden Vektors x0
*	@param[in] x0 Vektor des feineren Gitters
*	@param[in] x1 Zielvektor des groeberen Gitters
*	@param[in] sizeofLx Speicherort der verschiedenen Gittergroessen entsprechend der Level
*/
void rest( const int l, const double* restrict x0, double* restrict x1, const int* restrict sizeofLx );


/**
*	@brief V-Zyklus, iterativ
*
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] glaetter 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
*	@param[in] v1 Iterationsschritte beim Vorglätten
*	@param[in] v2 Iterationsschritte beim Nachglätten
* 	@param[in] sizeofL Speicherort der verschiedenen Gittergroessen entsprechend der Level
* 	@param[in] L Level
* 	@param[in] Lx Array welches auf die bisherigen Approximationen der verschiedenen Leveln zeigt
* 	@param[in] Lb Array welches auf die rechte Seite f(x,y) der verschiedenen Leveln zeigt
* 	@param[in] Lr Array welches auf die Residuen der verschiedenen Leveln zeigt
* 	@param[in] Lp Array welches auf die prolongierten Approximationen der verschiedenen Leveln zeigt
* 	@param[in] LU sparse LU-Zerlegung fuer Gitter in Level 0
*/
void V_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU);

/**
*	@brief W-Zyklus, rekursiv
*
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] glaetter 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
*	@param[in] v1 Iterationsschritte beim Vorglätten
*	@param[in] v2 Iterationsschritte beim Nachglätten
* 	@param[in] sizeofL Speicherort der verschiedenen Gittergroessen entsprechend der Level
* 	@param[in] L Level
* 	@param[in] Lx Array welches auf die bisherigen Approximationen der verschiedenen Leveln zeigt
* 	@param[in] Lb Array welches auf die rechte Seite f(x,y) der verschiedenen Leveln zeigt
* 	@param[in] Lr Array welches auf die Residuen der verschiedenen Leveln zeigt
* 	@param[in] Lp Array welches auf die prolongierten Approximationen der verschiedenen Leveln zeigt
* 	@param[in] LU sparse LU-Zerlegung fuer Gitter in Level 0
*/
void W_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU);

/**
*	@brief F-Zyklus, rekursiv
*
*	Eine von zwei Rekursionsstellen nutzt die Methode V_Zyklus
*
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] glaetter 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
*	@param[in] v1 Iterationsschritte beim Vorglätten
*	@param[in] v2 Iterationsschritte beim Nachglätten
* 	@param[in] sizeofL Speicherort der verschiedenen Gittergroessen entsprechend der Level
* 	@param[in] L Level
* 	@param[in] Lx Array welches auf die bisherigen Approximationen der verschiedenen Leveln zeigt
* 	@param[in] Lb Array welches auf die rechte Seite f(x,y) der verschiedenen Leveln zeigt
* 	@param[in] Lr Array welches auf die Residuen der verschiedenen Leveln zeigt
* 	@param[in] Lp Array welches auf die prolongierten Approximationen der verschiedenen Leveln zeigt
* 	@param[in] LU sparse LU-Zerlegung fuer Gitter in Level 0
*/
void F_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU);

/**
*	@brief Mehrgitterverfahren
*
*	@param[in] eps Zu erreichende Genauigkeit bis Abbruch der Iterration
*	@param[in] S Stencil. Wenn S = 5 oder 9, dann 5 bzw. 9-Punkte-Stern
*	@param[in] zyklus 0 = F-Zyklus; 1 = V-Zyklus; 2 = W-Zyklus
*	@param[in] glaetter 0 = Jacobi, 1 = Gauss-Seidel, 2 = SOR
*	@param[in] v1 Iterationsschritte beim Vorglätten
*	@param[in] v2 Iterationsschritte beim Nachglätten
* 	@param[in] sizeofL Speicherort der verschiedenen Gittergroessen entsprechend der Level
* 	@param[in] L Level
* 	@param[in] Lx Array welches auf die bisherigen Approximationen der verschiedenen Leveln zeigt
* 	@param[in] Lb Array welches auf die rechte Seite f(x,y) der verschiedenen Leveln zeigt
* 	@param[in] Lr Array welches auf die Residuen der verschiedenen Leveln zeigt
* 	@param[in] Lp Array welches auf die prolongierten Approximationen der verschiedenen Leveln zeigt
* 	@param[in] LU sparse LU-Zerlegung fuer Gitter in Level 0
*/
void mg( const double eps, const int S, const int zyklus, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU);

/**
*	@brief Haupteintritspunkt des Programms
*
*	Durch die Befehlszeilenargumente werden Parameter uebergeben. In der main wird überprüft ob die geforderten Werte gültig sind.
*	Anschliessend werden die Parameter Ausgegeben und das Mehrgitterverfahren wird angestossen.
*
*	@param[in] argc Anzahl der Befehlszeilenargumente. Es muessen genau 10 sein.
*	@param[in] argv 0 Auszufuehrende \n
*					1 N Gittergröße \n
*					2 L Anzahl Level \n
*					3 S Stencil 5 oder 9 \n
*					4 zyklus 0 = F-Zyklus; 1 = V-Zyklus; 2 = W-Zyklus \n
*					5 glaetter 0 = Jacobi; 1 = Gauss-Seidel; 2 = SOR \n
*					6 v1 Iterationsschritte beim Vorglätten \n
*					7 v2 Iterationsschritte beim Nachglätten \n
*					8 err, Fehlergenauigkeit eps = 1e-"err", mit err in {1,2,3,...} \n
*					9 Ausgabe der Loesung 0 = Nein; 1 = Ja
*/
int main ( int argc , char **argv )
{
	int N = atoi(argv[1]);
	int L = atoi(argv[2]);
	int S = atoi(argv[3]);
	int zyklus = atoi(argv[4]);
	int glaetter = atoi(argv[5]);
	int v1 = atoi(argv[6]);
	int v2 = atoi(argv[7]);

	if ( argc < 10)
	{
		printf("Die Anzahl der uebergebenen Argumente ist nicht ausreichend.\n");
		return 0;
	}

	if ( argc > 10)
	{
		printf("Die Anzahl der uebergebenen Argumente ist zu gross.\n");
		return 0;
	}

	if ( L < 2 || N < 2 )
	{
		printf("L und N muessen groesser 1 gewaehlt werden.\n");
		return 0;
	}

	int n = N;
	for (int i = 0; i < L; i++)
	{
		if ( n%2 == 0){
			printf("N muss auf allen Leveln ungerade sein.\n");
			return 0;
		}else{
			n = (n-1)/2;
		}
	}

	if ( S != 5 && S != 9){
		printf("Waehle 9 oder 5 Punkte Stern.\n");
		return 0;
	}

	//          0=F            V=1            W=2
 	if ( zyklus != 0 && zyklus != 1 && zyklus != 2 )
	{
		printf("Ein V-, W-, oder F-Zyklus muss gewaehlt werden.\nBitte gib entsprechend 1, 2 oder 0 ein.\n");
		return 0;
	}
	//       0=Jacobi      1=gaussSeide       2=SOR
	if ( glaetter != 0 && glaetter != 1 && glaetter != 2)
	{
		printf("Ein Jacobi, Gauss Seidel, oder SOR Glaetter muss gewaehlt werden.\nBitte gib entsprechend 0, 1 oder 2 ein.\n");
		return 0;
	}

	if ( v1 < 0 || v2 < 0)
	{
		printf("Die Anzahl der jeweiligen Vor- und Nachglaettungen duerfen nicht negativ sein.\n");
		return 0;
	}

	if ( atoi(argv[8]) < 1 )
	{
		printf("Der Fehler muss groesser 1 und damit kleiner 1e-1 gewaehlt werden.\n");
		return 0;
	}


	double eps = 1;
	for (int i = 0; i < atoi(argv[8]); i++ )
		eps = eps/10.0;

	double** Lx  = malloc( L * sizeof(double*) );
	double** Lb  = malloc( L * sizeof(double*) );
	double** Lr  = malloc( L * sizeof(double*) );
	double** Lp  = malloc( L * sizeof(double*) );
	int* sizeofL = malloc ( L * sizeof(int) );
	int size = N;
	for (int l = L-1; l >= 0; l--)
	{
		Lx[ l ]       = malloc( size*size * sizeof(double) );
		Lb[ l ]       = malloc( size*size * sizeof(double) );
		Lr[ l ]       = malloc( size*size * sizeof(double) );
		Lp[ l ]       = malloc( size*size * sizeof(double) );
		sizeofL[ l ]  = size;
		size          = (size - 1) / 2;
	}


	// start vector
	for (int i = 0; i < N*N; i++)
		Lx[ L-1 ][i] = 1;
	for (int l = L-2; l >= 0; l--)
		for (int i = 0; i < sizeofL[l]*sizeofL[l] ; i++)
			Lx[ l ][ i ] = 0;


	double h  = 1./(double)(N+1);
	double hh;
	if ( S == 5 )
		hh = h*h;
	if ( S == 9 )
		hh = 6*h*h;

	// f(x,y) = 2 * pi^2 * sin(pi * x) sin(pi * y)
	// b_{iN+j} = hh*f_ji = hh*f(jh,ih)
	for ( int i = 0; i < N; i++  )
		for ( int j = 0; j < N ; j++ )
			Lb[L-1][i*N + j] = hh*2*M_PI*M_PI * sin( M_PI * (j+1)*h ) * sin( M_PI * (i+1)*h  );


	int sizeLU = 0;
	if ( S == 5 )
		sizeLU = sizeofL[0]*sizeofL[0]*2*sizeofL[0];
	if ( S == 9 )
		sizeLU = sizeofL[0]*sizeofL[0]*(2*sizeofL[0]+3);
	double* LU = malloc( sizeLU * sizeof(double) );


	printf("Parameter\n");
	printf("-------------------------------------------------------------\n");
	if ( S == 5)
		printf("Stern:                          5-Punkte-Stern\n");
	else
		printf("Stern:                          9-Punkte-Stern\n");
	printf("Gittergroesse:                  %i x %i\n", N, N);
	printf("Level:                          %i\n", L);
	if ( zyklus == 0)
		printf("Zyklus:                         F-Zyklus\n");
	else if ( zyklus == 1)
		printf("Zyklus:                         V-Zyklus\n");
	else
		printf("Zyklus:                         W-Zyklus\n");
	if ( glaetter == 0)
		printf("Glaetter:                       Jacobi\n");
	else if ( glaetter == 1)
		printf("Glaetter:                       Gauss Seidel\n");
	else
		printf("Glaetter:                       SOR\n");
	printf("Vorglaettungen:                 %i\n", v1);
	printf("Nachglaettungen:                %i\n", v2);
	printf("Angestrebte Genauigkeit:        %.2e\n", eps);

	printf("\nErgebnis\n");
	printf("-------------------------------------------------------------\n");
	printf("LU-Faktorisierung\n");
	printf("Groebste Gittergroesse:         %i x %i\n", sizeofL[0], sizeofL[0]);
	clock_t gauss_start = clock();
		gauss(sizeofL[0], S, LU);
	clock_t gauss_end = clock();
	double gauss_seconds = (double)(gauss_end - gauss_start) / CLOCKS_PER_SEC;
	printf("Laufzeit:                       %f Sekunden\n", gauss_seconds );

	printf("\nMehrgitterverfahren\n");
	clock_t start = clock();
		mg( eps, S, zyklus, glaetter, v1, v2, sizeofL, L, Lx, Lb, Lr, Lp, LU);
	clock_t end = clock();
	double seconds = (double)(end - start) / CLOCKS_PER_SEC;

	printf("Laufzeit:                       %f Sekunden\n\n", seconds );
	printf("Laufzeit insgesamt:             %f Sekunden\n", (seconds + gauss_seconds) );

	if( atoi(argv[9]) )
	{
		printf("Loesung\n");
		printf("approx.  &&  analytisch\n");
		for ( int i = 0; i < N; i++  )
			for ( int j = 0; j < N ; j++  )
			{
				printf("%.4lf     ", Lx[L-1][i*N +j] );
				printf("%.4lf\n", sin( M_PI*(j+1)*h ) * sin( M_PI*(i+1)*h  ) ); //u_ij = u(ih, jh)
			}
	}

	return 0;

}

void gauss( const int N, const int S, double* restrict LU )
{
    // Nicht gespeichert:
    // L hat 1 auf der Diagonalen
    // U hat -1 auf (N+1)ter Nebendiagonale
    // Ersten N gehoeren zu L und die letzten N zu U

    int i, ii, j, k, Bk, bk;
    double a_ik;

	if ( S == 5 )
	{
		// set startvalue LU = sparse(A)
		for ( i = 0; i < N*N; i++)
		{
			LU[ i*2*N       ] = -1; // L
			for ( j = 1; j < N-1; j++)
				LU[ i*2*N + j ] = 0; // L
			LU[ i*2*N + N-1 ] = -1; // L
			LU[ i*2*N + N   ] =  4; // U
			LU[ i*2*N + N+1 ] = -1; // U
			for ( j = N+2; j < 2*N; j++)
				LU[ i*2*N + j ] =  0; // U
			//LU[ i*N + 2*N ] = -1; // U
		}

		for ( i = 0; i < N; i++)
			for ( j = 0; j < N-i; j++)
				LU[ i*2*N + j ] = 0; // L
		LU[ (N*N-1)*2*N + N + 1 ] = 0; // U

		for ( i = 0; i < N*N; i+=N)
		{
			LU[    i   *2*N + N-1 ] = 0; // L
			LU[ (i+N-1)*2*N + N+1 ] = 0; // U
		}


		// sparse Gauss kij

		for ( Bk = 0; Bk < N-1; Bk++) // entlang der Bloecke
		{

			for ( bk = 0; bk < N; bk++) // innerhalb des Blockes Bk
			{
				k = Bk*N + bk; // globaler = lokaler Index

				for ( i = k+1; i < k+1+N; i++)
				{

					ii = i - k; // Abweichung nach Links von der Diagonalen (Spalte N) in LU

					a_ik = LU[i*2*N + N-ii] / LU[k*2*N + N];
					LU[i*2*N + N-ii] = a_ik;


					for ( j = N; j < 2*N-1; j++)
					{
						LU[i*2*N + j-ii+1]  = LU[i*2*N + j-ii+1] -  a_ik * LU[ k*2*N + j+1] ;
					}
					LU[i*2*N + 2*N-ii]  = LU[i*2*N + 2*N-ii] + a_ik; // LU[ k*2*N + 2*N] existiert nicht immer -1


				}
			}

		}
		// letzter Block
			for ( bk = 0; bk < N-1; bk++)
			{
				k = (N-1)*N + bk;

				for ( i = k+1; i < k-bk+N; i++)
				{

					ii = i - k;

					a_ik = LU[i*2*N + N-ii] / LU[k*2*N + N];
					LU[i*2*N + N-ii] = a_ik;

					for ( j = N; j < 2*N-1-bk; j++)
					{
						LU[i*2*N + j-ii+1]  = LU[i*2*N + j-ii+1] -  a_ik * LU[ k*2*N + j+1] ;
					}

				}
			}
		return;
	}

	if (S == 9)
    {

		for( i = 0; i < N*N; i++)
        {
            LU[i*(2*N+3)]           = -1;
            LU[i*(2*N+3)+1]         = -4;
            LU[i*(2*N+3)+2]         = -1;
            for( j = 3; j < N; j++)
                LU[i*(2*N+3)+j] = 0;
            LU[i*(2*N+3)+N]         = -4;
            LU[i*(2*N+3)+N+1]       = 20;
            LU[i*(2*N+3)+N+2]       = -4;
            for( j = N+3; j < 2*N; j++)
                LU[i*(2*N+3)+N+3+j] = 0;
            LU[i*(2*N+3)+2*N]       = -1;
            LU[i*(2*N+3)+2*N+1]     = -4;
            LU[i*(2*N+3)+2*N+2]     = -1;
        }

        //Oberer linker Block Null setzen
        for( i = 0; i < N; i++)
            for( j = 0; j < N; j++)
                LU[i*(2*N+3)+j] = 0;

        //Unterer rechter Block Null setzen
        for( i = (N-1)*N; i < N*N; i++)
            for( j = 2*N; j < 2*N+3; j++)
                LU[i*(2*N+3)+j] = 0;

        for( i = 0; i < N*N; i += N)
        {
            //Nullen im linken Blockrand
            LU[i*(2*N+3)] = 0;
            LU[(i+N-1)*(2*N+3) +2] = 0;

            //Nullen im rechten Blockrand
            LU[(i+1)*(2*N+3) - 3] = 0;
            LU[(i+N)*(2*N+3) - 1] = 0;

            //Nullen im mittleren Block
            LU[i*(2*N+3) + N] = 0;
            LU[(i+N)*(2*N+3) - N -1 ] = 0;
        }

        for ( Bk = 0; Bk < N-2; Bk++) // entlang der Bloecke
        {
            for ( bk = 0; bk < N; bk++) // innerhalb des Blockes Bk
            {
                k = Bk*N + bk; // globaler = lokaler Index
				// 4
                for ( i = k+1; i < k+2+N; i++) // 5
                {
					// 1
                    ii = i - k; // Abweichung nach Links von der Diagonalen in LU

                    a_ik = LU[i*(2*N+3) + N+1-ii] / LU[k*(2*N+3) + N+1];
                    LU[i*(2*N+3) + N+1-ii] = a_ik;

					// Betrachte nur Eintraege rechts neben der Diagonalen in Zeile k
                    for ( j = N+1; j < 2*N+2; j++)
                    {
                        LU[i*(2*N+3) + j-ii+1]  = LU[i*(2*N+3) + j-ii+1] -  a_ik * LU[ k*(2*N+3) + j+1] ;
                    }


                }
            }

        }

		// vorletzter Block
		for ( bk = 0; bk < N-1; bk++) // innerhalb des Blockes Bk
		{
			k = (N-2)*N + bk; // globaler = lokaler Index
			// 4
			for ( i = k+1; i < k+2+N; i++) // 5
			{
				// 1
				ii = i - k; // Abweichung nach Links von der Diagonalen in LU

				a_ik = LU[i*(2*N+3) + N+1-ii] / LU[k*(2*N+3) + N+1];
				LU[i*(2*N+3) + N+1-ii] = a_ik;

				// Betrachte nur Eintraege rechts neben der Diagonalen in Zeile k
				for ( j = N+1; j < 2*N+2; j++)
				{
					LU[i*(2*N+3) + j-ii+1]  = LU[i*(2*N+3) + j-ii+1] -  a_ik * LU[ k*(2*N+3) + j+1] ;
				}

			}
		}
		// vorletzter Block letzte Zeile
			k = (N-2)*N + N-1; // globaler = lokaler Index
			for ( i = k+1; i < k+1+N; i++)
			{
				ii = i - k;

				a_ik = LU[i*(2*N+3) + N+1-ii] / LU[k*(2*N+3) + N+1];
				LU[i*(2*N+3) + N+1-ii] = a_ik;

				for ( j = N+1; j < 2*N+1; j++)
				{
					LU[i*(2*N+3) + j-ii+1]  = LU[i*(2*N+3) + j-ii+1] -  a_ik * LU[ k*(2*N+3) + j+1] ;
				}
			}


		// letzter Block

			for ( bk = 0; bk < N; bk++) // innerhalb des Blockes Bk
			{
				k = (N-1)*N + bk; // globaler = lokaler Index

				for ( i = k+1; i < k-bk+N; i++)
				{
					ii = i - k;

					a_ik = LU[i*(2*N+3) + N+1-ii] / LU[k*(2*N+3) + N+1];
					LU[i*(2*N+3) + N+1-ii] = a_ik;

					for ( j = N+1; j < 2*N-bk; j++)
					{
						LU[i*(2*N+3) + j-ii+1]  = LU[i*(2*N+3) + j-ii+1] -  a_ik * LU[ k*(2*N+3) + j+1] ;
					}


				}
			}

    }



}

void substitutionen( const int N, const int S, const double* restrict LU, double* restrict x, const double* restrict b )
{
	//substitionen
	int i, ii, j;
    double* z = malloc(N*N*sizeof(double));
    double tmp;

    // LU x = b
	if( S == 5)
	{

		// L z = b
		z[0] = b[0];
		for ( i = 1; i < N; i++)
		{
			tmp = 0;
			for ( j = 1; j < i+1; j++) // Abweichung nach Links von der Diagonalen (Spalte N) in LU
				tmp += LU[ i*2*N + N-j ] * z[i-j] ;

			z[i] = b[i] - tmp  ;
		}
		for ( i = N; i < N*N; i++)
		{
			tmp = 0;
			for (int j = 1; j < N+1; j++) // lokal
				tmp += LU[ i*2*N + N-j ] * z[i-j] ;

			z[i] = b[i] - tmp;
		}

		// U x = z
		x[N*N-1] = z[N*N-1] / LU[ (N*N-1)*2*N + N ];

		for ( i = N*N-2; i >= N*N-N; i--)
		{
			tmp = 0;
			ii = N*N - i;
			for ( j = 1; j < ii; j++)
				tmp += LU[ i*2*N + N+j ] * x[i+j] ;

			x[i] = (z[i] - tmp ) / LU[ i*2*N + N ];

		}

		for ( i = N*N-N-1; i >= 0; i--)
		{
			tmp = 0;
			for ( j = 1; j < N; j++)
				tmp += LU[ i*2*N + N+j ] * x[i+j] ;

			x[i] = (z[i] - tmp + x[i+N]) / LU[ i*2*N + N ];

		}

		return;
	}


	if( S == 9)
	{

		// L z = b
		z[0] = b[0];
		for ( i = 1; i < N+1; i++)
		{
			tmp = 0;
			for ( j = 1; j < i+1; j++) // Abweichung nach Links von der Diagonalen (Spalte N+1) in LU
				tmp += LU[ i*(2*N+3) + N+1-j ] * z[i-j] ;

			z[i] = b[i] - tmp  ;
		}
		for ( i = N+1; i < N*N; i++)
		{
			tmp = 0;
			for (int j = 1; j < N+2; j++) // lokal
				tmp += LU[ i*(2*N+3) + N+1-j ] * z[i-j] ;

			z[i] = b[i] - tmp;
		}


		// U x = z
		x[N*N-1] = z[N*N-1] / LU[ (N*N-1)*(2*N+3) + N+1 ];

		for ( i = N*N-2; i >= N*N-N-1; i--)
		{
			tmp = 0;
			ii = N*N - i;
			for ( j = 1; j < ii+1; j++)
				tmp += LU[ i*(2*N+3) + N+1+j ] * x[i+j] ;

			x[i] = (z[i] - tmp ) / LU[ i*(2*N+3) + N+1 ];

		}

		for ( i = N*N-N-2; i >= 0; i--)
		{
			tmp = 0;
			for ( j = 1; j < N+2; j++)
				tmp += LU[ i*(2*N+3) + N+1+j ] * x[i+j] ;

			x[i] = (z[i] - tmp) / LU[ i*(2*N+3) + N+1 ];

		}

	}

}

void mfMult( const int N, const double* restrict r, double* restrict y, const int S )
{
	int i, j;
	// A * r = y
	if ( S == 5 )
	{

	// core
	for (  i = 1; i < N-1; i++ )
		for (  j = 1; j < N-1; j++ )
			y[i*N + j] = 4*r[i*N + j] - r[i*N + j-1] - r[i*N + j+1] - r[(i-1)*N + j] - r[(i+1)*N + j];

	// margins
	for (  i = 1; i < N-1; i++ )
	{
		y[i]           = 4*r[i]           - r[i-1]           - r[i+1]           - r[N + i];           // lower row
		y[(N-1)*N + i] = 4*r[(N-1)*N + i] - r[(N-1)*N + i-1] - r[(N-1)*N + i+1] - r[(N-2)*N + i];     // upper row
		y[i*N]         = 4*r[i*N]         - r[i*N + 1]       - r[(i-1)*N]       - r[(i+1)*N];         // left column
		y[i*N + N-1]   = 4*r[i*N + N-1]   - r[i*N + N-2]     - r[(i-1)*N + N-1] - r[(i+1)*N + N-1];   // right column
	}
	// corners
	y[0]		 = 4*r[0]              	  - r[N]             - r[1];  				// bottom left
	y[N-1]		 = 4*r[N-1]          	  - r[2*N-1]         - r[N - 2];			// bottom right
	y[(N-1)*N] 	 = 4*r[(N-1)*N]     	  - r[(N-1)*N + 1]   - r[(N-2)*N];			// top left
	y[(N-1)*N + N-1] = 4*r[(N-1)*N + N-1] - r[(N-1)*N + N-2] - r[(N-2)*N + N-1];	// top right


	return;
	}

	if ( S == 9 )
	{
		// core
		for (  i = 1; i < N-1; i++ )
			for (  j = 1; j < N-1; j++ )
				y[i*N + j] = 20*r[i*N + j] - 4*( r[i*N + j-1] + r[i*N + j+1] + r[(i-1)*N + j] + r[(i+1)*N + j] )
										- r[(i+1)*N + j-1] - r[(i+1)*N + j+1] - r[(i-1)*N + j-1] - r[(i-1)*N + j+1];

		// margins
		for (  i = 1; i < N-1; i++ )
		{
			y[i]           = 20*r[i]     	   - 4*( r[i-1]          + r[i+1]           + r[N + i]) 	    - r[N + i-1]       - r[N + i+1];           // lower row
			y[(N-1)*N + i] = 20*r[(N-1)*N + i] - 4*(r[(N-1)*N + i-1] + r[(N-1)*N + i+1] + r[(N-2)*N + i])   - r[(N-2)*N + i-1] - r[(N-2)*N + i+1];     // upper row
			y[i*N]         = 20*r[i*N]         - 4*(r[i*N + 1]       + r[(i-1)*N]       + r[(i+1)*N]) 	    - r[(i+1)*N + 1]   - r[(i-1)*N + 1] ;      // left column
			y[i*N + N-1]   = 20*r[i*N + N-1]   - 4*(r[i*N + N-2]     + r[(i-1)*N + N-1] + r[(i+1)*N + N-1]) - r[(i+1)*N + N-2] - r[(i-1)*N + N-2];     // right column
		}

		// corners
		y[0]		 	 = 20*r[0]             - 4*(r[N]             + r[1])			 - r[N+1];  			// bottom left
		y[N-1]		     = 20*r[N-1]           - 4*(r[2*N-1]         + r[N - 2])		 - r[2*N-2];			// bottom right
		y[(N-1)*N] 	 	 = 20*r[(N-1)*N]       - 4*(r[(N-1)*N + 1]   + r[(N-2)*N])		 - r[(N-2)*N + 1];		// top left
		y[(N-1)*N + N-1] = 20*r[(N-1)*N + N-1] - 4*(r[(N-1)*N + N-2] + r[(N-2)*N + N-1]) - r[(N-2)*N + N-2];	// top right
	}


}

double norm( const int N, const double* restrict x )
{
	double norm = 0;
	for ( int i = 0; i < N*N; i++ )
		norm += x[i]*x[i];
	norm = sqrt(norm);
	return norm;
}

void jacobi( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S )
{

	double* y  = malloc(N*N*sizeof(double));
	double w = 0.6;

	if ( S == 5 )
		w = w * 0.25;
	if ( S == 9 )
		w = w * 0.05;

	mfMult( N, x, y, S);
	for (int i = 0; i < N*N; i++)
		y[i] = b[i] - y[i];

	for (int iter = 1; iter < maxIter; iter++)
	{

		// x = x + w*D^(-1)*( b - Ax )
		// x = x + w/4 * ( b - Ax ) // 5 Punkte Stern
		for (int i = 0; i < N*N; i++)
			x[i] = x[i] + w * y[i];

		mfMult( N, x, y, S); // fuer naechste Iteration
		for (int i = 0; i < N*N; i++)
			y[i] = b[i] - y[i];

	}

	free(y);

}

void gaussSeidel( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S )
{

	double* x0  = malloc(N*N*sizeof(double));

	if (S == 5)
	{
		const double w1 = 0.25;

		for (int iter = 1; iter < maxIter; iter++)
		{

			for (int i = 0; i < N*N; i++)
				x0[i] = x[i];

			// i = 0
			x[0] = w1 * ( b[0] + x[1]+x[N] );
			for (int i = 1; i < N-1; i++)
				x[i] = w1 * ( b[i]  + x[i-1] + x0[i+1]+x0[i+N] );
			// i = N-1
			x[N-1] = w1 * ( b[N-1] + x[N-2]  + x0[2*N-1] );

			for (int i = 1; i < N-1; i++)
			{
				// j=0
				x[i*N] = w1 * ( b[i*N] + x[i*N-N] + x0[i*N+1]+x0[i*N+N] );

				for (int j = 1; j < N-1; j++)
					x[i*N + j] = w1 * ( b[i*N + j] + x[i*N + j-1]+x[i*N + j-N] + x0[i*N + j+1]+x0[i*N + j+N] );

				//j=N-1
				x[i*N + N-1] = w1 * ( b[i*N + N-1] + x[i*N + N-2]+x[i*N -1] + x0[(i+2)*N-1] );
			}

			// i = N*N-N
			x[(N-1)*N] = w1 * ( b[(N-1)*N] + x[(N-2)*N] + x0[(N-1)*N+1] );
			for (int i = (N-1)*N+1; i < N*N-1; i++)
				x[i] = w1 * ( b[i]  + x[i-1]+x[i-N] + x0[i+1] );
			// i = N*N-1
			x[N*N-1] = w1 * ( b[N*N-1] + x[N*N-2]+ x[(N-1)*N-1] );


		}
		return;
	}

	if (S == 9)
    {
        const double w1 = 0.05;

        for (int iter = 1; iter < maxIter; iter++)
            {
            // i = 0
            x[0] = w1 * (b[0] + 4*(x[1] + x[N]) + x[N+1]);
            // i = 1,...,N-2
            for (int i = 1; i < N-1; i++)
                x[i] = w1 * (b[i] + 4*(x[i-1] + x[i+1] + x[i+N]) + x[i+N-1] + x[i+N+1]);
            // i = N-2
            x[N-1] = w1 * ( b[N-1] + 4*(x[N-2]  + x[2*N-1]) + x[2*N-2] );

            for (int i = 1; i < N-1; i++)
            {
                // j=0
                x[i*N] = w1 * ( b[i*N] + 4*(x[(i-1)*N] + x[i*N+1]+x[(i+1)*N]) + x[(i-1)*N+1] + x[(i+1)*N +1] );

                //Core
                for (int j = 1; j < N-1; j++)
                    x[i*N + j] = w1 * ( b[i*N + j] + 4*(x[i*N + j-1] + x[(i-1)*N + j] + x[i*N + j+1]+x[(i+1)*N + j]) + x[(i-1)*N + j-1] + x[(i-1)*N + j+1] + x[(i+1)*N + j-1] + x[(i+1)*N + j+1] );

                //j=N-1
                x[(i+1)*N-1] = w1 * ( b[(i+1)*N-1] + 4*(x[(i+1)*N-2]+ x[i*N -1] + x[(i+2)*N-1]) + x[i*N-2] + x[(i+2)*N-2]);
            }

            // i = N*N-N
            x[(N-1)*N] = w1 * ( b[(N-1)*N] + 4*(x[(N-2)*N] + x[(N-1)*N+1]) + x[(N-2)*N + 1] );
            // i = (N-1)*N+1,..., N*N - 2
            for (int i = (N-1)*N+1; i < N*N-1; i++)
                x[i] = w1 * ( b[i]  + 4*(x[i-1] + x[i-N] + x[i+1]) + x[i-N-1] + x[i-N+1] );
            // i = N*N-1
            x[N*N-1] = w1 * ( b[N*N-1] + 4*(x[N*N-2]+ x[(N-1)*N-1]) + x[(N-1)*N-2] );
        }
    }



}

void SOR( const int N, double* restrict x, const double* restrict b, const int maxIter, const int S )
{

	if (S == 5)
	{
		// omega = 1.9
		const double w1 = 1.9 * 0.25;
		const double w2 = 1 - 1.9;

		for (int iter = 1; iter < maxIter; iter++)
		{

			// i = 0
			x[0] = w1 * ( b[0] + x[1]+x[N] ) + w2*x[0];
			for (int i = 1; i < N-1; i++)
				x[i] = w1 * ( b[i]  + x[i-1] + x[i+1]+x[i+N] ) + w2*x[i];
			// i = N-1
			x[N-1] = w1 * ( b[N-1] + x[N-2]  + x[2*N-1] ) + w2*x[N-1];

			for (int i = 1; i < N-1; i++)
			{
				// j=0
				x[i*N] = w1 * ( b[i*N] + x[i*N-N] + x[i*N+1]+x[i*N+N] ) + w2*x[i*N];

				for (int j = 1; j < N-1; j++)
					x[i*N + j] = w1 * ( b[i*N + j] + x[i*N + j-1]+x[i*N + j-N] + x[i*N + j+1]+x[i*N + j+N] ) + w2*x[i*N + j];

				//j=N-1
				x[i*N + N-1] = w1 * ( b[i*N + N-1] + x[i*N + N-2]+x[i*N -1] + x[(i+2)*N-1] ) + w2*x[i*N + N-1];
			}

			// i = N*N-N
			x[(N-1)*N] = w1 * ( b[(N-1)*N] + x[(N-2)*N] + x[(N-1)*N+1] ) + w2*x[(N-1)*N];
			for (int i = (N-1)*N+1; i < N*N-1; i++)
				x[i] = w1 * ( b[i]  + x[i-1]+x[i-N] + x[i+1] ) + w2*x[i];
			// i = N*N-1
			x[N*N-1] = w1 * ( b[N*N-1] + x[N*N-2]+ x[(N-1)*N-1] ) + w2*x[N*N-1];


		}
		return;
	}

	if (S == 9)
    {
		// omega = 1.9
        const double w1 = 1.9 * 0.05;
        const double w2 = 1 - 1.9;

        for (int iter = 1; iter < maxIter; iter++)
            {
            // i = 0
            x[0] = w1 * (b[0] + 4*(x[1] + x[N]) + x[N+1]) + w2*x[0];
            // i = 1,...,N-2
            for (int i = 1; i < N-1; i++)
                x[i] = w1 * (b[i] + 4*(x[i-1] + x[i+1] + x[i+N]) + x[i+N-1] + x[i+N+1]) + w2 * x[i];
            // i = N-2
            x[N-1] = w1 * ( b[N-1] + 4*(x[N-2]  + x[2*N-1]) + x[2*N-2] ) + w2*x[N-1];

            for (int i = 1; i < N-1; i++)
            {
                // j=0
                x[i*N] = w1 * ( b[i*N] + 4*(x[(i-1)*N] + x[i*N+1]+x[(i+1)*N]) + x[(i-1)*N+1] + x[(i+1)*N +1] ) + w2*x[i*N];

                //Core
                for (int j = 1; j < N-1; j++)
                    x[i*N + j] = w1 * ( b[i*N + j] + 4*(x[i*N + j-1] + x[(i-1)*N + j] + x[i*N + j+1]+x[(i+1)*N + j]) + x[(i-1)*N + j-1] + x[(i-1)*N + j+1] + x[(i+1)*N + j-1] + x[(i+1)*N + j+1] ) + w2*x[i*N + j];

                //j=N-1
                x[(i+1)*N-1] = w1 * ( b[(i+1)*N-1] + 4*(x[(i+1)*N-2]+ x[i*N -1] + x[(i+2)*N-1]) + x[i*N-2] + x[(i+2)*N-2]) + w2*x[i*N + N-1];
            }

            // i = N*N-N
            x[(N-1)*N] = w1 * ( b[(N-1)*N] + 4*(x[(N-2)*N] + x[(N-1)*N+1]) + x[(N-2)*N + 1] ) + w2*x[(N-1)*N];
            // i = (N-1)*N+1,..., N*N - 2
            for (int i = (N-1)*N+1; i < N*N-1; i++)
                x[i] = w1 * ( b[i]  + 4*(x[i-1] + x[i-N] + x[i+1]) + x[i-N-1] + x[i-N+1] ) + w2*x[i];
            // i = N*N-1
            x[N*N-1] = w1 * ( b[N*N-1] + 4*(x[N*N-2]+ x[(N-1)*N-1]) + x[(N-1)*N-2] ) + w2*x[N*N-1];
        }
    }



}

void rest( const int l, const double* restrict x0, double* restrict x1, const int* restrict sizeofLx )
{
	int N0 = sizeofLx[l];   // size of x0
	int N1 = sizeofLx[l-1]; // size of x1

	for (int i = 0; i < N1; i++)
	{
		for (int j = 0; j < N1; j++)
		{

			x1[ i*N1 + j ] = 0.0625 * ( 4 *   x0[ i*2*N0+N0      + 2*j+1 ] +
										2 * ( x0[ i*2*N0+N0      + 2*j   ] + x0[ i*2*N0+N0   + 2*j+2 ]
											+ x0[ i*2*N0         + 2*j+1 ] + x0[ i*2*N0+2*N0 + 2*j+1 ]) +
											  x0[ i*2*N0+2*N0    + 2*j   ] + x0[ i*2*N0+2*N0 + 2*j+2 ]
											+ x0[ i*2*N0         + 2*j   ] + x0[ i*2*N0      + 2*j+2 ]);

			//x1[ i*N1 + j ] = x0[ i*2*N0+N0      + 2*j+1 ]; // einfache

		}
	}

}

void prol( const int l, const double* restrict x, double* restrict p, const int* restrict sizeofL )
{
	const int N0 = sizeofL[l];   // size of x
	const int N1 = sizeofL[l+1]; // size of p

	// Teil 1
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
			p[ i*2*N1+N1 + 2*j+1 ]  = x[ i*N0 + j ];

	// Teil 2
	for (int i = 0; i < N0; i++)
	{
		// j=0
			p[ i*2*N1+N1        ]  = 0.5 *                      x[ i*N0 ] ;
		// j= N0
			p[ i*2*N1+N1 + 2*N0 ]  = 0.5 *   x[ i*N0 + N0-1 ];

		for (int j = 1; j < N0; j++)
			p[ i*2*N1+N1 + 2 *j ]  = 0.5 * ( x[ i*N0 + j-1 ] + x[ i*N0 + j ] );
	}

	// Teil 3
	// untere Reihe
		for (int j = 0; j < N0; j++)
			p[            2*j+1 ]  = 0.5 *   x[              j ] ;
	// obere Reihe
		for (int j = 0; j < N0; j++)
			p[ N1*N1-N1 + 2*j+1 ]  = 0.5 *   x[ N0*N0 - N0 + j ] ;

	for (int i = 1; i < N0; i++)
		for (int j = 0; j < N0; j++)
			p[ i*2*N1      + 2*j+1 ]  = 0.5 * ( x[ i*N0-N0 + j ] + x[ i*N0 + j ] );


	// Teil 4
	//Ecken
	p[ 0 ]         = 0.25 * x[ 0 ];
	p[ N1-1 ]      = 0.25 * x[ N0-1 ];
	p[ N1*N1-N1 ]  = 0.25 * x[ N0*N0-N0 ];
	p[ N1*N1- 1 ]  = 0.25 * x[ N0*N0- 1 ];
	//linke Spalte
	for (int i = 0; i < N0-1; i++)
		p[ i*2*N1+2*N1  ]       = 0.25 * ( x[ i*N0        ] + x[ i*N0+N0        ] );
	//rechte Spalte
	for (int i = 0; i < N0-1; i++)
		p[ i*2*N1+2*N1 + N1-1 ] = 0.25 * ( x[ i*N0 + N0-1 ] + x[ i*N0+N0 + N0-1 ] );
	//untere reihe
	for (int j = 0; j < N0-1; j++)
		p[            2*j+2 ]   = 0.25 * ( x[            j ] + x[            j+1 ] );
	//obere reihe
	for (int j = 0; j < N0-1; j++)
		p[ N1*N1-N1 + 2*j+2 ]   = 0.25 * ( x[ N0*N0-N0 + j ] + x[ N0*N0-N0 + j+1 ] );
	// Kern
	for (int i = 0; i < N0-1; i++)
		for (int j = 0; j < N0-1; j++)
			p[ i*2*N1+2*N1 + 2*j+2 ]  = 0.25* ( x[ i*N0 + j ] + x[ i*N0 + j+1 ] + x[ i*N0 + N0 + j ] + x[ i*N0+N0 + j+1 ] );

}

void V_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU )
{

	for (int l = L-1; l > 0; l--)
	{
		// vorglaetten
		if ( glaetter == 0 )
			jacobi(sizeofL[l], Lx[l], Lb[l], v1, S);
		else if ( glaetter == 1)
			gaussSeidel(sizeofL[l], Lx[l], Lb[l], v1, S);
		else
			SOR(sizeofL[l], Lx[l], Lb[l], v1, S);

		// Bilde Residuum
		mfMult(sizeofL[l], Lx[l], Lr[l], S);
		for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
			Lr[l][i] = Lb[l][i] - Lr[l][i];

		// Restringiere
		rest( l, Lr[l], Lr[l-1], sizeofL);
		for (int i = 0; i < sizeofL[l-1] * sizeofL[l-1]; i++)
			Lb[l-1][i] = Lr[l-1][i];
	}

	// Loese genau
	substitutionen(sizeofL[0], S, LU, Lx[0], Lb[0]);

	for (int l = 1; l < L; l++)
	{
		// Prologieren
		prol(l-1, Lx[l-1], Lp[l], sizeofL);
		for (int i = 0; i < sizeofL[l]*sizeofL[l]; i++)
			Lx[l][i] += Lp[l][i];

		// nachglaetten
		if ( glaetter == 0 )
			jacobi(sizeofL[l], Lx[l], Lb[l], v2, S);
		else if ( glaetter == 1)
			gaussSeidel(sizeofL[l], Lx[l], Lb[l], v2, S);
		else
			SOR(sizeofL[l], Lx[l], Lb[l], v2, S);

	}

}

void W_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU )
{
	int l = L-1;
	// vorglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v1, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v1, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v1, S);

	// Bilde Residuum
	mfMult(sizeofL[l], Lx[l], Lr[l], S);
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lr[l][i] = Lb[l][i] - Lr[l][i];

	// Restringiere
	rest( l, Lr[l], Lr[l-1], sizeofL);
	l = l-1;
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lb[l][i] = Lr[l][i];

	// Loese genau
	if( l == 0 ){
		substitutionen(sizeofL[0], S, LU, Lx[0], Lb[0]);
	}else{
		W_Zyklus( S, glaetter, v1, v2, sizeofL, l+1, Lx, Lb, Lr, Lp, LU); // Rekursion
	}

	// Prologieren
	prol(l, Lx[l], Lp[l+1], sizeofL);
	l = l+1;
	for (int i = 0; i < sizeofL[l]*sizeofL[l]; i++)
		Lx[l][i] += Lp[l][i];

	// nachglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v2, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v2, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v2, S);

	// Bilde Residuum
	mfMult(sizeofL[l], Lx[l], Lr[l], S);
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lr[l][i] = Lb[l][i] - Lr[l][i];

	// Restringiere
	rest( l, Lr[l], Lr[l-1], sizeofL);
	l = l-1;
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lb[l][i] = Lr[l][i];

	if( l == 0 ){
		substitutionen(sizeofL[0], S, LU, Lx[0], Lb[0]);
	}else
	{
		W_Zyklus( S, glaetter, v1, v2, sizeofL, l+1, Lx, Lb, Lr, Lp, LU); // Rekursion
	}

	// Prologieren
	prol(l, Lx[l], Lp[l+1], sizeofL);
	l = l+1;
	for (int i = 0; i < sizeofL[l]*sizeofL[l]; i++)
		Lx[l][i] += Lp[l][i];

	// nachglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v2, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v2, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v2, S);

}

void F_Zyklus( const int S, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU )
{
	int l = L-1;
	// vorglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v1, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v1, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v1, S);

	// Bilde Residuum
	mfMult(sizeofL[l], Lx[l], Lr[l], S);
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lr[l][i] = Lb[l][i] - Lr[l][i];

	// Restringiere
	rest( l, Lr[l], Lr[l-1], sizeofL);
	l = l-1;
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lb[l][i] = Lr[l][i];

	// Loese genau
	if( l == 0 ){
		substitutionen(sizeofL[0], S, LU, Lx[0], Lb[0]);
	}else
	{
		F_Zyklus( S, glaetter, v1, v2, sizeofL, l+1, Lx, Lb, Lr, Lp, LU); // Rekursion
	}

	// Prologieren
	prol(l, Lx[l], Lp[l+1], sizeofL);
	l = l+1;
	for (int i = 0; i < sizeofL[l]*sizeofL[l]; i++)
		Lx[l][i] += Lp[l][i];

	// nachglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v2, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v2, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v2, S);

	// Bilde Residuum
	mfMult(sizeofL[l], Lx[l], Lr[l], S);
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lr[l][i] = Lb[l][i] - Lr[l][i];

	// Restringiere
	rest( l, Lr[l], Lr[l-1], sizeofL);
	l = l-1;
	for (int i = 0; i < sizeofL[l] * sizeofL[l]; i++)
		Lb[l][i] = Lr[l][i];

	if( l == 0 ){
		substitutionen(sizeofL[0], S, LU, Lx[0], Lb[0]);
	}else
	{
		V_Zyklus( S, glaetter, v1, v2, sizeofL, l+1, Lx, Lb, Lr, Lp, LU); // Rekursion
	}

	// Prologieren
	prol(l, Lx[l], Lp[l+1], sizeofL);
	l = l+1;
	for (int i = 0; i < sizeofL[l]*sizeofL[l]; i++)
		Lx[l][i] += Lp[l][i];

	// nachglaetten
	if ( glaetter == 0 )
		jacobi(sizeofL[l], Lx[l], Lb[l], v2, S);
	else if ( glaetter == 1)
		gaussSeidel(sizeofL[l], Lx[l], Lb[l], v2, S);
	else
		SOR(sizeofL[l], Lx[l], Lb[l], v2, S);

}

void mg( const double eps, const int S, const int zyklus, const int glaetter, const int v1, const int v2, const int* restrict sizeofL, const int L, double** Lx, double** Lb, double** Lr, double** Lp, const double* restrict LU )
{

	int N = sizeofL[L-1];

	for (int iter = 1; iter <= 10000; iter++)
	{

		if( zyklus == 0)
			F_Zyklus( S, glaetter, v1, v2, sizeofL, L, Lx, Lb, Lr, Lp, LU);
		else if( zyklus == 2 )
			W_Zyklus( S, glaetter, v1, v2, sizeofL, L, Lx, Lb, Lr, Lp, LU);
		else if ( zyklus == 1 )
			V_Zyklus( S, glaetter, v1, v2, sizeofL, L, Lx, Lb, Lr, Lp, LU);

		mfMult( N, Lx[L-1], Lr[L-1], S);
		for (int i = 0; i < N*N; i++)
			Lr[L-1][i] = Lb[L-1][i] - Lr[L-1][i];

		if( norm( N, Lr[L-1] ) < eps && iter < 10000)
		{
			printf("Benoetigte Iterationen:         %i \n", iter );
			printf("Norm des Residuums:             %.2e\n", norm( N, Lr[L-1] ) );
			break;
		}else if ( iter == 10000)
		{
			printf("***************************Achtung****************************\n");
			printf("Das Verfahren ist fruehzeitig nach %i Iterationen Abgebrochen.\n", iter );
			printf("Die dabei erreichte Norm des Residuums ist %.2e.\n", norm( N, Lr[L-1] ) );
			printf("**************************************************************\n");
		}
		//printf("Fehler %i %.2e\n", iter, norm( N, Lr[L-1] ) );

		for (int l = L-2; l >= 0; l--)
			for (int i = 0; i < sizeofL[l]*sizeofL[l] ; i++)
				Lx[ l ][ i ] = 0;

	}

}
