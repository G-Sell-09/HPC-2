/* 		Moegliche Abbruchbedingungen:
 * 			1) max_iter > 0, eps = -1 => Abbruch nach fester Anzahl von Iterationen
 * 			2) max_iter = -1, eps > 0 => Abbruch bei genuegend kleiner Residuumsnorm
 * 			3) max_iter > 0, eps > 0 => Abbruch nach fester Anzahl von Iterationen oder genuegend kleiner Residuumsnorm
 *
 * \brief SOR-Verfahren
 *
 * @param[in] x0 Startvektor
 * @param[in] b Rechte Seite
 * @param[in] omega Daempfungsparameter
 * @param[in] N Anzahl innerer Gitterpunkte in einer Dimension
 * @param[in] max_iter Maximale Iterationszahl
 * @param[in] eps Toleranzgrenze
 * @param[out] k Anzahl benoetigter Iterationen
 *
 * \author Torsten Foehr
 * \author Daniel Schuster
 *
 * \version 1.0
 * \copyright Team TODA
 *
*/
int sor (double *x0, double *b, double omega, int N, int max_iter, double eps)
{
	// 1) Vorbereitungen

	// Bestimmung der Abbruchbedingung des Verfahrens

	int break_condition;
	if (eps == -1)
	{
		// Nur max_iter entscheidet ueber Abbruch
		// => Abbruch nach fester Anzahl von Iterationen
		break_condition = 1;
	}
	else if (max_iter == -1)
	{
		// Nur eps entscheidet ueber Abbruch
		// => Abbruch bei genuegend kleiner Residuumsnorm
		break_condition = 2;
	}
	else
	{
		// max_iter und eps entscheiden uber Abbruch
		// => Abbruch nach fester Anzahl von Iterationen oder genuegend kleiner Residuumsnorm
		break_condition = 3;
	}

	// Residuum im Fall der Abbruchbedingungen 2 und 3
	double* res = (double *)malloc(N*N*sizeof(double));

	// Vorbereitung der ersten Iteration
	int k=0;

	double *xkm1 = (double *)malloc(N*N*sizeof(double));
	vec_copy(x0,xkm1,N*N);
	double *xk = x0;


	// 2) SOR - Schleife

	while(1)
	{
		// Angleich Iterationszahl
		k=k+1;


		// SOR - Updates

		/* Im Gegensatz zum Jacobi-Verfahren haengen die zu berechnenden Eintraege der neuen Iterierten nicht nur
		 * von Eintraegen der alten Iterierten ab, sondern werden zusaetzlich von bereits berechneten Eintraegen
		 * der neuen Iterierten beeinflusst. Aus diesem Grund koennen wir die Updates nicht fallweise wie im
		 * Jacobi-Verfahren berechnen. Weil die Reihenfolge der Vektoreintraege auf einer zeilenweisen Nummerierung
		 * der Gitterpunkte basiert, koennen wir das Gitter zum Updaten aber systematisch zeilenweise von unten
		 * nach oben durchlaufen. */


		// Schritt 1: Unterste innere Gitterzeile (j = 1)

		// Linke untere Ecke (j = 1, i = 1)
		xk[0] = omega*(b[0] + xkm1[1] + xkm1[N])/4 + (1-omega)*xkm1[0];

		// Punkte im Inneren der Zeile
		for (int i = 2; i <= N-1; i++) // Gitterspalte
		{
			xk[i-1] = omega*(b[i-1] + xk[i-2] + xkm1[i] + xkm1[N+(i-1)])/4 + (1-omega)*xkm1[i-1];
		}

		// Rechte untere Ecke (j = 1, i = N)
		xk[N-1] = omega*(b[N-1] + xk[N-2] + xkm1[N+(N-1)])/4 + (1-omega)*xkm1[N-1];


		// Schritt 2: Mittlere innere Gitterzeilen

		for (int j = 2; j <= N-1; j++) // Gitterzeile
		{
			// Punkt am linken Zeilenrand (i = 1)
			xk[(j-1)N] = omega*(b[(j-1)*N] + xk[(j-2)*N] + xkm1[(j-1)*N+1] + xkm1[j*N])/4 + (1-omega)*xkm1[(j-1)*N];

			// Punkte im Inneren der Zeile
			for (int i = 2; i <= N-1; i++) // Gitterspalte
			{
				xk[(j-1)N+(i-1)] = omega*(b[(j-1)*N+(i-1)] + xk[(j-2)*N+(i-1)] + xk[(j-1)*N+(i-2)] + xkm1[(j-1)*N+i] + xkm1[j*N+(i-1)])/4 + (1-omega)*xkm1[(j-1)*N+(i-1)];
			}

			// Punkt am rechten Zeilenrand (i = N)
			xk[(j-1)N+(N-1)] = omega*(b[(j-1)*N+(N-1)] + xk[(j-2)*N+(N-1)] + xk[(j-1)*N+(N-2)] + xkm1[j*N+(N-1)])/4 + (1-omega)*xkm1[(j-1)*N+(N-1)];
		}


		// Schritt 3: Oberste innere Gitterzeile (j = N)

		// Linke obere Ecke (j = N, i = 1)
		xk[(N-1)*N] = omega*(b[(N-1)*N] + xk[(N-2)*N] + xkm1[(N-1)*N+1])/4 + (1-omega)*xkm1[(N-1)*N];

		// Punkte im Inneren der Zeile
		for (int i = 2; i <= N-1; i++) // Gitterspalte
		{
			xk[(N-1)N+(i-1)] = omega*(b[(N-1)*N+(i-1)] + xk[(N-2)*N+(i-1)] + xk[(N-1)*N+(i-2)] + xkm1[(N-1)*N+i])/4 + (1-omega)*xkm1[(N-1)*N+(i-1)];
		}

		// Rechte obere Ecke (j = N, i = N)
		xk[(N-1)N+(N-1)] = omega*(b[(N-1)*N+(N-1)] + xk[(N-2)*N+(N-1)] + xk[(N-1)*N+(N-2)])/4 + (1-omega)*xkm1[(N-1)*N+(N-1)];


		// Abbruchbedingung
		switch(break_condition)
		{
			case 1:
			// Abbruch nach fester Anzahl von Iterationen
				if (k >= max_iter)
				{
					free(res);
					free (xkm1);
					return k;
				}
				break;
			case 2:
			// Abbruch bei genuegend kleiner Residuumsnorm
				get_residual(res,xk,b,N);
				if (norm(res,N*N) < eps)
				{
					free(res);
					free(xkm1);
					return k;
				}
				break;
			case 3:
			// Abbruch nach fester Anzahl von Iterationen oder genuegend kleiner Residuumsnorm
				get_residual(res,xk,b,N);
				if (k >= max_iter || norm(res, N*N) < eps)
				{
					free(res);
					free(xkm1);
					return k;
				}
				break;
		}


		// Vorbereitung fuer die naechste Iteration
		vec_copy(xk,xkm1,N*N);
	}
}
