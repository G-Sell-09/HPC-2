/*!
Skript zum berechnen des Skalarprodukts zweier Vektoren.

\brief Skalarprodukt

\author Robin Sell, 6071120
\author Neil Vetter, 6021336

\version 1.0

*/
#include <stdio.h>
#include <stdlib.h>

// Funktion zur Berechnung des Skalarprodukts
/// Funktion die das Skalrprodukt berechnet
int dot(int x[], int y[], int dim) {
  int vproduct = 0;
  for(int i=0;i<dim;i++){
    vproduct += x[i] * y[i];
  }
  return vproduct;
}

// Hauptroutine

/// Main
int main(void) {

  // int dim = atoi(argv[1]);
  int dim = 2;

  int *x = malloc(dim * sizeof(int));
  x[0] = 2;
  x[1] = 3;

  int *y = malloc(dim * sizeof(int));
  y[0] = 5;
  y[1] = 2;

  int lsg = dot(x,y,dim);

  printf("Das Skalarprodukt beträgt: %d\n", lsg);

  return 0;
}
