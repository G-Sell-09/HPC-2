// Hier wird jetzt programmiert

#include <stdio.h>

// Funktion zur Berechnung des Skalarprodukts
int dot(int x[], int y[], int dim) {
  int vproduct = 0;
  for(int i=0;i<dim;i++){
    int product = vproduct + x[i] * y[i];
  }
  return vproduct;
}

// Hauptroutine
//int main(int argc, char **argv[]) {
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
