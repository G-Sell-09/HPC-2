
void prol(int N, double* v_2h, double* v_h)
{

  // New bigger grid size
  int NN = 2*N+1

  // Boundary values
  // double bound = 0;

// Zunächst alle neuen Punkte die keinen Kontakt zum Rand haben
  // Punkte aus altem Gitter alle übernehmen
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // v_h[2i,2j] = v_2h[i,j]
      v_h[(2*i*NN) + (2*j)] = v_2h[i*N+j];
    }
  }

  // Im Skript falsch farbig markiert. Hier sind die Punkte die immer über und unter den alten Punkten liegen
  for (int i = 1; i < N; i++) {
    for (int j = 0; j < N; j++) {
      //v_h[2i+1,2j] = 1/2 * (v_2h[i,j]+v_2h[i+1,j])
      v_h[(2*(i+1)*NN) + (2*j)] = 0.5 * (v_2h[i*N+j] + v_2h[(i+1)*N+j]);
    }
  }

  // Im Skript falsch farbig markiert. Hier sind die Punkte die immer links und rechts neben den alten Punkten liegen
  for (int i = 0; i < N; i++) {
    for (int j = 1; j < N; j++) {
      // v_h[2i,2j+1] = 1/2 * (v_2h[i,j]+v_2h[i,j+1])
      v_h[(2*i*NN) + (2*j+NN)] = 0.5 * (v_2h[i*N+j] + v_2h[i*N+(j+1)]);
    }
  }

  // Punkte in der Mitte der alten Punkte
  for (int i = 1; i < N; i++) {
    for (int j = 1; j < N; j++) {
      // v_h[2i+1,2j+1] = 1/4 * (v_2h[i,j] + v_2h[i,j+1] + v_2h[i+1,j] + v_2h[i+1,j+1])
      v_h[(2*(i+1)*NN) + (2*j+NN)] = 0.25 * (v_2h[i*N+j] + v_2h[(i+1)*N+j] + v_2h[i*N+(j+1)] + v_2h[(i+1)*N+(j+1)]);
    }
  }

// Jetzt alle Punkte die Kontakt zum Rand haben
  // Alle Punkte die direkt neben altem Punkt und Randpunkt liegen
  for (int i = 0; i < N; i++) {

    // Punkte unten am Rand
    v_h[(2*i)+1]              = 0.5 * v_2h[i];
    // Punkte links am Rand
    v_h[NN + 2*i*NN]          = 0.5 * v_2h[i*N];
    // Punkte rechts am Rand
    v_h[(2*NN)-1 + 2*i*NN]    = 0.5 * v_2h[(N-1)+i*N];
    // Punkte oben am Rand
    v_h[(NN*(NN-1)+1) + 2*i]  = 0.5 * v_2h[(N*(N-1)) + i];

  }

  // Alle Punkte am unteren Rand ohne die Ecken
  for (int i = 1; i < N; i++) {

    // Punkte unten am Rand
    v_h[(2*i)]            = 0.25 * (v_2h[i-1] + v_2h[i]);
    // Punkte links am Rand
    v_h[2*i*NN]           = 0.25 * (v_2h[(i-1)*N] + v_2h[i*N];
    // Punkte rechts am Rand
    v_h[(2*NN*i)-1 + NN]  = 0.25 * (v_2h[(N-1) + (i-1)*N] + v_2h[(N-1) + (i)*N]);
    // Punkte oben am Rand
    v_h[NN*(NN-1) + 2*i]  = 0.25 * (v_2h[(N*(N-1)) + (i-1)] + v_2h[(N*(N-1)) + (i)];

  }

  // Alle Eckpunkte betrachten
  // unten links
  v_h[0]          = 0.25 * v_2h[0];
  // unten rechts
  v_h[NN-1]       = 0.25 * v_2h[N-1];
  //oben links
  v_h[NN*(NN-1)]  = 0.25 * v_2h[N*(N-1)];
  //oben rechts
  v_h[NN*NN-1]    = 0.25 * v_2h[N*N-1];


  return v_h;
}
