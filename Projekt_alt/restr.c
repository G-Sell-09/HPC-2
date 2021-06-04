
void restr(int N, double* v_h)
{

  // New smaller grid size
  int NN = (N-1)/2;

  // New grid values in bigger vector
  double v_2h = malloc(NN*NN*sizeof(double));

  // Full-weighted operator
  for (int i = 0; i <= NN; i++) {
    for (int j = 0; j <= NN; j++) {

      // v_h[2i,2j] = v_2h[i,j]
      v_2h[(2*i*NN)    + (2*j)]  =         v_2h[i*N+j];

      // Im Skript falsch farbig markiert. Hier sind die Punkte die immer über und unter den alten Punkten liegen
      //v_h[2i+1,2j] = 1/2 * (v_2h[i,j]+v_2h[i+1,j])
      v_2h[(2*(i+1)*NN)+ (2*j)]  = 0.5 *   (v_2h[i*N+j]  + v_2h[(i+1)*N+j]);

      // Im Skript falsch farbig markiert. Hier sind die Punkte die immer links und rechts neben den alten Punkten liegen
      // v_h[2i,2j+1] = 1/2 * (v_2h[i,j]+v_2h[i,j+1])
      v_2h[(2*i*NN)    + (2*j+NN)]= 0.5 *  (v_2h[i*N+j]  + v_2h[i*N+(j+1)]);

      // v_h[2i+1,2j+1] = 1/4 * (v_2h[i,j] + v_2h[i,j+1] + v_2h[i+1,j] + v_2h[i+1,j+1])


      v_2h[(i*NN)+j]= (1/16) * (4 * v_h[(2*i*NN) + (2*j)] \
      + 2 * ( v_h[(i+1)*N+j] + v_h[i*N+(j+1)] + v_h[(i+1)*N+(j+1)]);
    }
  }
  return v_h;
}
