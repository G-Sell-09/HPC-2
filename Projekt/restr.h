
void restr(double* v_h, double* v_2h, int N)
{

  // New smaller grid size
  int NN = (N-1)/2;

  // Boundary values
  // double bound = 0;

  // Full-weighted operator
  for (int i = 0; i < NN; i++) {
    for (int j = 0; j < NN; j++) {

      // 4*selben + 2* unten,links,rechts,oben (Cache-optimiert) + 1* unten links,rechts & oben links, rechts
      v_2h[(i*NN)+j]= (0.0625) * (4 * v_h[(N + 2*i*N +1) + (2*j)] \
      + 2 * ( v_h[((2*i*N)+1) + (2*j)] + v_h[(2*i*N) + (2*j+N)] + v_h[(2*(i)*N) + (2*(j+1)+N)] + v_h[((2*(i+1)*N)+1)] ) \
      + 1 * ( v_h[(2*i*N) + (2*j)] + v_h[(2*i*N) + (2*(j+1))] + v_h[(2*(i+1)*N) + (2*j)] + v_h[(2*(i+1)*N) + (2*(j+1))] ) );
    }
  }
}
