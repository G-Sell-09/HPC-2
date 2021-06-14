
void prol(double* v_2h, double* v_h, int N)
{

  // New finer grid size for boundary values = 0
  int NN = 2*N+1;

  // Buffer for v_h
  double *buffer = malloc(NN*NN*sizeof(double));
  duplicate(buffer, v_h, NN*NN);

// First look at all nodes that have no contact with the boundary values

  // First all old node values into the knew knode
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // v_h[2i,2j] = v_2h[i,j]
      v_h[(NN + 2*i*NN +1) + (2*j)] = v_2h[i*N+j];
    }
  }

  // new values for knodes that are above and below the old nodes
  for (int i = 1; i < N; i++) {
    for (int j = 0; j < N; j++) {
      v_h[((2*i*NN)+1) + (2*j)] = 0.5 * (v_2h[(i-1)*N+j] + v_2h[i*N+j]);
    }
  }

  // new values for knodes that are left and right of the old nodes
  for (int i = 0; i < N; i++) {
    for (int j = 1; j < N; j++) {
      v_h[(2*i*NN) + (2*j+NN)] = 0.5 * (v_2h[i*N+(j-1)] + v_2h[i*N+(j)]);
    }
  }

  // new values for knodes in the middle of the old nodes
  for (int i = 1; i < N; i++) {
    for (int j = 1; j < N; j++) {
      v_h[(2*i*NN) + (2*j)] = 0.25 * (v_2h[(i-1)*N+(j-1)] + v_2h[(i-1)*N+(j)] + v_2h[(i)*N+(j-1)] + v_2h[(i)*N+(j)]);
    }
  }

// Now look at all nodes that have contact with the boundary values

  // new values for knodes that are directly next to old node and the boundary node
  for (int i = 0; i < N; i++) {

    // nodes at the lower edge
    v_h[(2*i)+1]              = 0.5 * v_2h[i];
    // nodes at the left edge
    v_h[NN + 2*i*NN]          = 0.5 * v_2h[i*N];
    // nodes at the right edge
    v_h[(2*NN)-1 + 2*i*NN]    = 0.5 * v_2h[(N-1)+i*N];
    // nodes at the upper edge
    v_h[(NN*(NN-1)+1) + 2*i]  = 0.5 * v_2h[(N*(N-1)) + i];

  }

  // new values for knodes that are in between two old nodes and two boundary nodes
  for (int i = 1; i < N; i++) {

    // nodes at the lower edge
    v_h[(2*i)]            = 0.25 * (v_2h[i-1] + v_2h[i]);
    // nodes at the left edge
    v_h[2*i*NN]           = 0.25 * (v_2h[(i-1)*N] + v_2h[i*N]);
    // nodes at the right edge
    v_h[(2*NN*i)-1 + NN]  = 0.25 * (v_2h[(N-1) + (i-1)*N] + v_2h[(N-1) + (i)*N]);
    // nodes at the upper edge
    v_h[NN*(NN-1) + 2*i]  = 0.25 * (v_2h[(N*(N-1)) + (i-1)] + v_2h[(N*(N-1)) + (i)]);

  }

  // new values for knodes that are in the corner
  // lower left
  v_h[0]          = 0.25 * v_2h[0];
  // lower right
  v_h[NN-1]       = 0.25 * v_2h[N-1];
  // upper left
  v_h[NN*(NN-1)]  = 0.25 * v_2h[N*(N-1)];
  // upper right
  v_h[NN*NN-1]    = 0.25 * v_2h[N*N-1];

  // v_h = buffer + v_h
  vec2(v_h, buffer, 1,NN*NN);

}
