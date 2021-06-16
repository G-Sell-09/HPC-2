/*
 * \brief Restriction
 *
 * @param[in] v_h     vector with values of the finer grid
 * @param[in] v_2h    vector for the values of the rougher grid
 * @param[in] N       size of the finer grid
 *
 * @param[out] v_2h   values on the new rougher grid
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 *
*/
void restr(double* v_h, double* v_2h, int N)
{

  // new gridsize of the rougher grid
  int NN = (N-1)/2;
  int i,j;

  #pragma omp parallel
  {
    // full-weighted operator
    for (int i = 0; i < NN; i++) {
      for (int j = 0; j < NN; j++) {

        // cache-optimized order of old nodes to weight the new nodes
        v_2h[(i*NN)+j]= (0.0625) * (4 * v_h[(N + 2*i*N +1) + (2*j)] \
        + 2 * ( v_h[((2*i*N)+1) + (2*j)] + v_h[(2*i*N) + (2*j+N)] + v_h[(2*(i)*N) + (2*(j+1)+N)] + v_h[((2*(i+1)*N)+1)] ) \
        + 1 * ( v_h[(2*i*N) + (2*j)] + v_h[(2*i*N) + (2*(j+1))] + v_h[(2*(i+1)*N) + (2*j)] + v_h[(2*(i+1)*N) + (2*(j+1))] ) );
      }
    }
  }
}
