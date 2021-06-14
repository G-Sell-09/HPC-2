void cycle(int N,double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int mu, double nu1, double nu2, double w, int step, int *N_vec,int *N_pot, \
  double **lower, double **upper, myFuncDef smooth)
{

  int i;
  double tol = 0.0;

  for (i=1;i<=mu;i++)
  {
    // prolong and correct
    prol(ArrOfx[L-i], ArrOfx[L-1-i], N_vec[L-i]);
  }

  for (i=L-1-mu;i<L-1;i++)
  {
    // Pre-Smoothing with Jacobi, Gauss or SOR
    smooth(ArrOfx[i], ArrOfb[i], ArrOfr[i], w, nu1, tol, N_vec[i]);

    //restriction
    restr(ArrOfr[i], ArrOfr[i+1], N_vec[i]);

    duplicate(ArrOfb[i+1], ArrOfr[i+1], N_pot[i+1]);

    zeroes(ArrOfx[i+1], N_pot[i+1]);
  }
}



void w_cycle(int N,double **ArrOfx, double **ArrOfb, double **ArrOfr, int L, int mu, double nu1, double nu2, double w, int step, int *N_vec,int *N_pot, \
  double **lower, double **upper, myFuncDef smooth)
{

  int i,j;
  double S,T,U,V;
  double tol = 0.0;

  for (i=0;i<L-1;i++)
  {
    // Pre-Smoothing with Jacobi, Gauss or SOR
    smooth(ArrOfx[i], ArrOfb[i], ArrOfr[i], w, nu1, tol, N_vec[i]);

    //restriction
    restr(ArrOfr[i], ArrOfr[i+1], N_vec[i]);

    duplicate(ArrOfb[i+1], ArrOfr[i+1], N_pot[i+1]);

    zeroes(ArrOfx[i+1], N_pot[i+1]);
  }

  // lu decomposition for the roughest grid
  lu_Decomp(ArrOfx[L-1],ArrOfb[L-1],N_pot[L-1], lower, upper);


  for (i=1; i<2*mu; i++)
  {
    if (i<=mu){
      cycle(N, ArrOfx, ArrOfb, ArrOfr, L, i, nu1, nu2, w, step, N_vec, N_pot, lower, upper, smooth);
    }
    else{
      cycle(N, ArrOfx, ArrOfb, ArrOfr, L, 2*mu-i, nu1, nu2, w, step, N_vec, N_pot, lower, upper, smooth);
    }
    lu_Decomp(ArrOfx[L-1],ArrOfb[L-1],N_pot[L-1], lower, upper);
  }

  for (i=2; i<L+1; i++)
  {
    // prolong and correct
    prol(ArrOfx[L-i+1], ArrOfx[L-i], N_vec[L-i+1]);

    //smooth(ArrOfx[L-i], ArrOfb[L-i], ArrOfr[L-i], w, nu2, tol, N_vec[L-i]);
  }

  printf("Times in W step %d: for1: %f, LU: %f, for2: %f\n", step, T-S, U-T, V-U);
}
