void createArrInArr(double* ArrOfArr, int N, int Ind, double val)
{

  double* q = malloc(N*sizeof(double));
  //TODO is this for useless?
  for (int i=0;i<N;i++)
  {
    q[i] = val;
  }

  ArrOfArr[ind] = q;
}

void v_cycle(N,b,x0,L,nu1,nu2)
{

  int i,j;
  double w = 0.6;

  double *NN = malloc(L*sizeof(double));
  double *ArrOfx = malloc(L*sizeof(double));
  double *ArrOfb = malloc(L*sizeof(double));
  double *ArrOfr = malloc(L*sizeof(double));


  NN[0] = N;

  ArrOfx[0] = x0;
  ArrOfb[0] = b;

  for (i=1; i<L; i++)
  {
    NN[i]=(NN[i-1]-1)/2;
//TODO Abfrage ob der Wert Gerade
  }


  for (i=1;i<L;i++)
  {
    createArrInArr(ArrOfx, NN[i], i, 0.0);
    createArrInArr(ArrOfb, NN[i], i, 1.0);
    createArrInArr(ArrOfr, NN[i], i, 0.0);
  }



  for i...
      double *b = ArrOfb[i];


        rechne mit b







free(b);

}
