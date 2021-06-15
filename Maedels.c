void luZerlegung(double** LU, int gp, int s) {

        int n = gp*gp;                                          // n: Groesse der Matrix A

        #pragma omp parallel for
  for(int i=0;i<n;i++) {
                for(int j=0;j<n;j++) {
                        LU[i][j] = 0.0;
                }
        }

        // A fuellen

        if(s==1) {                                                              // 5-Punkte-Stern

                #pragma omp parallel
                {

                        #pragma omp for nowait
                        for(int i=0;i<n-gp;i++)
                        {
                                LU[i][i]    = 4.0;
                                LU[i][i+gp] = -1.0;
                                LU[i+gp][i] = -1.0;
                        }

                        #pragma omp for nowait
                        for(int i=n-gp;i<n;i++)
                        {
                                LU[i][i] = 4.0;
                        }

                        #pragma omp for collapse(2) nowait
                        for(int b=0;b<gp;b++)
                        {
                                for(int i=0;i<gp-1;i++)
                                {
                                        LU[b*gp+i][b*gp+i+1] = -1.0;
                                        LU[b*gp+i+1][b*gp+i] = -1.0;
                                }
                        }

                }
              }       else if(s==2) {                                 // 9-Punkte-Stern

                            #pragma omp parallel
                            {

                                    #pragma omp for nowait
                                    for(int i=0;i<n-gp;i++)
                                    {
                                            LU[i][i]    = 20.0;
                                            LU[i][i+gp] = -4.0;
                                            LU[i+gp][i] = -4.0;
                                    }

                                    #pragma omp for nowait
                                    for(int i=n-gp;i<n;i++)
                                    {
                                            LU[i][i] = 20.0;
                                    }

                                    #pragma omp for collapse(2) nowait
                                    for(int b=0;b<gp;b++)
                                    {
                                            for(int i=0;i<gp-1;i++)
                                            {
                                                    LU[b*gp+i][b*gp+i+1] = -4.0;
                                                    LU[b*gp+i+1][b*gp+i] = -4.0;
                                            }
                                    }

                                    #pragma omp for collapse(2) nowait
                                    for(int b=0;b<gp-1;b++)
                                    {
                                            for(int i=0;i<gp-1;i++)
                                            {
                                                    LU[b*gp+i+1][b*gp+i+gp] = -1.0;
                                                    LU[b*gp+i+gp][b*gp+i+1] = -1.0;
                                                    LU[b*gp+i][b*gp+i+gp+1] = -1.0;
                                                    LU[b*gp+i+gp+1][b*gp+i] = -1.0;
                                            }
                                    }

                            }

                    }
                    // LU-Zerlegung

                            for(int k=0;k<n-1;k++)
                            {
                                    for(int i=k+1;i<n;i++)
                                    {
                                            LU[i][k] = LU[i][k]/LU[k][k];

                                            for(int j=k+1;j<n;j++)
                                            {
                                                    LU[i][j] = LU[i][j]-LU[i][k]*LU[k][j];
                                            }
                                    }
                            }
