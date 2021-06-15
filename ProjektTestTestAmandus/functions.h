/*************************************************************************
 * \brief Contains all help functions
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 ***********************************************************************/

typedef void (*myMult)(double*, double*, int); ///< typedef function handle for choice of matrixfree mat-vec multiplication
typedef void (*mySmoother)(double*, double*, double* , double, double, double, int, myMult); ///< typedef function handle for choice of smoothing operator
typedef void (*myInitLU)(double**, double**, double**, int, int); ///< typedef function handle for choice of LU decomposition
typedef void (*myCycle)(double**, double**, double**, int, int, double, double, int, double, double, int*, int*, double**, double**, mySmoother, myMult); ///< typedef function handle for choice of cycle procedure
typedef void (*myInitXB)(double*, double*, double*, int , double); ///< typedef function handle for initilizing right side and initial values

/*!
 * \brief calculate the minimum of 2 values
 *
 * @param[in] x  value of x
 * @param[in] y  value of y
 *
 * @param[out] x value of x
 * @param[out] y value of y
*/
double min(double x, double y)
{
  if (x < y)
	 return x;
	else
	 return y;
}

/*!
 * \brief a = b + ska * c
 *
 * @param[in] a    vector a
 * @param[in] b    vector b
 * @param[in] c    vector c
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector a
*/
void vec1(double * restrict a,
              double * restrict b,
              double * restrict c,
              double ska,
              int N
              )
{
	int i;
	//#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = b[i] + ska * c[i];
  }
}

/*!
 * \brief a = a + ska * b
 *
 * @param[in] a    vector a
 * @param[in] b    vector b
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector a
*/
void vec2(double * restrict a,
              double * restrict b,
              double ska,
              int N
              )
{
	int i;
	//#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] += ska * b[i];
  }
}

/*!
 * \brief a = ska * a
 *
 * @param[in] a    vector a
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector a
*/
void vec3(double * restrict a,
          double ska,
          int N
          )
{
	int i;
	//#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = ska * a[i];
  }
}

/*!
 * \brief skalarproduct = a * b
 *
 * @param[in] a   vector a
 * @param[in] b   vector b
 * @param[in] N   size of the vectors
 *
 * @param[out] a  vector a
*/
double dot(double * restrict x,
        double * restrict y,
        int N)
{
  double vproduct = 0.0; ///< Skalarprodukt
	int i;
  #pragma omp parallel for schedule(static) reduction(+:vproduct)
  for(int i=0;i<N;i++)
  {
    vproduct += x[i] * y[i];
  }
  return vproduct;
}

/*!
 * \brief 2-norm of a vector x
 *
 * @param[in] x      vector x
 * @param[in] N      size of the vectors
 *
 * @param[out] norm  norm value
*/
double norm2(double * restrict x, int N)
{
  double norm = sqrt(dot(x,x,N));
  return norm;
}

/*!
 * \brief f = 2*pi²*sin(x*y)*sin(x*y)
 *
 * @param[in] x   vector x
 * @param[in] y   vector y
 *
 * @param[out] x  vector x
*/
double f(double x, double y)
{
  x = 2 * pow(M_PI, 2) * sin(M_PI*x) * sin(M_PI*y);
  return x;
}

/*!
 * \brief u = sin(x*y)*sin(x*y)
 *
 * @param[in] x   vector x
 * @param[in] y   vector y
 *
 * @param[out] x  vector x
*/
double u(double x, double y)
{
  x = sin(M_PI*x) * sin(M_PI*y);
  return x;
}

/*!
 * \brief gives current time stamp
 *
 * @param[out] ts   time stemp
*/
double getTimeStamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9;
}

/*!
 * \brief calculates residual of Ax-b for given matrixfree mat-vec multiplication
 *
 * @param[in] res     residual
 * @param[in] x       vector x
 * @param[in] b       vector b
 * @param[in] N       vector size
 * @param[in] mfMult  function handle to the 5-point or 9-point matrixfree mat-vec multiiplication
 *
 * @param[out] res    residual
*/
void residual(double* res, double* x, double *b, int N, myMult mfMult)
{
  mfMult(x, res, N);
  vec1(res, b, res, -1, N*N);
}

/*!
 * \brief a = b duplication
 *
 * @param[in] a   vector a
 * @param[in] b   vector b
 * @param[in] N   vector size
 *
 * @param[out] a  vector a
*/
void duplicate(double *a, double *b, int N)
{
  //#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = b[i];
  }
}

/*!
 * \brief a = 0
 *
 * @param[in] a   vector a
 * @param[in] N   vector size
 *
 * @param[out] a  vector a
*/
void zeroes(double *a, int N)
{
  //#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i]=0.0;
  }
}

/*!
 * \brief calculates the Level L until grid "small enough"
 *
 * @param[in] N   vector size
 * @param[in] L   level
 *
 * @param[out] L  level
*/
int calculateL(int N, int L)
{
	int N_help=N;
  L = 1;
	for (int i = 0; i < 100; i++)
	{
		if (N_help > 55)
		{
			L += 1;
			N_help = (N_help-1)/2;
		}
    else
    break;
	}
  return L;
}

/*!
 * \brief calculates the grid size for each level L
 *
 * @param[in]   N_vec array with values of N on each level
 * @param[in]   N_pot array with values of N squared on each level
 * @param[in]   L level
 * @param[in]   N initial grid size
 *
 * @param[out]  N_vec array with values of N on each level
*/
int initN(int *N_vec,int *N_pot, int L, int N)
{
	N_vec[0] = N;
  N_pot[0] = N*N;
  for (int i=1; i<L; i++)
  {
    N_vec[i]=(N_vec[i-1]-1)/2;
    N_pot[i] = N_vec[i]*N_vec[i];
    if ((N_vec[i]%2)==0 && i!=L-1)
    {
      printf("The grid size can not be an even number on any level except the last one. Please check the entered values for N and L! \n");
    	return 1;
    }
  }
	return 0;
}

/*!
 * \brief x = sin(k*pi*x*h)*sin(l*pi*y*h)
 *
 * @param[in] x   vector x
 * @param[in] y   vector y
 * @param[in] k   factor k
 * @param[in] l   factor l
 * @param[in] h   step size
 *
 * @param[out] x  vector x
*/
double initial_val(double x, double y, int k, int l, double h)
{
	x = sin(k*M_PI*x*h) + sin(l*M_PI*y*h);
}

/*!
 * \brief initialize starting values, right hand side and the analytical solution for 5-point stencil
 *
 * @param[in] x       vecotr x
 * @param[in] b       vector b
 * @param[in] N       initial grid size
 * @param[in] h       step size
 *
 * @param[out] x      starting vector x
 * @param[out] x_sol  vector with analytical solution
 * @param[out] b      right hand side
*/
void initXB5(double *x, double *x_sol, double *b, int N, double h)
{
	// int k,l = 1;
  //#pragma omp parallel
  {
    for (int i = 0; i < N; i++  )
    {
      for (int j = 0; j < N ; j++ )
      {
        x[i*N + j] = 1.0; //initial_val(i+1, j+1, k, l, h);
        x_sol[i*N + j] = u((i+1)*h, (j+1)*h);
        b[i*N + j] = h*h*f((i+1)*h, (j+1)*h);
      }
    }
  }
}

/*!
 * \brief initialize starting values, right hand side and the analytical solution for 9-point stencil
 *
 * @param[in] x       vecotr x
 * @param[in] b       vector b
 * @param[in] N       initial grid size
 * @param[in] h       step size
 *
 * @param[out] x      starting vector x
 * @param[out] x_sol  vector with analytical solution
 * @param[out] b      right hand side
*/
void initXB9(double *x, double *x_sol, double *b, int N, double h)
{
	// int k,l = 1;
  //#pragma omp parallel
  {
    for (int i = 0; i < N; i++  )
    {
      for (int j = 0; j < N ; j++ )
      {
        x[i*N + j] = 0.0; //initial_val(i+1, j+1, k, l, h);
        x_sol[i*N + j] = u((i+1)*h, (j+1)*h);
        b[i*N + j] = 6*h*h*f((i+1)*h, (j+1)*h);
      }
    }
  }
}

/*!
 * \brief checks inputs from main for incorrect parameters
 *
 * @param[in] L         level of the multigrid method
 * @param[in] N         gridsize
 * @param[in] GL        name of chosen smoothing operator
 * @param[in] stencil   chosen stencil
 * @param[in] cyc       chosen cycle procedure
 * @param[in] argc      number of input arguments from main()
 * @param[in] argv      input arguments from main()
 * @param[in] mu        number of iterations for mu_cycle()
 *
 * @param[out] {0,1}    0 if the inputs are valid, 1 otherwise
*/
int inputs(int L, int N, char *GL, int stencil, char *cyc, int argc, char **argv, int *mu)
{
  if( argc < 7 )
  {
    *mu = 0;
  }
  else
  {
    *mu = atoi(argv[6]);
  }

  if ((L < 0) || (N<1) )
  {
    printf("N and L have to be positive integers\n");
    return 1;
  }

  if (L-*mu < 1)
  {
    printf("mu: %d is not valid for level: %d. Please choose different values\n", *mu, L);
    return 1;
  }

  if (strcmp(cyc,"v") != 0 && strcmp(cyc,"w") != 0 && strcmp(cyc, "mu") != 0)
  {
    printf("No valid pattern given. Please choose v, w or mu\n");
    return 1;
  }
  else if (strcmp(GL, "jac") != 0 && strcmp(GL,"gau") != 0 && strcmp(GL,"sor") != 0)
  {
    printf("No valid smoothing operator. Please choose jac, gau or sor\n");
    return 1;
  }
  else if (stencil != 5 && stencil != 9)
  {
    printf("No valid stencil for input: %d. Please choose 5 or 9.\n", stencil);
    return 1;
  }
  else if ((strcmp(cyc, "mu") == 0) && (*mu < 1))
  {
    printf("No cycle procedure for input: %d. Please choose 1, 2 or 3.\n", stencil);
    return 1;
  }
  else
  {
    return 0;
  }
}
