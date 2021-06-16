/**************************************************************//***********
 * \brief Contains all help functions
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
 ***********************************************************************/

typedef void (*myMult)(double*, double*, int);   ///< typedef function handle for choice of matrixfree mat-vec multiplication
typedef void (*mySmoother)(double*, double*, double* , double, double, double, int, myMult);   ///< typedef function handle for choice of smoothing operator
typedef void (*myInitLU)(double**, int, int);   ///< typedef function handle for choice of LU decomposition
typedef void (*myCycle)(double**, double**, double**, int, int, double, double, int, double, double, int*, int*, double**, mySmoother, myMult);   ///< typedef function handle for choice of cycle procedure
typedef void (*myInitXB)(double*, double*, double*, int , double);   ///< typedef function handle for initilizing right side and initial values


/*!
 * \brief a = b + ska * c
 *
 * @param[in] a    vector
 * @param[in] b    vector
 * @param[in] c    vector
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector b + ska * c
*/
void vec1(double * restrict a,
              double * restrict b,
              double * restrict c,
              double ska,
              int N
              )
{
	int i;
	#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = b[i] + ska * c[i];
  }
}

/*!
 * \brief a = a + ska * b
 *
 * @param[in] a    vector
 * @param[in] b    vector
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector  a + ska * b
*/
void vec2(double * restrict a,
              double * restrict b,
              double ska,
              int N
              )
{
	int i;
	#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] += ska * b[i];
  }
}

/*!
 * \brief a = ska * a
 *
 * @param[in] a    vector
 * @param[in] ska  skalar value
 * @param[in] N    size of the vectors
 *
 * @param[out] a   vector ska * a
*/
void vec3(double * restrict a,
          double ska,
          int N
          )
{
	int i;
	#pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = ska * a[i];
  }
}

/*!
 * \brief skalarproduct = a' * b
 *
 * @param[in] a   vector
 * @param[in] b   vector
 * @param[in] N   size of the vectors
 *
 * @param[out] a  vector a' * b
*/
double dot(double * restrict x,
        double * restrict y,
        int N)
{
  double vproduct = 0.0;
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
 * @param[in] x      vector
 * @param[in] N      size of the vectors
 *
 * @param[out] norm  value of the 2 norm of vector x
*/
double norm2(double * restrict x, int N)
{
  double norm = sqrt(dot(x,x,N));
  return norm;
}

/*!
 * \brief calculate the minimum of 2 values
 *
 * @param[in] x  skalar
 * @param[in] y  skalar
 *
 * @param[out] x vskalar
 * @param[out] y skalar
*/
double min(double x, double y)
{
  if (x < y)
	 return x;
	else
	 return y;
}

/*!
 * \brief f = 2*pi²*sin(x*y)*sin(x*y)
 *
 * @param[in] x   vector
 * @param[in] y   vector
 *
 * @param[out] x  vector with values of function f
*/
double f(double x, double y)
{
  x = 2 * pow(M_PI, 2) * sin(M_PI*x) * sin(M_PI*y);
  return x;
}

/*!
 * \brief u = sin(x*y)*sin(x*y)
 *
 * @param[in] x   vector
 * @param[in] y   vector
 *
 * @param[out] x  vector with values of function u
*/
double u(double x, double y)
{
  x = sin(M_PI*x) * sin(M_PI*y);
  return x;
}

/*!
 * \brief reports current time stamp
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
 * @param[in] res     vector
 * @param[in] x       vector
 * @param[in] b       vector
 * @param[in] N       size of vectors
 * @param[in] mfMult  function handle to the used matrixfree matvec-multiiplication
 *
 * @param[out] res    residual Ax-b
*/
void residual(double* res, double* x, double *b, int N, myMult mfMult)
{
  mfMult(x, res, N);
  vec1(res, b, res, -1, N*N);
}

/*!
 * \brief a = b duplication
 *
 * @param[in] a   vector
 * @param[in] b   vector
 * @param[in] N   size of vectors
 *
 * @param[out] a  vector b
*/
void duplicate(double *a, double *b, int N)
{
  #pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = b[i];
  }
}

/*!
 * \brief a = 0
 *
 * @param[in] a   vector
 * @param[in] N   size of vector
 *
 * @param[out] a  vector with all values = 0
*/
void zeroes(double *a, int N)
{
  #pragma omp parallel for
  for (int i=0;i<N;i++)
  {
    a[i] = 0.0;
  }
}

/*!
 * \brief calculates the Level L until grid "rough enough" to be solved directly
 *
 * @param[in] N   gridsize
 * @param[in] L   level = 0
 *
 * @param[out] L  optimized level
*/
int calculateL(int N, int L)
{
	int N_help=N;
  L = 1;
  while(N_help>50)
  {
    L += 1;
    N_help = (N_help-1)/2;
  }
  return L;
}

/*!
 * \brief calculates the grid size for each level L and catches wrong inputs for N
 *
 * @param[in]   N_vec array with values of N on each level
 * @param[in]   N_pot array with values of N squared on each level
 * @param[in]   L level
 * @param[in]   N initial grid size
 *
 * @param[out]  N_vec array with values of N on each level
 * @param[out]  N_pot array with values of N squared on each level
*/
int initN(int *N_vec,int *N_pot, int L, int N)
{
	N_vec[0] = N;
  N_pot[0] = N*N;
  for (int i=1; i<L; i++)
  {
    N_vec[i]=(N_vec[i-1]-1)/2;
    N_pot[i] = N_vec[i]*N_vec[i];
    if (((N_vec[i]%2)==0 && i!=L-1) || (N_vec[i]==1 && i!=L-1))
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
 * @param[in] x   vector
 * @param[in] y   vector
 * @param[in] k   factor in sin
 * @param[in] l   factor in sin
 * @param[in] h   step size
 *
 * @param[out] x  vector with values of initial values
*/
double initial_val(double x, double y, int k, int l, double h)
{
	x = sin(k*M_PI*x*h) + sin(l*M_PI*y*h);
}

/*!
 * \brief initialize starting values, right hand side and the analytical solution for 5-point stencil
 *
 * @param[in] x       vecotr
 * @param[in] b       vector
 * @param[in] N       initial grid size
 * @param[in] h       step size
 *
 * @param[out] x      initial values
 * @param[out] x_sol  analytical solution
 * @param[out] b      right hand side
*/
void initXB5(double *x, double *x_sol, double *b, int N, double h)
{
	// int k,l = 1;
  #pragma omp parallel
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
 * @param[in] x       vector
 * @param[in] b       vector
 * @param[in] N       initial grid size
 * @param[in] h       step size
 *
 * @param[out] x      initial values
 * @param[out] x_sol  vector with analytical solution
 * @param[out] b      right hand side
*/
void initXB9(double *x, double *x_sol, double *b, int N, double h)
{
	// int k,l = 1;
  #pragma omp parallel
  {
    for (int i = 0; i < N; i++  )
    {
      for (int j = 0; j < N ; j++ )
      {
        x[i*N + j] = 1.0; //initial_val(i+1, j+1, k, l, h);
        x_sol[i*N + j] = u((i+1)*h, (j+1)*h);
        b[i*N + j] = 6*h*h*f((i+1)*h, (j+1)*h);
      }
    }
  }
}

/*!
 * \brief catches input errors from main
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
    *mu = 1;
  }
  else
  {
    *mu = atoi(argv[6]);
  }

  if ((L < 0) || (N<1) )
  {
    printf("N has to be > 1 and L > 0\n");
    return 1;
  }

  if (L-*mu < 1 && L!=1)
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
