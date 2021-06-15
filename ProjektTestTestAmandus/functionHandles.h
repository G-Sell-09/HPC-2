/*!
 * \brief Picks the correct smoothing operator for chosen stencil and initializes the relaxation factor.
 *
 * @param[in] GL       {"jac", "gau", "sor"} selects the smoothing operator
 * @param[in] w        initializes the appropriate relaxation factor
 * @param[in] stencil  {5, 9} accounts for the choice of stencil
 *
 * @param[out] smooth  {&jac, &gau, &sor} typdef mySmoother as described in functions.h
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
mySmoother smoothHandle(char* GL, double *w, int stencil) {

  mySmoother smooth;
  if (strcmp(GL, "jac") == 0 && stencil == 5)
	{
		smooth = &jac5;
    *w = 0.6;
	}
  else if (strcmp(GL, "jac") == 0 && stencil == 9)
	{
		smooth = &jac9;
    *w = 0.6;
	}
	else if (strcmp(GL, "gau") == 0  && stencil == 5)
	{
		smooth = &sor5;
    *w = 1.0;
	}
	else if (strcmp(GL, "sor") == 0  && stencil == 5)
	{
		smooth = &sor5;
    *w = 1.95;
	}
  else if (strcmp(GL, "gau") == 0  && stencil == 9)
	{
		smooth = &sor9;
    *w = 1.0;
	}
  else if (strcmp(GL, "sor") == 0  && stencil == 9)
	{
		smooth = &sor9;
    *w = 1.95;
	}
	return smooth;
}

/*!
 * \brief Picks a function handle to initialize the correct LU decomposition.
 *
 * @param[in] stencil    {5, 9} accounts for the choice of stencil
 *
 * @param[out] initLU    {&initLU5, &initLU9} typdef myinitLU as described in functions.h
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
myInitLU initLUHandle(int stencil)
{
  myInitLU initLU;
  if (stencil==5)
  {
    initLU = &initLU5;
  }
  else if (stencil==9)
  {
    initLU = &initLU9;
  }
  return initLU;
}

/*!
 * \brief Picks the correct function handle for, matrix-free Mat-Vector multiplication to a given stencil.
 *
 * @param[in] stencil   {5, 9} accounts for the choice of stencil
 *
 * @param[out] mfMult   {&mfMult5, &mfMult9} typdef myMult as described in functions.h
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
myMult multHandle(int stencil)
{
  myMult mfMult;
  if (stencil==5)
  {
    mfMult = &mfMult5;
  }
  else if (stencil==9)
  {
    mfMult = &mfMult9;
  }
  return mfMult;
}

/*!
 * \brief Provides a function handle for either v_cycle or w_cylce used in the multigrid method.
 *
 * @param[in] stencil   {"v", "w"} v-cycle or w-cycle available
 *
 * @param[out] cycle    {&v_cycle_rec, &w_cycle_rec} typedef myCycle as described in functions.h
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
myCycle cycleHandle(char* cyc, int *mu)
{
  myCycle cycle;
  if (strcmp(cyc,"v") == 0)
  {
    cycle = &v_cycle;
    *mu = 0;
  }
  else if (strcmp(cyc,"w") == 0)
  {
    cycle = &mu_cycle;
    *mu = 2;
  }
  else if (strcmp(cyc, "mu") == 0)
  {
    cycle = &mu_cycle;
  }
  return cycle;
}

/*!
 * \brief Given a stencil picks the correct methode to initialize the right side and initial values.
 *
 * @param[in] stencil   {5, 9} accounts for the choice of stencil
 *
 * @param[out] initXB   {&initXB5, &initXB9} typdef myInitXB as described in functions.h
 *
 * \author Robin Sell
 * \author Neil Vetter
 *
 * \version 1.0
 * \copyright HPC-2 Team RSNV
*/
myInitXB initXBHandle(int stencil)
{
  myInitXB initXB;
  if (stencil == 5)
  {
    initXB = &initXB5;
  }
  else if (stencil == 9)
  {
    initXB = &initXB9;
  }
  return initXB;
}
