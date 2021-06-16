----Multigrid method----

author: Robin Sell
author: Neil Vetter

version:   1.0
copyright: HPC-2 Team RSNV

compile with: "gcc -O3 -o grid  ./main.c -lm"          for serial mode
compile with: "gcc -03 -fopenmp -o gridP ./main.c -lm"  for openMP parallel mode
Note: Parallizing boosts about 15-25% depending on the used smoothing operator and the grid size, but it does not work well on Amandus.

call with    : ./grid [N] [L] [GL] [stencil] [cyc] [mu]
call with    : ./gridP [N] [L] [GL] [stencil] [cyc] [mu]

Alternativly : bash run5.sh              // Tests for 5 point stencil
		           bash run9.sh             // Tests for 9 point stencil
		           bash run5_parallel.sh    // Tests for 5 point stencil with openMP
		           bash run9_parallel.sh    // Tests for 9 point stencil with openMP

Doxygen with : Doxygen config


Input parameter:

 + [N]       : Gridsize of the initial grid.
 	        			- Can be any positive integer that produces uneven gridsizes on every level but the last.

 + [L]       : Represents the level of the inital grid. From here the method resticts the grid in L-1 steps down to level 1.
 								- Can be any positive integer, for which the resulting thinnest grid is equal or greater than 1x1.
 								- Or can be 0, in which case an optimal L is calculated (pre choosen value is a roughest grid with N<50).

 + [GL]      : Is the smoothing operator used in the pre-/post-smoothing process.
 								- Can be 	 "jac"  -> Jacobi       with relaxation factor : w = 0.6.
				 				- "gau"  -> Gauß-Seidel  with relaxation factor : w = 1.0.
 								- "sor"  -> SOR		 with relaxation factor : w = 1.95.

 + [stencil] : Is the stencil used in the matrixfree matvec-multiplikation and the sor smoothing operator.
 								- Gives a choice between 5 and 9. Represents either the classic 5 point stencil or the compact 9 point stencil.

 + [cyc]     : Cycle procedure used in the multigrid method.
 								- Can be either "v"     ->  v cycle
 				 												"w"     ->  w cycle
 				 												"mu"    ->  mu_cycle

 + [mu]      : Optional parameter in mu_cycle that controls the number of rekursions on each level.
 								- Can be any positive integer up to L-1


Static test parameters:

	nu1,nu2 = 30        	# of iterations per pre-/post-smoothing
	error_tol = 1e-8	   breaking tolerance for the multigrid
	maxIter = 1000		   maximum number of iterations
