#! /bin/bash


gcc -O3 -o  grid ./main.c -lm
echo
echo "9-point stencil v-cycle"
echo
echo "First start with Jacobi and 255x255 grid"
echo
./grid 255 4 jac 9 v
echo
echo "Now the same 255x255 grid with Gauss-Seidel"
echo
./grid 255 4 gau 9 v
echo
echo "And now a 255x255 grid with SOR"
echo
./grid 255 4 sor 9 v
echo
echo "Use the best smoothing operator SOR for a 511x511 grid now"
echo
./grid 511 5 sor 9 v
echo
echo "And now SOR for a 751x751 grid"
echo
./grid 751 5 sor 9 v


echo
echo "We now use the 9-point stencil with the w-cycle"
echo
echo "Again first with Jacobi and a 255x255 grid"
echo
./grid 255 4 jac 9 w
echo
echo "Now the same 255x255 grid with Gauss-Seidel"
echo
./grid 255 4 gau 9 w
echo
echo "And now a 255x255 grid with SOR"
echo
./grid 255 4 sor 9 w
echo
echo "Again a 511x511 grid with SOR"
echo
./grid 511 5 sor 9 w
echo
echo "And SOR for a 751x751 grid"
echo
./grid 751 5 sor 9 w


echo
echo "Now at the end of the 9-point stencil we also try the mu3-cycle"
echo
echo "Jacobi: 255x255 grid"
echo
./grid 255 4 jac 9 mu 3
echo
echo "Gauss-Seidel: 255x255 grid"
echo
./grid 255 4 gau 9 mu 3
echo
echo "sor: 255x255 grid"
echo
./grid 255 4 sor 9 mu 3
echo
echo "sor: 511x511 grid"
echo
./grid 511 5 sor 9 mu 3
echo
echo "sor: 751x751 grid"
echo
./grid 751 5 sor 9 mu 3
echo
echo "End"
echo



exit 0
