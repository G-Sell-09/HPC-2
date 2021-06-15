#! /bin/bash

export OMP_NUM_THREADS=1

echo
echo "5-point stencil v-cycle"
echo
gcc -O3 -o  grid ./main.c -lm
./grid 255 4 jac 5 v
./grid 255 4 gau 5 v
./grid 255 4 sor 5 v

echo
echo "5-point stencil w-cycle"
echo
./grid 255 4 jac 5 w
./grid 255 4 gau 5 w
./grid 255 4 sor 5 w

echo
echo "5-point stencil mu3-cycle"
echo
./grid 255 4 jac 5 mu 3
./grid 255 4 gau 5 mu 3
./grid 255 4 sor 5 mu 3


echo
echo "9-point stencil v-cycle"
echo

./grid 255 4 jac 9 v
./grid 255 4 gau 9 v
./grid 255 4 sor 9 v

echo
echo "9-point stencil w-cycle"
echo
./grid 255 4 jac 9 w
./grid 255 4 gau 9 w
./grid 255 4 sor 9 w

echo
echo "5-point stencil mu3-cycle"
echo
./grid 255 4 jac 9 mu 3
./grid 255 4 gau 9 mu 3
./grid 255 4 sor 9 mu 3





















exit 0
