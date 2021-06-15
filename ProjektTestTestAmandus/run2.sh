#! /bin/bash

echo
echo "5-point stencil v-cycle"
echo
gcc -O3 -o  grid ./main.c -lm
./grid 1023 5 jac 5 v
./grid 1023 5 gau 5 v
./grid 1023 5 sor 5 v

echo
echo "5-point stencil w-cycle"
echo
./grid 1023 5 jac 5 w
./grid 1023 5 gau 5 w
./grid 1023 5 sor 5 w


echo
echo
echo
echo
echo "9-point stencil v-cycle"
echo

./grid 1023 5 jac 9 v
./grid 1023 5 gau 9 v
./grid 1023 5 sor 9 v

echo
echo "9-point stencil w-cycle"
echo
./grid 1023 5 jac 9 w
./grid 1023 5 gau 9 w
./grid 1023 5 sor 9 w



echo
echo "Jetzt wird parallelisiert"
echo
echo "5-point stencil v-cycle"
echo
gcc -O3 -o -fopenmp  grid ./main.c -lm
./grid 1023 5 jac 5 v
./grid 1023 5 gau 5 v
./grid 1023 5 sor 5 v

echo
echo "5-point stencil w-cycle"
echo
./grid 1023 5 jac 5 w
./grid 1023 5 gau 5 w
./grid 1023 5 sor 5 w


echo
echo
echo
echo
echo "9-point stencil v-cycle"
echo

./grid 1023 5 jac 9 v
./grid 1023 5 gau 9 v
./grid 1023 5 sor 9 v

echo
echo "9-point stencil w-cycle"
echo
./grid 1023 5 jac 9 w
./grid 1023 5 gau 9 w
./grid 1023 5 sor 9 w






















exit 0
