#include <stdio.h>
#include <stdlib.h>
#include <time.h>



__global__ void residual(float* x_now, float* x_next, float* b, int N)
{
    int j;
    float sigma = 0.0;
    int row = blockIdx.x*blockDim.x+threadIdx.x;


    if (row < N)
    {
      printf(" %d ", row);
      if (row == 0)
      {
        sigma += 4*x_now[0] - x_now[N] - x_now[1]; //bottom left

      }
      else if (row == N-1)
      {
        sigma += 4*x_now[(N-1)*N] - x_now[(N-1)*N + 1] - x_now[(N-2)*N]; //top left

        for (j=1; j<N-1;j++)
          sigma += 4*x_now[(N-1)*N + j] - x_now[(N-1)*N + j-1] - x_now[(N-1)*N + j+1] - x_now[(N-2)*N + j]; //top edge

        sigma += 4*x_now[(N-1)*N + N-1] - x_now[(N-1)*N + N-2] - x_now[(N-2)*N + N-1]; //top right
      }
      else
      {
        sigma += 4*x_now[row*N] - x_now[row*N + 1] - x_now[(row-1)*N] - x_now[(row+1)*N]; //left edge (i = row)

        for (j=1; j<N-1;j++)
          sigma += 4*x_now[row*N + j] - x_now[row*N + j-1] - x_now[row*N + j+1] - x_now[(row-1)*N + j] - x_now[(row+1)*N + j]; // inner points (i=row)

        sigma += 4*x_now[row*N + N-1]   - x_now[row*N + N-2]     - x_now[(row-1)*N + N-1] - x_now[(row+1)*N + N-1]; // right edge
      }


    x_next[row] = sigma;
  }
}


// Kernel Funktion fuer die Durchfuehrung einer Jacobi-Iteration
__global__ void jacobiOnDevice(float* x_now, float* x_next, float* b, int N)
{
    int j;
    float sigma = 0.0;
    float omega = 0.6;
    float aDiag = 0.25;
    int row = blockIdx.y*blockDim.y+threadIdx.y;


    if (row < N)
    {
      if (row == 0)
      {
        sigma += 4*x_now[0] - x_now[N] - x_now[1]; //bottom left

        for (j=1; j<N-1;j++)
          sigma += 4*x_now[j] - x_now[j-1] - x_now[j+1] - x_now[N + j];  // bottom edge

        sigma += 4*x_now[N-1] - x_now[2*N-1] - x_now[N - 2];	//bottom right
      }
      else if (row == N-1)
      {
        sigma += 4*x_now[(N-1)*N] - x_now[(N-1)*N + 1] - x_now[(N-2)*N]; //top left

        for (j=1; j<N-1;j++)
          sigma += 4*x_now[(N-1)*N + j] - x_now[(N-1)*N + j-1] - x_now[(N-1)*N + j+1] - x_now[(N-2)*N + j]; //top edge

        sigma += 4*x_now[(N-1)*N + N-1] - x_now[(N-1)*N + N-2] - x_now[(N-2)*N + N-1]; //top right
      }
      else
      {
        sigma += 4*x_now[row*N] - x_now[row*N + 1] - x_now[(row-1)*N] - x_now[(row+1)*N]; //left edge (i = row)

        for (j=1; j<N-1;j++)
          sigma += 4*x_now[row*N + j] - x_now[row*N + j-1] - x_now[row*N + j+1] - x_now[(row-1)*N + j] - x_now[(row+1)*N + j]; // inner points (i=row)

        sigma += 4*x_now[row*N + N-1]   - x_now[row*N + N-2]     - x_now[(row-1)*N + N-1] - x_now[(row+1)*N + N-1]; // right edge
      }


    x_next[row] = aDiag * omega *(b[row] - sigma);
  }
}



int main(int argc, char* argv[]){


    int N = atoi(argv[1]);
    int iter = 10;
    int k;
    clock_t before = clock();

    float *x_now = (float*)malloc(N*sizeof(float));
    float *x_next = (float*)malloc(N*sizeof(float));
    float *b = (float*)malloc(N*sizeof(float));
    float res, dev_res;

    float *dev_x_now,*dev_x_next,*dev_b;


    //Allokiere Speicher im globalen Speicher der GPU
    cudaMalloc((void**)&dev_x_now,N*sizeof(float));
    cudaMalloc((void**)&dev_x_next,N*sizeof(float));
    cudaMalloc((void**)&dev_b,N*sizeof(float));

    //Füllen der Arrays auf der CPU
    for (int i=0;i<N;i++)
    {
      x_now[i] = 1.0;
      x_next[i] = 0.0;
      b[i] = 1.0;
    }

    //Kopiere Daten auf GPU in globalen Speicher
    cudaMemcpy(dev_x_now,x_now,N*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_x_next,x_next,N*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b,b,N*sizeof(float),cudaMemcpyHostToDevice);

    //Baue 2D Gitter von Blocks der Größe 32x32 Threads
    int nblocks = (N+32)/32;
    //dim3 gridDim(nblocks,nblocks);
    //dim3 blockDim(32,32);

    //Aufruf des Jacobi Verfahrens auf der GPU
    // for (k=0; k<iter; k++)
    // {
    //     if (k%2)
    //         jacobiOnDevice<<<gridDim,blockDim>>>(dev_x_next, dev_x_now, dev_b, N);
    //     else
    //         jacobiOnDevice<<<gridDim,blockDim>>>(dev_x_now, dev_x_next, dev_b, N);
    // }

    residual<<<nblocks,32>>>(dev_x_now, dev_x_next, dev_b, N);

    //Ergebnis zurück auf den Host kopieren
    // if (k%2)
    //   cudaMemcpy(x_next,dev_x_now,N*sizeof(float),cudaMemcpyDeviceToHost);
    // else
    //   cudaMemcpy(x_next,dev_x_next,N*sizeof(float),cudaMemcpyDeviceToHost);

    cudaMemcpy(x_next,dev_x_next, N*sizeof(float),cudaMemcpyDeviceToHost);


    clock_t after = clock();
    clock_t difference = clock() - before;
    int msec = difference * 1000 / CLOCKS_PER_SEC;

    printf("\nTime taken %d.%d seconds \n\n",msec/1000,msec%1000);

    printf("Result: C= ");

    for (int i=0;i<N;i++){
        printf(" %f",x_next[i]);
    }

    cudaFree(dev_x_now);
    cudaFree(dev_x_next);
    cudaFree(dev_b);

    free(x_now);
    free(x_next);
    free(b);

    return 0;
}
