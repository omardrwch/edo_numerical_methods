#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h> 

int func(double t, double* y, double *dy, int M)
{
  dy[0] = -sin(t);
  return 0;
}

int volterra(double t, double* y, double* dy, int M)
{
  // parameters
  double a = 1.0;
  double b = 2.0;
  double c = 3.0;
  double d = 1.4;
  
  dy[0] = a*y[0] - b*y[0]*y[1];
  dy[1] = c*y[0]*y[1] - d*y[1];

  return 0;  
}

int multiplyAddRK4 (double** Y, double* out, int col, double* k, double constant, int M)
{
  int i;
  for (i=0;i<M; i++)
  {
    out[i] = Y[i+1][col] + constant*k[i];
  }    
  return 0;
}

int rk4(int (*f)(double, double*, double*, int), double** Y, double* yi, double ti, double tf, int N, int M)
{
  int i, j;
  double t = ti;
  double h = (tf - ti)/(N-1);
  double* yn = (double*) malloc(sizeof(double*)*M);
  double* aux = (double*) malloc(sizeof(double*)*M);
  double* k1 = (double*) malloc(sizeof(double*)*M);
  double* k2 = (double*) malloc(sizeof(double*)*M);
  double* k3 = (double*) malloc(sizeof(double*)*M);
  double* k4 = (double*) malloc(sizeof(double*)*M);

  Y[0][0] = ti; 
  
  for(j=0; j<M; j++)
    Y[j+1][0] = yi[j];
  
  i = 1;
  while (i < N)
  {
    // Save value of Y[j+1][i-1] in yn[j]
    for(j=0; j<M; j++)
      yn[j] = Y[j+1][i-1];
    
    // k1 <- f(tn, yn)
    (*f)(t, yn, k1, M);
    
    // aux <- yn + (h/2)*k1
    multiplyAddRK4(Y, aux, i-1, k1, h/2, M);
    
    // k2 <- f(tn + h/2, aux)
    (*f)(t + h/2, aux, k2, M);
    
    //aux <= yn + (h/2)*k2
    multiplyAddRK4(Y, aux, i-1, k2, h/2, M);
    // k3 <- f(tn +h/2, aux)
    (*f)(t + h/2, aux, k3, M);
    
    //aux <= yn + h*k3
    multiplyAddRK4(Y, aux, i-1, k3, h, M);    
    // k4 <- f(tn + h, aux)
    (*f)(t + h, aux, k4, M);
    
    for(j=0; j<M; j++)
      Y[j+1][i] = Y[j+1][i-1] + (h/6)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
    
    t += h;
    Y[0][i] = t;
    i++;
  }

  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(yn);
  free(aux);
  return 0;
}

void saveDoubleVectorData(double** vec, int length, char* filename)
{
  int i;
  FILE* fp;
  fp = fopen(filename, "w");
  
  for(i=0; i<length; i++)
  {
    fprintf(fp, "%f %f %f\n", vec[0][i], vec[1][i], vec[2][i]);
  }
  fclose(fp);
}

int main(int argc, char* argv[])
{
  // Time evaluation
  struct timeval t1, t2;
  double elapsedTime;
  
  int i;
  int N = 10000;
  int M = 2;
  
  // Time interval
  double ti = 0;
  double tf = 20;
  // Initial condition
  double yi[] = {1.0, 1.0};
  
  // Solution matrix
  double** out = (double**) malloc(sizeof(double**)*(M+1));
  for(i = 0; i<=M; i++)
    out[i] = (double*) malloc(sizeof(double*)*N);
  
  // start timer
  gettimeofday(&t1, NULL);
  printf("\nSolving equation ...\n\n");
  
  // solve equation
  rk4(&volterra, out,  yi, ti, tf, N, M);
  
  // stop timer
  gettimeofday(&t2, NULL);
  
  // compute and print the elapsed time in ms
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  printf("Elapsed time: %f ms\n\n", elapsedTime);
  
  printf("Saving data ...\n");
  saveDoubleVectorData(out, N, "rk4.dat");
  printf("Finished.\n\n");
 
  for(i = 0; i<=M; i++)
    free(out[i]);
  free(out);
  return 0;
}