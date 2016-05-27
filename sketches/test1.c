#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double func(double t, double y)
{
  return t;
}

double* eulerEDO(double (*f)(double, double), double ti, double yi, double tf, int N)
{
  int i;
  double t = ti;
  double h = (tf - ti)/(N-1);
  double* y;
  y = (double*) malloc(sizeof(double)*N);
  
  y[0] = yi;
  i = 1;
  while (i < N)
  {
    y[i] = y[i-1] + h*(*f)(t, y[i-1]);
    t += h;
    i++;
  }
  return y;
}

void printDoubleVector(double* vec, int length)
{
  int i;
  for(i=0; i<length; i++)
    printf("%f ", vec[i]);
  printf("\n");
}

void saveDoubleVectorData(double* vec, int length, char* filename)
{
  int i;
  FILE* fp;
  fp = fopen(filename, "w");
  
  for(i=0; i<length; i++)
  {
    fprintf(fp, "%f \n", vec[i]);
  }
  fclose(fp);
}

int main(int argc, char* argv[])
{
  int N = 1000;
  double ti = 0;
  double tf = 2*3.141592;
  double yi = 1;
  double* y = eulerEDO(&func, ti, yi, tf, N);
  
  saveDoubleVectorData(y, N, "test1.dat");
 
  free(y);
  return 0;
}