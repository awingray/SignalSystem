#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.141592653589793238462643383279502884;

void *safeMalloc(int sz) {
  void *p = calloc(sz, 1);
  if (p == NULL) {
    fprintf(stderr, "Fatal error: safeMalloc(%d) failed.\n", sz);
    exit(EXIT_FAILURE);
  }
  return p;
}

int *makeIntArray(int n) {
  /* allocates dynamic int array of size/length n */
  return safeMalloc(n*sizeof(int));
}

double *makeDoubleArray(int n) {
  /* allocates dynamic double array of size/length n */
  return safeMalloc(n*sizeof(double));
}

void destroyArray(void *p) {
  free(p);
}

void printIntArray(int length, int *arr) {
  printf("[");
  if (length > 0) {
    printf("%d", arr[0]);
    for (int i=1; i<length; i++) {
      printf(",%d", arr[i]);
    }
  }
  printf("]\n");
}

void printDoubleArray(int length, double *arr) {
  printf("[");
  if (length > 0) {
    printf("%.2lf", arr[0]);
    for (int i=1; i<length; i++) {
      printf(",%.2lf", arr[i]);
    }
  }
  printf("]\n");
}

void demo() {
  int n = 10;
  int *intArray;
  double *dblArray;

  intArray = makeIntArray(n);
  dblArray = makeDoubleArray(n);
  for (int i=0; i < n; i++) {
    intArray[i] = i;
    dblArray[i] = 1.0/(i+1);
  }
  printf("intArray="); printIntArray(n, intArray);
  printf("dblArray="); printDoubleArray(n, dblArray);
  destroyArray(dblArray);
  destroyArray(intArray);
}

int main(int argc, char *argv[]) {
  /* demo();   replace by your own code */
  
  int x, y;
  double r, t;
  scanf("%d %d", &x, &y);
  r = sqrt(x*x+y*y);
  t = atan(y/x);
  printf("%f %f", r, t);
  return 0;
}
