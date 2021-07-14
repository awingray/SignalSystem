#include <stdio.h>
#include <stdlib.h>

int *readSignal(int *len) {
  int *x;
  char c;
  scanf("%d:", len);
  x = calloc(*len, sizeof(int));
  do c = getchar(); while (c != '[');
  if (len > 0) {
    scanf("%d", &x[0]);
    for (int i=1; i < *len; i++) scanf(",%d", &x[i]);
  }
  do c = getchar(); while (c != ']');
  return x;
}

void printSignal(int len, int *x) {
  printf("%d: [", len);
  if (len > 0) {
    printf("%d", x[0]);
    for (int i=1; i < len; i++) printf(",%d", x[i]);
  }
  printf("]\n");
}

int *correlator1D(int lenH, int *h, int lenX, int *x, int *lenY) {
  /* assumption: 0 < lenH <= lenX */
  int len = lenX - lenH + 1;  /* length steady state=lenX-lenH+1 */
  int *y = malloc(len*sizeof(int));
  for (int i=0; i<len; i++) {
    y[i] = 0;
    for (int j=0; j<lenH; j++) {
      y[i] += h[j]*x[i+j];
    }
  }
  *lenY = len;
  return y;
}

int main(int argc, char *argv[]) {
  int lenX, *x, lenH, *h, lenY, *y;
  /* the first line of the input contains the template h */
  h = readSignal(&lenH);
  /* the second line of the input contains the signal x */
  x = readSignal(&lenX);
  /* compute the correlation y = corr(h,x) */
  y = correlator1D(lenH, h, lenX, x, &lenY);
  /* print correlation vector */
  printSignal(lenY, y);
  /* clean up */
  free(y);
  free(x);
  free(h);
  return 0;
}
