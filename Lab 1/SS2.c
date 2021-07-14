#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.141592653589793238462643383279502884;

int main(int argc, char *argv[]) {
  
  double x, y;
  double r, theta;
  scanf("%lf %lf", &r, &theta);
  x = r*cos(theta);
  y = r*sin(theta);
  printf("%.2lf %.2lf\n", x, y);
  return 0;
}
