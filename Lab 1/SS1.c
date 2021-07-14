#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.141592653589793238462643383279502884;

int main(int argc, char *argv[]) {
  
  int x, y;
  double r, theta;
  scanf("%d %d", &x, &y);
  r = sqrt(x*x + y*y);
  theta = atan2(y, x);
	if(theta >= pi) {
		theta = theta - 2*pi;
	}
	if(theta < -1*pi) {
		theta = theta + 2*pi;
	}
  printf("%.2lf %.2lf\n", r, theta);
  return 0;
}
