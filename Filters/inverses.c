#include <stdio.h>
#include <stdlib.h>

typedef unsigned int uint;

int main(int argc, char *argv[]) {
  uint prime = 40961;
  uint omega = 243;
  for (int n=8192; n > 1; n/=2) {
    /* find inverse m of n, i.e. (m*n)%prime == 1 */
    int m;
    for (m=1; (m*n)%prime != 1; m++);
    printf("%u: omega=%u, inverse=%u\n", n, omega, m);
    /* Reduction property: if omega is a primitive root of unity of 
     * order 2N, then omega*omega is a primitive root of order N.
     */
    omega = (omega*omega)%prime;
  }
  return 0;
}
