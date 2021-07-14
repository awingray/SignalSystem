#include <stdio.h>
#include <stdlib.h>

typedef unsigned int uint;

/* The following table has the following format (for each index i):
 *   N = table[i][0] (a power of two), 
 *   p = table[i][1] (a prime),
 *   w = table[i][2] (a primitive root of unity)
 */
uint table[14][3] = {
  {8192,40961,243},
  {4096,61441,39003},
  {2048,61441,16290},
  {1024,64513,57751},
  {512,64513,49440},
  {256,64513,45056},
  {128,64513,12565},
  {64,65089,53832},
  {32,65089,56855},
  {16,65521,61640},
  {8,65521,57852},
  {4,65521,41224},
  {2,65521,65520},
  {1,65521,1}
};


uint powmod(uint base, uint exponent, uint prime) {
  /* This function computes: base raised to the power exponent modulus prime
   * in math notation: (base^exponent) mod prime
   */
  uint pm = 1;
  base = base%prime;
  while (exponent > 0) {
    if (exponent%2 == 1) { /* exponent is odd */
      pm = (pm*base)%prime;
    }
    exponent /= 2;
    base = (base*base)%prime;
  } 
  return pm;
}

int main(int argc, char *argv[]) {
  int i;
  for (i=0; i < 14; i++) {
    /* note the uint (not int) */
    uint N=table[i][0], p=table[i][1], w=table[i][2];
    printf("N=%4u: p=%5u, w=%5u, w^N=%u\n", N, p, w, powmod(w, N, p)); 
  };
  
  return EXIT_SUCCESS;
}
