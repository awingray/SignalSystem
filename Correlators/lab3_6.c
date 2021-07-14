#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned int uint;

uint prime = 40961;
uint primroot = 243;  /* primitive 8192-th root of unity */

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
uint powmod(uint base, uint exponent, uint prime) {

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

uint egcd(uint a, uint b, int *p, int *q) {

  int x=0, lx=1, y=1, ly=0;
  while (b > 0) {
    int q = a/b;
    int h = a%b;
    a = b;
    b = h;
    
    h = x;
    x = lx - q*x;
    lx = h;

    h = y;
    y = ly - q*y;
    ly = h;
    }
    *p = lx;
    *q = ly;
  return a;
}

uint inverse(uint a, uint prime) {

  int x, y;
  egcd(a, prime, &x, &y);
  if (x < 0) {
    x += (1 + (-x)/prime)*prime; 
  }
  return (uint)x%prime;
}

int findRoot(int length, uint *omega) {
  int len;
  *omega = primroot;
  for (len=8192; len/2 >= length; len/=2) {
    *omega = ((*omega)*(*omega)) % prime;
  }
  return len;
}


void recNTT(uint length, uint root, uint prime, uint *a) {

  uint i, j, len, omega, sqroot, *odd, *even;
  if (length == 1)
    return;

  len = length/2;
  even = malloc(len*sizeof(uint));
  odd = malloc(len*sizeof(uint));
  for (i=j=0; i<len; i++) {
    even[i] = a[j++];
    odd[i] = a[j++];
  }

  sqroot = (root*root)%prime;
  recNTT(len, sqroot, prime, even);
  recNTT(len, sqroot, prime, odd);

  omega = 1;
  for (i=0; i<len; i++) {
    uint h = (omega*odd[i])%prime;
    a[i] = (even[i] + h)%prime;
    a[i+len] = (even[i] + prime - h)%prime;
    omega = (omega*root)%prime;
  }

  free(even);
  free(odd);
}

void ntt(uint prime, uint omega, uint length, uint *x) {
  recNTT(length, omega, prime, x);
}

void intt(uint prime, uint omega, uint length, uint *x) {
  uint i, inv, invomega;
  invomega = powmod(omega, length-1, prime);
  recNTT(length, invomega, prime, x);
  inv = inverse(length, prime);
  for (i=0; i<length; i++) x[i] = (x[i]*inv)%prime;
}

uint *conv(uint lenx, int *x, uint leny, int *y, uint *lenxy) {
  uint i, omega, lenz, *tmp, *z;
  uint cnt, idx;

  *lenxy = lenx + leny - 1;  
  lenz = findRoot(*lenxy, &omega);
  
  z = calloc(lenz, sizeof(uint));
  tmp = calloc(lenz, sizeof(uint));
  memcpy(tmp, x, lenx*sizeof(uint));
  memcpy(z, y, leny*sizeof(uint));
	
  ntt(prime, omega, lenz, tmp);
  ntt(prime, omega, lenz, z);
  
  for (i=0; i < lenz; i++) z[i] = (z[i]*tmp[i]) % prime;
 
  free(tmp);
  intt(prime, omega, lenz, z);
  return (uint*)z;
}

int main(int argc, char *argv[]) {
  uint *x, *h, *y;
  uint lenX, lenH, lenY;
  
  h = readSignal(&lenH);
  x = readSignal(&lenX);
  
  y = conv(lenH, h, lenX, x, &lenY);
  
  printSignal(lenY, y);
  
  free(x);
  free(h);
  free(y);
  return 0;
}
