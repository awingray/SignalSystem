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

void printSignal(uint len, double *x) {
  printf("%d: [", len);
  if (len > 0) {
    printf("%.5lf", x[0]);
    for (int i=1; i < len; i++) printf(",%.5lf", x[i]);
  }
  printf("]\n");
}

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

uint egcd(uint a, uint b, int *p, int *q) {
  /* Extended Euclidean algorithm:
   * computes p, q such that a*p + b*q = gcd(x,y) 
   * and returns gcd(a,b).
   */
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
  /* returns b such that (a*b)%prime==1 */
  int x, y;
  egcd(a, prime, &x, &y);
  if (x < 0) {
    x += (1 + (-x)/prime)*prime;  /* note: -x to make the number unsigned! */
  }
  return (uint)x%prime;
}

int findRoot(int length, uint *omega) {
  int len;
  *omega = primroot;
  for (len=8192; len/2 >= length; len/=2) {
    *omega = ((*omega)*(*omega))%prime;
  }
  return len;
}


void recNTT(uint length, uint root, uint prime, uint *a) {
  /* precondition: length must be a power of two */
  uint i, j, len, omega, sqroot, *odd, *even;
  if (length == 1) {
    return;
  }
  /* odd/even part */
  len = length/2;
  even = malloc(len*sizeof(uint));
  odd = malloc(len*sizeof(uint));
  for (i=j=0; i<len; i++) {
    even[i] = a[j++];
    odd[i] = a[j++];
  }
  /* divide: recursive NTT */
  sqroot = (root*root)%prime;
  recNTT(len, sqroot, prime, even);
  recNTT(len, sqroot, prime, odd);
  /* conquer: merge results */
  omega = 1;
  for (i=0; i<len; i++) {
    /* invariant: omega = (root^i)%prime */
    uint h = (omega*odd[i])%prime;
    a[i] = (even[i] + h)%prime;
    a[i+len] = (even[i] + prime - h)%prime;
    omega = (omega*root)%prime;
  }
  /* clean up */
  free(even);
  free(odd);
}

void ntt(uint prime, uint omega, uint length, uint *x) {
  recNTT(length, omega, prime, x);
}

void intt(uint prime, uint omega, uint length, uint *x) {
  uint i, inv, invomega = powmod(omega, length-1, prime);
  recNTT(length, invomega, prime, x);
  inv = inverse(length, prime);
  for (i=0; i<length; i++) {
    x[i] = (x[i]*inv)%prime;
  }
}

uint *correlator1D(uint lenx, int *x, uint leny, int *y, uint *lenxy) {
  uint i, omega, lenz, *tmp, *z, *_z, _len;
	uint cnt, idx;

  *lenxy = lenx + leny - 1;  
  lenz = findRoot(*lenxy, &omega);
  
  z = calloc(lenz, sizeof(uint));
	_z = calloc(lenz, sizeof(uint));
  tmp = calloc(lenz, sizeof(uint));
  memcpy(tmp, x, lenx*sizeof(uint));
  memcpy(z, y, leny*sizeof(uint));
	
	_len = lenz - 1;
	
	for (idx=0; idx < _len; idx++, _len--) {
		cnt = tmp[idx];
		tmp[idx] = tmp[_len];
		tmp[_len] = cnt;
	}
  ntt(prime, omega, lenz, tmp);
  ntt(prime, omega, lenz, z);

  for (i=0; i<lenz; i++) z[i] = (z[i]*tmp[i])%prime;
  
  free(tmp);
  intt(prime, omega, lenz, z);
	
	_z[0] = z[lenz-1];
	for (i=1; i<lenz; i++) _z[i] = z[i-1];

	free(z);
  return (uint*)_z;
}

void sqSum(uint len, int *x, uint len0, uint len1, double *sqsum, double *mean) {
	for(int i=len0; i < len1 && i < len; i++) {
		*sqsum += pow(x[i], 2);
		*mean += x[i];
	}
	*mean /= len1 - len0;
}

double *pearsonCorrelator1D(uint lenH, int *h, uint lenX, int *x, uint *lenY) {
	uint lenz, *z;
	int i, j;
	z = correlator1D(lenH, h, lenX, x, &lenz);
	
	uint _len = lenX - lenH + 1;
	double *p = calloc(_len, sizeof(double));
	
	double meanX, meanY, sqsum1, sqsum2;
	double denom, numer;

	sqSum(lenX, h, 0, lenH, &sqsum1, &meanX);
	for(i=0; i < _len; i++) {
		meanY = 0;
		sqsum2 = 0;
		sqSum(lenX, x, i, lenH+i, &sqsum2, &meanY);
		numer = z[i]-lenH*meanX*meanY;
		denom = sqrt(sqsum1 - lenH*pow(meanX, 2))*sqrt(sqsum2 - lenH*pow(meanY, 2));
		p[i] = roundf(numer/denom*100000)/100000.;
	}
	
	return p;
}

int main(int argc, char *argv[]) {
  int *x, *h;
	double *y;
	uint lenX, lenH, lenY;
  /* the first line of the input contains the template h */
  h = readSignal(&lenH);
  /* the second line of the input contains the signal x */
  x = readSignal(&lenX);
  /* compute the correlation y = corr(h,x) */
  y = pearsonCorrelator1D(lenH, h, lenX, x, &lenY);
  /* print correlation vector */
  printSignal(lenX - lenH + 1, y);
  /* clean up */
  free(y);
  free(x);
  free(h);
  return 0;
}