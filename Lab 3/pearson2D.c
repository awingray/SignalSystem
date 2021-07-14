#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#define FORWARD  1
#define REVERSE -1

/* The following threshold should be tweaked for the application */
#define THRESHOLD 0.9

/** Dynamic 2D arrays *******************************/

int **allocIntArr2D(int width, int height) {
  int **arr;
  arr = malloc(height*sizeof(int *));
  if (arr) {
    arr[0] = calloc(width*height, sizeof(int));
    if (arr[0] == NULL) {
      free(arr);
      arr = NULL;
    } else {
      int i;
      for (i=1; i<height; ++i) {
        arr[i] = arr[i-1] + width;
      }
    }
  }
  return arr;
}

void freeIntArr2D(int **arr) {
  free(arr[0]);
  free(arr);
}

double **allocDoubleArr2D(int width, int height) {
  double **arr;
  arr = malloc(height*sizeof(double *));
  if (arr) {
    arr[0] = calloc(width*height, sizeof(double));
    if (arr[0] == NULL) {
      free(arr);
      arr = NULL;
    } else {
      int i;
      for (i=1; i<height; ++i) {
        arr[i] = arr[i-1] + width;
      }
    }
  }
  return arr;
}

void freeDoubleArr2D(double **arr) {
  free(arr[0]);
  free(arr);
}

/*********** images *********/

struct GrayImage {
  int width, height;
  int **im;
};

typedef struct GrayImage GrayImage;

int allocImage(int width, int height, GrayImage *image) {
  image->width = width;
  image->height = height;
  image->im = allocIntArr2D(width, height);
  if (image->im == NULL) {
    return -1;
  }
  return 0;
}

void freeImage(GrayImage image) {
  freeIntArr2D(image.im);
}

int readPGM(char *filename, GrayImage *image) {
  FILE *f;
  int w, h, r, c, k;
  unsigned char *buf;
  int **im;

  f = fopen(filename, "rb");
  fscanf (f, "P5");
  do {
    c = fgetc(f);
  } while (c != '\n');

  c = fgetc(f);
  while (c == '#') {
    do {
      c = fgetc(f);
    } while (c != '\n');
    c = fgetc(f);
  }
  ungetc (c, f);

  fscanf (f, "%d %d", &w, &h);
  do {
    c = fgetc(f);
  } while (c != '\n');

  fscanf (f, "255");
  do {
    c = fgetc(f);
  } while (c != '\n');

  if (allocImage(w, h, image) == -1) {
    fclose(f);
    return -1;
  }

  im = image->im;
  buf = malloc(w*h*sizeof(unsigned char));
  if (buf) {
    fread(buf, sizeof(unsigned char), w*h, f);
    k = 0;
    for (r=0; r<h; ++r) {
      for (c=0; c<w; ++c) {
        im[r][c] = buf[k++];
      }
    }
    free(buf);
  } else {
    for (r=0; r<h; ++r) {
      for (c=0; c<w; ++c) {
        im[r][c] = fgetc(f);
      }
    }
  }

  fclose(f);
  return 0;
}

int writePGM(char *filename, GrayImage image) {
  FILE *f;
  int w, h, **im;
  unsigned char *buf;
  int r, c, k;

  f = fopen(filename, "wb");
  if (!f) return -1;
  fprintf (f, "P5\n");
  w = image.width;
  h = image.height;
  fprintf (f, "%d %d\n255\n", w, h);

  im = image.im;
  buf = malloc(w*h*sizeof(unsigned char));
  if (buf) {
    k = 0;
    for (r=0; r<h; ++r) {
      for (c=0; c<w; ++c) {
        buf[k++] = im[r][c];
      }
    }
    fwrite(buf, sizeof(unsigned char), w*h, f);
    free(buf);
  } else {
    for (r=0; r<h; ++r) {
      for (c=0; c<w; ++c) {
        fputc(im[r][c], f);
      }
    }
  }
  fclose(f);
  return 0;
}

/*************************************/

void mirror(GrayImage image) {
  int w = image.width;
  int h = image.height;
  int r, c, tmp;
  int **im = image.im;

  /* mirror/reflect each row */
  for (r=0; r<h; r++) {
    for (c=0; c<w/2; c++) {
      tmp = im[r][c];
      im[r][c] = im[r][w-1-c];
      im[r][w-1-c] = tmp;
    }
  }

  /* mirror/reflect rows */
  for (r=0; r<h/2; r++) {
    for (c=0; c<w; c++) {
      tmp = im[r][c];
      im[r][c] = im[h-1-r][c];
      im[h-1-r][c] = tmp;
    }
  }  
}

int powerOfTwo(int n) {
  int p2;
  p2 = 1;
  while (p2 < n) {
    p2 *= 2;
  }
  return p2;
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  FORWARD gives forward transform
   dir =  REVERSE gives reverse transform
*/

void fft1D(int direction, int len, double *re, double *im) {
   long i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;
   int bits;
   /* compute the number of bits */
   bits = 1;
   l = 2;
   while (l < len) {
     l = 2*l;
     bits++;
   }
   /* Do the bit reversal */
   i2 = len >> 1;
   j = 0;
   for (i=0; i<len-1; i++) {
     if (i < j) {
       tx = re[i];
       ty = im[i];
       re[i] = re[j];
       im[i] = im[j];
       re[j] = tx;
       im[j] = ty;
     }
     k = i2;
     while (k <= j) {
       j -= k;
       k >>= 1;
     }
     j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0; l<bits; l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<len;i+=l2) {
            i1 = i + l1;
            t1 = u1 * re[i1] - u2 * im[i1];
            t2 = u1 * im[i1] + u2 * re[i1];
            re[i1] = re[i] - t1; 
            im[i1] = im[i] - t2;
            re[i] += t1;
            im[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (direction == FORWARD)
        c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for reverse transform */
   if (direction == REVERSE) {
      for (i=0;i<len;i++) {
         re[i] /= len;
         im[i] /= len;
      }
   }
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir can be FORWARD or REVERSE.
*/

void fft2D(int direction, int width, int height, double **re, double **im) {
  /* width and height must both be powers of two ! */
  /* Note that the parameters re and im denotes the REal part and 
   * the IMaginary part of the 2D image.
   */
  printf("fft2D: YOU HAVE TO IMPLEMENT THE BODY OF THIS FUNCTION YOURSELF\n");
  printf("MAKE USE OF THE fft1D() FUNCTION\n");
  exit(0);
}

void fftCorrelator(GrayImage image, GrayImage mask,
		   int corrWidth, int corrHeight, double **corr) {
  double **imreal, **imimag;
  double **maskreal, **maskimag;
  double realPart, imagPart;
  int r, c;

  int width = powerOfTwo(image.width + mask.width - 1);
  int height = powerOfTwo(image.height + mask.height - 1);

  /* allocate real/imaginary parts */
  imreal   = allocDoubleArr2D(width, height);
  imimag   = allocDoubleArr2D(width, height);
  maskreal = allocDoubleArr2D(width, height);
  maskimag = allocDoubleArr2D(width, height);

  /* copy image into real part */
  for (r=0; r<image.height; ++r) {
    for (c=0; c<image.width; ++c) {
      imreal[r][c] = image.im[r][c];
    }
  }

  /* copy mirror of mask into real part */
  for (r=0; r<mask.height; ++r) {
    for (c=0; c<mask.width; ++c) {
      maskreal[r][c] = mask.im[mask.height-1-r][mask.width-1-c];
    }
  }

  /* now perform FFT on both complex images */
  fft2D(FORWARD, width, height, imreal, imimag);
  fft2D(FORWARD, width, height, maskreal, maskimag);
  
  /* pairwise complex multiplication: im=im*mask */
  for (r=0; r<height; ++r) {
    for (c=0; c<width; ++c) {
      /* (a,b)*(c,d)=(ac-bd, ad+bc) */
      realPart = imreal[r][c]*maskreal[r][c] -
                 imimag[r][c]*maskimag[r][c];
      imagPart = imreal[r][c]*maskimag[r][c] +
                 imimag[r][c]*maskreal[r][c];      
      imreal[r][c] = realPart;
      imimag[r][c] = imagPart;
    }
  }

  /* inverse fft */
  fft2D(REVERSE, width, height, imreal, imimag);
  
  /* copy result into corr */
  for (r=0; r<corrHeight; ++r) {
    for (c=0; c<corrWidth; ++c) {
      corr[r][c] = imreal[r][c];
    }
  }
  /* deallocate memory */
  freeDoubleArr2D(maskimag);
  freeDoubleArr2D(maskreal);
  freeDoubleArr2D(imimag);
  freeDoubleArr2D(imreal);
}

/***********************************************/

int **prefSum(int width, int height, int **im) {
  int **ps = allocIntArr2D(width, height);
  int r, c;
  /*  first row */
  ps[0][0] = im[0][0];
  for (c=1; c<width; ++c) {
    ps[0][c] = im[0][c] + ps[0][c-1];
  }
  /* other rows */
  for (r=1; r<height; ++r) {
    /* first column */
    ps[r][0] = im[r][0] + ps[r-1][0];
    /* remaining columns */
    for (c=1; c<width; ++c) {
      ps[r][c] = ps[r-1][c] + ps[r][c-1] - ps[r-1][c-1] + im[r][c];
    }
  }
  return ps;
}

int **prefSquaredSum(int width, int height, int **im) {
  int **ps = allocIntArr2D(width, height);
  int r, c;
  /*  first row */
  ps[0][0] = im[0][0]*im[0][0];
  for (c=1; c<width; ++c) {
    ps[0][c] = im[0][c]*im[0][c] + ps[0][c-1];
  }
  /* other rows */
  for (r=1; r<height; ++r) {
    /* first column */
    ps[r][0] = im[r][0]*im[r][0] + ps[r-1][0];
    /* remaining columns */
    for (c=1; c<width; ++c) {
      ps[r][c] = ps[r-1][c] + ps[r][c-1] - ps[r-1][c-1] + im[r][c]*im[r][c];
    }
  }
  return ps;
}

int sum(int row0, int col0, int row1, int col1, int **prefsum) {
  int s=0;
  
  row0--;
  col0--;
  row1--;
  col1--;

  if (0<=row1) {
    if (0<=col1) s += prefsum[row1][col1];
    if (0<=col0) s -= prefsum[row1][col0];
  }

  if (0<=row0) {
    if (0<=col1) s -= prefsum[row0][col1];
    if (0<=col0) s += prefsum[row0][col0];
  }

  return s;
}

int sum0(int row0, int col0, int row1, int col1, int **im) {
  int r, c;
  int s = 0;
  for (r=row0; r<row1; ++r) {
    for (c=col0; c<col1; ++c) {
      s += im[r][c];
    }
  }
  return s;
}

int sqsum(int row0, int col0, int row1, int col1, int **im) {
  int r, c;
  int s = 0;
  for (r=row0; r<row1; ++r) {
    for (c=col0; c<col1; ++c) {
      s += im[r][c]*im[r][c];
    }
  }
  return s;
}

int pearsonCorrelator(GrayImage image, GrayImage mask,
		        int corrWidth, int corrHeight, double **corr) {
  int **im  = image.im;
  int **ms  = mask.im;

  int r, c;
  int dr, dc; /* delay in r and c */
  int ir0, ic0, ir1, ic1;
  int mr0, mc0, mr1, mc1;

  int npoints, sx, sx2, sy, sy2;
  double pc, denom, mx, my;
  int **imprefsum;
  int **msprefsum;
  int **imprefsqsum;
  int **msprefsqsum;

  /* compute standard correlation (no pearson yet!) */
  fftCorrelator(image, mask, corrWidth, corrHeight, corr);

  /* compute prefix sum of image (a.k.a. summed are table) */
  imprefsum = prefSum (image.width, image.height, im);
  msprefsum = prefSum (mask.width, mask.height, ms);
  imprefsqsum = prefSquaredSum (image.width, image.height, im);
  msprefsqsum = prefSquaredSum (mask.width, mask.height, ms);

  for (r=0; r<corrHeight; ++r) {
    dr = r - mask.height + 1; 
    for (c=0; c<corrWidth; ++c) {
      dc = c - mask.width + 1;

      /* compute bounding box of overlap in image */
      ir0 = MAX(dr, 0);
      ic0 = MAX(dc, 0);
      ir1 = MIN(image.height, dr+mask.height);
      ic1 = MIN(image.width, dc+mask.width);

      /* compute bounding box of overlap in mask */
      mr0 = MAX(-dr, 0);
      mc0 = MAX(-dc, 0);
      mr1 = MIN(mask.height, image.height-dr);
      mc1 = MIN(mask.width, image.width-dc);

      /* compute number of points in overlap */
      npoints = (ir1-ir0)*(ic1-ic0);

      /* compute pearson coefficient for delay (dr,dc) */
      sx = sum(ir0, ic0, ir1, ic1, imprefsum);
      sx2 = sum(ir0, ic0, ir1, ic1, imprefsqsum);
      
      sy = sum(mr0, mc0, mr1, mc1, msprefsum);
      sy2 = sum(mr0, mc0, mr1, mc1, msprefsqsum);
      
      mx = (double)sx/npoints;
      my = (double)sy/npoints;

      if ((sx2 == mx*sx) || (sy2 == my*sy)) {
        pc = 0;
      } else {
        denom = (sx2-mx*sx)*(sy2-my*sy);
	denom = sqrt(denom);
        pc = corr[r][c] - (((double)sx)/npoints)*sy;
        pc /= denom;
      }
      corr[r][c] = pc;
    }
  }

  freeIntArr2D(msprefsqsum);
  freeIntArr2D(imprefsqsum);
  freeIntArr2D(msprefsum);
  freeIntArr2D(imprefsum);
}

int drawbox(int y0, int x0, int y1, int x1, int**im) {
  int x, y, newbox = 1;
  for (x=x0; x<x1; x++) {
    if (im[y0][x] == 255 || im[y1][x] == 255) {
      newbox = 0;
    }
    im[y0][x] = im[y1][x] = 255;
  }
  for (y=y0+1; y<y1-1; y++) {
    if (im[y][x0] == 255 || im[y][x1] == 255) {
      newbox = 0;
    }    
    im[y][x0] = im[y][x1] = 255;
  }
  return newbox;
}

int match(double p, int mw, int mh,
	  double **pcorr, GrayImage image) {
  int **im  = image.im;
  int r, c, cnt = 0;
  for (r=mh; r<image.height; ++r) {
    for (c=mw; c<image.width; ++c) {
      if (pcorr[r][c] >= THRESHOLD) {
        cnt += drawbox(r-mh, c-mw, r, c, im);
      }
    }
  }
  return cnt;
}

int main (int argc, char **argv) {
  GrayImage image;
  GrayImage mask;
  double **corr;

  /* read input image */
  if (readPGM("emmius.pgm", &image) == -1) {
    fprintf (stderr, "Error: opening image file 'emmius.pgm' failed\n");
    return EXIT_FAILURE;
  }

  /* read mask image */  
  if (readPGM("M.pgm", &mask) == -1) {
    fprintf (stderr, "Error: opening mask image file 'M.pgm' failed\n");
    return EXIT_FAILURE;
  }

  /* create correlation image */
  int corrWidth = image.width + mask.width - 1;
  int corrHeight = image.height + mask.height - 1;
  corr = allocDoubleArr2D(corrWidth, corrHeight);

  /* correlate */
  pearsonCorrelator(image, mask, corrWidth, corrHeight, corr);
  
  /* output image and print number of matches */
  int cnt = match(THRESHOLD, mask.width, mask.height, corr, image);
  printf("%d\n", cnt);
  writePGM("match.pgm", image);

  /* deallocate images */
  free(corr[0]);
  free(corr);
  freeImage(mask);
  freeImage(image);

  return EXIT_SUCCESS;
}
