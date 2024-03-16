
/*
 *  Scales an array
 */
  void
  dscal(int n,double sa,double sx[], int incx){
    int i,j;

    if (incx == 1) {
      for (i=0; i<n; i++)
        sx[i] *= sa;
    } else {
      j = 0;
      for (i=0; i<n; i++) {
        sx[j] *= sa;
        j += incx; // sequential: j depends on previous j
        // write j = i * incx or for (j = 0: j < n * incx; j += incx)
      }
    }
  }
