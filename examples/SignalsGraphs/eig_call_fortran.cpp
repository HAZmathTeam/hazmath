#if defined (__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

  void  dsyev_(char *chv, char *chu, int *ni, double *a, int *nj, double *w, double *work, int *n5, int *info );
  
  void drvdsyev(int ni, double *a, int nj, double *w, int *info ){
    // char *chv=(char *)malloc(sizeof(char));
    // char *chu=(char *)malloc(sizeof(char));
    // chv="V";
    // chu="U";
    int n5=5*ni;
    char chu='U',chv='V';
    double *work=(double *)calloc(n5,sizeof(double));
    dsyev_(&chv, &chu, &ni, a, &nj, w, work, &n5, info );
    if(work) free(work);
    return;
  }
}
