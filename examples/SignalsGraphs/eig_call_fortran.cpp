#if defined (__cplusplus)
extern "C" {
#endif

#include <stdlib.h>

  void  dsyev_(char *chv, char *chu, int *ni, double *a, int *nj, double *w, double *work, int *n5, int *info );
  
  void drvdsyev(int ni, double *a, int nj, double *w, int *info ){
    int n5=5*ni;
    char chu='U',chv='V';
    double *work=(double *)calloc(n5,sizeof(double));
    dsyev_(&chv, &chu, &ni, a, &nj, w, work, &n5, info );
    if(work) free(work);
    return;
  }
}
