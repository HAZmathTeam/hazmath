#include "hazmath.h"
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
typedef struct {
  INT num;// how many elements in a name
  REAL *val;
  INT  id; //for the permutation
} names;
INT compare0(names *a1,  names *a2) {
// compares two names with nc components each:
  INT k,lt=0,nc1=a1->num,nc2=a2->num;
  if(nc1>nc2) nc1=nc2;
  for(k=0;k<nc1;k++) {
    if(a1->val[k]<a2->val[k]){lt=-1;break;}
    else if(a1->val[k]>a2->val[k]){lt=1;break;}
    else{continue;}
  }
  return lt;
}
void lexsort(REAL *x, INT nr, INT nc, INT *p, INT *invp) {
  INT k;
  names *tosort = (names *)calloc(nr,sizeof(names));
  REAL *y = (REAL *)calloc(nr*nc,sizeof(REAL));
  memcpy(y,x,nr*nc*sizeof(REAL));
  REAL *ptr=y;
  for (k=0; k<nr;k++) {
    tosort[k].num=nc;
    tosort[k].val=ptr;
    tosort[k].id=k;
    ptr+=nc;
  }
  //Sort now.
  qsort((void *) tosort,nr,sizeof(names),(testit )compare0 );
  for (k=0; k<nr;k++){
    p[k]=tosort[k].id;
    invp[p[k]]=k;
//    fprintf(stdout,"\nk=%i,p[k]=%i",k,p[k]);
  }
  if(y) free(y);
  if(tosort) free(tosort);
  return;
}
