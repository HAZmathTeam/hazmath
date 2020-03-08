#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
#ifndef DIM
#define DIM 3
#endif

INT reverse(void *arr,INT length, size_t elsize)
{
  /* permutes a[0],...a_[length-1] to a[length-1],...a[0]. */
  INT i,j,k,nnn=(INT)(length/2);
  void *swap=(void *)malloc(elsize);
  //  reverses ordering in an INT array;
  void *arrk=arr+elsize*(length-1);
  void *arri=arr;
  for(i=0;i<nnn;i++){
    memcpy(swap,arri,elsize);
    memcpy(arri,arrk,elsize);
    memcpy(arrk,swap,elsize);
    arri+=elsize;
    arrk-=elsize;
  }
  if(swap)free(swap);
  return 0;
}
INT main(INT argc, char **argv)
{
  INT i=-1,j=-1,k=-1,kperm,dim=DIM;
  INT iii[10]={2,8,1,9,3,6,4,5,7,10};
  REAL rrr[10]={2.,8.,1.,9.,3.,6.,4.,5.,7.,10.};
  const char *ch[10]={"z2.","cav8.","vv1.","9.yy","3.xx","cc6.","rr4.","jj5.","ff7.","fu10."};
  fprintf(stdout,"\nint [");
  for(i=0;i<10;i++)fprintf(stdout," %d ",iii[i]);
  fprintf(stdout,"]\n");
  reverse(iii,10,sizeof(iii[0]));
  fprintf(stdout,"\nint2 [");
  for(i=0;i<10;i++)fprintf(stdout," %d ",iii[i]);
  fprintf(stdout,"]\n");
  //
  fprintf(stdout,"\nreal [");
  for(i=0;i<10;i++)fprintf(stdout," %3.1f ",rrr[i]);
  fprintf(stdout,"]\n");
  reverse(rrr,10,sizeof(rrr[0]));
  fprintf(stdout,"\nreal2 [");
  for(i=0;i<10;i++)fprintf(stdout," %3.1f ",rrr[i]);
  fprintf(stdout,"]\n");
  //
  fprintf(stdout,"\nchar [");
  for(i=0;i<10;i++)fprintf(stdout," %s ",ch[i]);
  fprintf(stdout,"]\n");
  reverse(ch,10,sizeof(ch[0]));
  fprintf(stdout,"\nchar2 [");
  for(i=0;i<10;i++)fprintf(stdout," %s ",ch[i]);
  fprintf(stdout,"]\n");
  return 0;
}
