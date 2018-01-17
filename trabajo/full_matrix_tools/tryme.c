#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hazmath.h"
INT main(int argc, char **argv){
  INT m=4,n=11,p=7;
  INT i,j,k;
  REAL *xm=calloc(2*n+2*m,sizeof(REAL));
  REAL *xn=xm+m;
  REAL *ym=xn+n;
  REAL *yn=ym+m;
  //
  REAL *a=calloc(m*n+n*p+m*p,sizeof(REAL));
  REAL *b=a+m*n;
  REAL *c = b+n*p;
  FILE *fp=fopen("abc.input","r");
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      k=fscanf(fp,"%lg", (a+i*n+j));    
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      k=fscanf(fp,"%lg", (b+i*p+j));    
    }
  }
  for (i=0;i<m;i++){
    for (j=0;j<p;j++){
      k=fscanf(fp,"%lg", (c+i*p+j));    
    }
  }
  for (i=0;i<m;i++){
    k=fscanf(fp,"%lg", (xm+i));    
  }
  for (i=0;i<n;i++){
    k=fscanf(fp,"%lg", (xn+i));    
  }
  for (i=0;i<m;i++){
    k=fscanf(fp,"%lg", (ym+i));    
  }
  for (i=0;i<n;i++){
    k=fscanf(fp,"%lg", (yn+i));    
  }
  fprintf(stdout,"\n");
  print_full_mat(m,n,a,"aa");
  print_full_mat(n,p,b,"bb");
  print_full_mat(m,1,xm,"xxm");
  print_full_mat(n,1,xn,"xxn");
  abybfull(m,p,c,a,b,n);  
  print_full_mat(m,p,c,"cc");
  abyvfull(m,ym,a,xn,n);
  print_full_mat(m,1,ym,"yym");
  fprintf(stdout,"\n");
  atbyvfull(m,yn,a,xm,n);
  print_full_mat(n,1,yn,"yyn");
  fprintf(stdout,"\n");
  free(a);
  free(xm);
  return 0;
}

