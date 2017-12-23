#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hazmath.h"
#include "grid_defs.h"
/**********************************************************************
 * Routines to save to file the meshin different formats.  Outputs in
 * a "HAZMAT" format or VTK format. 
 *
 *************************************************************************/

/* read */
void rveci(FILE *fp, INT *vec, INT *nn)       
  /* reads a vector of integers of size nn from a file fp*/
{

  INT n,dummy;
  INT *vec_end;
  n = *nn;
  vec_end  =  vec + n;
  for ( ; vec < vec_end; ++vec)
    dummy=fscanf(fp,"%i",vec);
  //fprintf(stdout,"Read %d INTEGERS", n);
  return;
}
void rvecd(FILE *fp,  REAL *vec, INT *nn)
  /* reads a vector of REALS of size nn from a file fp*/
{
  INT n,dummy;
  REAL *vec_end;  
  n= *nn;
  vec_end =  vec + n;
  for ( ; vec < vec_end; ++vec)
    dummy=fscanf(fp,"%lg",vec);
  //fprintf(stdout,"Read %d REALS", n);
  return;
}
void prtmat(const  INT n, const INT m, REAL *A,const char *varname)
{
  /*prints a matrix A with (n) rows and (m) columns in matlab format */
  INT i,j,n1=n-1;
  if(varname==NULL){
    fprintf(stdout,"\nA=[");
  }else{
    fprintf(stdout,"\n%s=[",varname);
  }
  for (i = 0; i<n;i++){
    for(j=0;j<m;j++){
      fprintf(stdout,"%23.16e ", A[m*i+j]);
    }
    if(i!=n1){
      fprintf(stdout,";");
    }else{
      fprintf(stdout,"];\n");
    }
  }
}

FILE *HAZ_fopen1( const char *fname, const char *mode )
{
  FILE   *fp;
  fp = fopen(fname,mode);
  if ( fp == NULL ) {
    fprintf(stderr,"ERROR: Cannot open %s\n",fname);
    exit(2);
  }
  return fp;
} 
