/*
 *  io.c
 *
 *  Created by James Adler and Xiaozhe Hu on 3/6/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

void iarray_print(INT *vec, INT n   )
{
    /* prints a vector of integers of size nn */
    INT *vec_end;
    
    vec_end  =  vec + n;
    
    fprintf(stdout,"\n");

    for ( ; vec < vec_end; ++vec)
        fprintf(stdout, "%i\n  ",*vec);
    
    fprintf(stdout,"\n");
    
    return;
}

void array_print(REAL *vec, INT n   )
{
    /* prints a vector of integers of size nn */
    REAL *vec_end;
    
    vec_end  =  vec + n;
    
    fprintf(stdout,"\n");
    
    for ( ; vec < vec_end; ++vec)
        fprintf(stdout, "%e\n  ",*vec);
    
    fprintf(stdout,"\n");
    
    return;
}

void csr_print_matlab(FILE* fid,dCSRmat *A)
{
  /* prints a csr matrix in matlab output*/    
  INT i,j1,j2,j; /* Loop Indices */

  for(i=0;i<A->row;i++) {
    j1 = A->IA[i]-1;
    j2 = A->IA[i+1]-1;
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%d\t%d\t%25.16e\n",i+1,A->JA[j],A->val[j]);
    }
  }		
  return;
}

void dvector_print(FILE* fid,dvector *b)
{
  /* prints a csr matrix in matlab output*/    
  INT i; /* Loop Indices */

  for(i=0;i<b->row;i++) {
    fprintf(fid,"%25.16e\n",b->val[i]);
  }
  return;
}
