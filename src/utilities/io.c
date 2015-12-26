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


void dcsr_write_dcoo (const char *filename,
                      dCSRmat *A)
{
    /**
     * \fn void dcsr_write_dcoo (const char *filename, dCSRmat *A)
     *
     * \brief Write a matrix to disk file in IJ format (coordinate format)
     *
     * \param A         pointer to the dCSRmat matrix
     * \param filename  char for vector file name
     *
     * \note
     *
     *      The routine writes the specified REAL vector in COO format.
     *      Refer to the reading subroutine \ref fasp_dcoo_read.
     *
     * \note File format:
     *   - The first line of the file gives the number of rows, the
     *   number of columns, and the number of nonzeros.
     *   - Then gives nonzero values in i j a(i,j) format.
     *
     */
    
    const INT m = A->row, n = A->col;
    INT i, j;
    
    FILE *fp = fopen(filename, "w");
    
    if ( fp == NULL ) {
        printf("### ERROR: Cannot open %s!\n", filename);
        chkerr(ERROR_OPEN_FILE, __FUNCTION__);
    }
    
    printf("%s: writing to file %s...\n", __FUNCTION__, filename);
    
    fprintf(fp,"%d  %d  %d\n",m,n,A->nnz);
    for ( i = 0; i < m; ++i ) {
        for ( j = A->IA[i]; j < A->IA[i+1]; j++ )
            fprintf(fp,"%d  %d  %0.15e\n",i,A->JA[j],A->val[j]);
    }
    
    fclose(fp);
}

