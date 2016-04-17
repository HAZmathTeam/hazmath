/*
 *  io.c
 *
 *  Created by James Adler and Xiaozhe Hu on 3/6/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 *  Routines for reading and writing to file or screen.
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

void icsr_print_matlab(FILE* fid,iCSRmat *A)
{
  /* prints a csr matrix in matlab output*/    
  INT i,j1,j2,j; /* Loop Indices */

  for(i=0;i<A->row;i++) {
    j1 = A->IA[i]-1;
    j2 = A->IA[i+1]-1;
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%d\t%d\n",i+1,A->JA[j]);
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


void dvec_write (const char *filename,
                      dvector *vec)
{
    /**
     * \fn void dvec_write (const char *filename, dvector *vec)
     *
     * \brief Write a dvector to disk file
     *
     * \param vec       Pointer to the dvector
     * \param filename  File name
     *
     * \author Xiaozhe Hu
     * \date   03/02/2016
     */
    
    INT m = vec->row, i;
    
    FILE *fp = fopen(filename,"w");
    
    if ( fp == NULL ) {
        printf("### ERROR: Cannot open %s!\n", filename);
        chkerr(ERROR_OPEN_FILE, __FUNCTION__);
    }
    
    printf("%s: writing to file %s...\n", __FUNCTION__, filename);
    
    fprintf(fp,"%d\n",m);
    
    for ( i = 0; i < m; ++i ) fprintf(fp,"%0.15e\n",vec->val[i]);
    
    fclose(fp);
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

/*** Auxillary Files (some from Ludmil) *******************************************************/

/****************************************************************************************/
void rveci_(FILE *fp, INT *vec, INT *nn)       
/* reads a vector of integers of size nn from a file fp*/
{
	
  INT n;
  INT *vec_end;
  n = *nn;
  vec_end  =  vec + n;
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%i",vec);
  //fprintf(stdout,"Read %d INTEGERS", n);
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void rvecd_(FILE *fp,  REAL *vec, INT *nn)
/* reads a vector of REALS of size nn from a file fp*/
{
  INT n;
  REAL *vec_end;  
  n= *nn;
  vec_end =  vec + n;
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%lg",vec);
  //fprintf(stdout,"Read %d REALS", n);
  return;
}
/****************************************************************************************/

/****************************************************************************************/
FILE* HAZ_fopen( char *fname, char *mode )
{
  /* ..............................................................
     . A graceful version of fopen(). It checks if the file has .
     . been successfully opened.  If  that is  not  the case  a .
     . message is printed and the program is exited.            .
     .............................................................. */

  FILE   *fp;

  fp = fopen(fname,mode);
  if ( fp == NULL ) {
    fprintf(stderr,"Cannot open %s  -Exiting\n",fname);
    exit(255);
  }
  return fp;
}
/****************************************************************************************/
