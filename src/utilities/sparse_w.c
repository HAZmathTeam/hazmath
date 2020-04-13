/*! \file src/utilities/alloc.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 5/13/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  memory allocating functions using void arrays for structures.  
 *
 *  \note: created by ludmil zikatanov on 20200412
 *  
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn dCSRmat *dcsr_create_w(const INT m, const INT n, const INT nnz)
 *
 * \brief Create a dCSRmat sparse matrix. Uses void array for the
 * whole matrix. the void array contains in first position the whole struct. 
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A the dCSRmat matrix. All of this matrix can be freed by
 *         free((void *)A) or even free(A).
 *
 *  \note: created by ludmil zikatanov on 20200412
 */
dCSRmat *dcsr_create_w (const INT m,		\
			const INT n,		\
			const INT nnz)
{
  dCSRmat *A;
  size_t structby=sizeof(dCSRmat);// size of the struct
  size_t realby=sizeof(REAL),intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby; //at least space for structure. 
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 ) 
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*realby;
  void *w=(void *)calloc(total/size(char),sizeof(char));
  dCSRmat *A=(dCSRmat *)w;
  w+=1*structby; 
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (REAL *)calloc(nnz, sizeof(REAL));
    w+=nnz*realby;
  }
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;// end of it
  return A;
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_alloc_w(const INT m, const INT n, const INT nnz, dCSRmat *A)
 *
 * \brief Allocate dCSRmat sparse matrix memory space using one void array. uses dcsr_create_w.
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 */
void dcsr_alloc_w(const INT m,
		  const INT n,
		  const INT nnz,
		  dCSRmat *A)
{
  if(A != NULL) {free(A);A=NULL};
  A=dcsr_create_w(m,n,nnz);
  return;
}
/***********************************************************************************************/
/**
 * \fn iCSRmat icsr_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create iCSRmat sparse matrix using a void array; the
 *        structure itself is part of the void array.
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   a pointer to an iCSRmat matrix
 *
 */
iCSRmat *icsr_create_w(const INT m,		\
		       const INT n,		\
		       const INT nnz)
{
  iCSRmat *A;
  size_t structby=sizeof(iCSRmat);// size of the struct
  size_t intby=sizeof(INT);// size of ints
  size_t total=1*structby; //space for the structure. 
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 ) 
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*intby;
  void *w=(void *)calloc(total/size(char),sizeof(char));
  iCSRmat *A=(iCSRmat *)w;
  w+=1*structby; 
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (INT *)calloc(nnz, sizeof(INT));
    w+=nnz*intby;
  }
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;// end of it
  return A;
}

/**
 * \fn dCOOmat *dcoo_create_w(INT m, INT n, INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space using one
 * contguous void array for all data including the structure itself 
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   A pointer to a dCOOmat matrix
 *
 */
dCOOmat *dcoo_create_w(INT m,			\
		       INT n,			\
		       INT nnz)
{
  size_t structby=sizeof(dCOOmat),realby=sizeof(REAL),intby=sizeof(INT);
  size_t total=(1*structby+(3+2*nnz)*intby+nnz*realby)/sizeof(char);
  void *w=(void *)calloc(total, sizeof(char));
  //sturture
  dCOOmat *A=(dCOOmat *)w;
  w+=1*structby;
  // arrays;
  A->rowind = (INT *)w;
  A->colind = A->rowind + nnz;
  w+=2*nnz*intby;
  A->val    = (REAL *)w;
  w+=nnz*realby; // end of it....
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  return A;
}

/*EOF*/
