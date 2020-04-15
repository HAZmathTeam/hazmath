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
 * \return A a pointer to a dCSRmat matrix. All of this matrix can be
 *         freed by free((void *)A) or even free(A).
 *
 *  \note: created by ludmil zikatanov on 20200412
 */
dCSRmat *dcsr_create_w (const INT m,		\
			const INT n,		\
			const INT nnz)
{
  dCSRmat *A=NULL;
  size_t structby=sizeof(dCSRmat);// size of the struct
  size_t realby=sizeof(REAL),intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+3*intby; //at least space for structure. 
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 ) 
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*realby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  A=(dCSRmat *)w;
  w+=1*structby; 
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (REAL *)w;
    w+=nnz*realby;// end of it. 
  }
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
  // this should never be used, use dcsr_create_w
  dCSRmat *Atmp=dcsr_create_w(m,n,nnz);
  memcpy(A,Atmp,1*sizeof(dCSRmat));
  return;
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
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  // arrays;
  A->rowind = (INT *)w;
  w+=nnz*intby;
  A->colind = (INT *)w;
  w+=nnz*intby;
  A->val    = (REAL *)w;
  w+=nnz*realby; // end of it....
  return A;
}
/*************************************************************************/
/***********************************************************************************************/
/*!
 * \fn dvector dvec_create_w (const INT m)
 *
 * \brief Create a dvector of given length: use void array to pack the
 * structure in.
 *
 * \param m    length of the dvector
 *
 * \return pointer to u   The new dvector
 *
 */
dvector *dvec_create_w(const INT m)
{
  dvector *u=NULL;
  size_t structby=sizeof(dvector);// size of the struct
  size_t realby=sizeof(REAL),intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+1*intby; //space for structure and size. 
  if (m > 0 )
    total+=m*realby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  u=(dvector *)w;
  w+=1*structby; 
  INT *mm=(INT *)w;
  u->row = mm[0]=m;
  w+=1*intby;
  u->val = NULL;
  if ( m > 0 ) {
    u->val = (REAL *)w;
    w+=m*realby;//end
  }
  //  fprintf(stdout,"\nVVVVVV=%d, %e\n",u->row,u->val[0]);fflush(stdout);
  return u;
}
/***********************************************************************************************/
/*!
 * \fn void dvec_alloc_w(const INT m, dvector *u)
 *
 * \brief Allocate a dvector of given length in a void array. 
 *
 * \param m    length of the dvector
 * \param u    Pointer to dvector (OUTPUT)
 *
 */
void dvec_alloc_w(const INT m,
                 dvector *u)
{
  //  if(u!=NULL){
  //    free(u);// u->val should not be defined here. 
  //    u=NULL;
  //  }
  dvector *utmp=dvec_create_w(m);
  memcpy(u,utmp,sizeof(dvector));
  //  fprintf(stdout,"\nUUUUU=%d, %e\n",u->row,u->val[0]);fflush(stdout);
  return;
}

/***********************************************************************************************/
/***********************************************************************************************/
/**
 * \fn iCSRmat icsr_create_w (const INT m, const INT n, const INT nnz)
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
  iCSRmat *A=NULL;
  size_t structby=sizeof(iCSRmat);// size of the struct
  size_t intby=sizeof(INT);// size of ints
  size_t total=1*structby; //space for the structure. 
  if ( m > 0 )
    total+=(m+1)*intby;
  if ( n > 0 ) 
    total+=nnz*intby;
  if ( nnz > 0 )
    total+=nnz*intby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  A=(iCSRmat *)w;
  w+=1*structby; 
  A->IA = NULL;
  A->JA = NULL;
  A->val = NULL;
  INT *mn_nnz=(INT *)w;  
  A->row=mn_nnz[0]=m; A->col=mn_nnz[1]=n; A->nnz=mn_nnz[2]=nnz;
  w+=3*intby;
  if ( m > 0 ) {
    A->IA = (INT *)w;
    w+=(m+1)*intby;
  }
  if ( n > 0 ) {
    A->JA = (INT *)w;
    w+=nnz*intby;
  }
  if ( nnz > 0 ) {
    A->val = (INT *)w;
    w+=nnz*intby; // end of it
  }
  return A;
}
/***********************************************************************************************/
/*!
 * \fn void ivec_alloc_w (const INT m, ivector *u)
 *
 * \brief Allocate an ivector of given length and pack it in a void array;
 *
 * \param m   length of the ivector
 * \param u   Pointer to ivector (OUTPUT)
 *
 */
/************************************************************/
/*!
 * \fn ivector *ivec_create_w (const INT m)
 *
 * \brief Create an ivector of given length
 *
 * \param m   length of the ivector
 *
 * \return u  The new ivector
 *
 */
ivector *ivec_create_w(const INT m)
{
  ivector *u=NULL;
  size_t structby=sizeof(ivector);// size of the struct
  size_t intby=sizeof(INT);// size of ints and reals
  size_t total=1*structby+1*intby; //space for structure. 
  if (m > 0 )
    total+=m*intby;
  void *w=(void *)calloc(total/sizeof(char),sizeof(char));
  u=(ivector *)w;
  w+=1*structby; 
  INT *mm=(INT *)w;
  u->row = mm[0]=m;
  w+=1*intby;
  u->val = NULL;
  if ( m > 0 ) {
    u->val = (INT *)w;
    w+=m*intby;
  }
  return u;
}
/**********************************************************************/
void ivec_alloc_w (const INT m,			\
                 ivector *u)
{
  /* if(u!=NULL){ */
  /*   free(u);// u->val should not be defined here.  */
  /*   u=NULL; */
  /* } */
  ivector *utmp=ivec_create_w(m);
  memcpy(u,utmp,sizeof(ivector));
  return;
}
/*EOF*/
