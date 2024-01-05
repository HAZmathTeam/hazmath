/*! \file src/utilities/io.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  Routines for reading and writing to file or screen.
 *
 * \note: modified by Xiaozhe Hu on 10/31/2016
 * \note: modified 20171119 (ltz) added output using scomplex
 * \note: modifided by jha on 02/22/2019 for 0-1 fix.
 */

#include "hazmath.h"

/*
 * \fn chkn(INT n, const INT nmin, const INT nmax)
 *
 * \brief checks if n is in the closed interval [nmin, nmax] and if
 *    not, changes the value of n to be the closest end of the interval
 *
 * \param n     integer to check
 * \param nmin  min allowed value
 * \param nmax  max allowed value
 *
 */
INT chkn(INT n, const INT nmin, const INT nmax)
{
  INT nnew=n;
  if(n > nmax) {
    fprintf(stdout, "\n Input value too large: %lld ; Changing to the max allowed: %lld\n",(long long )n,(long long )nmax);
    nnew=nmax;
  } else if(n < nmin) {
    fprintf(stdout, "\n Input value too small: %lld ; Changing to the min allowed: %lld\n",(long long )n,(long long )nmin);
    nnew=nmin;
  }
  return nnew;
}
/***********************************************************************************************/
/*!
 * \fn void iarray_print(INT *vec, INT n)
 *
 * \brief print an integer array on screen
 *
 * \param vec   Pointer to the INT array
 * \param n     Length of the array
 *
 */
void iarray_print(INT *vec,
                  INT n)
{
  // local variable
  INT *vec_end;
  vec_end  =  vec + n;

  fprintf(stdout,"\n");

  // main loop
  for ( ; vec < vec_end; ++vec)
    fprintf(stdout, "%lld  ",(long long )(*vec));

  fprintf(stdout,"\n");

  return;
}

/***********************************************************************************************/
/*!
 * \fn void array_print(REAL *vec, INT n)
 *
 * \brief print a REAL array on screen
 *
 * \param vec   Pointer to the REAL array
 * \param n     Length of the array
 *
 */
void array_print(REAL *vec,
                 INT n)
{
  // local variable
  REAL *vec_end;
  vec_end  =  vec + n;

  fprintf(stdout,"\n");

  for ( ; vec < vec_end; ++vec)
    fprintf(stdout, "%e  ",*vec);

  fprintf(stdout,"\n");

  return;

}

/***********************************************************************************************/
/**
 * \fn void ivector_write (const char *filename, ivector *vec)
 *
 * \brief Write a ivector to disk file in coordinate format
 *
 * \param vec       Pointer to the dvector
 * \param filename  File name
 *
 */
void ivector_write(const char *filename,
                   ivector *vec)
{
  INT m = vec->row, i;

  FILE *fp = fopen(filename,"w");

  if ( fp == NULL ) {
    check_error(ERROR_OPEN_FILE, __FUNCTION__);
  }

  fprintf(stdout,"%%%%%s: writing to file %s...\n", __FUNCTION__, filename);

  // write number of nonzeros
  fprintf(fp,"%lld %d\n",(long long )m,0);

  // write index and value each line
  for ( i = 0; i < m; ++i ) fprintf(fp,"%lld %lld\n",(long long )i,(long long )vec->val[i]);

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void dvector_print(FILE* fid,dvector *b)
 *
 * \brief print a dvector to a file
 *
 * \param fid  Pointer to the file
 * \param b    Pointer to the dvector
 *
 */
void dvector_print(FILE* fid,
                   dvector *b)
{
  // local vairable
  INT i;

  // main loop
  for(i=0;i<b->row;i++) {
    fprintf(fid,"%25.16e\n",b->val[i]);
  }
  return;
}

/***********************************************************************************************/
/**
 * \fn void dvector_write (const char *filename, dvector *vec)
 *
 * \brief Write a dvector to disk file
 *
 * \param vec       Pointer to the dvector
 * \param filename  File name
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 */
void dvector_write (const char *filename,
                    dvector *vec)
{
  INT m = vec->row, i;

  FILE *fp = fopen(filename,"w");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is writing to file %s...\n", __FUNCTION__, filename);

  fprintf(fp,"%lld\n",(long long )m);

  for ( i = 0; i < m; ++i ) fprintf(fp,"%0.15e\n",vec->val[i]);

  fclose(fp);
}


/***********************************************************************************************/
/**
 * \fn void dvector_read (const char *filename, dvector *b)
 *
 * \brief Read b from a disk file in array format
 *
 * \param filename  File name for vector b
 * \param b         Pointer to the dvector b (output)
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 *
 */
void dvector_read (const char *filename,
                   dvector *b)
{

  INT  i, n;
  long long i_in;
  REAL value;

  FILE *fp = fopen(filename,"r");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is reading file %s...\n", __FUNCTION__, filename);

  fscanf(fp,"%lld",&i_in); n=(INT )i_in;

  dvec_alloc(n,b);

  for ( i = 0; i < n; ++i ) {

    fscanf(fp, "%le", &value);
    b->val[i] = value;

    if ( value > BIGREAL ) {
      fprintf(stderr,"### ERROR: Wrong value = %lf\n", value);
      dvec_free(b);
      fclose(fp);
      exit(ERROR_INPUT_PAR);
    }

  }

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void ivector_print(FILE* fid,ivector *b)
 *
 * \brief print a ivector to a file
 *
 * \param fid  Pointer to the file
 * \param b    Pointer to the Integer dvector
 *
 */
void ivector_print(FILE* fid,
                   ivector *b)
{
  // local vairable
  INT i;

  // main loop
  for(i=0;i<b->row;i++) {
    fprintf(stdout,"i = %lld\n",(long long )i);
    fprintf(fid,"%lld\n",(long long )b->val[i]);
  }
  return;
}


/*******************************************************************/
/*  \fn void print_full_mat(const  INT n, const INT m, REAL *A,const char *varname)
 *
 *
 *  \note: prints a matrix A with (n) rows and (m) columns in matlab
 *         format e.g. 2 x 2 identity is printed as I2=[1. 0.;0. 1.];
 *         if th varname= "I2"
 *
 */
void print_full_mat(const  INT n, const INT m, REAL *A,const char *varname)
{
  INT nprt=16400,mprt=1025;
  if( (n<1) || (m<1) ) return;
  INT i,j,n1=n-1;
  if(n<=nprt) nprt=n;
  if(m<=mprt) mprt=m;
  if(varname==NULL){
    fprintf(stdout,"\nA=[");
  }else{
    fprintf(stdout,"\n%s=[",varname);
  }
  for (i = 0; i<nprt;i++){
    for(j=0;j<mprt;j++){
      fprintf(stdout,"%23.16e ", A[m*i+j]);
    }
    if(i!=n1){
      fprintf(stdout,";");
    }else{
      fprintf(stdout,"];\n");
    }
  }
  return;
}
/*******************************************************************/
/*  \fn void print_full_mat_l(const  INT n, const INT m, REAL16 *A,const char *varname)
 *
 *
 *  \note: prints a long double matrix A with (n) rows and (m) columns
 *         in matlab format e.g. 2 x 2 identity is printed as
 *         I2=[1. 0.;0. 1.]; if the varname= "I2"
 *
 */
void print_full_mat_l(const  INT n, const INT m, REAL16 *A,const char *varname)
{
  INT nprt=16400,mprt=1025;
  if( (n<1) || (m<1) ) return;
  INT i,j,n1=n-1;
  if(n<=nprt) nprt=n;
  if(m<=mprt) mprt=m;
  if(varname==NULL){
    fprintf(stdout,"\nA=[");
  }else{
    fprintf(stdout,"\n%s=[",varname);
  }
  for (i = 0; i<nprt;i++){
    for(j=0;j<mprt;j++){
      fprintf(stdout,"%.18Le ", A[m*i+j]);
    }
    if(i!=n1){
      fprintf(stdout,";");
    }else{
      fprintf(stdout,"];\n");
    }
  }
  return;
}
/*******************************************************************/
/*  \fn void print_full_mat_int(const  INT n, const INT m, INT *A,const char *varname)
 *
 *
 *  \note: prints an INTEGER matrix A with (n) rows and (m) columns in
 *         matlab format e.g. 2 x 2 integer matrix is printed as
 *         X=[1 2; 3 4]; if the varname= "X"
 *
 */
void print_full_mat_int(const  INT n, const INT m, INT *A,const char *varname)
{
  INT nprt=16400,mprt=1025;
  if( (n<1) || (m<1) ) return;
  INT i,j,n1=n-1;
  if(n<=nprt) nprt=n;
  if(m<=mprt) mprt=m;
  if(varname==NULL){
    fprintf(stdout,"\nintA=[");
  }else{
    fprintf(stdout,"\n%s=[",varname);
  }
  for (i = 0; i<nprt;i++){
    for(j=0;j<mprt;j++){
      fprintf(stdout,"%16lld ", (long long )A[m*i+j]);
    }
    if(i!=n1){
      fprintf(stdout,";");
    }else{
      fprintf(stdout,"];\n");
    }
  }
  return;
}
/***********************************************************************************************/
/*!
 * \fn void csr_print_matlab(FILE* fid,dCSRmat *A)
 *
 * \brief print a dCSRmat format sparse matrix to a file
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the dCSRmat format sparse matrix
 *
 */
void csr_print_matlab(FILE* fid,
                      dCSRmat *A)
{
  // local variables
  INT i,j1,j2,j;

  // main loop
  for(i=0;i<A->row;i++) {
    j1 = A->IA[i];
    j2 = A->IA[i+1];
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%lld %lld %25.16e\n",(long long )(i+1),(long long )A->JA[j]+1,A->val[j]);
    }
  }

  return;
}
/***********************************************************************************************/
/*!
 * \fn void bdcsr_print_matlab(FILE* fid,dCSRmat *A)
 *
 * \brief print a dCSRmat format sparse matrix to a file
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the dCSRmat format sparse matrix
 *
 */
void bdcsr_print_matlab(FILE* fid,
                        block_dCSRmat *A)
{
  // local variables
  dCSRmat Amerge = bdcsr_2_dcsr(A);
  //fprintf(stdout,"row = %d, col = %d, nnz = %d\n", Amerge.row, Amerge. col, Amerge.nnz);
  //dcsr_write_dcoo("Amerge.dat", &Amerge);
  csr_print_matlab(fid,&Amerge);
  dcsr_free(&Amerge);

  return;
}
void scomplex_print_matlab(const char *fname,scomplex *sc)
{
  INT ns_max=100000;
  FILE *fp=fopen(fname,"w");
  if(sc->ns>ns_max){
    fprintf(fp,"\n%%%% Too large:elements=%lld>%lld\n",(long long )sc->ns,(long long )ns_max);
    fclose(fp);
    return;
  }
  INT i, j,n1=sc->n+1,dim=sc->n,nv=sc->nv,ns=sc->ns;
  fprintf(fp,"\nt=[");
  for(j=0;j<n1;j++){
    for(i=0;i<ns-1;++i){
      fprintf(fp,"%lld,",(long long )sc->nodes[i*n1+j]);
    }
    fprintf(fp,"%lld;\n",(long long )sc->nodes[(ns-1)*n1+j]);
  }
  fprintf(fp,"];t=t';\n");
  fprintf(fp,"\nnbr=[");
  for(j=0;j<n1;j++){
    for(i=0;i<ns-1;i++){
      fprintf(fp,"%lld,",(long long )sc->nbr[i*n1+j]);
    }
    fprintf(fp,"%lld;\n",(long long )sc->nbr[(ns-1)*n1+j]);
  }
  fprintf(fp,"];nbr=nbr';\n");
  //
  fprintf(fp,"\nxp=[");
  for(j=0;j<dim;j++){
    for(i=0;i<nv-1;i++){
      fprintf(fp,"%.16e,",sc->x[i*dim+j]);
    }
    fprintf(fp,"%.16e;\n",sc->x[(nv-1)*dim+j]);
  }
  fprintf(fp,"];xp=xp';\n");
  fprintf(fp,"\nib=[");
  for(i=0;i<nv-1;i++){
    fprintf(fp,"%lld,",(long long )sc->bndry[i]);
  }
  fprintf(fp,"%lld];ib=ib';\n",(long long )sc->bndry[nv-1]);
  fclose(fp);
  return;
}
/***********************************************************************************************/
/*!
 * \fn void csr_print_native(FILE* fid,dCSRmat *A)
 *
 * \brief print a dCSRmat format sparse matrix to a file in the
 * follwoing way: number of rows, number of columns, number of
 * nonzeroes; next line: ia[k] k=0:(nrows+1); last line ja[k],k=1:nnz
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the dCSRmat format sparse matrix
 *
 */
void csr_print_native(FILE* fid,
                      dCSRmat *A, dvector *rhs)
{
  // local variables
  INT i,shift;

  if(A->IA[0]==0) shift=0; else  shift=-1;
  // shift -1 if it was fortran, starting from 1.
  fprintf(fid,"%lld %lld %lld\n",(long long )A->row, (long long )A->col, (long long )A->nnz);
  for(i=0;i<(A->row+1);i++)
    fprintf(fid,"%lld ",((long long )A->IA[i]+(long long )shift));
  fprintf(fid,"\n");
  for(i=0;i<A->nnz;i++)
    fprintf(fid,"%lld ",((long long )A->JA[i]+(long long )shift));
  fprintf(fid,"\n");
  for(i=0;i<A->nnz;i++)
    fprintf(fid,"%23.16e ",A->val[i]);
  fprintf(fid,"\n");
  if(rhs)
    for(i=0;i<rhs->row;i++)
      fprintf(fid,"%23.16e ",rhs->val[i]);
  return;
}

/***********************************************************************************************/
/*!
 * \fn void icsr_print_matlab(FILE* fid,dCSRmat *A)
 *
 * \brief print a iCSRmat format sparse matrix to a file
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the iCSRmat format sparse matrix
 *
 * \todo add shift flag -- Xiaozhe Hu
 *
 */
void icsr_print_matlab(FILE* fid,
                       iCSRmat *A)
{
  // local variables
  INT i,j1,j2,j;

  // main loop; comma separated
  for(i=0;i<A->row;i++) {
    j1 = A->IA[i];
    j2 = A->IA[i+1];
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%lld,%lld\n",(long long )(i+1),(long long )(A->JA[j]+1));
    }
  }
  return;
}
/***********************************************************************************************/
/*!
 * \fn void icsr_print_rows(FILE* fid,iCSRmat *A,const char *sname)
 *
 * \brief print a iCSRmat format sparse matrix to a file with VALUES
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the iCSRmat format sparse matrix
 * \param sname just a string so the printout starts with "sname*"
 *
 * \todo
 *
 */
void icsr_print_rows(FILE* fid,
		     iCSRmat *A,const char *sname)
{
  // local variables
  INT i,j;
  fprintf(fid,"\n%%%% %s:",sname);
  for(i=0;i<A->row;++i){
    if((A->IA[i+1]-A->IA[i])<=0) continue;
      fprintf(fid,"\nrow[%lld]=[ ",(long long )i);
      for(j=A->IA[i];j<A->IA[i+1];++j){
	fprintf(fid,"%lld ",(long long )A->JA[j]);
      }
      fprintf(fid,"]");
  }
  fflush(fid);
  return;
}
/***********************************************************************************************/
/*!
 * \fn void icsr_print_matlab_val(FILE* fid,iCSRmat *A)
 *
 * \brief print a iCSRmat format sparse matrix to a file with VALUES
 *
 * \param fid  Pointer to the file
 * \param A    Pointer to the iCSRmat format sparse matrix
 *
 * \todo
 *
 */
void icsr_print_matlab_val(FILE* fid,
			   iCSRmat *A)
{
  // local variables
  INT i,j1,j2,j;

  // main loop; comma separated
  for(i=0;i<A->row;i++) {
    j1 = A->IA[i];
    j2 = A->IA[i+1];
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%lld,%lld,%lld\n",(long long )(i+1),(long long )(A->JA[j]+1),(long long )(A->val[j]));
    }
  }
  return;
}

/***********************************************************************************************/
/*!
 * \fn void dvec_write (const char *filename, dvector *vec)
 *
 * \brief Write a dvector to disk file
 *
 * \param vec       Pointer to the dvector
 * \param filename  File name
 *
 */
void dvec_write (const char *filename,
                 dvector *vec)
{
  // local variables
  INT m = vec->row, i;

  FILE *fp = fopen(filename,"w");

  if ( fp == NULL )
    check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is writing to file %s...\n", __FUNCTION__, filename);

  fprintf(fp,"%lld\n",(long long )m);

  //main loop
  for ( i = 0; i < m; ++i ) fprintf(fp,"%0.15e\n",vec->val[i]);

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void ivec_write (const char *filename, ivector *vec)
 *
 * \brief Write an ivector to disk file
 *
 * \param vec       Pointer to the ivector
 * \param filename  File name
 *
 */
void ivec_write (const char *filename,
                 ivector *vec)
{
  // local variables
  INT m = vec->row, i;

  FILE *fp = fopen(filename,"w");

  if ( fp == NULL )
    check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: ivector written to file: %s...\n", __FUNCTION__, filename);

  fprintf(fp,"%lld\n",(long long )m);

  //main loop
  for ( i = 0; i < m; ++i ) fprintf(fp,"%10lld\n",(long long )vec->val[i]);

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void ddense_write(const char *filename, dDENSEmat *A)
 *
 * \brief Write a dDENSEmat matrix to disk file in row-wise format
 *
 * \param A         pointer to the dDENSEmat matrix
 * \param filename  char for file name
 *
 */
void ddense_write(const char *filename,
                  dDENSEmat *A)
{
  // local variables
  const INT n = A->row, m = A->col;
  const INT nnz = n*m;
  INT i;

  FILE *fp = fopen(filename, "w");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is writing to file %s...\n", __FUNCTION__, filename);

  // main loop
  fprintf(fp,"%lld  %lld\n",(long long )n,(long long )m);
  for (i = 0; i < nnz; ++i) {
    fprintf(fp,"%0.15e\n", A->val[i]);
  }

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_write_dcoo (const char *filename, dCSRmat *A)
 *
 * \brief Write a dCSRmat matrix to disk file in IJ format (coordinate format)
 *
 * \param A         pointer to the dCSRmat matrix
 * \param filename  char for vector file name
 *
 */
void dcsr_write_dcoo (const char *filename,
                      dCSRmat *A)
{
  // local variables
  const INT m = A->row, n = A->col;
  INT i, j;

  FILE *fp = fopen(filename, "w");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is writing to file %s...\n", __FUNCTION__, filename);

  // main loop
  fprintf(fp,"%lld  %lld  %lld\n",(long long )m,(long long )n,(long long )A->nnz);
  for ( i = 0; i < m; ++i ) {
    for ( j = A->IA[i]; j < A->IA[i+1]; j++ )
      fprintf(fp,"%lld  %lld  %0.15e\n",(long long )i,(long long )A->JA[j],A->val[j]);
  }

  fclose(fp);
}

/***********************************************************************************************/
/*!
 * \fn void icsr_write_icoo (const char *filename, dCSRmat *A)
 *
 * \brief Write an iCSRmat matrix to disk file in IJ format (coordinate format)
 *
 * \param A         pointer to the iCSRmat matrix
 * \param filename  char for vector file name
 *
 */
void icsr_write_icoo (const char *filename,
                      iCSRmat *A)
{
  // local variables
  const INT m = A->row, n = A->col;
  INT i, j;

  FILE *fp = fopen(filename, "w");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: iCSRmat is written to file: %s...\n", __FUNCTION__, filename);
  // main loop
  fprintf(fp,"%lld  %lld  %lld\n",(long long )m,(long long )n,(long long )A->nnz);
  if(A->val!=NULL){
    for ( i = 0; i < m; ++i ) {
      for ( j = A->IA[i]; j < A->IA[i+1]; j++ )
	fprintf(fp,"%lld  %lld  %lld\n",(long long )i,(long long )A->JA[j],(long long )A->val[j]);
    }
  }else{
    for ( i = 0; i < m; ++i ) {
      for ( j = A->IA[i]; j < A->IA[i+1]; j++ )
	fprintf(fp,"%lld  %lld  %lld\n",(long long )i,(long long )A->JA[j],(long long )1);
    }
  }

  fclose(fp);
  return;
}

/***********************************************************************************************/
/*!
 * \fn void bdcsr_write_dcoo(FILE* fid,dCSRmat *A)
 *
 * \brief print a dCSRmat format sparse matrix to a file
 *
 * \param filename  char for vector file name
 * \param A    Pointer to the dCSRmat format sparse matrix
 *
 */
void bdcsr_write_dcoo(const char *filename,
                      block_dCSRmat *A)
{
  // local variables
  dCSRmat Amerge = bdcsr_2_dcsr(A);

  // write
  dcsr_write_dcoo(filename, &Amerge);

  // free
  dcsr_free(&Amerge);

  return;
}


/***********************************************************************************************/
/**
 * \fn void dcoo_read_dcsr(const char *filename, dCSRmat *A)
 *
 * \brief Read A from matrix disk file in IJ format and convert to CSR format
 *
 * \param filename  File name for matrix
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 */
void dcoo_read_dcsr (const char *filename,
                     dCSRmat *A)
{
  INT  i,j,k,m,n,nnz;
  REAL value;

  FILE *fp = fopen(filename,"r");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is reading file %s...\n", __FUNCTION__, filename);

  fscanf(fp,"%lld %lld %lld",(long long *)&m,(long long *)&n,(long long *)&nnz);

  dCOOmat Atmp=dcoo_create(m,n,nnz);

  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%lld %lld %le", (long long *)&i, (long long *)&j,&value) != EOF ) {
      Atmp.rowind[k]=i; Atmp.colind[k]=j; Atmp.val[k] = value;
    }
    else {
      check_error(ERROR_WRONG_FILE, "dcoo_read_dcsr");
    }
  }

  fprintf(stdout,"%%%%%s: HAZMATH reading file %s is DONE. \n", __FUNCTION__, filename); fflush(stdout);
  fclose(fp);

  dcoo_2_dcsr(&Atmp,A);
  dcoo_free(&Atmp);
}

/**
 * \fn void dbsr_read (const char *filename, dBSRmat *A)
 *
 * \brief Read A from a disk file in dBSRmat format
 *
 * \param filename   File name for matrix A
 * \param A          Pointer to the dBSRmat A
 *
 * \note
 *   This routine reads a dBSRmat matrix from a disk file in the following format:
 *
 * \note File format:
 *   - ROW, COL, NNZ
 *   - nb: size of each block
 *   - storage_manner: storage manner of each block
 *   - ROW+1: length of IA
 *   - IA(i), i=0:ROW
 *   - NNZ: length of JA
 *   - JA(i), i=0:NNZ-1
 *   - NNZ*nb*nb: length of val
 *   - val(i), i=0:NNZ*nb*nb-1
 *
 * \author Xiaozhe Hu
 * \date   10/29/2010
 */
void dbsr_read (const char  *filename,
                dBSRmat     *A)
{
    INT     ROW, COL, NNZ, nb, storage_manner;
    INT     i, n;
    INT     index;
    REAL    value;
    size_t  status;

    FILE *fp = fopen(filename,"r");

    if ( fp == NULL ) check_error(ERROR_OPEN_FILE, filename);

    printf("%%%%%s: HAZMATH is reading file %s...\n", __FUNCTION__, filename);

    status = fscanf(fp, "%lld %lld %lld", (long long *)&ROW,(long long *)&COL,(long long *)&NNZ); // read dimension of the problem
    check_error(status, filename);
    A->ROW = ROW; A->COL = COL; A->NNZ = NNZ;

    status = fscanf(fp, "%lld", (long long *)&nb); // read the size of each block
    check_error(status, filename);
    A->nb = nb;

    status = fscanf(fp, "%lld", (long long *)&storage_manner); // read the storage_manner
    check_error(status, filename);
    A->storage_manner = storage_manner;

    // allocate memory space
    dbsr_alloc(ROW, COL, NNZ, nb, storage_manner, A);

    // read IA
    status = fscanf(fp, "%lld", (long long *)&n);
    check_error(status, filename);
    for ( i = 0; i < n; ++i ) {
      status = fscanf(fp, "%lld", (long long *)&index);
        check_error(status, filename);
        A->IA[i] = index;
    }

    // read JA
    status = fscanf(fp, "%lld", (long long *)&n);
    check_error(status, filename);
    for ( i = 0; i < n; ++i ){
      status = fscanf(fp, "%lld", (long long *)&index);
        check_error(status, filename);
        A->JA[i] = index;
    }

    // read val
    status = fscanf(fp, "%lld", (long long *)&n);
    check_error(status, filename);
    for ( i = 0; i < n; ++i ) {
        status = fscanf(fp, "%le", &value);
        check_error(status, filename);
        A->val[i] = value;
    }

    fclose(fp);
}


/*** Auxillary Files *********************************************************************/
/*
 * \fn void read_eof(const char *fname, void **data, INT *ndata, const char dtype,const SHORT print_level)
 *
 * \brief Reads data until EOF until eof. skips any control chars at
 * the beginning of the file (ASCII < 0).
 *
 * The result is stored in data[0] which is allocated here (REAL* for
 *  double/float or INT* for ints).
 *
 * dtype is either 'R', or 'r' for REAL or 'I', 'i' for INT.
 *
*/
void read_eof(const char *fname, void **data, INT *ndata, const char dtype,const SHORT print_level)
/*************************************************************************************/
{
  SHORT data_type=-1;
  char *type_name=NULL;
  if(dtype=='i' || dtype=='I'){
    data_type=TRUE;// integers
    type_name=strdup("INTs");
  } else if(dtype=='r' || dtype=='R'){
    data_type=FALSE;
    type_name=strdup("REALs");
  } else{
    fprintf(stderr,"*** Wrong type in %s; must be \'I\' for integer type  or \'R\' for double/float",__FUNCTION__);
    exit(9);
  }
  INT k,count=0,ichk=-1;
  long long lint0;
  long double dl;
  char ch;
  FILE *fp=fopen(fname,"r");
  if (!fp) {
    perror("fopen");
    return;
  }
  if(print_level>10)
    fprintf(stdout,"\n%%%%Reading %s from %s\n",type_name,fname);
  ch=fgetc(fp); while((INT )ch < 0){ch=fgetc(fp);count++;}
  if(count){
    if(print_level>10)
      fprintf(stdout,"%%Read: %lld control chars...", (long long int)count);
    fseek(fp,count*sizeof(char),SEEK_SET);
  } else rewind(fp);
  if(data_type){
    k=0;
    while(1){
      if(feof(fp)) break;
      ichk=fscanf(fp,"%lld", &lint0);
      if(ichk<0) break;
      k++;
    }
  } else {
    k=0;
    while(1){
      if(feof(fp)) break;
      ichk=fscanf(fp,"%Lg", &dl);
      if(ichk<0) break;
      k++;
    }
  }
  if(count){
    fseek(fp,count*sizeof(char),SEEK_SET);
  } else {
    rewind(fp);
  }
  /* we allocate k reals array*/
  ndata[0]=k;
  if(data_type){
    data[0]=(void *)calloc(ndata[0],sizeof(INT));
    rveci_(fp,(INT *)data[0],ndata);
  }else{
    data[0]=(void *)calloc(ndata[0],sizeof(REAL));
    rvecd_(fp,(REAL *)data[0],ndata);
  }
  if(print_level>10)
    fprintf(stdout,"\n%%%%Read %lld %s\n",(long long int)ndata[0],type_name);
  return;
}
/****************************************************************************************/

/*!
 * \fn void rveci_(FILE *fp, INT *vec, INT *nn)
 *
 * \brief Reads a vector of integers of size nn from a file fp
 *
 * \param fp        FILE ID
 * \param vec       Where to store vector
 * \param nn        Size of Vector
 *
 */
void rveci_(FILE *fp,INT *vec,INT *nn)
{
  // local variables
  INT n;
  INT *vec_end;
  n = *nn;
  vec_end  =  vec + n;
  long long readint;
  // main loop
  for ( ; vec < vec_end; ++vec) {
  //  fscanf(fp,"%lld",(long long *)vec);
    fscanf(fp,"%lld",&readint);
    *vec = readint;
  }
  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
 * \fn void rvecd_(FILE *fp, REAL *vec, INT *nn)
 *
 * \brief Reads a vector of doubles of size nn from a file fp
 *
 * \param fp        FILE ID
 * \param vec       Where to store vector
 * \param nn        Size of Vector
 *
 */
void rvecd_(FILE *fp,
            REAL *vec,
            INT *nn)
{
  // local variables
  INT n;
  REAL *vec_end;
  n= *nn;
  vec_end =  vec + n;

  // main loop
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%lg",vec);
  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
 * \fn FILE* HAZ_fopen( char *fname, char *mode )
 *
 * \brief A graceful version of fopen(). It checks if the file has
 *     been successfully opened.  If  that is  not  the case  a
 *     message is printed and the program is exited.
 *
 * \param fname     Filename
 * \param mode      read or write
 *
 */
FILE* HAZ_fopen(char *fname, char *mode )
{
  // local variable
  FILE   *fp;

  fp = fopen(fname,mode);
  if ( fp == NULL ) {
    fprintf(stderr,"Cannot open %s  -Exiting\n",fname);
    exit(255);
  }
  return fp;
}
/****************************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_sol_onV_vtk(char *namevtk,mesh_struct *mesh,REAL *sol,INT ncomp)
 *
 * \brief Dumps solution data to vtk format only on vertices
 *
 * \param namevtk  Filename
 * \param mesh     Mesh struct to dump
 * \param sol      solution vector to dump
 * \param ncomp:   Number of components to the solution
 *
 */
void dump_sol_onV_vtk(char *namevtk,
                      mesh_struct *mesh,
                      REAL *sol,
                      INT ncomp)
{
  // Basic Quantities
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;

  // VTK needed Quantities
  INT tcell=-10;
  INT k=-10,j=-10,kndl=-10;
  char *tfloat=strndup("Float64",8);
  char *tinto=strndup("Int64",6);
  char *endian=strndup("LittleEndian",13);
  /*
    What endian?:

    Intel x86; OS=MAC OS X: little-endian
    Intel x86; OS=Windows: little-endian
    Intel x86; OS=Linux: little-endian
    Intel x86; OS=Solaris: little-endian
    Dec Alpha; OS=Digital Unix: little-endian
    Dec Alpha; OS=VMS: little-endian
    Hewlett Packard PA-RISC; OS=HP-UX: big-endian
    IBM RS/6000; OS=AIX: big-endian
    Motorola PowerPC; OS=Mac OS X:  big-endian
    SGI R4000 and up; OS=IRIX: big-endian
    Sun SPARC; OS=Solaris: big-endian
  */

  /*
    Types of cells for VTK

    VTK_VERTEX (=1)
    VTK_POLY_VERTEX (=2)
    VTK_LINE (=3)
    VTK_POLY_LINE (=4)
    VTK_TRIANGLE(=5)
    VTK_TRIANGLE_STRIP (=6)
    VTK_POLYGON (=7)
    VTK_PIXEL (=8)
    VTK_QUAD (=9)
    VTK_TETRA (=10)
    VTK_VOXEL (=11)
    VTK_HEXAHEDRON (=12)
    VTK_WEDGE (=13)
    VTK_PYRAMID (=14)
  */

  const INT LINE=3;
  const INT TRI=5;
  const INT TET=10;

  if(dim==1) {
    tcell=LINE; /* line */
  } else if(dim==2) {
    tcell=TRI; /* triangle */
  } else {
    tcell=TET; /* tet */
  }
  // Open File for Writing
  FILE* fvtk = HAZ_fopen(namevtk,"w");

  // Write Headers
  fprintf(fvtk, \
	  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n", \
	  endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(long long )nv,(long long )nelm);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">", \
  	  tfloat);

  // Dump coordinates
  if(dim == 1) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],0e0,0e0);
    }
  } else if(dim == 2) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k],0e0);
    }
  } else {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k], \
  	      mesh->cv->z[k]);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");

  // Dump solution Data on Vertices of mesh
  fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
  INT i=0;
  for(i=0;i<ncomp;i++) {
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %lld\" Format=\"ascii\">",tfloat,(long long )i);
    for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[i*nv+k]);
    fprintf(fvtk,"</DataArray>\n");
  }
  fprintf(fvtk,"</PointData>\n");

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %lld ",(long long )mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %lld ",(long long )mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %lld ",(long long )tcell);

  // Put in remaining headers
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Cells>\n");
  fprintf(fvtk,"</Piece>\n");
  fprintf(fvtk,"</UnstructuredGrid>\n");
  fprintf(fvtk,"</VTKFile>\n");

  fclose(fvtk);
  free(tfloat);
  free(tinto);
  free(endian);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_sol_vtk(char *namevtk,char *varname,mesh_struct *mesh,fespace *FE,REAL *sol)
 *
 * \brief Dumps solution data to vtk format.  Tries to do best interpolation for given FE space.
 *
 * \param namevtk  Filename
 * \param varname  String for variable name
 * \param mesh     Mesh struct to dump
 * \param FE       FE space of solution
 * \param sol      solution vector to dump
 *
 */
void dump_sol_vtk(char *namevtk,char *varname,mesh_struct *mesh,fespace *FE,REAL *sol)
{
  // Basic Quantities
  INT i;
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;

  // VTK needed Quantities
  INT tcell=-10;
  INT k=-10,j=-10,kndl=-10;
  char *tfloat="Float64", *tinto="Int64", *endian="LittleEndian";

  /*
    What endian?:

    Intel x86; OS=MAC OS X: little-endian
    Intel x86; OS=Windows: little-endian
    Intel x86; OS=Linux: little-endian
    Intel x86; OS=Solaris: little-endian
    Dec Alpha; OS=Digital Unix: little-endian
    Dec Alpha; OS=VMS: little-endian
    Hewlett Packard PA-RISC; OS=HP-UX: big-endian
    IBM RS/6000; OS=AIX: big-endian
    Motorola PowerPC; OS=Mac OS X:  big-endian
    SGI R4000 and up; OS=IRIX: big-endian
    Sun SPARC; OS=Solaris: big-endian
  */

  /*
    Types of cells for VTK

    VTK_VERTEX (=1)
    VTK_POLY_VERTEX (=2)
    VTK_LINE (=3)
    VTK_POLY_LINE (=4)
    VTK_TRIANGLE(=5)
    VTK_TRIANGLE_STRIP (=6)
    VTK_POLYGON (=7)
    VTK_PIXEL (=8)
    VTK_QUAD (=9)
    VTK_TETRA (=10)
    VTK_VOXEL (=11)
    VTK_HEXAHEDRON (=12)
    VTK_WEDGE (=13)
    VTK_PYRAMID (=14)
  */

  const INT LINE=3;
  const INT TRI=5;
  const INT TET=10;

  if(dim==1) {
    tcell=LINE; /* line */
  } else if(dim==2) {
    tcell=TRI; /* triangle */
  } else {
    tcell=TET; /* tet */
  }
  // Open File for Writing
  FILE* fvtk = HAZ_fopen(namevtk,"w");

  // Write Headers
  fprintf(fvtk,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(long long )nv,(long long )nelm);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">",tfloat);

  // Dump vertex coordinates and solution
  if(dim == 1) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],0e0,0e0);
    }
  } else if(dim == 2) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k],0e0);
    }
  } else {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k], \
              mesh->cv->z[k]);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");

  // Dump Solution
  // Depending on the FE space, we will dump things differently
  // since we need to project to the vertices
  REAL* sol_on_V=NULL;
  if(FE->FEtype==0) { // P0 - only have cell data
    fprintf(fvtk,"<CellData Scalars=\"scalars\">\n");
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component - %s\" Format=\"ascii\">",tfloat,varname);
    for(k=0;k<nelm;k++) fprintf(fvtk," %23.16e ",sol[k]);
    fprintf(fvtk,"</DataArray>\n");
    fprintf(fvtk,"</CellData>\n");
  } else if(FE->FEtype>0 && FE->FEtype<20) { // PX elements (assume sol at vertices comes first)
    fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component - %s\" Format=\"ascii\">",tfloat,varname);
    for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[k]);
    fprintf(fvtk,"</DataArray>\n");
    fprintf(fvtk,"</PointData>\n");
  } else { // Vector Elements
    sol_on_V = (REAL *) calloc(dim*mesh->nv,sizeof(REAL));
    Project_to_Vertices(sol_on_V,sol,FE,mesh);
    fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
    for(i=0;i<dim;i++) {
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component - %s%lld\" Format=\"ascii\">",tfloat,varname,(long long )i);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol_on_V[i*nv+k]);
      fprintf(fvtk,"</DataArray>\n");
    }
    fprintf(fvtk,"</PointData>\n");
  }

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %lld ",(long long )mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %lld ",(long long )mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %lld ",(long long )tcell);

  // Put in remaining headers
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Cells>\n");
  fprintf(fvtk,"</Piece>\n");
  fprintf(fvtk,"</UnstructuredGrid>\n");
  fprintf(fvtk,"</VTKFile>\n");

  fclose(fvtk);
  if(sol_on_V) free(sol_on_V);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_blocksol_vtk(char *namevtk,char *varname,mesh_struct *mesh,block_fespace *FE,REAL *sol)
 *
 * \brief Dumps solution data to vtk format.  Tries to do best interpolation for given FE space.
 *
 * \param namevtk  Filename
 * \param varname  String for variable names
 * \param mesh     Mesh struct to dump
 * \param FE       Block FE space of solution
 * \param sol      solution vector to dump
 *
 */
void dump_blocksol_vtk(char *namevtk,char **varname,mesh_struct *mesh,block_fespace *FE,REAL *sol)
{
  // Basic Quantities
  INT i,nsp;
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;

  // VTK needed Quantities
  INT tcell=-10;
  INT k=-10,j=-10,kndl=-10;
  char *tfloat="Float64", *tinto="Int64", *endian="LittleEndian";

  /*
    What endian?:

    Intel x86; OS=MAC OS X: little-endian
    Intel x86; OS=Windows: little-endian
    Intel x86; OS=Linux: little-endian
    Intel x86; OS=Solaris: little-endian
    Dec Alpha; OS=Digital Unix: little-endian
    Dec Alpha; OS=VMS: little-endian
    Hewlett Packard PA-RISC; OS=HP-UX: big-endian
    IBM RS/6000; OS=AIX: big-endian
    Motorola PowerPC; OS=Mac OS X:  big-endian
    SGI R4000 and up; OS=IRIX: big-endian
    Sun SPARC; OS=Solaris: big-endian
  */

  /*
    Types of cells for VTK

    VTK_VERTEX (=1)
    VTK_POLY_VERTEX (=2)
    VTK_LINE (=3)
    VTK_POLY_LINE (=4)
    VTK_TRIANGLE(=5)
    VTK_TRIANGLE_STRIP (=6)
    VTK_POLYGON (=7)
    VTK_PIXEL (=8)
    VTK_QUAD (=9)
    VTK_TETRA (=10)
    VTK_VOXEL (=11)
    VTK_HEXAHEDRON (=12)
    VTK_WEDGE (=13)
    VTK_PYRAMID (=14)
  */

  const INT LINE=3;
  const INT TRI=5;
  const INT TET=10;

  if(dim==1) {
    tcell=LINE; /* line */
  } else if(dim==2) {
    tcell=TRI; /* triangle */
  } else {
    tcell=TET; /* tet */
  }
  // Open File for Writing
  FILE* fvtk = HAZ_fopen(namevtk,"w");

  // Write Headers
  fprintf(fvtk,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(long long )nv,(long long )nelm);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">",tfloat);

  // Dump vertex coordinates and solution
  if(dim == 1) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],0e0,0e0);
    }
  } else if(dim == 2) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k],0e0);
    }
  } else {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k], \
              mesh->cv->z[k]);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");

  // Dump Solution for each FE space
  REAL* sol_on_V=NULL;
  REAL* solptr=NULL;
  INT* P0cntr = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT anyP0=0;
  INT spcntr = 0;
  fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
  for(nsp=0;nsp<FE->nspaces;nsp++) {
    // Depending on the FE space, we will dump things differently
    // since we need to project to the vertices
    if(FE->var_spaces[nsp]->FEtype==0) { // P0 - only have cell data
      // We need to save this for later so mark
      anyP0=1;
      P0cntr[nsp] = 1;
    } else if((FE->var_spaces[nsp]->FEtype>0 && FE->var_spaces[nsp]->FEtype<20) || FE->var_spaces[nsp]->FEtype==103) { // PX elements (assume sol at vertices comes first)
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %lld - %s\" Format=\"ascii\">",tfloat,(long long )nsp,varname[nsp]);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[spcntr + k]);
      fprintf(fvtk,"</DataArray>\n");
    } else if(FE->var_spaces[nsp]->FEtype==99) { // Single DoF constraint element (just plot single value everywhere)
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %lld - %s\" Format=\"ascii\">",tfloat,(long long )nsp,varname[nsp]);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[spcntr]);
      fprintf(fvtk,"</DataArray>\n");
    } else { // Vector Elements
      sol_on_V = (REAL *) calloc(dim*mesh->nv,sizeof(REAL));
      solptr = sol+spcntr;
      Project_to_Vertices(sol_on_V,solptr,FE->var_spaces[nsp],mesh);
      for(i=0;i<dim;i++) {
        fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %lld - %s%lld\" Format=\"ascii\">",tfloat,(long long )nsp,varname[nsp],(long long )i);
        for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol_on_V[i*nv+k]);
        fprintf(fvtk,"</DataArray>\n");
      }
      if(sol_on_V) free(sol_on_V);
    }
    spcntr += FE->var_spaces[nsp]->ndof;
  }
  fprintf(fvtk,"</PointData>\n");

  // Now go back and dump P0 stuff
  if(anyP0) {
    spcntr=0;
    fprintf(fvtk,"<CellData Scalars=\"scalars\">\n");
    for(nsp=0;nsp<FE->nspaces;nsp++) {
      if(P0cntr[nsp]==1) {
        fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %lld - %s\" Format=\"ascii\">",tfloat,(long long )nsp,varname[nsp]);
        for(k=0;k<nelm;k++) fprintf(fvtk," %23.16e ",sol[spcntr+k]);
        fprintf(fvtk,"</DataArray>\n");
      }
      spcntr += FE->var_spaces[nsp]->ndof;
    }
    fprintf(fvtk,"</CellData>\n");
  }
  if(P0cntr) free(P0cntr);

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %lld ",(long long )mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %lld ",(long long )mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %lld ",(long long )tcell);

  // Put in remaining headers
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Cells>\n");
  fprintf(fvtk,"</Piece>\n");
  fprintf(fvtk,"</UnstructuredGrid>\n");
  fprintf(fvtk,"</VTKFile>\n");

  fclose(fvtk);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void create_pvd(char *namepvd, INT nfiles,char *vtkfilename,char *filetype)
 *
 * \brief Dumps solution data in vtk format to a single file.  Useful for timestepping
 * \note  File names of vtk file must have same structure
 *
 * \param namepvd      Filename
 * \param nfiles       Number of files to store (i.e. timesteps)
 * \param vtkfilename  Filename structure of vtu files.
 * \param filetype     Name for types of files (i.e. "timestep")
 *
 */
void create_pvd(char *namepvd, INT nfiles,char *vtkfilename,char *filetype)
{
  // VTK needed Quantities
  //  What endian?:
  //    Intel x86; OS=MAC OS X: little-endian
  //    Intel x86; OS=Windows: little-endian
  //    Intel x86; OS=Linux: little-endian
  //    Intel x86; OS=Solaris: little-endian
  //    Dec Alpha; OS=Digital Unix: little-endian
  //    Dec Alpha; OS=VMS: little-endian
  //    Hewlett Packard PA-RISC; OS=HP-UX: big-endian
  //    IBM RS/6000; OS=AIX: big-endian
  //    Motorola PowerPC; OS=Mac OS X:  big-endian
  //    SGI R4000 and up; OS=IRIX: big-endian
  //    Sun SPARC; OS=Solaris: big-endian

  char *endian="LittleEndian";
  INT i;

  // Open File for Writing
  FILE* fvtk = HAZ_fopen(namepvd,"w");

  // Write Headers
  fprintf(fvtk,"<?xml version=\"1.0\"?>\n");
  fprintf(fvtk,
          "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"%s\" compressor=\"vtkZLibDataCompressor\">\n", \
          endian);
  fprintf(fvtk,"<Collection>\n");
  char filecounter[40];
  for(i=0;i<nfiles;i++) {
    sprintf(filecounter,"%s%03lld.vtu",vtkfilename,(long long )i);
    fprintf(fvtk,"<DataSet %s=\"%lld\" group=\"\" part=\"0\" file=\"%s\"/>\n",filetype,(long long )i,filecounter);
  }

  // Put in remaining headers
  fprintf(fvtk,"</Collection>\n");
  fprintf(fvtk,"</VTKFile>\n");

  fclose(fvtk);

  return;
}
/******************************************************************************/

/***********************************************************************************************/
/*!
 * \fn void debug_print(char *string,INT kill)
 *
 * \brief print debug message
 *
 * \param string  Pointer to the INT array
 * \param kill    if 1 kill the program
 *
 */
void debug_print(char* string, INT kill)
{
  fprintf(stdout,"%%%%%s\n",string);fflush(stdout);

  if(kill) exit(0);

  return;
}

/****************************************************************/
scomplex *hazr(char *namein)
{
  // READING a mesh in hazmath format.
  FILE *fmeshin;
  INT nv,ns,dim,nholes=0;
  INT k=-10,j=-10,m=-10;
  fmeshin=HAZ_fopen(namein,"r");
  /* *******************************************
     read a hazmath mesh file.
     *******************************************    */
  fscanf(fmeshin,"%lld %lld %lld %lld\n",(long long *)&ns,(long long *)&nv,(long long *)&dim,(long long *)&nholes);
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv,dim); // cannot read different dimensions.
  INT dim1=sc->n+1;
  for (j=0;j<dim1;j++) {
    for (k=0;k<ns;k++){
      m=dim1*k+j;
      /* shift if needed */
      fscanf(fmeshin,"%lld", (long long *)(sc->nodes+m));
      //      fprintf(stdout,"\n%lld",sc->nodes[m]);
    }
  }
  for (k=0;k<sc->ns;k++){
    fscanf(fmeshin,"%lld", (long long *)(sc->flags+k));
    //      fprintf(stdout,"\n%d",sc->flags[k]);
  }
  for(j=0;j<dim;j++){
    for(k=0;k<nv;k++){
      fscanf(fmeshin,"%lg",sc->x+k*dim+j);
    }
  }
  for(k=0;k<nv;k++){
    fscanf(fmeshin,"%lld", (long long *)(sc->bndry+k));
  }
  /* //NOT USED: a function value */
  /* if(sc->fval){ */
  /*   for(k=0;k<n;k++){ */
  /*     fscanf(fmeshin," %23.16g ", sc->fval+k); */
  /*   } */
  /* } */
  fclose(fmeshin);
  //  haz_scomplex_print(sc,0,"HERE");fflush(stdout);
  fprintf(stdout,"\n%%INPUT mesh: (hazmath) read from:%s\n",namein);
  return sc;
}
/********************************************************************************/
/*
 * Routines to save to file the mesh in different formats. uses the
 * simplicial complex data structure (scomplex *sc)
*/
/********************************************************************************/
/*!
 * \fn void hazw(char *nameout,scomplex *sc, const INT shift)
 *
 * \brief Write a simplicial complex to a file in a "hazmath" format.
 *
 * \param nameout   File name
 * \param sc        Pointer to a simplicial complex
 * \param shift     integer added to the elements of arrays such as sc->nodes.
 *
 * \note The data is organized as follows:
 *
 * 0. num_simplices,num_vertices,dimension,connected_components(bndry)-1;
 *
 * 1. As every simplex has (dim+1) vertices, the next
 *    (dim+1)*num_simplices integers are: 1st vertex for each simplex
 *    (num_simplices integers); 2nd vertex for every simplex, etc.
 *
 * 2. num_simplices integers are the flags (tags) associated with
 *    every element: could be something characterizing the element.
 *
 * 3. (dim)*num_vertices REALs (coordinates of the vertices) e.g.:
 *    num_vertices REALs with 1st coordinate, num_vertices REALs with
 *    2nd coordinate,...
 *
 * 4. num_vertices integers with tages (boundary codes) for every vertex.
 *
 */
/********************************************************************************/
void hazw(char *nameout,scomplex *sc, const INT shift)
{
  // WRITING in HAZMATH format.
  FILE *fmesh;
  INT n=sc->nv,ns=sc->ns, dim=sc->n,ndl=sc->n+1;
  INT *je = sc->nodes, *ib=sc->bndry;
  REAL *x = sc->x;
  INT k=-10,j=-10,kndl=-10;
  fmesh=HAZ_fopen(nameout,"w");
  /* *******************************************
     HAZMAT way of writing mesh file. sc->bndry_cc is the number of
     connected components on the boundary. sc->cc is the number of
     connected components domains.
     * TODO add sc->cc to the reading.
     */
  fprintf(fmesh,"%lld %lld %lld %lld\n",(long long )ns,(long long )n,(long long )dim,(long long )(sc->bndry_cc-1)); /* this is the
							       number of
							       holes;*/
  /* fprintf(stdout,"%lld %lld %lld\n",(long long )n,(long long )ns,(long long )sizeof(ib)/sizeof(INT)); */
  for (j=0;j<ndl;j++) {
    for (k=0;k<ns;k++){
      kndl=ndl*k+j;
      /* shift if needed */
      fprintf(fmesh," %lld ",(long long )(je[kndl]+shift));
    }
    fprintf(fmesh,"\n");
  }
  for (k=0;k<ns;k++){
    fprintf(fmesh," %lld ", (long long )sc->flags[k]);
  }
  fprintf(fmesh,"\n");
  for(j=0;j<dim;j++){
    for(k=0;k<n;k++){
      fprintf(fmesh," %23.16g ",x[k*dim+j]);
      /*      fprintf(stdout," (%i,%i) %23.16g ",k,j,x[k*dim+j]); */
    }
    fprintf(fmesh,"\n");
  }
  for(k=0;k<n;k++){
    fprintf(fmesh," %lld ", (long long )ib[k]);
  }
  fprintf(fmesh,"\n");
  //NOT USED: write a function value
  /* if(sc->fval){ */
  /*   for(k=0;k<n;k++){ */
  /*     fprintf(fmesh," %23.16g ", sc->fval[k]); */
  /*   } */
  /*   fprintf(fmesh,"\n"); */
  /* } */
  fprintf(stdout,"\n%%Output (hazmath) written on:%s\n",nameout);
  fclose(fmesh);
  return;
}
/********************************************************************************/
/*!
 * \fn void mshw(char *namemsh,scomplex *sc, const INT shift0)
 *
 * \brief Write a simplicial complex to a file in a ".msh" format. No
 *        boundary information is written. The data is organized as
 *        described at
 *        https://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
 *
 * \param namemsh   File name
 * \param sc        Pointer to a simplicial complex
 * \param shift0    integer added to the elements of arrays (here always=1).
 *
 */
/********************************************************************************/
void mshw(char *namemsh,scomplex *sc, const INT shift0)
{
  // WRITING in .msh format.
  INT shift=shift0;// this is fake because shift must be 1 below, nno zero node nnumbers:
  shift=1;
  INT n=sc->nv,ns=sc->ns, dimbig=sc->nbig,dim=sc->n,ndl=sc->n+1;
  INT *je = sc->nodes, *material=sc->flags;//*ib=sc->bndry;
  REAL *x = sc->x;
  INT k=-10,j=-10,kndl=-10;
  FILE *fmesh;
  fmesh=HAZ_fopen(namemsh,"w");
  /*
     MSH way of writing mesh file.
     fmt=0 is ASCII -- only this is supported now.
  */
  int fmt=0,num_el_tags=1;
  size_t data_size=sizeof(double);
  float ver=2.0;
  int el_type;
  switch(dim){
  case 1:
    el_type=1;//interval;
    break;
  case 2:
    el_type=2;//triangle 3 is the quad
    break;
  default:
    el_type=4;//tet
    break;
  }
  // writing:
  fprintf(fmesh,"%s\n","$MeshFormat");
  fprintf(fmesh,"%.1f %d %ld\n",ver,fmt,data_size);
  fprintf(fmesh,"%s\n","$EndMeshFormat");
  fprintf(fmesh,"%s\n","$Nodes");
  fprintf(fmesh,"%lld\n",(long long )n);
  for(k=0;k<n;k++){
    fprintf(fmesh,"%lld",(long long )(k+shift));
    for(j=0;j<dimbig;j++){
      fprintf(fmesh," %23.16e",x[k*dimbig+j]);
    }
    fprintf(fmesh,"\n");
  }
  fprintf(fmesh,"%s\n","$EndNodes");
  fprintf(fmesh,"%s\n","$Elements");
  fprintf(fmesh,"%lld\n",(long long )ns);
  for (k=0;k<ns;k++){
    // element number, 1 tag=material property of the element;
    fprintf(fmesh,"%lld %lld %lld %lld", (long long )(k+shift),\
	    (long long )el_type,			   \
	    (long long )num_el_tags,			   \
	    (long long )material[k]);
    for (j=0;j<ndl;j++) {
      kndl=ndl*k+j;
      /* shift if needed */
      fprintf(fmesh," %lld",(long long )(je[kndl]+shift));
    }
    fprintf(fmesh,"\n");
  }
  fprintf(fmesh,"%s\n","$EndElements");
  /* for(k=0;k<n;k++){ */
  /*   fprintf(fmesh," %lld ", (long long )ib[k]); */
  /* } */
  fprintf(fmesh,"\n");
  fprintf(stdout,"\n%%Output (MSH) written on:%s\n",namemsh);
  fclose(fmesh);
  return;
}
/*******************************************************************************/
/*******************************************************************************/
void matlw(scomplex *sc, char *namematl)
{
  FILE *fp;
  INT ns=sc->ns,nv=sc->nv,n=sc->n,n1=n+1,j=-10,k=-10;
  INT *nodes=sc->nodes;
  REAL *x=sc->x;
  fp=HAZ_fopen(namematl,"w");
  //  if(!fp) fp=stdout;
  if(!fp) fp=stdout;
  fprintf(stdout,"\n%lld %lld %lld\n",(long long )ns,(long long )nv,(long long )n);
  fflush(stdout);
  fprintf(fp,"\nt=[");
  fflush(fp);
  for (j=0;j<ns;j++){
    /*  for (j=0;j<ns;j++){*/
    for (k=0;k<n1;k++) {
      fprintf(fp,"%lld ",(long long )nodes[j*n1+k]+(long long )1);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"];\n");
  fprintf(fp,"\nx=[");
  for (j=0;j<nv;j++){
    for (k=0;k<n;k++) {
      fprintf(fp,"%23.16e ",x[j*n+k]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"];\n");
  if(n<3) {
    fprintf(fp,"\ntriplot(t,x(:,1),x(:,2));hold on");
  }else {
    fprintf(fp,"\ntetramesh(t,x(:,1),x(:,2),x(:,3));hold on");
  }
  fprintf(fp,"\nplot(x(:,1),x(:,2),'ro');");
  fprintf(fp,"\nplot(x(:,1),x(:,2),'r*'); hold off\n");
  fprintf(stdout,"%%Output (matlab) written on:%s\n",namematl);
  fclose(fp);
  return;
}
/*****************************************************************************************/
/*!
 * \fn void print_matlab_vector_field(dvector* ux, dvector* uy, dvector* uz, fespace* FE )
 *
 * \brief print a vector field of functions in a PX space in Matlab format.
 *
 * \param ux   Point to x component of vector field.
 * \param uy   Point to y component of vector field.
 * \param ux   Point to z component of vector field.
 * \param FE   FE space of each component (must be the same)
 *
 * \note Only designed for 3D for now.
 */
void print_matlab_vector_field(dvector* ux, dvector* uy, dvector* uz, fespace* FE )
{
  FILE *fid = fopen("output/usol_vfield.mat","w");
  INT i;
  for(i=0; i<ux->row; i++) {
    fprintf(fid,"%f\t%f\t%f\t%f\t%f\t%f\n",FE->cdof->x[i],FE->cdof->y[i],FE->cdof->z[i],ux->val[i],uy->val[i],uz->val[i]);
  }
  fclose(fid);
  return;
}


/************************* io functions returning pointers ***********/
/**
 * \fn void dcoo_read_eof_dcsr_p(const char *fname,const INT *size,const unsigned char format)
 *
 * \brief Read A from matrix disk file in I,J,V format and convert to CSR format
 *        first reads (as reals) until the end of file to count the number of nonzeroes
 *        then reqinds the file and reads all nonzeroes.
 *        The number of rows is the max of all integers from I[] and size[0]
 *        The number of columns is the max of all integers from J[] and size[1]
 *        If size is NULL is the same as size[0]=size[1]=0
 *
 * \param filename  File name for matrix
 * \param size      Pointer to an array with two integers or NULL.
 * \param format    if format='b' or 'B' reads npy format. if format='a' or 'A' it reads ASCII
 *
 * \note File format (ASCII):
 *   - i  j  a_ij        % i, j a_ij
 *
 * \return A         Pointer to the matrix in CSR format
 *
 */
dCSRmat *dcoo_read_eof_dcsr_p(const char *fname,const INT *size,const unsigned char format)
{
  SHORT ascii;
  if(format=='a' || format=='A'){
    ascii=TRUE;// ascii format
  } else if(format=='b' || format=='B'){
    ascii=FALSE; // binary npy
  } else{
    fprintf(stderr,"*** Wrong type in %s; must be \'A\' for ascii format type  or \'B\' for npy format",__FUNCTION__);
    exit(10);
  }
  INT i,i1,i2,nnz;
  dCOOmat acoo=dcoo_create(0,0,0);
  if(ascii){
    read_eof(fname,(void **)&acoo.val,&nnz,'R',0);
  } else {
    // READ NPY:
    read_npy_eof(fname,(void **)&acoo.val,&nnz,'R',0);
  }
  // check if we read anything...
  if(!nnz) check_error(ERROR_WRONG_FILE, __FUNCTION__);
  //
  nnz=(INT )nnz/3; // this is for coo matrix;
  acoo.rowind = (INT *)calloc(nnz, sizeof(INT));
  acoo.colind = (INT *)calloc(nnz, sizeof(INT));
  acoo.row=-1;acoo.col=-1;
  acoo.nnz=nnz;
  acoo.nnz=0;
  for(i=0;i<nnz;++i){
    i1=(INT )rint(acoo.val[3*i]);
    i2=(INT )rint(acoo.val[3*i+1]);
    if(i1>acoo.row)acoo.row=i1;
    if(i2>acoo.col)acoo.col=i2;
    acoo.rowind[i]=i1;
    acoo.colind[i]=i2;
    acoo.val[acoo.nnz]=acoo.val[3*i+2];
    acoo.nnz++;
  }
  acoo.row++;acoo.col++;
  acoo.val=realloc(acoo.val,acoo.nnz*sizeof(REAL));
  if(size){
    if(acoo.row<size[0])
      acoo.row=size[0];
    if(acoo.col<size[1])
      acoo.col=size[1];
  }
  dCSRmat *a=malloc(1*sizeof(dCSRmat));
  dcoo_2_dcsr(&acoo,a);
  //  dCSRmat *a=dcoo_2_dcsr_p(&acoo);
  dcoo_free(&acoo); // just one free if it was allocated with "something_p"
  return a;
}


/***********************************************************************************************/
/**
 * \fn dvector *dvector_read_eof_p(const char *fname, const unsigned char format)
 *
 * \brief Read b from a disk file until EOF. B is just a column of REALS numbers
 *
 * \param  filename  File name for vector b
 * \param format    if format='b' or 'B' reads npy format. if format='a' or 'A' it reads ASCII
 *
 * \note File Format:
 *   - val_j, j=0:nrow-1
 *
 * \return b         Pointer to the dvector b (output)
 *
 */
dvector *dvector_read_eof_p(const char *fname, const unsigned char format)
{
  dvector *b=malloc(1*sizeof(dvector));
  b->row=0;
  b->val=NULL;
  if(format=='a' || format=='A'){
    read_eof(fname,(void **)&b->val,&b->row,'R',0);
  } else if(format=='b' || format=='B'){
    // READ NPY:
    read_npy_eof(fname,(void **)&b->val,&b->row,'R',0);
  } else {
    fprintf(stderr,"*** Wrong type in %s; must be \'A\' for ascii format type  or \'B\' for npy format",__FUNCTION__);
    exit(11);
  }
  if(!b->row)
    check_error(ERROR_WRONG_FILE, __FUNCTION__);

  return b;
}

/*
 * \fn ivector *ivector_read_eof_p(const char *fname, const unsigned char format)
 *
 * \brief Read b from a disk file until EOF. B is just a column of INT32 numbers
 *
 * \param  filename  File name for vector b
 * \param format    if format='b' or 'B' reads npy format. if format='a' or 'A' it reads ASCII
 *
 * \note File Format:
 *   - val_j, j=0:nrow-1
 *
 * \return b         Pointer to the ivector b (output)
 *
 */
ivector *ivector_read_eof_p(const char *fname, const unsigned char format)
{
  ivector *b=malloc(1*sizeof(ivector));
  b->row=0;
  b->val=NULL;
  if(format=='a' || format=='A'){
    read_eof(fname,(void **)&b->val,&b->row,'I',0);
  } else if(format=='b' || format=='B'){
    // READ NPY:
    read_npy_eof(fname,(void **)&b->val,&b->row,'I',0);
  } else {
    fprintf(stderr,"*** Wrong type in %s; must be \'A\' for ascii format type  or \'B\' for npy format",__FUNCTION__);
    exit(11);
  }
  if(!b->row)
    check_error(ERROR_WRONG_FILE, __FUNCTION__);

  return b;
}
/*********************************************************************/
/**
 * \fn dCSRmat *dcoo_read_dcsr_p(FILE *fp)
 *
 * \brief Read A from a file in IJ format and convert to CSR format
 *
 * \param fp        File descriptor for the input file
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 */
dCSRmat *dcoo_read_dcsr_p(FILE *fp)
{
  INT  i,j,k,m,n,nnz;
  //INT  val;
  REAL value;
  INT offset;

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  char *buffer=malloc(512*sizeof(char));
  SHORT    status = SUCCESS;

  while ( status == SUCCESS ) {

      offset = ftell(fp);
      //val = fscanf(fp,"%s",buffer);
      fscanf(fp,"%s",buffer);
      if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
          fgets(buffer,512,fp); // skip rest of line
          continue;
      }
      else {
          break;
      }

  }
  free(buffer);
  // move back to the beginning of the current line
  fseek(fp, offset, SEEK_SET);

  // now start to read
  fscanf(fp,"%lld %lld %lld",(long long *)&m,(long long *)&n,(long long *)&nnz);

  dCOOmat *Atmp=dcoo_create_p(m,n,nnz);

  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%lld %lld %le", (long long *)&i, (long long *)&j, &value) != EOF ) {
      Atmp->rowind[k]=i; Atmp->colind[k]=j; Atmp->val[k] = value;
    }
    else {
      check_error(ERROR_WRONG_FILE, "dcoo_read_dcsr_p");
    }
  }
  dCSRmat *A=dcoo_2_dcsr_p(Atmp);
  free((void *)Atmp);
  return A;
}

/*********************************************************************/
/**
 * \fn dCSRmat *dcoo_read_dcsr_p(FILE *fp)
 *
 * \brief Read A from a file in IJ format and convert to CSR format (index starts at 1)
 *
 * \param fp        File descriptor for the input file
 * \param A         Pointer to the CSR matrix
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 * \note   This is for the case the index starts at 1 insread of 0
 *
 */
dCSRmat *dcoo_read_dcsr_p_1(FILE *fp)
{
  INT  i,j,k,m,n,nnz;
  //INT  val;
  REAL value;
  INT offset;

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  char *buffer=malloc(512*sizeof(char));
  SHORT    status = SUCCESS;

  while ( status == SUCCESS ) {

      offset = ftell(fp);
      //val =
      fscanf(fp,"%s",buffer);
      if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
          fgets(buffer,512,fp); // skip rest of line
          continue;
      }
      else {
          break;
      }

  }
  free(buffer);
  // move back to the beginning of the current line
  fseek(fp, offset, SEEK_SET);

  // now start to read
  fscanf(fp,"%lld %lld %lld",(long long *)&m,(long long *)&n,(long long *)&nnz);

  dCOOmat *Atmp=dcoo_create_p(m,n,nnz);

  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%lld %lld %le", (long long *)&i, (long long *)&j, &value) != EOF ) {
      Atmp->rowind[k]=i-1; Atmp->colind[k]=j-1; Atmp->val[k] = value;
    }
    else {
      check_error(ERROR_WRONG_FILE, "dcoo_read_dcsr_p");
    }
  }
  dCSRmat *A=dcoo_2_dcsr_p(Atmp);
  free((void *)Atmp);
  return A;
}


/*****************************************************************/
/**
 * \fn dvector *dvector_read_p(FILE *fp,dvector *b)
 *
 * \brief Read b from a disk file in array format
 *
 * \param filename  File name for vector b
 * \param b         Pointer to the dvector b (output)
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 *
 */
dvector *dvector_read_p(FILE *fp)
{
  dvector *b=NULL;
  INT  i, n;
  REAL value;
  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);
  fscanf(fp,"%lld",(long long *)&n);
  b=dvec_create_p(n);
  for ( i = 0; i < n; ++i ) {
    fscanf(fp, "%lg", &value);
    b->val[i] = value;
    if ( value > BIGREAL ) {
      fprintf(stderr,"### ERROR: Wrong value = %le\n", value);
      free((void *)b);
      exit(ERROR_INPUT_PAR);
    }
  }
  return b;
}
/********************************************************************/
/********************************************************************************/
/*!
 * \fn INT features_r(features *feat,scomplex *sc, const INT do_map)
 *
 * \brief Reads file with features (coordinates of points to refine
 *        around them). Allocates memory for the components of struct
 *        features (amr.h).
 *
 * \param sc    simplicial complex defining the FE grid.
 *
 * \param feat  feature structure:
 *
 * \param do_map if not zero, then the simplicial complex coordinates
 *               are mapped so that they are in the cube defined by
 *               the min and max coords of the data in "feat".
 *
 * \note Ludmil (20221020)
 */
/********************************************************************************/
INT features_r(features *feat,scomplex *sc, const INT do_map, const REAL scale)
{
  /*
     dimbig is the spatial dimension; dim can be either dimbig or
     something smaller as which will indicate the features are on a plane for which all coords x[dim:dimbig-1]=vfill
     Upon return feat->x[] here is dimbig by f->n array
  */
  INT dim=feat->n, dimbig=feat->nbig;
  /***************************************************/
  INT ichk,k=0,count=0,i,j,m;
  /* INT *scnjn=NULL; */
  /* REAL *xstar0=NULL; */
  REAL xtmp;
  char ch;
  if(!feat->fpf) {
    feat->nf=0;
    feat->x=NULL;
    fprintf(stderr,"****ERROR:Could not open file for reading in %s; no features read",__FUNCTION__);
    return 3;
  }
  /*
     Read csv file. When exported from excel often a file has 3-4 control
     chars at the beginning we need to skip these.
  */
  /*read special chars if any*/
  ch=fgetc(feat->fpf); while((INT )ch < 0){ch=fgetc(feat->fpf);count++;}
  if(count){
    fprintf(stdout,"%%Read: %lld control chars...", (long long int)count);
    fseek(feat->fpf,count*sizeof(char),SEEK_SET);
  } else rewind(feat->fpf);
  k=0;
  while(1){
    //    if(feof(feat->fpf) || k>40000) break;
    if(feof(feat->fpf)) break;
    for (j=0;j<dim;j++){
      ichk=fscanf(feat->fpf,"%lg", &xtmp);
      if(ichk<0) break;
      k++;
    }
  }
  k=k/dim;
  /* read now k numbers from the file */
  if(count){
    fseek(feat->fpf,count*sizeof(char),SEEK_SET);
  } else {
    rewind(feat->fpf);
  }
  feat->nf=k;
  /* we allocate always the max of dim or dimbig */
  if(dimbig>dim) {
    feat->x=(REAL *)calloc(dimbig*(feat->nf),sizeof(REAL));
    //    fprintf(stdout,"\ndimbig=%lld, dim=%lld",(long long int)dimbig,(long long int )dim);fflush(stdout);
  }else{
    feat->x=(REAL *)calloc(dim*(feat->nf),sizeof(REAL));
  }
  for (i=0;i<k;i++){
    for(j=0;j<dim;j++){
      count = fscanf(feat->fpf,"%lg", (feat->x+dim*i+j));
    }
  }
  //  fprintf(stdout,"%%%%Read %lld coords\n",(long long int)k);
  fclose(feat->fpf);
  // if dimbig>dim we fill the remaining coordinates with feat->fill.
  if(dimbig > dim) {
    k=feat->nf-1;
    for(i=k;i>0;i--){
      for(m=dim-1;m>=0;m--){
	feat->x[dimbig*i+m]=feat->x[dim*i+m];
	//	fprintf(stdout,"(%lld,%lld): %lld  %lld\n",(long long int)i,(long long int)m,(long long int)(dimbig*i+m),(long long int)(dim*i+m));  fflush(stdout);
      }
      for(m=dim;m<dimbig;m++){
	feat->x[dimbig*i+m]=feat->fill;
      }
    }
  }
  if(sc!=NULL && do_map){
    cube2simp *c2s=cube2simplex(dim);//now we have the vertices of the unit cube in bits
    REAL *vc=calloc(4*dim*c2s->nvcube,sizeof(REAL));
    REAL *xmintmp=vc;// maps to [0...0]
    REAL *xmaxtmp=xmintmp+dim*(c2s->nvcube-1); // last vertex
    REAL *xmin=xmaxtmp+dim*(c2s->nvcube-1); // last vertex
    REAL *xmax=xmin+dim*(c2s->nvcube-1); // last vertex
    INT kdimi;
    for(i=0;i<dim;i++){
      xmax[i]=feat->x[i];
      xmin[i]=xmax[i];
      kdimi=dim+i;
      for (k=1;k<feat->nf;k++){
	if(feat->x[kdimi] > xmax[i]){
	  xmax[i]=feat->x[kdimi];
	}
	if(feat->x[kdimi] < xmin[i]){
	  xmin[i]=feat->x[kdimi];
	}
	kdimi+=dim;
      }
    }
    // compute the diagonal length
    /* REAL diag0=0e0; */
    /* for(i=0;i<dim;i++){ */
    /*   diag0+=(xmax[i]-xmin[i])*(xmax[i]-xmin[i]); */
    /* } */
    /* diag0=pow(0.5,((REAL )feat->n))*(sqrt(diag0)/pow((REAL )feat->nf,1./((REAL )feat->n))); */
    /* //  fprintf(stdout,"\nNew:::diag-length=%.16e\n",diag0); */
    /* for(i=0;i<dim;i++){ */
    /*   xmax[i]+=diag0*scale; */
    /*   xmin[i]-=diag0*scale; */
    /* } */
    /* print_full_mat(1,dim,xmin,"xmin"); */
    /* print_full_mat(1,dim,xmax,"xmax"); */
    for(i=0;i<dim;i++){
      xmaxtmp[i]=xmax[i]+(scale-1e0)*(xmax[i]-xmin[i]);
      xmintmp[i]=xmin[i]-(scale-1e0)*(xmax[i]-xmin[i]);
    }
    /* print_full_mat(1,dim,xmintmp,"xmintmp"); */
    /* print_full_mat(1,dim,xmaxtmp,"xmaxtmp"); */
    for(j=1;j<c2s->nvcube-1;j++){
      for(i=0;i<dim;i++){
	vc[j*dim+i]=xmintmp[i]+(xmaxtmp[i]-xmintmp[i])*(c2s->bits[dim*j+i]);
      }
    }
    mapit(sc,vc);
    free(vc);
    cube2simp_free(c2s);
  }
  return 0;
}

/************************************************************************************************/
/*
 * \fn void read_npy_eof(const char *fname, void **data, INT *ndata, const char dtype,const SHORT print_level)
 *
 * \brief Reads data from npy file until EOF. skips any header, so it is insensitive to the shape.
 * The result is stored in data[0] which is allocated here (REAL* for
 *  double/float or INT* for ints).
 *
 * dtype is either 'R', or 'r' for REAL or 'I', 'i' for INT.
 *
 *
 *   NPY format description as found on GitHub:
 *   https://github.com/numpy/numpy/blob/main/doc/neps/nep-0001-npy-format.rst
 *
 * \note The first 6 bytes are a magic string: exactly "x93NUMPY". The
 *       next 1 byte is an unsigned byte: the major version number of
 *       the file format, e.g. x01.
 *
 *       The next 1 byte is an unsigned byte: the minor version number
 *       of the file format, e.g. x00. Note: the version of the file
 *       format is not tied to the version of the numpy package.
 *
 *       The next 2 bytes form a little-endian unsigned short int: the
 *       length of the header data HEADER_LEN.
 *
 *       The next HEADER_LEN bytes form the header data describing the
 *       array's format. It is an ASCII string which contains a Python
 *       literal expression of a dictionary. It is terminated by a
 *       newline ('n') and padded with spaces ('x20') to make the
 *       total length of the magic string + 4 + HEADER_LEN be evenly
 *       divisible by 16 for alignment purposes.
 *
 *       The dictionary contains three keys:
 *
 *        "descr" : dtype.descr: An object that can be passed as an
 *                  argument to the numpy.dtype() constructor to
 *                  create the array's dtype.
 *
 * "fortran_order": bool: Whether the array data is
 *                        Fortran-contiguous or not. Since
 *                        Fortran-contiguous arrays are a common form
 *                        of non-C-contiguity, we allow them to be
 *                        written directly to disk for efficiency.
 *
 * "shape" : tuple of int The shape of the array. For repeatability
 *           and readability, this dictionary is formatted using
 *           pprint.pformat() so the keys are in alphabetic order.
 *
 *         Following the header comes the array data. If the dtype
 *         contains Python objects (i.e. dtype.hasobject is True),
 *         then the data is a Python pickle of the array. Otherwise
 *         the data is the contiguous (either C- or Fortran-,
 *         depending on fortran_order) bytes of the array. Consumers
 *         can figure out the number of bytes by multiplying the
 *         number of elements given by the shape (noting that shape=()
 *         means there is 1 element) by dtype.itemsize.
 *
 * Format Specification: Version 2.0:
 *
 *         The version 1.0 format only * allowed the array header to
 *         have a total size of 65535 * bytes. This can be exceeded by
 *         structured arrays with a large number of columns. The
 *         version 2.0 format extends the header size to 4
 *         GiB. numpy.save will automatically save in 2.0 format if
 *         the data requires it, else it will always use the more
 *         compatible 1.0 format.
 *
 *         The description of the fourth element of the header
 *         therefore has become:
 *
 *         The next 4 bytes form a little-endian unsigned int: the
 *         length of the header data HEADER_LEN.
 *
 * 20221115 (ltz)
 *
 */
/************************************************************************************************/
void read_npy_eof(const char *fname, void **data, INT *ndata, const char dtype,const SHORT print_level)
/*************************************************************************************/
{
  SHORT data_type=-1;
  char *type_name=NULL;
  if(dtype=='i' || dtype=='I'){
    data_type=TRUE;// integers
    type_name=strdup("INTs");
  } else if(dtype=='r' || dtype=='R'){
    data_type=FALSE;
    type_name=strdup("REALs");
  } else{
    fprintf(stderr,"*** Wrong type in %s; must be \'I\' for integer type  or \'R\' for double/float",__FUNCTION__);
    exit(9);
  }
  if(print_level>10)
    fprintf(stdout,"\n%%%%Reading %s from %s\n",type_name,fname);
  FILE *fp = fopen(fname, "rb");
  if (!fp) {
    perror("fopen");
    return;
  }
  unsigned char buffer[8];// 8 bytes first
  size_t ret[4];
  unsigned short int nb = 8;
  ret[0]= fread(buffer,sizeof(unsigned char),nb, fp);
  if(print_level>10)
    fprintf(stdout,"%%%%magic=%#2x; header++=%s\n", buffer[0], &buffer[1]);
  ret[0]*=sizeof(unsigned char);
  unsigned short int vmajor=(unsigned short int )buffer[6];
  unsigned short int vminor=(unsigned short int )buffer[7];
  if(print_level>10)
    fprintf(stdout,"%%%%NPY format version=%1d.%1d\n", vmajor,vminor);
  unsigned short int buffint;
  unsigned int buffint2;
  size_t hl;
  if(vmajor<2){
    ret[1] = fread(&buffint,sizeof(unsigned short int),1, fp);
    ret[1]*=sizeof(unsigned short int);
    hl=(size_t )buffint;
  } else {
    ret[1] = fread(&buffint2,sizeof(unsigned int),1, fp);
    ret[1]*=sizeof(unsigned int);
    hl=(size_t )buffint2;
  }
  if(print_level>10)
    fprintf(stdout,"\n%%%%HEADER LENGTH (in chars)=%ld\n",hl);
  char *header=calloc(hl,sizeof(char));
  ret[2] = fread(header,sizeof(char),hl, fp);
  ret[2]*=sizeof(char);
  if(print_level>10)
    fprintf(stdout,"%%%%HEADER=%s\n",header);
  //
  size_t hb=ret[0]+ret[1]+ret[2];
  //
  if(print_level>10)
    fprintf(stdout,"\n%%%%Bytes_to_skip=%ld\n",hb);
  free(header);
  ndata[0]=0;
  size_t read_size=0;
  if(data_type)
    read_size=sizeof(INT);
  else
    read_size=sizeof(REAL);
  data[0]=calloc(1,read_size);
  while(1){
    if(feof(fp)) break;
    ret[3]=fread((REAL *)data[0],read_size,1,fp);
    if(ret[3]<0) break;
    ndata[0]++;
  }
  rewind(fp);
  fseek(fp,hb,SEEK_SET);
  data[0]=realloc(data[0],ndata[0]*read_size);
  ret[3]=fread((REAL *)data[0],read_size,ndata[0],fp);
  fclose(fp);
  ndata[0]=ret[3];
  if(print_level>10)
    fprintf(stdout,"\n%%%%Read %lld %s\n",(long long int)ndata[0],type_name);
  free(type_name);
  return;
}

/*********************************************************************/

/*EOF*/
