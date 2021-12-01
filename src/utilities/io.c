/*! \file src/utilities/io.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  Routines for reading and writing to file or screen.
 *
 * \note: modified by Xiaozhe Hu on 10/31/2016
 * \note: modified 20171119 (ltz) added output using scomplex
 * \note: modifided by jha on 02/22/2019 for 0-1 fix.  need somone to check last three functions by ltz
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
    fprintf(stdout, "\n Input value too large: %d ; Changing to the max allowed: %d\n",n,nmax);
    nnew=nmax;
  } else if(n < nmin) {
    fprintf(stdout, "\n Input value too small: %d ; Changing to the min allowed: %d\n",n,nmin);
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
    fprintf(stdout, "%i\n  ",*vec);

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
    fprintf(stdout, "%e\n  ",*vec);

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
  fprintf(fp,"%d\n",m);

  // write index and value each line
  for ( i = 0; i < m; ++i ) fprintf(fp,"%d %d\n",i,vec->val[i]);

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

  fprintf(fp,"%d\n",m);

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

  int  i, n;
  REAL value;

  FILE *fp = fopen(filename,"r");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is reading file %s...\n", __FUNCTION__, filename);

  fscanf(fp,"%d",&n);

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
    fprintf(stdout,"i = %d\n", i);
    fprintf(fid,"%d\n",b->val[i]);
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
  INT nprt=1025,mprt=1025;
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
  INT nprt=1025,mprt=1025;
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
  INT nprt=1025,mprt=1025;
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
      fprintf(stdout,"%16i ", A[m*i+j]);
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
      fprintf(fid,"%d\t%d\t%25.16e\n",i+1,A->JA[j]+1,A->val[j]);
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
  fprintf(fid,"%d %d %d\n",A->row, A->col, A->nnz);
  for(i=0;i<(A->row+1);i++)
    fprintf(fid,"%d ",A->IA[i]+shift);
  fprintf(fid,"\n");
  for(i=0;i<A->nnz;i++)
    fprintf(fid,"%d ",A->JA[i]+shift);
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
      fprintf(fid,"%d,%d\n",i+1,A->JA[j]+1);
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
      fprintf(fid,"\nrow[%d]=[ ",i);
      for(j=A->IA[i];j<A->IA[i+1];++j){
	fprintf(fid,"%d ",A->JA[j]);
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
      fprintf(fid,"%d,%d,%d\n",i+1,A->JA[j]+1,A->val[j]);
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

  fprintf(fp,"%d\n",m);

  //main loop
  for ( i = 0; i < m; ++i ) fprintf(fp,"%0.15e\n",vec->val[i]);

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
  fprintf(fp,"%d  %d\n", n, m);
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
  fprintf(fp,"%d  %d  %d\n",m,n,A->nnz);
  for ( i = 0; i < m; ++i ) {
    for ( j = A->IA[i]; j < A->IA[i+1]; j++ )
      fprintf(fp,"%d  %d  %0.15e\n",i,A->JA[j],A->val[j]);
  }

  fclose(fp);
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
  int  i,j,k,m,n,nnz;
  REAL value;

  FILE *fp = fopen(filename,"r");

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fprintf(stdout,"%%%%%s: HAZMATH is reading file %s...\n", __FUNCTION__, filename);

  fscanf(fp,"%d %d %d",&m,&n,&nnz);

  dCOOmat Atmp=dcoo_create(m,n,nnz);

  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%d %d %le", &i, &j, &value) != EOF ) {
      Atmp.rowind[k]=i; Atmp.colind[k]=j; Atmp.val[k] = value;
    }
    else {
      check_error(ERROR_WRONG_FILE, "dcoo_read_dcsr");
    }
  }

  fclose(fp);

  dcoo_2_dcsr(&Atmp,A);
  dcoo_free(&Atmp);
}

/*** Auxillary Files *********************************************************************/
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

  // main loop
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%i",vec);
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
  fprintf(fvtk, \
	  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n", \
	  endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",nv,nelm);
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
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %i\" Format=\"ascii\">",tfloat,i);
    for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[i*nv+k]);
    fprintf(fvtk,"</DataArray>\n");
  }
  fprintf(fvtk,"</PointData>\n");

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %i ",mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %i ",mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %i ",tcell);

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
  fprintf(fvtk,"<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",nv,nelm);
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
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component - %s%i\" Format=\"ascii\">",tfloat,varname,i);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol_on_V[i*nv+k]);
      fprintf(fvtk,"</DataArray>\n");
    }
    fprintf(fvtk,"</PointData>\n");
  }

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %i ",mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %i ",mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %i ",tcell);

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
  fprintf(fvtk,"<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",nv,nelm);
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
    } else if(FE->var_spaces[nsp]->FEtype>0 && FE->var_spaces[nsp]->FEtype<20) { // PX elements (assume sol at vertices comes first)
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %i - %s\" Format=\"ascii\">",tfloat,nsp,varname[nsp]);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[spcntr + k]);
      fprintf(fvtk,"</DataArray>\n");
    } else if(FE->var_spaces[nsp]->FEtype==99) { // Single DoF constraint element (just plot single value everywhere)
      fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %i - %s\" Format=\"ascii\">",tfloat,nsp,varname[nsp]);
      for(k=0;k<nv;k++) fprintf(fvtk," %23.16e ",sol[spcntr]);
      fprintf(fvtk,"</DataArray>\n");
    } else { // Vector Elements
      sol_on_V = (REAL *) calloc(dim*mesh->nv,sizeof(REAL));
      solptr = sol+spcntr;
      Project_to_Vertices(sol_on_V,solptr,FE->var_spaces[nsp],mesh);
      for(i=0;i<dim;i++) {
        fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %i - %s%i\" Format=\"ascii\">",tfloat,nsp,varname[nsp],i);
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
        fprintf(fvtk,"<DataArray type=\"%s\" Name=\"Solution Component %i - %s\" Format=\"ascii\">",tfloat,nsp,varname[nsp]);
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
  for(k=1;k<=nelm;k++) fprintf(fvtk," %i ",mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %i ",mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %i ",tcell);

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
    sprintf(filecounter,"%s%03d.vtu",vtkfilename,i);
    fprintf(fvtk,"<DataSet %s=\"%d\" group=\"\" part=\"0\" file=\"%s\"/>\n",filetype,i,filecounter);
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
  fscanf(fmeshin,"%i %i %i %i\n",&ns,&nv,&dim,&nholes);
  scomplex *sc = (scomplex *)haz_scomplex_init(dim,ns,nv,dim); // cannot read different dimensions.
  INT dim1=sc->n+1;
  for (j=0;j<dim1;j++) {
    for (k=0;k<ns;k++){
      m=dim1*k+j;
      /* shift if needed */
      fscanf(fmeshin,"%i", sc->nodes+m);
      //      fprintf(stdout,"\n%d",sc->nodes[m]);
    }
  }
  for (k=0;k<sc->ns;k++){
    fscanf(fmeshin,"%d", sc->flags+k);
    //      fprintf(stdout,"\n%d",sc->flags[k]);
  }
  for(j=0;j<dim;j++){
    for(k=0;k<nv;k++){
      fscanf(fmeshin,"%lg",sc->x+k*dim+j);
    }
  }
  for(k=0;k<nv;k++){
    fscanf(fmeshin,"%d", sc->bndry+k);
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

/**********************************************************************
 * Routines to save to file the mesh in different formats.
 * uses the simplicial complex data structure (scomplex *sc).
 *
 *************************************************************************/
void hazw(char *nameout,scomplex *sc, const int shift)
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
  fprintf(fmesh,"%i %i %i %i\n",ns,n,dim,(sc->bndry_cc-1)); /* this is the
							       number of
							       holes;*/
  /* fprintf(stdout,"%i %i %li\n",n,ns,sizeof(ib)/sizeof(INT)); */
  for (j=0;j<ndl;j++) {
    for (k=0;k<ns;k++){
      kndl=ndl*k+j;
      /* shift if needed */
      fprintf(fmesh," %d ", je[kndl]+shift);
    }
    fprintf(fmesh,"\n");
  }
  for (k=0;k<ns;k++){
    fprintf(fmesh," %d ", sc->flags[k]);
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
    fprintf(fmesh," %i ", ib[k]);
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
/* WRITE mesh on VTU file*/
void vtkw(char *namevtk, scomplex *sc, const INT shift, const REAL zscale)
{
  if((sc->n!=2)&&(sc->n!=3))
    fprintf(stderr,"\n*** ERR(%s; dim=%d): NO vtk files for dim .eq. 1 or (dim .gt. 3).\n",__FUNCTION__,sc->n);
  FILE *fvtk;
  INT nv=sc->nv,ns=sc->ns, n=sc->n,n1=n+1;
  INT *nodes = sc->nodes, *ib=sc->bndry;
  REAL *x = sc->x;
  INT tcell=-10;
  INT k=-10,j=-10;
  char *tfloat="Float64", *tinto="Int64", *endian="LittleEndian";
  /*
    what endian?:

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
  if(n==1)
    tcell=LINE; /* triangle */
  else if(n==2)
    tcell=TRI; /* triangle */
  else
    tcell=TET; /* tet */

  /* VTK format writing the mesh for plot */
  fvtk=HAZ_fopen(namevtk,"w");
  fprintf(fvtk,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",nv,ns);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">",tfloat);
  for (j=0;j<nv;j++){
    for (k=0;k<n;k++) {
      fprintf(fvtk,"%.8f ",x[j*n+k]);
    }
    for (k=0;k<(3-n);k++) {
      fprintf(fvtk,"%.8f ",0.);
    }
  }
  /* if(n == 2){ */
  /*   for (j=0;j<nv;j++){ */
  /*     for (k=0;k<n;k++) { */
  /* 	fprintf(fvtk,"%.8f ",x[j*n+k]); */
  /*     } */
  /*     fprintf(fvtk,"%.8f ",0.); */
  /*   } */
  /* }  else { */
  /*   for (j=0;j<nv;j++){ */
  /*     for (k=0;k<n-1;k++) { */
  /* 	fprintf(fvtk,"%.8f ",x[j*n+k]); */
  /*     } */
  /*     fprintf(fvtk,"%.8f ",x[j*n+n-1]*zscale); */
  /*   } */
  /* } */
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");
  fprintf(fvtk,"<CellData Scalars=\"scalars\">\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"%s (layer)\" Format=\"ascii\">",tfloat,"L");
  for(k=0;k<ns;k++) fprintf(fvtk," %g ",(REAL )sc->flags[k]);//element flags
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</CellData>\n");
  // Dump v_bdry Data to indicate if vertices are boundaries
  fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"v_bdry\" Format=\"ascii\">",tinto);
  for(k=0;k<nv;k++) fprintf(fvtk," %i ",ib[k]);
  fprintf(fvtk,"</DataArray>\n");
  /*NOT USED: if(sc->fval){ */
  /*   fprintf(fvtk,"<DataArray type=\"%s\" Name=\"ele\" Format=\"ascii\">",tfloat); */
  /*   for(k=0;k<nv;k++) fprintf(fvtk," %e ",sc->fval[k]); */
  /*   fprintf(fvtk,"</DataArray>\n"); */
  /* } */
  // Dump information about connected components.  For now only assume
  // 1 connected region and at most 2 connected boundaries.  Positive
  // integers indicate connected components of a domain Negative
  // integers indicate connected components of the boundaries Example:
  // A cube (1 connected domain and 1 connected boundary) would be 1
  // on the interior and -1 on points on the boundary A cube with a
  // hole (1 connected domain and 2 connected boundaries) would have 1
  // on the points in the interior and -1 on points on the outer
  // boundary and -2 on the inner boundary If NULL, then one connected
  // region and boundary.
  if(sc->bndry_cc>1) {
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectedcomponents\" Format=\"ascii\">",tinto);
    for(k=0;k<nv;k++) {
      if(ib[k]==0) {
	fprintf(fvtk," %i ",1);
      } else if(ib[k]==1) {
	fprintf(fvtk," %i ",-1);
      } else if(ib[k]==-1) {
	fprintf(fvtk," %i ",-2);
      } else {
	fprintf(fvtk," %i ",ib[k]);
      }
    }
    fprintf(fvtk,"</DataArray>\n");
  }
  fprintf(fvtk,"</PointData>\n");
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=ns;k++) fprintf(fvtk," %i ",k*n1);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  /* for(k=0;k<ns;k++){ */
  /*   kn1=k*n1; */
  /*   for(j=0;j<n1;j++) fprintf(fvtk," %i ",nodes[kn1 + j]); */
  /* } */
  for (j=0;j<ns;j++){
    /*  for (j=0;j<ns;j++){*/
    for (k=0;k<n1;k++) {
      fprintf(fvtk,"%d ",nodes[j*n1+k]+shift);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=ns;k++)
    fprintf(fvtk," %i ",tcell);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Cells>\n");
  fprintf(fvtk,"</Piece>\n");
  fprintf(fvtk,"</UnstructuredGrid>\n");
  fprintf(fvtk,"</VTKFile>\n");
  fprintf(stdout,"%%Output (vtk) written on:%s\n",namevtk);
  fclose(fvtk);
  return;
}
void matlw(scomplex *sc, char *namematl)
{
  FILE *fp;
  INT ns=sc->ns,nv=sc->nv,n=sc->n,n1=n+1,j=-10,k=-10;
  INT *nodes=sc->nodes;
  REAL *x=sc->x;
  fp=HAZ_fopen(namematl,"w");
  //  if(!fp) fp=stdout;
  if(!fp) fp=stdout;
  fprintf(stdout,"\n%i %i %i\n",ns,nv,n);
  fflush(stdout);
  fprintf(fp,"\nt=[");
  fflush(fp);
  for (j=0;j<ns;j++){
    /*  for (j=0;j<ns;j++){*/
    for (k=0;k<n1;k++) {
      fprintf(fp,"%i ",nodes[j*n1+k]+1);
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
 * \fn void dcoo_read_eof_dcsr_p(FILE *fp,INT *size)
 *
 * \brief Read A from matrix disk file in I,J,V format and convert to CSR format
 *        first reads until the end of file to count the number of nonzeroes
 *        then reqinds the file and reads all nonzeroes.
 *        The number of rows is the max of all integers from I[] and size[0]
 *        The number of columns is the max of all integers from J[] and size[1]
 *        If size is NULL is the same as size[0]=size[1]=0
 *
 * \param filename  File name for matrix
 * \param size      Pointer to an array with two integers or NULL.
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 * \return A         Pointer to the matrix in CSR format
 *
 */
dCSRmat *dcoo_read_eof_dcsr_p (FILE *fp,INT *size)
{
  int  i,j,k,m,n,nnz,ichk;
  REAL value;
  dCSRmat *A=NULL;
  if(size){
    m=size[0];
    n=size[1];
  } else {
    m=-1;
    n=-1;
  }
  nnz=0;
  while(!feof(fp)){
    ichk=fscanf(fp,"%d %d %lg", &i,&j, &value);
    if(ichk<0) break;
    if(i>m) m=i;
    if(j>n) n=j;
    nnz++;
  }
  fprintf(stdout,"%s FOUND %d records.\n",__FUNCTION__,nnz);
  // check if we read anything...
  if(!nnz) check_error(ERROR_WRONG_FILE, __FUNCTION__);
  //
  m++;n++; //bc indices strat from zero so row and coll are with one more.
  dCOOmat *Atmp=dcoo_create_p(m,n,nnz);
  //    fprintf(stdout,"\nCheck structure: %d %d %d\n\n",Atmp->row,Atmp->col,Atmp->nnz);fflush(stdout);
  rewind(fp);
  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%d %d %lg", &i, &j, &value) != EOF ) {
      Atmp->rowind[k]=i; Atmp->colind[k]=j; Atmp->val[k] = value;
    } else {
      check_error(ERROR_WRONG_FILE, __FUNCTION__);
    }
  }
  A=dcoo_2_dcsr_p(Atmp);
  free(Atmp); // just one free if it was allocated with "something_p"
  return A;
}
/***********************************************************************************************/
/**
 * \fn dvector *dvector_read_eof_p(FILE *fp)
 *
 * \brief Read b from a disk file until EOF in array format
 *
 * \param  filename  File name for vector b
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 *
 * \return b         Pointer to the dvector b (output)
 *
 */
dvector *dvector_read_eof_p(FILE *fp)
{
  int  i, n,ichk;
  REAL value;
  dvector *b=NULL;
  n=0;
  while(!feof(fp)){
    ichk=fscanf(fp,"%lg", &value);
    if(ichk<0) break;
    n++;
  }
  fprintf(stdout,"%s: FOUND %d records.\n",__FUNCTION__,n);
  if(!n)
    check_error(ERROR_WRONG_FILE, __FUNCTION__);
  b=dvec_create_p(n);
  rewind(fp);
  for ( i = 0; i < n; ++i ) {
    fscanf(fp, "%lg", &value);
    b->val[i] = value;
    if ( value > BIGREAL ) {
      fprintf(stderr,"### ERROR: Wrong value = %lf\n", value);
      fclose(fp);
      if(b) free(b);
      exit(ERROR_INPUT_PAR);
    }
  }
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
  int  i,j,k,m,n,nnz;
  REAL value;

  if ( fp == NULL ) check_error(ERROR_OPEN_FILE, __FUNCTION__);

  fscanf(fp,"%d %d %d",&m,&n,&nnz);

  dCOOmat *Atmp=dcoo_create_p(m,n,nnz);

  for ( k = 0; k < nnz; k++ ) {
    if ( fscanf(fp, "%d %d %le", &i, &j, &value) != EOF ) {
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
  fscanf(fp,"%d",&n);
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

/*EOF*/
