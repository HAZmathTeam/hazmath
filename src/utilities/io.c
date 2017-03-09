/*! \file src/utilities/io.c
 *
 *  Created by James Adler and Xiaozhe Hu on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  Routines for reading and writing to file or screen.
 *
 * \note: modified by Xiaozhe Hu on 10/31/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
void iarray_print(INT *vec,
                  INT n)
{

    /*!
     * \fn void iarray_print(INT *vec, INT n)
     *
     * \brief print an integer array on screen
     *
     * \param vec   Pointer to the INT array
     * \param n     Length of the array
     *
     */

    /* prints a vector of integers of size n */
    INT *vec_end;
    
    vec_end  =  vec + n;
    
    fprintf(stdout,"\n");

    for ( ; vec < vec_end; ++vec)
        fprintf(stdout, "%i\n  ",*vec);
    
    fprintf(stdout,"\n");
    
    return;
}

/***********************************************************************************************/
void array_print(REAL *vec,
                 INT n)
{

    /*!
     * \fn void array_print(REAL *vec, INT n)
     *
     * \brief print a REAL array on screen
     *
     * \param vec   Pointer to the REAL array
     * \param n     Length of the array
     *
     */

    /* prints a vector of integers of size nn */
    REAL *vec_end;
    
    vec_end  =  vec + n;
    
    fprintf(stdout,"\n");
    
    for ( ; vec < vec_end; ++vec)
        fprintf(stdout, "%e\n  ",*vec);
    
    fprintf(stdout,"\n");
    
    return;

}

/***********************************************************************************************/
void dvector_print(FILE* fid,
                   dvector *b)
{
    /*!
     * \fn void dvector_print(FILE* fid,dvector *b)
     *
     * \brief print a dvector to a file
     *
     * \param fid  Pointer to the file
     * \param b    Pointer to the dvector
     *
     */

  /* prints a dvector in matlab output*/
  INT i; /* Loop Indices */

  for(i=0;i<b->row;i++) {
    fprintf(fid,"%25.16e\n",b->val[i]);
  }
  return;
}

/***********************************************************************************************/
void csr_print_matlab(FILE* fid,
                      dCSRmat *A)
{

    /*!
     * \fn void csr_print_matlab(FILE* fid,dCSRmat *A)
     *
     * \brief print a dCSRmat format sparse matrix to a file
     *
     * \param fid  Pointer to the file
     * \param A    Pointer to the dCSRmat format sparse matrix
     *
     */

  /* prints a csr matrix in matlab output*/
  INT i,j1,j2,j; /* Loop Indices */
  INT shift_flag = 0; /* Check if Indexing starts at 0 or 1 */

  if(A->IA[0]==0) {
    printf("hello\n\n");
    dcsr_shift(A, 1);  // shift A
    shift_flag = 1;
  }

  for(i=0;i<A->row;i++) {
    j1 = A->IA[i]-1;
    j2 = A->IA[i+1]-1;
    for(j=j1;j<j2;j++) {
      fprintf(fid,"%d\t%d\t%25.16e\n",i+1,A->JA[j],A->val[j]);
    }
  }

  if(shift_flag==1) {
    dcsr_shift(A, -1);  // shift A back
  }

  return;
}

/***********************************************************************************************/
void icsr_print_matlab(FILE* fid,iCSRmat *A)
{
    /*!
     * \fn void icsr_print_matlab(FILE* fid,dCSRmat *A)
     *
     * \brief print a iCSRmat format sparse matrix to a file
     *
     * \param fid  Pointer to the file
     * \param A    Pointer to the iCSRmat format sparse matrix
     *
     */

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

/***********************************************************************************************/
void dvec_write (const char *filename,
                 dvector *vec)
{
    /*!
     * \fn void dvec_write (const char *filename, dvector *vec)
     *
     * \brief Write a dvector to disk file
     *
     * \param vec       Pointer to the dvector
     * \param filename  File name
     *
     */

    INT m = vec->row, i;

    FILE *fp = fopen(filename,"w");

    if ( fp == NULL ) {
        printf("### ERROR: Cannot open %s!\n", filename);
        check_error(ERROR_OPEN_FILE, __FUNCTION__);
    }

    printf("%s: writing to file %s...\n", __FUNCTION__, filename);

    fprintf(fp,"%d\n",m);

    for ( i = 0; i < m; ++i ) fprintf(fp,"%0.15e\n",vec->val[i]);

    fclose(fp);
}

/***********************************************************************************************/
void dcsr_write_dcoo (const char *filename,
                      dCSRmat *A)
{
    /*!
     * \fn void dcsr_write_dcoo (const char *filename, dCSRmat *A)
     *
     * \brief Write a dCSRmat matrix to disk file in IJ format (coordinate format)
     *
     * \param A         pointer to the dCSRmat matrix
     * \param filename  char for vector file name
     *
     */
    
    const INT m = A->row, n = A->col;
    INT i, j;
    
    FILE *fp = fopen(filename, "w");
    
    if ( fp == NULL ) {
        printf("### ERROR: Cannot open %s!\n", filename);
        check_error(ERROR_OPEN_FILE, __FUNCTION__);
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
/* */
{
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
{
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
void dump_sol_onV_vtk(char *namevtk,trimesh *mesh,REAL *sol,INT ncomp)
{
  /*!
   * \fn void dump_sol_onV_vtk(char *namevtk,trimesh *mesh,REAL *sol,INT ncomp)
   *
   * \brief Dumps solution data to vtk format
   *
   * \param namevtk  Filename
   * \param mesh     Mesh struct to dump
   * \param sol      solution vector to dump
   * \param ncomp:   Number of components to the solution
   *
   */

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
  for(k=1;k<=nelm;k++) fprintf(fvtk," %i ",mesh->el_v->IA[k]-1);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %i ",mesh->el_v->JA[kndl + j]-1);
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
void create_pvd(char *namepvd,INT nfiles,char *vtkfilename,char *filetype)
{

  /*!
   * \fn void create_pvd(char *namevtk,trimesh *mesh,REAL *sol,INT ncomp)
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
