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

/******************************************************************************/
void dump_sol_onV_vtk(char *namevtk,trimesh *mesh,REAL *sol,INT ncomp)
{

  /* Dumps solution data to vtk format 
  *
  * Input:
  *   mesh:     Mesh struct to dump
  *    sol:     solution vector to dump
  *  ncomp:     Number of components to the solution
  * 
  * Output:
  *  namevtk  File name of vtk file
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
  const INT TRI=5;  
  const INT TET=10;
  
  if(dim==2) 
    tcell=TRI; /* triangle */
  else 
    tcell=TET; /* tet */

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
  if(dim == 2) {
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
