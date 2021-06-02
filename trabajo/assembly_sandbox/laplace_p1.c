/*! \file src/amr/laplace_p1.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note Solving Laplace equation with linear finite elements on sequence of grids. 
 *
 *  \note: modified by ltz on 20210531
 *
 */
/*********************************************************************/
#include "hazmath.h"
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT dim=2;
  char *meshfile=NULL;
  switch(dim){
  case 3:
    meshfile=strdup("../../examples/grids/3D/unitCUBE_n9.haz");
    break;
  default:
    meshfile=strdup("../../examples/grids/2D/unitSQ_n3.haz");
  }
  FILE *fp=fopen(meshfile,"r");
  fprintf(stdout,"\n%%%%%%Reading mesh file(dim=%d)...",dim);fflush(stdout);
  scomplex *sc=haz_scomplex_read(fp,0);// 0 is the print level; reads the mesh
  fprintf(stdout,"done;\n");fflush(stdout);
  dCSRmat *A=malloc(sizeof(dCSRmat));
  dCSRmat *M=malloc(sizeof(dCSRmat));
  assemble_p1(sc,A,M);
  INT print_level=6;
  if(print_level > 5)
    csr_print_matlab(stdout,M);
  dcsr_free(A);  
  dcsr_free(M);  
  haz_scomplex_free(sc);  
  return 0;
}
