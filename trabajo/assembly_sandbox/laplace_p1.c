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
#include "meshes_inline.h"
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT dim=2;
  scomplex *sc;
  switch(dim){
  case 3:
    sc=mesh3d();
    break;
  default:
    sc=mesh2d();
  }
  dCSRmat *A=malloc(sizeof(dCSRmat));
  dCSRmat *M=malloc(sizeof(dCSRmat));
  assemble_p1(sc,A,M);
  INT print_level=0;
  if(print_level > 5)
    csr_print_matlab(stdout,M);
  dcsr_free(A);  
  dcsr_free(M);  
  haz_scomplex_free(sc);  
  return 0;
}


