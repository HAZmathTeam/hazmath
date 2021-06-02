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
  INT dim=3;
  INT j,k;
  scomplex *sc;
  switch(dim){
  case 3:
    sc=mesh3d();    
    break;
  default:
    sc=mesh2d();
  }
  scomplex *sctop=NULL;
  INT ref_levels=6;
  //
  ivector *marked=malloc(1*sizeof(ivector));
  marked->row=0;
  marked->val=NULL;
  // end intialization 
  for(j=0;j<ref_levels;j++){
    /* choose the finest grid */
    sctop=scfinest(sc);
    /* mark everything */
    marked->row=sctop->ns;
    marked->val=realloc(marked->val,sctop->ns*sizeof(INT));
    for(k=0;k<sctop->ns;k++) marked->val[k]=TRUE;
    // now we refine.
    refine(1,sc,marked);
    /* free */
    haz_scomplex_free(sctop);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc);
  /* write the output mesh file:    */
  //  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtkw("out.vtu",sc,0,1.);
  //
  dCSRmat *A=malloc(sizeof(dCSRmat));
  dCSRmat *M=malloc(sizeof(dCSRmat));
  assemble_p1(sc,A,M);
  INT print_level=0;
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  if(print_level > 5)
    csr_print_matlab(stdout,M);  
  dcsr_free(A);  
  dcsr_free(M);  
  haz_scomplex_free(sc);  
  return 0;
}


