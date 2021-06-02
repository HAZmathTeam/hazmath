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
  INT i,j,k,ij;
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
    marked->val=realloc(marked->val,marked->row*sizeof(INT));
    for(k=0;k<marked->row;k++) marked->val[k]=TRUE;
    // now we refine.
    refine(1,sc,marked);
    /* free */
    haz_scomplex_free(sctop);
    /*  MAKE sc to be the finest grid only */
  }
  scfinalize(sc);
  /* write the output mesh file:    */
  //  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview: can be viewed with
     "paraview out.vtu" */
   vtkw("out.vtu",sc,0,1.); 
  //
  dCSRmat *A=malloc(sizeof(dCSRmat));
  dCSRmat *M=malloc(sizeof(dCSRmat));
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  assemble_p1(sc,A,M);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n\n%%%%%%elapsed CPU time for assembling the global mass and stiffness matrices = %.3f seconds.\n\n",
         (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  // Now we solve the Neumann problem (A+M) u = f;
  for(i=0;i<A->nnz;++i){
    A->val[i]+=M->val[i];
  }
  // RHS set up: We take as a right hand side the function f(x)=1; so
  //  our right hand side is just M*ones(n,1);
  dvector f=dvec_create(A->row);
  dvec_set(f.row,&f,1e0);// this we can have as function value at the vertices in general. 
  dvector rhs=dvec_create(M->row);// this is the "algebraic" right hand side 
  for(i=0;i<M->row;++i){
    rhs.val[i]=0e0;
    for(ij=M->IA[i];ij<M->IA[i+1];++ij){
      rhs.val[i]+=M->val[ij]*f.val[j];
    }
  }
  // since A and M have the same sparsity pattern we 
  INT print_level=0;
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  if(print_level > 5)
    csr_print_matlab(stdout,M);  
  dcsr_free(A);  
  dcsr_free(M);  
  haz_scomplex_free(sc);  
  return 0;
}


