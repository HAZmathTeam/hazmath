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
static dvector fe_set(scomplex *sc,				\
		      dCSRmat *A,dCSRmat *M,			\
		      const double alpha,			\
		      const double gamma)
{
  // gamma is the reaction coefficient; alpha is the diffusion coefficient
  REAL fi;
  INT i,ij,j;
  // Now we solve the Neumann problem (A+M) u = f; this also changes A->val 
  for(i=0;i<A->nnz;++i)
    A->val[i]=alpha*A->val[i] + gamma*M->val[i];
  /**/
  // RHS set up: We take as a right hand side the function f(x)=1;
  // the right hand side is just M*ones(n,1);
  dvector f=dvec_create(A->row);
  dvec_set(f.row,&f,1e0);// this we can have as function value at the vertices in general. 
  dvector rhs=dvec_create(M->row);// this is the "algebraic" right hand side 
  for(i=0;i<M->row;++i){
    fi=0e0;
    for(ij=M->IA[i];ij<M->IA[i+1];++ij){
      j=M->JA[ij];
      fi+=M->val[ij]*f.val[j];
    }
    rhs.val[i]=fi;
  }
  dvec_free(&f);
  return rhs;
}
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT dim=2;
  INT jlevel,k;
  scomplex *sc;
  switch(dim){
  case 3:
    sc=mesh3d();    
    break;
  default:
    sc=mesh2d();
  }
  scomplex *sctop=NULL;
  INT ref_levels=10;
  //
  ivector *marked=malloc(1*sizeof(ivector));
  marked->row=0;
  marked->val=NULL;
  dCSRmat A,M;
  dvector rhs;
  // end intialization 
  for(jlevel=0;jlevel<ref_levels;jlevel++){
    /* choose the finest grid */
    sctop=scfinest(sc);
    clock_t clk_assembly_start = clock(); // begin assembly timing;
    assemble_p1(sctop,&A,&M);
    // fe_set(mesh,stiffness,mass,diffusion_coeff, reaction_coefficient)
    rhs=fe_set(sctop,&A,&M,1.0,1.0);
    clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
    fprintf(stdout,"\n%%%%%%level=%d; CPUtime(assembly) = %.3f sec",jlevel,
	    (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
    dcsr_free(&A);
    dcsr_free(&M);
    dvec_free(&rhs);
    /* mark everything */
    marked->row=sctop->ns;
    marked->val=realloc(marked->val,marked->row*sizeof(INT));
    for(k=0;k<marked->row;k++) marked->val[k]=TRUE;
    // now we refine.
    refine(1,sc,marked);
    //    haz_scomplex_print(sctop,0,"In assembly sandbox");
    /* free */    
    haz_scomplex_free(sctop);
    /*  MAKE sc to be the finest grid only */
  }
  scfinalize(sc);
  clock_t clk_assembly_start = clock(); // begin assembly timing last level
  assemble_p1(sc,&A,&M);
  rhs=fe_set(sc,&A,&M,1.0,1.0);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%level=%d; CPUtime(assembly) = %.3f sec",ref_levels, \
	  (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  /* write the output mesh file:    */
  //  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview: can be viewed with
     "paraview out.vtu" */
  vtkw("out.vtu",sc,0,1.); 
  // since A and M have the same sparsity pattern we 
  INT print_level=0;
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  if(print_level > 5)
    csr_print_matlab(stdout,&M);  
  dcsr_free(&A);
  dcsr_free(&M);
  dvec_free(&rhs);
  haz_scomplex_free(sc);  
  return 0;
}


