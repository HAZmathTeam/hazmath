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
static dvector fe_sol(scomplex *sc,				\
		      const REAL alpha,				\
		      const REAL gamma)
{
  // fe_sol(mesh,stiffness,mass,diffusion_coeff, reaction_coefficient)
  // gamma is the reaction coefficient; alpha is the diffusion coefficient
  REAL fi;
  INT i,ij,j,solver_flag=-10,print_level=0;
  dCSRmat A,M;
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  assemble_p1(sc,&A,&M);
  // Now we solve the Neumann problem (A+M) u = f;
  for(i=0;i<A.nnz;++i)
    A.val[i]=alpha*A.val[i] + gamma*M.val[i];
  /*We have Neumann problem, so no boundary codes*/
  memset(sc->bndry,0,sc->nv*sizeof(INT));
  // RHS set up: We take as a right hand side the function f(x)=1;
  // the right hand side is just M*ones(n,1);
  dvector f=dvec_create(A.row);
  dvector sol=f;
  dvec_set(f.row,&f,1e0);// this we can have as function value at the vertices in general. 
  dvector rhs=dvec_create(M.row);// this is the "algebraic" right hand side 
  for(i=0;i<M.row;++i){
    fi=0e0;
    for(ij=M.IA[i];ij<M.IA[i+1];++ij){
      j=M.JA[ij];
      fi+=M.val[ij]*f.val[j];
    }
    rhs.val[i]=fi;
  }
  //  dvector_print(stdout,&rhs);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(assembly) = %.3f sec",
	  (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  if(print_level > 5)
    csr_print_matlab(stdout,&M);  
  dcsr_free(&M); // not needed
  /*SOLVER SOLVER*/
  // use the same as f for the solution;
  linear_itsolver_param linear_itparam;
  linear_itparam.linear_precond_type = PREC_AMG;
  //  linear_itparam.linear_precond_type = PREC_DIAG;
  //  linear_itparam.linear_precond_type = SOLVER_AMG;
  linear_itparam.linear_itsolver_type     = SOLVER_CG;
  linear_itparam.linear_stop_type         = STOP_REL_RES;
  // Solver parameters
  linear_itparam.linear_tol      = 1e-8;
  linear_itparam.linear_maxit    = 500;
  linear_itparam.linear_print_level    = 10;
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_print(&amgparam);
  memset(sol.val,0,sol.row*sizeof(REAL));
  switch(linear_itparam.linear_precond_type){
  case SOLVER_AMG:
    // Use AMG as iterative solver
    solver_flag = linear_solver_amg(&A, &rhs, &sol, &amgparam);
    break;
  case PREC_AMG:
    // Use Krylov Iterative Solver with AMG
    solver_flag = linear_solver_dcsr_krylov_amg(&A, &rhs, &sol, &linear_itparam, &amgparam);
    break;
  case  PREC_DIAG:
    solver_flag = linear_solver_dcsr_krylov_diag(&A, &rhs, &sol, &linear_itparam);
    break;
  default:
    solver_flag = linear_solver_dcsr_krylov_amg(&A, &rhs, &sol, &linear_itparam, &amgparam);
    break;
  }
  dcsr_free(&A);
  dvec_free(&rhs);
  return sol;// return solution
}
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT dim=4;
  INT jlevel,k;
  scomplex *sc;
  switch(dim){
  case 4:
    sc=mesh4d();    
    break;
  case 3:
    sc=mesh3d();    
    break;
  default:
    sc=mesh2d();
  }
  scomplex *sctop=NULL;
  INT ref_levels=2;
  //
  ivector marked;
  marked.row=0;
  marked.val=NULL;
  dvector sol;
  // end intialization 
  for(jlevel=0;jlevel<ref_levels;jlevel++){
    /* choose the finest grid */
    sctop=scfinest(sc);
    sol=fe_sol(sctop,1.0,1.0);
    /* mark everything */
    marked.row=sctop->ns;
    marked.val=realloc(marked.val,marked.row*sizeof(INT));
    for(k=0;k<marked.row;k++) marked.val[k]=TRUE;
    // now we refine.
    refine(1,sc,&marked);
    //    haz_scomplex_print(sctop,0,"In assembly sandbox");
    /* free */
    dvec_free(&sol);
    haz_scomplex_free(sctop);
    /*  MAKE sc to be the finest grid only */
  }
  scfinalize(sc);
  vtkw("out.vtu",sc,0,1.); 
  sol=fe_sol(sc,1.0,1.0);
  /* write the output mesh file:    */
  //  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview: can be viewed with
     "paraview out.vtu" */
  // since A and M have the same sparsity pattern we 
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  dvec_free(&sol);
  ivec_free(&marked);
  haz_scomplex_free(sc);  
  return 0;
}


