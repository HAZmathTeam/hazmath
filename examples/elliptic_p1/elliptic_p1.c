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
#include "stereo_g.h"
/****************************************************************************/
// refinement type: 1 is uniform and 0 is newest vertex bisection
#ifndef UNIFORM_REFINEMENT
#define UNIFORM_REFINEMENT 1
#endif
/////////////////////////////////////////////////////////////////
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 10
#endif

#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 4
#endif

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
  //  dvector_print(stdout,&f);
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
  linear_itparam.linear_maxit    = 100;
  linear_itparam.linear_restart       = 100;
  linear_itparam.linear_print_level    = 10;
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);

  // adjust AMG parameters if needed
  // General AMG parameters
  amgparam.AMG_type             = UA_AMG;
  amgparam.print_level          = 3;
  amgparam.maxit                = 1;
  amgparam.max_levels           = 10;
  amgparam.coarse_dof           = 2000;
  amgparam.cycle_type           = W_CYCLE;
  amgparam.smoother             = SMOOTHER_GS;
  amgparam.presmooth_iter       = 1;
  amgparam.postsmooth_iter      = 1;
  amgparam.coarse_solver        = SOLVER_UMFPACK;
  //amgparam.relaxation           = 1.0;
  //amgparam.polynomial_degree    = 2;
  //amgparam.coarse_scaling       = ON;
  //amgparam.amli_degree          = 2;
  //amgparam.amli_coef            = NULL;
  //amgparam.nl_amli_krylov_type  = SOLVER_VFGMRES;

  // Aggregation AMG parameters
  amgparam.aggregation_type     = VMB;
  amgparam.strong_coupled       = 0.0;
  amgparam.max_aggregation      = 100;

  //amgparam.tentative_smooth     = 0.67;
  //amgparam.smooth_filter        = OFF;

  // print AMG paraemters
  //  param_amg_print(&amgparam);

  /* Actual solve */
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
  SHORT uniref=UNIFORM_REFINEMENT;
  INT ref_levels=REFINEMENT_LEVELS;
  INT dim=SPATIAL_DIMENSION;// 2d,3d,4d... example
  INT jlevel,k;
  scomplex *sc=NULL;
  switch(dim){
  case 5:
    sc=mesh5d();
    uniref=0;
    break;
  case 4:
    sc=mesh4d();    
    uniref=0;
    break;
  case 3:
    sc=mesh3d();    
    if( uniref ){
      for(jlevel=0;jlevel<=ref_levels;++jlevel){
	uniformrefine3d(sc);
	sc_vols(sc);
      }
    }
    break;
  default:
    sc=mesh2d();
    if( uniref ){
      for(jlevel=0;jlevel<=ref_levels;++jlevel){
	uniformrefine2d(sc);
	sc_vols(sc);
      }
    }
  }
  scomplex *sctop=NULL;
  ivector marked;
  marked.row=0;
  marked.val=NULL;
  dvector sol;
  // end intialization
  if(! (uniref) ){ 
    //    fprintf(stdout,"\nlevels=%d",ref_levels);fflush(stdout);
    for(jlevel=0;jlevel<ref_levels;jlevel++){
      /* choose the finest grid */
      sctop=scfinest(sc);
      // solve the FE
      sctop->vols=realloc(sctop->vols,sctop->ns*sizeof(REAL));
      sol=fe_sol(sctop,1.0,1.0);
      /* mark everything; or use an estimator */
      marked.row=sctop->ns;
      marked.val=realloc(marked.val,marked.row*sizeof(INT));
      /*mark everything*/
      for(k=0;k<marked.row;k++) marked.val[k]=TRUE;
      // now we refine.
      refine(1,sc,&marked);
      /* free */
      dvec_free(&sol);
      haz_scomplex_free(sctop);
      /*  MAKE sc to be the finest grid only */
    }
    scfinalize(sc);
    sc  ->vols=realloc(sc->vols,sc->ns*sizeof(REAL));
    sc_vols(sc);
  }
  //  icsr_print_matlab(stdout,sc->parent_v);
  //  haz_scomplex_print(sc,0,__FUNCTION__);
  sol=fe_sol(sc,1.0,1.0);
  // find the boundary simplicial complex:
  INT idsc;
  scomplex *dsc=sc_bndry(sc);
  if(dim==4 || dim ==3){
    // stereographic projection (only in 4D)
    idsc = stereo_g(dsc);
    /* write the output mesh file:    */
    //  hazw(grid_file,sc,0);
  } else{
    idsc=0;
  }
  /* WRITE THE OUTPUT vtu file for paraview: can be viewed with
     paraview */
  if(dim==2){
    vtkw("2d.vtu",sc,0,1.);
  } else if(dim==3){
    vtkw("3d.vtu",sc,0,1.);
    if(dsc) 
      vtkw("stereo3d.vtu",dsc,0,1.);
  } else if(dim==4){
    if(dsc)
      vtkw("stereo4d.vtu",dsc,0,1.);
  } else {
    fprintf(stdout,"\nNO PLOT: Dimension=%d is too large for plotting",dim);
  }
  fprintf(stdout,"\n\n");
  dvec_free(&sol);
  ivec_free(&marked);
  haz_scomplex_free(sc);  
  haz_scomplex_free(dsc);
  return 0;
}


