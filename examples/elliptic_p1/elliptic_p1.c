/*! \file examples/elliptic_p1/elliptic_p1.c
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
#include "functions_elliptic_p1.h"
#include "supporting_elliptic_p1.h"
/****************************************************************************/
/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 4
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 1
#endif
/**/
#ifndef REFINEMENT_TYPE
#define REFINEMENT_TYPE 11
#endif
/****************************************************************************/
int main(int argc, char *argv[])
{
  // Overall CPU Timing
  clock_t clk_all_start = clock();
  //
  INT ref_levels=REFINEMENT_LEVELS;
  INT dim=SPATIAL_DIMENSION;
  INT ref_type=REFINEMENT_TYPE; // >10 uniform refinement;
  //
  INT jlevel,k;
  scomplex **sc_all=NULL,*sc=NULL,*sctop=NULL;
  /**/
  INT ndiv=(1<<ref_levels);
  //  fprintf(stdout,"ndiv=%lld\n",(long long )ndiv);
  // Time the mesh generation and FE setup
  fprintf(stdout,"\nMeshing in dimension=%lld ...",(long long )dim);
  clock_t clk_mesh_start = clock();
  sc_all=mesh_cube_init(dim,ndiv,ref_type);
  fprintf(stdout,"\nElements = Simplexes = %12lld;\nDoF      = Vertices  = %12lld\n",(long long )sc_all[0]->ns,(long long )sc_all[0]->nv); fflush(stdout);
  /**/
  sc=sc_all[0];
  scfinalize(sc,(INT )1);  
  sc_vols(sc);
  clock_t clk_mesh_end = clock();
  /* Assemble */
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  dCSRmat A=dcsr_create(0,0,0);
  dvector rhs, sol;
  rhs.row=0; rhs.val=NULL;
  sol.row=0; sol.val=NULL;
  // Time the mesh generation and FE setup
  REAL alpha=1e0,gamma=0e0;
  /* call_assembly_w_dg(sc,			\ */
  /* 	 &A,					\ */
  /* 	 &rhs,					\ */
  /* 	 alpha,gamma); */
  /* dcsr_write_dcoo("Adg.dat",&A); */
  /********************************************************/
  call_assembly(sc,				\
		&A,				\
		&rhs,				\
		alpha,gamma);
  //dcsr_write_dcoo("A.dat",&A);
  short todraw=1;
  draw_grids(todraw, sc,&sol);
  /* write the output mesh file:    */
  /* hazw("output/mesh.haz",sc,0); */
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  /*Solve*/
  clock_t clk_solver_start = clock(); // End of timing for mesh and FE setup  
  solveit(&A,&rhs,&sol);
  clock_t clk_solver_end = clock(); // End of timing for mesh and FE setup
  // err computation |u_I-u_h|_{1} and use rhs as a working space
  REAL err1,err0,uinterp;
  INT i,ij,j;
  // rhs= A*(err)
  for(i=0;i<sc->nv;++i){
    err0=0e0;
    for(ij=A.IA[i];ij<A.IA[i+1];++ij){
      j=A.JA[ij];
      exactsol(&uinterp,&sc->x[j*dim],0e0,dim,NULL);
      err0+=A.val[ij]*(sol.val[j]-uinterp);
    }
    rhs.val[i]=err0;
  }
  // err0=err'*A*err.
  err0=0e0;
  for(i=0;i<sc->nv;++i){
    exactsol(&uinterp,&sc->x[i*dim],0e0,dim,NULL);
    err0+=rhs.val[i]*(sol.val[i]-uinterp);
  }
  fprintf(stdout,"\n|u_I-u_h|_1=%.6e\n",sqrt(err0));
  haz_scomplex_free(sc_all[0]);
  free(sc_all);  
  dvec_free(&sol);
  dvec_free(&rhs);
  dcsr_free(&A);
  clock_t clk_all_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(mesh)     = %10.3f sec\n%%%%%%CPUtime(assembly) = %10.3f sec\n%%%%%%CPUtime(solver)   = %10.3f sec\n%%%%%%CPUtime(all)      = %10.3f sec\n", \
	  (REAL ) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC,	\
	  (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC, \
	  (REAL ) (clk_solver_end - clk_solver_start)/CLOCKS_PER_SEC,	\
	  (REAL ) (clk_all_end - clk_all_start)/CLOCKS_PER_SEC);
  return 0;
}
