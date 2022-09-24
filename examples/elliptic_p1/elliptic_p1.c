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
#include "supporting_elliptic_p1.h"
/****************************************************************************/
/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 3
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 5
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
  clock_t clk_mesh_start = clock();
  fprintf(stdout,"\nMeshing in dimension=%lld ...",(long long )dim);
  sc_all=mesh_cube_init(dim,ndiv,ref_type);
  fprintf(stdout,"\nElements = Simplexes = %12lld;\nDoF      = Vertices  = %12lld\n",(long long )sc_all[0]->ns,(long long )sc_all[0]->nv); fflush(stdout);
  /**/
  sc=sc_all[0];
  scfinalize(sc,(INT )1);  
  sc_vols(sc);
  clock_t clk_mesh_end = clock();
  dCSRmat A=dcsr_create(0,0,0);
  dvector rhs, sol;
  rhs.row=0; rhs.val=NULL;
  sol.row=0; sol.val=NULL;
  // Time the mesh generation and FE setup
  /* fe_sol(sc,					\ */
  /* 	 &A,					\ */
  /* 	 &rhs,					\ */
  /* 	 &sol,					\ */
  /* 	 1e0,1e0); */
  fe_sol_no_dg(sc,					\
	 &A,					\
	 &rhs,					\
	 &sol,					\
	 1e0,1e0);
    //  short todraw=0;
    //  draw_grids(todraw, sc,&sol);
    /* write the output mesh file:    */
    /* hazw("output/mesh.haz",sc,0); */
  dvec_free(&sol);
  dvec_free(&rhs);
  dcsr_free(&A);
  haz_scomplex_free(sc_all[0]);  
  free(sc_all);  
  clock_t clk_all_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(mesh)     = %10.3f sec\n%%%%%%CPUtime(all)      = %10.3f sec\n",
	  (REAL ) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC,	\
	  (REAL ) (clk_all_end - clk_all_start)/CLOCKS_PER_SEC);
  return 0;
}
