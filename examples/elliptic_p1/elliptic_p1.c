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
/****************************************************************************/

/* 
 * refinement type: .gt. 10 is uniform refinement and .le. 10
 *                  (typically 0) is the newest vertex bisection
*/
#ifndef REFINEMENT_TYPE
#define REFINEMENT_TYPE 11
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 6
#endif
/**/
#ifndef SPATIAL_DIMENSION
#define SPATIAL_DIMENSION 3
#endif
/**/
#ifndef SET_BNDRY_CODES
#define SET_BNDRY_CODES 1
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
  dvector sol;
  INT ndiv=(1<<ref_levels);
  fprintf(stdout,"ndiv=%ld\n",(long )ndiv);
  //  exit(55);
  sol.row=0;  sol.val=NULL;
  // Time the mesh generation and FE setup
  clock_t clk_mesh_start = clock();
  sc_all=mesh_cube_init(dim,ndiv,ref_type);
  sc=sc_all[0];
  /**/
  scfinalize(sc,(INT )1);  
  sc_vols(sc);
  clock_t clk_mesh_end = clock();
  fprintf(stdout,"\n%%%%%%CPUtime(mesh) = %.3f sec\n",
	  (REAL ) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  // Time the mesh generation and FE setup
  // Time the mesh generation and FE setup
  sol=fe_sol_no_dg(sc,1.0,1.0);
  //  short todraw=0;
  //  draw_grids(todraw, sc,&sol);
  /* write the output mesh file:    */
  /* hazw("output/mesh.haz",sc,0); */
  dvec_free(&sol);
  haz_scomplex_free(sc_all[0]);  
  free(sc_all);  
  clock_t clk_all_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(all) = %.3f sec\n",
	  (REAL ) (clk_all_end - clk_all_start)/CLOCKS_PER_SEC);
  return 0;
}
