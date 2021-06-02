/*! \file examples/basic_elliptic/amr_grids.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note This example highlights some of the features of the simple
 * mesh generator included with HAZmath. It is only to illustrate how
 * to use the mesh refinement.
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/
#include "solve_estimate_mark.h"
/*********************************************************************/
#include "connected_components.h"
/*********************************************************************/
// an interior vertex is defined as the vertex such that all faces attached to it are interior.
//
INT main(INT   argc,   char *argv[])
{
  INT i;
  FILE *fp;
  // fp=stdin;
  //  fp=HAZ_fopen("inputs/2d_2L.input","r");
  // fp=HAZ_fopen("inputs/2d_2L.input","r");
  fp=HAZ_fopen("inputs/4d_cube.input","r");
  /*
    PARSE THE INPUT.
  */
  input_grid *g=parse_input_grid(fp);
  input_grid_print(g);
  fclose(fp);
  /*
    GENERATE INITIAL GRID AND DECLARE VARIABLES.
  */
  scomplex *sc=generate_initial_grid(g);
  scomplex *sctop=NULL;
  INT ref_levels=g->nref, amr_marking_type=g->mark_type,j;
  dvector *solfem=NULL,*estimator=NULL;
  ivector *marked=NULL;
  void *all=NULL;
  REAL *xstar=NULL;
  INT nstar,dim=sc->n;
  if(amr_marking_type==0){
    // refine ref_levels;
    refine(ref_levels,sc,NULL);
  } else if(amr_marking_type==33){

    REAL h = 1.0/128;  // step distance of points
    REAL threshold = h; // threshold for close to the points or not

    /*
    INT nstep = 0;
    nstar= 2 + nstep*4; // refining near several points: (even number )
    xstar=(REAL *)calloc(nstar*dim,sizeof(REAL));
    xstar[0*dim+0]=1.666667e-1;
    xstar[0*dim+1]=6.666667e-1;
    xstar[1*dim+0]=8.333333e-1;
    xstar[1*dim+1]=6.666667e-1;
    for (i=1; i<nstep; i++){
      // points near (1/6, 2/3)
      xstar[(2+(i-1)*4)*dim+0]=1.666667e-1 + i*h;
      xstar[(2+(i-1)*4)*dim+1]=6.666667e-1;
      xstar[(3+(i-1)*4)*dim+0]=1.666667e-1 - i*h;
      xstar[(3+(i-1)*4)*dim+1]=6.666667e-1;

      // points near (5/6, 2/3)
      xstar[(4+(i-1)*4)*dim+0]=8.333333e-1 + i*h;
      xstar[(4+(i-1)*4)*dim+1]=6.666667e-1;
      xstar[(5+(i-1)*4)*dim+0]=8.333333e-1 - i*h;
      xstar[(5+(i-1)*4)*dim+1]=6.666667e-1;
    }*/

    nstar = g->num_refine_points;
    xstar = g->data_refine_points;

    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /* MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=mark_near_points(sctop,nstar,xstar, threshold);
      refine(1,sc,marked);
      /* free */
      haz_scomplex_free(sctop);
    }
    ivec_free(marked);
    //free(xstar);
  } else {
    /*
      Use "all" here can pass data around. Below we make 4 dvectors
      and one ivector, just as example. A good example for using the
      array will be to pass the whole hierarchy, not only the fine
      grid via all, e.g. all=(coid *)sc
    */
    all=(void *)malloc(5*sizeof(dvector)+sizeof(ivector));
    /**/
    for(j=0;j<ref_levels;j++){
      /*
       * SELECT the finest grid:
       */
      sctop=scfinest(sc);
      /*
       * SOLVE on the finest (for now) grid
       */
      solfem=(dvector *)exmpl_solve(sctop,all);
      /*
       * ESTIMATE using the numerical solution and the data stored in *all
       */
      estimator=(dvector *)exmpl_estimate(sctop,solfem,all);
      /*
       * MARK: marked is an ivector with num.rows=the number of
       *       simplices; its componenets are nonzero if the simplex
       *       is marked for refinement and 0 if it is not marked for
       *       refinement.
       */
      marked=exmpl_mark(sctop,estimator,all);
      /*
       *  Refine the grid. this always refines 1 time, but since we
       *  are in a loop, it will refine ref_levels times;
       */
      haz_scomplex_print(sctop,0,__FUNCTION__);
      refine(1,sc,marked);
      /* free */
      haz_scomplex_free(sctop);
      dvec_free(solfem);
      dvec_free(estimator);
      ivec_free(marked);
    }
    free(all);
  }
  /*  MAKE sc to be the finest grid only */
  //haz_scomplex_print(sc,0,"ZZZ");fflush(stdout);
  scfinalize(sc);
  haz_scomplex_print(sc,0,"XXX");fflush(stdout);
  //find_cc_bndry_cc(sc);
  /* write the output mesh file:    */
  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview:    */
  if(dim <4)
    vtkw(g->fvtu,sc,0,1.);
  /*FREE*/
  input_grid_free(g);
  haz_scomplex_free(sc);
  return 0;
}
