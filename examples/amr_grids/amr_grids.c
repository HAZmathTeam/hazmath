/*! \file examples/basic_elliptic/basic_elliptic.c
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
INT main(INT   argc,   char *argv[])
{
  FILE *fp=stdin;     
  //  fp=HAZ_fopen("2d_square.input","r");
  /*
    PARSE THE INPUT.
  */
  input_grid *g=parse_input_grid(fp);
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
  if(amr_marking_type){
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
       * SOLVE on the finest (for now grid)
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
      refine(1,sc,marked); 
      /* free */
      haz_scomplex_free(sctop);
      dvec_free(solfem);
      dvec_free(estimator);
      ivec_free(marked);
    }
  } else {
    // refine ref_levels;
    refine(ref_levels,sc,NULL);
  }
  /*  MAKE sc to be the finest grid only */
  scfinalize(sc);
  /* write the output mesh file:    */
  hazw(g->fgrid,sc,0);
  /* WRITE THE OUTPUT vtu file for paraview:    */
  vtkw(g->fvtu,sc,0,1.);
  /*FREE*/
  input_grid_free(g);
  free(all);
  haz_scomplex_free(sc);
  return 0;
}
