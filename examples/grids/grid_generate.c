/*! \file examples/basic_elliptic/grid_generate.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2015/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simple grids formed by
 * macroelements. Each macroelement is assumed to be isomorphic to the
 * n-dimensional unit cube. 
 *
 * \note This example is just to generate grids for testing. In
 * general one can use any generator to generate such grids.
 *
 */
/*********** HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
/*********************************************************************/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
INT main(INT argc, char **argv)
{
  //  INT i=-1;
  FILE *fp=stdin;  /* OR:   fp=HAZ_fopen("grid.input","r"); */
  input_grid *g=parse_input_grid(fp);
  input_grid_print(g);
  scomplex *sc=generate_grid(g);
    fprintf(stdout,"Writing a vtk file...\n");
  vtkw("newmesh.vtu",sc,0,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  input_grid_free(g);  
  return 0;
}
/*********************EOF**********************************/
