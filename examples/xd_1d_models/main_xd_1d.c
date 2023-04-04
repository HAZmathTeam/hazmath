/*! \file examples/amr_grids/xd_1d.c
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
/**********************************************************************/
#ifndef MAX_NODES_PER_SIMPLEX
#define MAX_NODES_PER_SIMPLEX  1
#endif
/**/
#ifndef OUTER_SPATIAL_DIMENSION
#define OUTER_SPATIAL_DIMENSION 3
#endif
/**/
#ifndef REFINEMENT_LEVELS
#define REFINEMENT_LEVELS 150
#endif

/****************************************************************************************8*/
#include "definitions_xd_1d.h"
#include "supporting_xd_1d.h"
void xd_1d_lib(const INT dimbig, const INT max_nodes_in, const INT ref_levels_in, \
	       const char *idir, const char *odir);
/****************************************************************************************8*/
INT main(INT argc, char *argv[]) {
  INT dimbig=OUTER_SPATIAL_DIMENSION;
  char *idir = str_add_dim(dimbig,"./input/1d_nets_","d/");
  char *odir = strdup("./output/");
  xd_1d_lib(dimbig, MAX_NODES_PER_SIMPLEX, REFINEMENT_LEVELS,idir, odir);  
  return 0;
}
