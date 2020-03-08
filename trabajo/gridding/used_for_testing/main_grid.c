#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
#include "grid_defs.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
scomplex *macro_split(input_grid *g0,cube2simp *c2s);
/****************************************************************/
INT main(INT argc, char **argv)
{
  INT i=-1;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  scomplex *sc=macro_split(g,c2s);
  if(dim==2||dim==3) {
    fprintf(stdout,"Writing vtk file...\n");
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  /*FREE*/
  haz_scomplex_free(sc);
  input_grid_free(g); 
  cube2simp_free(c2s);
  fprintf(stdout,"\nDone.\n");
  return 0;
}
/*********************EOF**********************************/
