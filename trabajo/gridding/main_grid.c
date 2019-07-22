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
//#include "grid_defs.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
INT main(INT argc, char **argv)
{
  //  INT i=-1;
  input_grid *g=parse_input_grid("grid.input");
  if(g->print_level>3)
    input_grid_print(g);
  scomplex *sc=generate_grid(g);
  fprintf(stdout,"Writing vtk file...\n");
  vtkw("newmesh.vtu",sc,0,0,1.);
  /*FREE*/
  haz_scomplex_free(sc);
  input_grid_free(g);  
  return 0;
}
/*********************EOF**********************************/
