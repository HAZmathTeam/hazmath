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
/**********************************************************************/
/* void input_grid_free(input_grid *g); */
/* void input_grid_print(input_grid *g); */
/* void coo2csr(INT nrow,INT ncol,INT nnz,					\ */
/* 	     INT *row_idx,INT *col_idx, void *aval,			\ */
/* 	     INT *ia,INT *ja, void *bval,				\ */
/* 	     size_t elsize); */
/* char **splits(char *s, const char *d, INT *num); */
/* void read_data(char *data_coordsystems,		\ */
/* 	       char *data_vertices,		\ */
/* 	       char *data_edges,		\ */
/* 	       input_grid *g); */
/* void get_out(char *pattern, size_t le); */
/* char *make_string_from_file(FILE *the_file, size_t *length_string); */
/* char *get_substring(char *pattern,		\ */
/* 		    size_t *length_substring,	\ */
/* 		    char *the_string); */
/* input_grid *parse_input_grid(const char *input_file_grid); */
/********************************************************************/
/***************************************************************/
int main(int argc, char **argv){
  char input_grid_file[256]={"grid.input"};
  input_grid *g=parse_input_grid(input_grid_file);
  input_grid_print(g);  
  input_grid_free(g);
  return 0;
}
