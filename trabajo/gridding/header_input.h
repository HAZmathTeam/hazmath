#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#define PI 3.14159265358979323851281e+00
#define EXP_OF_1  2.71828182845904523542817e+00

#ifndef REAL
#define REAL double
#endif

#ifndef INT
#define INT int
#endif

#ifndef SHORT
#define SHORT short
#endif



#ifndef FILENAMESIZE
#define FILENAMESIZE  1024
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE  0
#endif

#ifndef nula
#define nula  -1
#endif

#ifndef uint
#define uint unsigned int
#endif
/*================================================================*/
typedef struct {
  char *title; // the title of the input
  INT dim; // the dimension of the problem.
  //----------------
  // output flags
  //----------------
  SHORT print_level;   /**< print level */
  //----------------
  // files
  //----------------
  char *fgrid;  /* grid file name */
  char *dgrid;   /* output directory */
  char *fvtu;  /* grid file name */
  char *dvtu;   /* output directory */
  INT ncsys; /* number of coordinate systems */
  REAL *ox; /* origins of the coordinate systems */
  INT *systypes; /** types for the coord. system */
  INT *syslabels; /** labels for the coord. system */
  INT nv; /* number of vertices in the graph describing the
	     computational domain */
  REAL *xv; /* coordinates for each vertex [nv][dim]*/
  INT *csysv; /* coordinate system labels for vertices [nv]*/
  INT *bcodesv; /* boundary codes for vertices [nv]*/
  INT ne; /* number of edges/segments */ 
  REAL *xe; /* coordinates for each midpoint of an edge [ne][dim]*/
  INT *seg;/* segments array of size ne by 3. For every edge:
	      (v1,v2,divisions) with v1<v2 */
  INT nmacro;/*number of macroelements*/
  INT *macroel; /* macroelements: macroelement label, vertices forming
		   a macro element, macroelement material */ 
  INT nmacrofaces;  /*number of macroelement faces that are marked
		      with codes; boundary or internal it does not
		      matter */
  INT *macrofaces;   /* faces and boundary codes of faces */ 
}input_grid; /** Input GRID parameters */
/**********************************************************************/
void ilexsort(const INT nr, const INT nc,INT *a,INT *p);
/* void dlexsort(const INT nr, const INT nc,REAL *a,INT *p); */
/* void isi_sort(INT n, INT *a); */
/* void isi_sortp(const INT n, INT *a, INT *p, INT *invp); */
/* void dsi_sortp(const INT n, REAL *a, INT *p, INT *invp); */
