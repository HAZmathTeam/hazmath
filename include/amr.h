//
//  amr.h
//
//
//  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#ifndef MAXFILENAMESIZE
#define MAXFILENAMESIZE 1024
#endif
#ifndef _amr_h
#define _amr_h
#endif
//#ifndef INPUT_GRID_DATA_
//#define INPUT_GRID_DATA_ "title{","file_grid{","file_vtu{","data_coordsystems{","data_vertices{","data_edges{","data_macroelements{","data_macrofaces{","dimension{","num_coordsystems{","num_vertices{","num_edges{","num_macroelements{","num_macrofaces{","num_refinements{","refinement_type{","amr_marking_type{","err_stop_refinement{","print_level{","num_refine_points{","data_refine_points{"
//#endif
#ifndef DEFAULT_GRID_DATA_
#define DEFAULT_GRID_DATA_ "title{Grid on a cube (-1,1)x(-1,1)x(-1,1);5x4x3 lattice}dimension{3}print_level{0}file_grid{mesh3d.haz}file_vtu{mesh3d.vtu} num_edges{3}data_edges{0 1 3  0 2 4 0 4 5}num_vertices{8} data_vertices{0 0 0. 0. 0. 1 0 0. 0. 1. 2 0 0. 1. 0. 3 0 0. 1. 1. 4 0 1. 0. 0. 5 0 1. 0. 1. 6 0 1. 1. 0. 7 0 1. 1. 1.}  num_macroelements{1}  data_macroelements{0 1 2 3 4 5 6 7 -1}num_macrofaces{6} data_macrofaces{0 1 2 3 1 0 4 1 5 1 4 7 5 6 1 2 6 3 7 1 0 4 2 6 1 1 5 7 3 1}num_coordsystems{1}data_coordsystems{0 0. 0. 0. 0}num_refinements{0}refinement_type{0}err_stop_refinement{-1.e-10}{amr_marking_type{0}\0"
#endif
/* defaults to unit cube in 3d and criss-cross grid 3x4x5 . */
//INT max_chars_input_grid_file=((1<<15) - 1); //maxcount=(1<<15-1);
/*******************************************************************/
typedef struct /* n-homogenous simplicial complex */
{
  INT nbig; /* the dimension of the space in which SC is embedded */
  INT n; /* the dimension of SC */
  INT nv; /* number of 0-dimensional simplices */
  INT ns; /* number of n-dimensional simplices */
  INT level; /* level of refinement */
  INT *marked; /*whether marked or not*/
  INT *gen; /* array to hold the simplex generation during refinement */
  INT *nbr; /* array (ns by n+1) to hold the neighbors */
  INT *parent; /*parent of each simplex*/
  INT *child0;
  INT *childn; /*children (0&n) of each simplex*/
  INT *bndry; /* nv boundary codes for vertices */
  INT *csys; /* nv coord system for every vertex point */
  INT *nodes; /*0-dim simplices opposite to (n-1) dimensional
		neighbors, i.e. the simplex2vertex incidence (ns by
		(n+1)) */
  iCSRmat *parent_v; /* for each vertex added in a refinement gives the
		       two vertices forming the edge where the
		       refinement vertex was put */
  INT *flags; /*flag of the simplex, e.g. for evaluation of
		   piece-wise defined functions*/
  REAL *x; /*(nv x n) array to hold the coordinates of vertices */
  REAL *vols; /* volumes of the simplices */
  REAL factorial; /*n factorial */
  INT *bcodesf; /* codes for boundary faces */
  INT *isbface; /* indicator if a face is on the boundary */
  iCSRmat *bfs; /* bfs structure for the simplices. works for more than one connected components */
  INT *etree; /* bfs tree (etree[k]=unique_ancestor_of_k in the BFS
		 tree) */
  INT ref_type; /*n connected components */
  INT cc; /*n connected components */
  INT bndry_cc; /*connected components on the boundary */
} scomplex;
/* /\*================================================================*\/ */
/* typedef struct /\* a macroelement (isomrphic to the hypercube */
/* 		  usually) *\/ */
/* { */
/*   INT type; /\* the type of the macro element: not used at the moment *\/ */
/*   INT *csys; /\* coordinate system for every vertex *\/ */
/*   REAL *xmac; /\* coordinates of the vertices in the corresponding */
/* 		 coordinate system ; *\/ */
/*   REAL *xemac; /\* coordinates of the midpoint of the edges if */
/* 		  needed. These are computed depending on the coord */
/* 		  system and are not an input *\/ */
/* } macroelement; */
/* /\*================================================================*\/ */
typedef struct /* a coordinate system */
{
  INT type; /* the type of the coordinate system: 0 is cartesian, 1 if
	       it is polar, 2 is cyllindical and so on */
  REAL *o; /* coordinates of the origin */
  scomplex *parent; /*parent complex */
} coordsystem;
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
  char *fvtu;  /* grid file name */
  INT ncsys; /* number of coordinate systems */
  REAL *ox; /* origins of the coordinate systems */
  INT *systypes; /** types for the coord. system */
  INT *syslabels; /** labels for the coord. system */
  INT nv; /* number of vertices in the graph describing the
	     computational domain */
  REAL *xv; /* coordinates for each vertex [nv][dim]*/
  INT *csysv; /* coordinate system labels for vertices [nv]*/
  INT *labelsv; /* coordinate system labels for vertices [nv]*/
  INT *bcodesv; /* boundary codes for vertices [nv]*/
  INT ne; /* number of edges/segments */
  REAL *xe; /* coordinates for each midpoint of an edge [ne][dim]*/
  INT *seg;/* segments array of size ne by 3. For every edge:
	      (v1,v2,divisions) with v1<v2 */
  INT nel;/*number of macroelements*/
  INT *mnodes; /* macroelements: macroelement label, vertices forming
		   a macro element, macroelement material */
  INT nf;  /*number of macroelement faces that are marked
	     with codes; boundary or internal it does not
	     matter */
  INT *mfaces;   /* faces and boundary codes of faces */
  INT nref;   /* number of refinements (AMR)*/

  INT ref_type;   /* INIT refinement type: If ref_type =
		     -2,-1,0,1,2,3,4, then newest vertex bisection in
		     any d.  the negative and positive choices
		     indicate where the edges should point in the
		     initial grid: left, right forth, back, etc */
  INT mark_type;   /* AMR marking type (0)uniform; nonzero: user defined.. */
  REAL err_stop;   /* stop tolerance for AMR */
  // ----- refine point ---- //
  INT  num_refine_points;
  REAL *data_refine_points;
}input_grid; /** Input GRID parameters */
/*************************************************************/
typedef struct /* n-homogenous simplicial SUBcomplex */
{
  INT nbig; /* spatial dimension in which the sub-SC is embedded*/
  INT n; /* dimension of subSC */
  INT ns; /* number of n-1 dimensional simplices forming the subcomplex */
  INT *elf; /* relation with the parent SC (element to face map)  */
  INT *nodes; /* simplex-vertex incidence (ns by n+1)) */
  INT *flags; /* flag for the the simplex, indicating boundary or not
		 and boundary code */
  REAL *normals; /* coordinates of the normal vectors to the simplices */
  REAL *areas; /* areas (could be called volumes) of the simplices from
		 subSC */
  scomplex *parent; /*parent complex */
} subscomplex;

typedef struct /* n-dimensional uniform grid */
{ INT n; /* spatial dimension of the grid */
  INT ugtype; /* type = 0 is uniform grid; type =1 dx and dy must be
	       arrays giving the distances between x_{i} and x_{i+1},
	       where x is the vector with coordinates. not implemented
	       yet for types .ne. 0*/
  INT nall; /*
	      the total number of VERTICES in the lattice (the
	      dimension of ug->data)
	    */
  INT *ndiv; /*number of divisions in each direction ndiv[dim].
	       NOTE: nall=(ndiv[dim-1]+1)*...*(ndiv[0]+1) */
  int nvcube; /* number of vertices on the unit cube in R^n=2^{n}.*/
  unsigned int *bits; /* the binary digits of all the integers from 0
		to 2^{n-1} as an array. These are also the coordinates
		of the vertices of the unit cube in R^n*/
  REAL *xo; /* coordinates of the origin xo[dim]*/
  REAL *xn; /* coordinates of the max corner(NE in 2D) xn[dim]*/
  REAL *dx; /*the step sizes in each dimension:
	      dx[i]=(xn[i]-xo[i])/ndiv[i], i=1:n*/
  REAL *data; /*data given over this uniform grid.*/
  INT dataid; /* integer to identify data id. if 0, data is ignored
		 (not freed at the end basicaly and not alocated) */
  FILE *fp; /* file from which the data is read */
} unigrid;


typedef struct /* structure to support splitting unit cube into simplices */
{ INT n; /* spatial dimension of the grid */
  INT nvcube; /* number of vertices on the unit cube in R^n=2^{n}.*/
  INT nvface; /* number of vertices on a face of the the unit cube in
		 R^n=2^{n-1}.*/
  INT ns; /* number of n dimensional simplices in the unit
	     cube(n_factorial of them) */
  INT ne; // number of edges in the cube.
  INT nf; // number of n-1 dimensional faces in the cube.
  unsigned int *bits; /* the binary digits of all the integers from 0
		to 2^{n-1} as an array. These are also the coordinates
		of the vertices of the unit cube in R^n*/
  INT *edges; /* the array containing the ends of edges of the unit cube. */
  INT *faces; /* the array containing faces (consistantly ordered)
		 unit cube. */
  INT *nodes; /* the array describing each of the n factorial simplices
		in the unit cube */
  INT *perms; /* the n by nvcube array describing the permutations
		which give different orderings of the cube's vertices
		so that we get a criss-cross grid in n
		dimensions. There are n such permutations reflecting
		the cube across the hyperplanes meeting at the vertex
		with coordinates [1,1,1,1,...,1] perms[nodes] will
		give consistent splitting of neighboring cubes.
	     */
} cube2simp;
/*=================================================================*/
typedef struct /* macroelement complex (isomorphic to
		  parallelepiped w 2^{n} vertices) */
{
  INT nel; /*number of macro elements*/
  INT nf; /*number of macro faces*/
  INT nfi; /*number of interior faces*/
  INT nfb; /*number of boundary faces*/
  INT **nd; /*number of divisions for every macro element*/
  INT **elneib; /*
		   element neighboring list where the position of the
		   neighbor is the same as the position of the face in
		   cube2simp->faces shared by the two el.
		*/
  INT **el2fnum; /* for every element this gives the local to global
		    face number map;   */
  INT **iindex; /* used to remove repeated vertices */
  iCSRmat *fullel2el; /* full element to element which has also the
			 number of common vertices as entries in the matrix; */
  INT *bcodesf; /* codes for boundary faces */
  INT *isbface; /* indicator if a face is on the boundary */
  INT *flags; /*flags (material) of the macro elements*/
  INT cc; /*n connected components */
  INT bndry_cc; /*connected components on the boundary */
  iCSRmat *bfs; /* bfs levels structure for the mesh el2el only if
		   they share a face; */
  INT *etree; /* bfs tree (el by el only if the elems share a face) */
} macrocomplex;
/*********************************************************************/
/*==================================================================*/
typedef struct /* features (to refine around these) */
{
  INT nbig; /* dimension in which this is "embedded", i.e. one
		 coordinate can be left free. for example if features
		 are given by (x,y) coordinates in 3D, then this will
		 be 3. and the z coordinate will take value
		 equal to features->fill" */
  INT n; /* real dimension, i.e. if features are given by (x,y)
	      coords in 3d, this will be 2*/
  INT nf; /* number of points in the feature */
  REAL *x; /* x[j*dim+i],j=0:nfeatures-1, i=0:dim. last few
	      coordinates are filled with the value in fill */
  REAL fill; /*
		mock value for the remaining coordinates from n to
		nbig
	     */
  FILE *fpf;
} features;
/* END OF STRUCTURE DEFINITIONS */
