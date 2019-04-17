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
#include <math.h>
#include <limits.h>
#include <time.h>

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
  INT *csys; /* nv coord system for a point */
  INT *nodes; /*0-dim simplices opposite to (n-1) dimensional
		neighbors, i.e. the simplex-vertex incidence (ns by
		(n+1)) */
  INT *flags; /*flag of the simplex, e.g. for evaluation of
		   piece-wise defined functions*/  
  REAL *x; /*(nv x n) array to hold the coordinates of vertices */
  REAL *vols; /* volumes of the simplices */
  REAL *fval; /* function values at vertices (could be ANY, but is
		 used for elevation in hydrolena) */
  REAL factorial; /*n factorial */
} scomplex;

typedef struct /* a coordinate system */
{
  INT type; /* the type of the coordinate system: 0 is cartesian, 1 if
	       it is polar, 2 is cyllindical and so on */
  REAL *o; /* coordinates of the origin */
  scomplex *parent; /*parent complex */
} coordsystem;

typedef struct /* a macroelement (isomrphic to the hypercube
		  usually) */
{
  INT type; /* the type of the macro element: not used at the moment */
  INT *csys; /* coordinate system for every vertex */
  REAL *xmac; /* coordinates of the vertices in the corresponding
		 coordinate system ; */
  REAL *xemac; /* coordinates of the midpoint of the edges if
		  needed. These are computed depending on the coord
		  system and are not an input */
} macroelement;

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
  INT ns; /* number of n dimensional simplices in the unit
	     cube(n_factorial of them) */
  INT ne; // number of edges in the cube.
  INT nf; // number of n-1 dimensional faces in the cube.
  unsigned int *bits; /* the binary digits of all the integers from 0
		to 2^{n-1} as an array. These are also the coordinates
		of the vertices of the unit cube in R^n*/
  INT *edges; /* the array containing the ends of edges of the unit cube. */
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

