/*
 *  basis.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 2/1/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
 *
 */

// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"

/* Compute the basis functions for triangles or tetrahedra
 *
 *  Discussion:  Typically this involves DOF defined on either the 
 *  vertices, edges, or faces.  In most cases, the basis elements are
 *  defined using the standard Lagrange finite-element basis functions
 *  using barycentric coordinates.  An illustration of the 2D triangle is
 *  here.
 *
 *  The physical element:
 *
 *    In this picture, we don't mean to suggest that the bottom of
 *    the physical triangle is horizontal.  However, we do assume that
 *    each of the sides is a straight line, and that the intermediate
 *    points are exactly halfway on each side.
 *
 *    |
 *    |
 *    |        3
 *    |       / \
 *    |      /   \
 *    Y    e31    e23
 *    |    /       \
 *    |   /         \
 *    |  1----e12-----2
 *    |
 *    +--------X-------->
 *
 *  Reference element T3:
 *
 *    In this picture of the reference element, we really do assume
 *    that one side is vertical, one horizontal, of length 1.
 *
 *    |
 *    |
 *    1  3
 *    |  |\
 *    |  | \
 *    S e31 e23
 *    |  |   \
 *    |  |    \
 *    0  1-e12-2
 *    |
 *    +--0--R--1-------->
 *
 *  Determine the (R,S) coordinates corresponding to (X,Y).
 *
 *  What is happening here is that we are solving the linear system:
 *
 *    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
 *    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
 *
 *  by computing the inverse of the coefficient matrix and multiplying
 *  it by the right hand side to get R and S.
 *
 *  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
 *  for R and S.
 *  
 *  For quadratic elements:
 *    |
 *    1  3
 *    |  |\
 *    |  | \
 *    S  5  6 
 *    |  |   \
 *    |  |    \
 *    0  1--4--2
 *    |
 *    +--0--R--1-------->
 */

/****************************************************************************************************************************/
/* Compute Standard Lagrange Finite Element Basis Functions (PX) at a particular point in 2 or 3D*/
/* For now, we only assume Linears or Quadratic Elements (P1 or P2) */
void PX_H1_basis(REAL *p,REAL *dpx,REAL *dpy,REAL *dpz,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh mesh) 
{
  /*
   *    INPUT:
   *          x,y,z               Coordinates on physical triangle where to compute basis
   *          mesh                Mesh Data
   *	      dof                 DOF for the given element (in this case vertices and their global numbering)
   *          porder              Order of elements (1 or 2 for now)
   *          cv                  Coordinates of Vertices
   *          v_per_elm           Number of vertices per element
   *          dim                 Dimension of Problem
   *    OUTPUT:
   *          p(v_per_elm)	  Basis functions at particular point (1 for each vertex)
   *          dpx,dpy(3)	  Derivatives of basis functions at each vertex evaluated at given point
   */
	

    
  REAL dp1r,dp2r,dp3r,dp4r,dp5r,dp6r,dp7r,dp8r,dp9r,dp10r;
  REAL dp1s,dp2s,dp3s,dp4s,dp5s,dp6s,dp7s,dp8s,dp9s,dp10s;
  REAL dp1t,dp2t,dp3t,dp4t,dp5t,dp6t,dp7t,dp8t,dp9t,dp10t;
  REAL onemrst;
  INT i;
    
    printf("hello-basis0\n");

    
  // Get Mesh Data
  INT v_per_elm = mesh.v_per_elm;
  INT dim = mesh.dim;
    
    printf("hello-basis1\n");
    
    printf("v_per_elm = %d\n", v_per_elm);
    printf("dim = %d\n", dim);


  REAL* xp = (REAL *) calloc(v_per_elm,sizeof(REAL));
    
    printf("hello-basis11\n");

    
  REAL* yp = (REAL *) calloc(v_per_elm,sizeof(REAL));
    
    printf("hello-basis12\n");

    
  REAL* zp = NULL;
    
    printf("hello-basis2\n");


  coordinates* cv = mesh.cv;


  
  // 2D and 3D is slightly different
  if(dim==2) {
  
    // Get Physical Coordinates of Vertices
    for (i=0; i<v_per_elm; i++) {
      xp[i] = cv->x[dof[i]-1];
      yp[i] = cv->y[dof[i]-1];
    }
		
    // Get coordinates on reference triangle
    REAL det = (xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0]);
	
    REAL r = ((yp[2]-yp[0])*(x-xp[0]) + (xp[0]-xp[2])*(y-yp[0]))/det;
	
    REAL drdx = (yp[2]-yp[0])/det;
    REAL drdy = (xp[0]-xp[2])/det;
	
    REAL s = ((yp[0]-yp[1])*(x-xp[0]) + (xp[1]-xp[0])*(y-yp[0]))/det;
	
    REAL dsdx = (yp[0]-yp[1])/det;
    REAL dsdy = (xp[1]-xp[0])/det;
	
    /*  Get the basis functions for linear elements on each node.
     *  The basis functions can now be evaluated in terms of the
     *  reference coordinates R and S.  It's also easy to determine
     *  the values of the derivatives with respect to R and S.
     */
	
    onemrst = 1 - r - s;
    if(porder==1) {
      p[0] = onemrst;
      p[1] = r;
      p[2] = s;
      dp1r = -1;
      dp2r = 1;
      dp3r = 0;
      dp1s = -1;
      dp2s = 0;
      dp3s = 1;
    } else if(porder==2) {
      p[0] = 2*onemrst*(onemrst-0.5);
      p[1] = 2*r*(r-0.5);
      p[2] = 2*s*(s-0.5);
      p[3] = 4*r*onemrst;
      p[4] = 4*s*onemrst;
      p[5] = 4*r*s;	
      dp1r = 4*r+4*s-3;
      dp2r = 4*r - 1;
      dp3r = 0;
      dp4r = 4-8*r-4*s;
      dp5r = -4*s;
      dp6r = 4*s;
      dp1s = dp1r;
      dp2s = 0;
      dp3s = 4*s-1;
      dp4s = -4*r;
      dp5s = 4-4*r-8*s;
      dp6s = 4*r;
    } else {
      printf("You are attempting P%d elements.  Code can only do P1 or P2 at the moment.",porder);
      exit(0);
    }
	
    /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
     *  to (X,Y) using the chain rule.
     */
	
    dpx[0] = dp1r * drdx + dp1s * dsdx;
    dpy[0] = dp1r * drdy + dp1s * dsdy;
    dpx[1] = dp2r * drdx + dp2s * dsdx;
    dpy[1] = dp2r * drdy + dp2s * dsdy;
    dpx[2] = dp3r * drdx + dp3s * dsdx;
    dpy[2] = dp3r * drdy + dp3s * dsdy;
    if(porder==2) {
      dpx[3] = dp4r * drdx + dp4s * dsdx;
      dpy[3] = dp4r * drdy + dp4s * dsdy;
      dpx[4] = dp5r * drdx + dp5s * dsdx;
      dpy[4] = dp5r * drdy + dp5s * dsdy;
      dpx[5] = dp6r * drdx + dp6s * dsdx;
      dpy[5] = dp6r * drdy + dp6s * dsdy;
    }
	
  } else if (dim==3) {
    zp = (REAL *) calloc(v_per_elm,sizeof(REAL));
    // Get Nodes and Physical Coordinates
    for (i=0; i<v_per_elm; i++) {
      xp[i] = cv->x[dof[i]-1];
      yp[i] = cv->y[dof[i]-1];
      zp[i] = cv->z[dof[i]-1];
    }
		
    // Get coordinates on reference triangle
    REAL det = (xp[3]-xp[0])*((yp[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(zp[1]-zp[0])) \
      - (xp[2]-xp[0])*((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0])) \
      + (xp[1]-xp[0])*((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0]));
	
    REAL r = ((xp[3]-xp[0])*((y-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(z-zp[0])) \
	      - (xp[2]-xp[0])*((y-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(z-zp[0])) \
	      + (x-xp[0])*((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0])))/det;
	
    REAL drdx = ((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0]))/det;
    REAL drdy = ((xp[3]-xp[0])*(zp[2]-zp[0]) - (xp[2]-xp[0])*(zp[3]-zp[0]))/det;
    REAL drdz = ((xp[2]-xp[0])*(yp[3]-yp[0]) - (xp[3]-xp[0])*(yp[2]-yp[0]))/det;
	
    REAL s = ((xp[3]-xp[0])*((yp[1]-yp[0])*(z-zp[0])-(y-yp[0])*(zp[1]-zp[0])) \
	      - (x-xp[0])*((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0])) \
	      + (xp[1]-xp[0])*((y-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(z-zp[0])))/det;
	
    REAL dsdx = -((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0]))/det;
    REAL dsdy = ((xp[1]-xp[0])*(zp[3]-zp[0]) - (xp[3]-xp[0])*(zp[1]-zp[0]))/det;
    REAL dsdz = ((xp[3]-xp[0])*(yp[1]-yp[0]) - (xp[1]-xp[0])*(yp[3]-yp[0]))/det;
	
    REAL t = ((x-xp[0])*((yp[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(zp[1]-zp[0])) \
	      - (xp[2]-xp[0])*((yp[1]-yp[0])*(z-zp[0])-(y-yp[0])*(zp[1]-zp[0])) \
	      + (xp[1]-xp[0])*((yp[2]-yp[0])*(z-zp[0])-(y-yp[0])*(zp[2]-zp[0])))/det;
	
    REAL dtdx = ((yp[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(zp[1]-zp[0]))/det;
    REAL dtdy = ((xp[2]-xp[0])*(zp[1]-zp[0]) - (xp[1]-xp[0])*(zp[2]-zp[0]))/det;
    REAL dtdz = ((xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0]))/det;
	
    /*  Get the basis functions for linear elements on each node.
     *  The basis functions can now be evaluated in terms of the
     *  reference coordinates R and S.  It's also easy to determine
     *  the values of the derivatives with respect to R and S.
     */
	
    onemrst = 1 - r - s - t;
    if(porder==1) {
      p[0] = onemrst;
      p[1] = r;
      p[2] = s;
      p[3] = t;
      dp1r = -1;
      dp2r = 1;
      dp3r = 0;
      dp4r = 0;
      dp1s = -1;
      dp2s = 0;
      dp3s = 1;
      dp4s = 0;
      dp1t = -1;
      dp2t = 0;
      dp3t = 0;
      dp4t = 1;
    } else if(porder==2) {
      p[0] = onemrst*(1 - 2*r - 2*s - 2*t);
      p[1] = 2*r*(r-0.5);
      p[2] = 2*s*(s-0.5);
      p[3] = 2*t*(t-0.5);
      p[4] = 4*r*onemrst;
      p[5] = 4*s*onemrst;
      p[6] = 4*t*onemrst;
      p[7] = 4*r*s;
      p[8] = 4*r*t;
      p[9] = 4*s*t;
	
      dp1r = 4*r+4*s+4*t-3;
      dp2r = 4*r-1;
      dp3r = 0;
      dp4r = 0;
      dp5r = 4*onemrst - 4*r;
      dp6r = -4*s;
      dp7r = -4*t;
      dp8r = 4*s;
      dp9r = 4*t;
      dp10r = 0;
      dp1s = 4*r+4*s+4*t-3;
      dp2s = 0;
      dp3s = 4*s-1;
      dp4s = 0;
      dp5s = -4*r;
      dp6s = 4*onemrst - 4*s;
      dp7s = -4*t;
      dp8s = 4*r;
      dp9s = 0;
      dp10s = 4*t;
      dp1t = 4*r+4*s+4*t-3;
      dp2t = 0;
      dp3t = 0;
      dp4t = 4*t-1;
      dp5t = -4*r;
      dp6t = -4*s;
      dp7t = 4*onemrst - 4*t;
      dp8t = 0;
      dp9t = 4*r;
      dp10t = 4*s;
    } else {
      printf("You are attempting P%d elements.  Code can only do P1 or P2 at the moment.",porder);
      exit(0);
    }
	
    /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
     *  to (X,Y) using the chain rule.
     */
	
    dpx[0] = dp1r * drdx + dp1s * dsdx + dp1t * dtdx;
    dpy[0] = dp1r * drdy + dp1s * dsdy + dp1t * dtdy;
    dpz[0] = dp1r * drdz + dp1s * dsdz + dp1t * dtdz;
    dpx[1] = dp2r * drdx + dp2s * dsdx + dp2t * dtdx;
    dpy[1] = dp2r * drdy + dp2s * dsdy + dp2t * dtdy;
    dpz[1] = dp2r * drdz + dp2s * dsdz + dp2t * dtdz;
    dpx[2] = dp3r * drdx + dp3s * dsdx + dp3t * dtdx;
    dpy[2] = dp3r * drdy + dp3s * dsdy + dp3t * dtdy;
    dpz[2] = dp3r * drdz + dp3s * dsdz + dp3t * dtdz;
    dpx[3] = dp4r * drdx + dp4s * dsdx + dp4t * dtdx;
    dpy[3] = dp4r * drdy + dp4s * dsdy + dp4t * dtdy;
    dpz[3] = dp4r * drdz + dp4s * dsdz + dp4t * dtdz;
    if(porder==2) {
      dpx[4] = dp5r * drdx + dp5s * dsdx + dp5t * dtdx;
      dpy[4] = dp5r * drdy + dp5s * dsdy + dp5t * dtdy;
      dpz[4] = dp5r * drdz + dp5s * dsdz + dp5t * dtdz;
      dpx[5] = dp6r * drdx + dp6s * dsdx + dp6t * dtdx;
      dpy[5] = dp6r * drdy + dp6s * dsdy + dp6t * dtdy;
      dpz[5] = dp6r * drdz + dp6s * dsdz + dp6t * dtdz;
      dpx[6] = dp7r * drdx + dp7s * dsdx + dp7t * dtdx;
      dpy[6] = dp7r * drdy + dp7s * dsdy + dp7t * dtdy;
      dpz[6] = dp7r * drdz + dp7s * dsdz + dp7t * dtdz;
      dpx[7] = dp8r * drdx + dp8s * dsdx + dp8t * dtdx;
      dpy[7] = dp8r * drdy + dp8s * dsdy + dp8t * dtdy;
      dpz[7] = dp8r * drdz + dp8s * dsdz + dp8t * dtdz;
      dpx[8] = dp9r * drdx + dp9s * dsdx + dp9t * dtdx;
      dpy[8] = dp9r * drdy + dp9s * dsdy + dp9t * dtdy;
      dpz[8] = dp9r * drdz + dp9s * dsdz + dp9t * dtdz;
      dpx[9] = dp10r * drdx + dp10s * dsdx + dp10t * dtdx;
      dpy[9] = dp10r * drdy + dp10s * dsdy + dp10t * dtdy;
      dpz[9] = dp10r * drdz + dp10s * dsdz + dp10t * dtdz;
    }
  } else {
    baddimension();
  }
    
    printf("hello-basis3\n");
	
  if(xp) free(xp);
  if(yp) free(yp);
  if(zp) free(zp);
	
  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/* Compute Standard Quadratic Finite Element Basis Functions (P2) at a particular point */
/* Also compute the 2nd derivatives for some reason... */
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh mesh) 
{
	
  /*
   *    INPUT:
   *          x,y,z               Coordinates on physical triangle where to compute basis
   *          mesh                Mesh Data
   *	      dof                 DOF for the given element (in this case vertices and their global numbering)
   *          porder              Order of elements (1 or 2 for now)
   *    OUTPUT:
   *          p(v_per_elm)	  Basis functions at particular point (1 for each vertex)
   *          dpx,dpy(3)	  Derivatives of basis functions at each vertex evaluated at given point
   *          dpxx,dpyy,dpxy      2nd Derivatives of basis functions
   */
	
	
  REAL dp1r,dp2r,dp3r,dp4r,dp5r,dp6r;
  REAL dp1s,dp2s,dp3s,dp4s,dp5s,dp6s;
  REAL onemrs;
  INT i;
  INT v_per_elm = mesh.v_per_elm;
  REAL* xp = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* yp = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL dp1rr,dp2rr,dp3rr,dp4rr,dp5rr,dp6rr;
  REAL dp1ss,dp2ss,dp3ss,dp4ss,dp5ss,dp6ss;
  REAL dp1rs,dp2rs,dp3rs,dp4rs,dp5rs,dp6rs;

  coordinates* cv = mesh.cv;
	
  // Get Nodes and Physical Coordinates of vertices only
  for (i=0; i<3; i++) {
    xp[i] = cv->x[dof[i]-1];
    yp[i] = cv->y[dof[i]-1];
  }
		
  // Get coordinates on reference triangle
  REAL det = (xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0]);
	
  REAL r = ((yp[2]-yp[0])*(x-xp[0]) + (xp[0]-xp[2])*(y-yp[0]))/det;
	
  REAL drdx = (yp[2]-yp[0])/det;
  REAL drdy = (xp[0]-xp[2])/det;
	
  REAL s = ((yp[0]-yp[1])*(x-xp[0]) + (xp[1]-xp[0])*(y-yp[0]))/det;
	
  REAL dsdx = (yp[0]-yp[1])/det;
  REAL dsdy = (xp[1]-xp[0])/det;
	
  /*  Get the basis functions for quadratic elements on each node and edge.
   *  The basis functions can now be evaluated in terms of the
   *  reference coordinates R and S.  It's also easy to determine
   *  the values of the derivatives with respect to R and S.
   */
	
  onemrs = 1 - r - s;
  p[0] = 2*onemrs*(onemrs-0.5);
  p[1] = 2*r*(r-0.5);
  p[2] = 2*s*(s-0.5);
  p[3] = 4*r*onemrs;
  p[4] = 4*s*onemrs;
  p[5] = 4*r*s;	
  dp1r = 4*r+4*s-3;
  dp1rr = 4.0;
  dp1rs = 4.0;
  dp2r = 4*r - 1;
  dp2rr = 4.0;
  dp2rs = 0.0;
  dp3r = 0.0;
  dp3rr = 0.0;
  dp3rs = 0.0;
  dp4r = 4.0-8*r-4*s;
  dp4rr = -8.0;
  dp4rs = -4.0;
  dp5r = -4.0*s;
  dp5rr = 0.0;
  dp5rs = -4.0;
  dp6r = 4.0*s;
  dp6rr = 0.0;
  dp6rs = 4.0;
  dp1s = dp1r;
  dp1ss = 4.0;
  dp2s = 0.0;
  dp2ss = 0.0;
  dp3s = 4.0*s-1.0;
  dp3ss = 4.0;
  dp4s = -4.0*r;
  dp4ss = 0;
  dp5s = 4.0-4*r-8*s;
  dp5ss = -8.0;
  dp6s = 4.0*r;
  dp6ss = 0.0;
	
	
  /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
   *  to (X,Y) using the chain rule.
   */
	
  dpx[0] = dp1r * drdx + dp1s * dsdx;
  dpy[0] = dp1r * drdy + dp1s * dsdy;
  dpx[1] = dp2r * drdx + dp2s * dsdx;
  dpy[1] = dp2r * drdy + dp2s * dsdy;
  dpx[2] = dp3r * drdx + dp3s * dsdx;
  dpy[2] = dp3r * drdy + dp3s * dsdy;
  dpx[3] = dp4r * drdx + dp4s * dsdx;
  dpy[3] = dp4r * drdy + dp4s * dsdy;
  dpx[4] = dp5r * drdx + dp5s * dsdx;
  dpy[4] = dp5r * drdy + dp5s * dsdy;
  dpx[5] = dp6r * drdx + dp6s * dsdx;
  dpy[5] = dp6r * drdy + dp6s * dsdy;
  dpxx[0] = dp1rr*drdx*drdx + 2*dp1rs*drdx*dsdx + dp1ss*dsdx*dsdx;
  dpxx[1] = dp2rr*drdx*drdx + 2*dp2rs*drdx*dsdx + dp2ss*dsdx*dsdx;
  dpxx[2] = dp3rr*drdx*drdx + 2*dp3rs*drdx*dsdx + dp3ss*dsdx*dsdx;
  dpxx[3] = dp4rr*drdx*drdx + 2*dp4rs*drdx*dsdx + dp4ss*dsdx*dsdx;
  dpxx[4] = dp5rr*drdx*drdx + 2*dp5rs*drdx*dsdx + dp5ss*dsdx*dsdx;
  dpxx[5] = dp6rr*drdx*drdx + 2*dp6rs*drdx*dsdx + dp6ss*dsdx*dsdx;
  
  dpyy[0] = dp1rr*drdy*drdy + 2*dp1rs*drdy*dsdy + dp1ss*dsdy*dsdy;
  dpyy[1] = dp2rr*drdy*drdy + 2*dp2rs*drdy*dsdy + dp2ss*dsdy*dsdy;
  dpyy[2] = dp3rr*drdy*drdy + 2*dp3rs*drdy*dsdy + dp3ss*dsdy*dsdy;
  dpyy[3] = dp4rr*drdy*drdy + 2*dp4rs*drdy*dsdy + dp4ss*dsdy*dsdy;
  dpyy[4] = dp5rr*drdy*drdy + 2*dp5rs*drdy*dsdy + dp5ss*dsdy*dsdy;
  dpyy[5] = dp6rr*drdy*drdy + 2*dp6rs*drdy*dsdy + dp6ss*dsdy*dsdy;
  
  dpxy[0] = dp1rr*drdy*drdx + dp1rs*drdy*dsdx + dp1rs*drdx*dsdy + dp1ss*dsdx*dsdy;
  dpxy[1] = dp2rr*drdy*drdx + dp2rs*drdy*dsdx + dp2rs*drdx*dsdy + dp2ss*dsdx*dsdy;
  dpxy[2] = dp3rr*drdy*drdx + dp3rs*drdy*dsdx + dp3rs*drdx*dsdy + dp3ss*dsdx*dsdy;
  dpxy[3] = dp4rr*drdy*drdx + dp4rs*drdy*dsdx + dp4rs*drdx*dsdy + dp4ss*dsdx*dsdy;
  dpxy[4] = dp5rr*drdy*drdx + dp5rs*drdy*dsdx + dp5rs*drdx*dsdy + dp5ss*dsdx*dsdy;
  dpxy[5] = dp6rr*drdy*drdx + dp6rs*drdy*dsdx + dp6rs*drdx*dsdy + dp6ss*dsdx*dsdy;
	
  if(xp) free(xp);
  if(yp) free(yp);
	
  return;
}
/****************************************************************************************************************************/

/* /\****************************************************************************************************************************\/ */
/* void ned_basis(REAL *phi,REAL *cphi,REAL x,REAL y,INT *jel_n,INT *jel_ed,INT *ied_n,INT *jed_n,REAL *xn,REAL *yn,INT element_order,INT mydim)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Build basis functions *************************************************\/ */
	
/*   /\* Compute the Nedelec Elements */
/*    *    INPUT: */
/*    *          x,y              Coordinates on physical triangle where to compute basis */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *			jel_ed			 Edges belonging to particular element */
/*    *			ied_n,jed_n		 Edge to Node Map */
/*    *			element_order	 Number of nodes per element */
/*    *			mydim			 Dimension of Problem */
/*    *          xn,yn            Coordinates of nodes */
/*    *    OUTPUT: */
/*    *          phi(3,2)         Basis functions at particular point (2 for */
/*    *                           each edge, 6 total) (from reference triangle) */
/*    *          cphi(3)          2D Curl of Basis functions at particular point  */
/*    *                           (1 for each edge) (from reference triangle) */
/*    *\/ */
	
/*   INT i,ica,n1,n2,z,ihi,ilo; */
/*   INT mark1 = -1; */
/*   INT mark2 = -1; */
/*   REAL* p; */
/*   REAL* dpx; */
/*   REAL* dpy; */
/*   //REAL p[element_order],dpx[element_order],dpy[element_order]; */
/*   INT* ip = calloc(element_order,sizeof(INT)); */
/*   //INT ip[element_order]; */
/*   INT nedge = 3*(mydim-1); */
/*   INT ie[nedge]; */
	
/*   /\* Get Linear Basis Functions for particular element *\/ */
/*   p = calloc(element_order,sizeof(REAL)); */
/*   dpx = calloc(element_order,sizeof(REAL)); */
/*   dpy = calloc(element_order,sizeof(REAL)); */
/*   lin_tri_2D(p,dpx,dpy,x,y,jel_n,xn,yn,element_order); */
	
	
/*   REAL elen; */
/*   REAL* tau_edge = calloc(mydim,sizeof(REAL)); */
/*   REAL* midpt = calloc(mydim,sizeof(REAL)); */
/*   //REAL tau_edge[mydim],midpt[mydim]; */
	
/*   // Get Nodes correspond to linear basis elements */
/*   for (i=0; i<element_order; i++) { */
/*     ip[i] = jel_n[i]; */
/*   } */
	
/* #define max(a,b) ((a) > (b) ? (a) : (b)) */
	
	
/*   /\* Now, with the linear basis functions p, dpx, and dpy, the nedelec elements are  */
/*    * phi_eij = |eij|*(p(i)grad(p(j)) - p(j)grad(p(i))) */
/*    * |eij| = sqrt(|xj-xi|^2)   */
/*    *\/ */
/*   //INT* ie = calloc((mydim-1)*3,sizeof(INT)); */
	
/*   // Go through each edge and get length and find the corresponding nodes */
/*   for (i=0; i<nedge; i++) { */
/*     ie[i] = jel_ed[i]; */
/*     ica = ied_n[ie[i]-1]; */
/*     n1 = jed_n[ica-1]; */
/*     n2 = jed_n[ica]; */
/*     edge_stats(&elen,tau_edge,midpt,ie[i],xn,yn,NULL,ied_n,jed_n,mydim); */
		
/*     // Find out which linear basis elements line up with nodes on this edge */
/*     for (z=0; z<element_order; z++) { */
/*       if (ip[z]==n1) { */
/* 	mark1=z; */
/*       } */
/*       if (ip[z]==n2) { */
/* 	mark2=z; */
/*       } */
/*     } */
/*     // Make sure orientation is correct always go from i->j if nj > ni */
/*     if (max(n1,n2)==n1) { */
/*       ihi = mark1; */
/*       ilo = mark2; */
/*     } else { */
/*       ihi = mark2; */
/*       ilo = mark1; */
/*     } */
		
/*     phi[i*mydim+0] = elen*(p[ilo]*dpx[ihi] - p[ihi]*dpx[ilo]); */
/*     phi[i*mydim+1] = elen*(p[ilo]*dpy[ihi] - p[ihi]*dpy[ilo]); */
		
/*     /\* Now compute Curls */
/*      * In 2D curl v = (-dy,dx)*(v1,v2)^T = (dx,dy)(0 1;-1 0)(v1,v2)^T = div (Jv) */
/*      * curl phi_eij = |eij|*(grad(p(i))*(J*grad(p(j)))-grad(p(j))*(J*grad(p(i))) */
/*      * This results from the fact that the p's are linear... */
/*      *\/ */
		
/*     cphi[i] = elen*2*(dpx[ilo]*dpy[ihi] - dpy[ilo]*dpx[ihi]); */
/*   } */
	
/*   //if(ie) free(ie); */
/*   if(p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(ip) free(ip); */
/*   if(tau_edge) free(tau_edge); */
/*   if(midpt) free(midpt); */
	
/*   return; */
/* } */
/* /\***************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void rt_basis(REAL *phi,REAL *dphi,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_face,INT *iface_n,INT *jface_n, \ */
/* 	      REAL *f_area,REAL *f_norm,INT face_order,REAL *xn,REAL *yn,REAL *zn,INT element_order,INT mydim) */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Build basis functions *************************************************\/ */
	
/*   /\* Compute the Raviart Thomas Elements */
/*    *    INPUT: */
/*    *          x,y,z                Coordinates on physical triangle where to compute basis */
/*    *          jel_n	         Nodes belonging to particular element */
/*    *	    jel_face             Faces belonging to particular element */
/*    *          el_face              Ordering of faces on element */
/*    *	    iface_n,jface_n	 Face to Node Map */
/*    *          fel_order            Nodes belonging to face and ordering */
/*    *          f_area               Area of faces */
/*    *          f_norm               Normal vector of faces */
/*    *          face_order           Number of faces per element */
/*    *	    element_order	 Number of nodes per element */
/*    *	    mydim		 Dimension of Problem */
/*    *          xn,yn,zn             Coordinates of nodes */
/*    *    OUTPUT: */
/*    *          phi(face_order,mydim)    Basis functions at particular point (face_order for */
/*    *                                   each face, 6 total in 2D 12 in 3D) (from reference triangle) */
/*    *          Dphi(face_order)         Div of Basis functions at particular point */
/*    *\/ */
	
/*   INT i,j,ica,icb,jcnt; */
/*   REAL* p; */
/*   REAL* dpx; */
/*   REAL* dpy; */
/*   REAL* dpz; */
/*   INT* ipf = calloc(mydim,sizeof(INT)); */
/*   INT myf; */
/*   REAL farea; */
/*   REAL* xf = calloc(mydim,sizeof(REAL)); */
/*   REAL* yf = calloc(mydim,sizeof(REAL)); */
/*   REAL* zf; */
/*   if(mydim==3) { zf = calloc(mydim,sizeof(REAL)); } */

/*   INT elnd,ef1,ef2,ef3; */
	
/*   /\* Get Linear Basis Functions for particular element *\/ */
/*   p = calloc(element_order,sizeof(REAL)); */
/*   dpx = calloc(element_order,sizeof(REAL)); */
/*   dpy = calloc(element_order,sizeof(REAL)); */
	
/*   if(mydim==2) { */
/*     lin_tri_2D(p,dpx,dpy,x,y,jel_n,xn,yn,element_order); */
/*   } else if(mydim==3) { */
/*     dpz = calloc(element_order,sizeof(REAL)); */
/*     lin_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_n,xn,yn,zn,element_order); */
/*   } else { */
/*     printf("You have now entered the twilight zone, in rt_basis()"); */
/*     exit(0); */
/*   } */
	
/*   // Go through each face and find the corresponding nodes */
/*   if(mydim==2) { */
/*     for (i=0; i<face_order; i++) { */
/*       myf = jel_face[i]; */
/*       ica = iface_n[myf-1]-1; */
/*       icb = iface_n[myf]-1; */
/*       jcnt=0; */
/*       for(j=ica;j<icb;j++) { */
/* 	ipf[jcnt] = jface_n[j]; */
/* 	xf[jcnt] = xn[ipf[jcnt]-1]; */
/* 	yf[jcnt] = yn[ipf[jcnt]-1]; */
/* 	jcnt++; */
/*       } */
    
/*       // Get the area and normal vector of the face */
/*       farea = f_area[myf-1]; */
/*       /\* nx = f_norm[(myf-1)*mydim]; *\/ */
/*       /\* ny = f_norm[(myf-1)*mydim+1]; *\/ */

/*       /\* // Determine proper orientation of basis vectors  Compute n^(\perp)*t.  If + use face_node ordering, if - switch sign *\/ */
/*       /\* tx = xf[1]-xf[0]; *\/ */
/*       /\* ty = yf[1]-yf[0]; *\/ */
/*       /\* mysign = -ny*tx + nx*ty; *\/ */
/*       /\* if(mysign<0) { *\/ */
/*       /\* 	nf1 = 1; *\/ */
/*       /\* 	nf2 = 0; *\/ */
/*       /\* } else { *\/ */
/*       /\* 	nf1 = 0; *\/ */
/*       /\* 	nf2 = 1; *\/ */
/*       /\* } *\/ */

/*       // Loop through Nodes on element to find corresponding nodes */
/*       for(j=0;j<element_order;j++) { */
/* 	elnd = jel_n[j]; */
/* 	if(ipf[0]==elnd) { */
/* 	  ef1 = j; */
/* 	} */
/* 	if(ipf[1]==elnd) { */
/* 	  ef2 = j; */
/* 	} */
/*       } */
     
/*       /\* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 2D */
/*        * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i))) */
/*        * |fij| = |eij| */
/*        *\/ */
/*       phi[i*mydim+0] = farea*(p[ef1]*dpy[ef2] - p[ef2]*dpy[ef1]); */
/*       phi[i*mydim+1] = farea*(-p[ef1]*dpx[ef2] + p[ef2]*dpx[ef1]); */
		
/*       // Compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i))) */
/*       dphi[i] = 2*farea*(dpx[ef1]*dpy[ef2] - dpx[ef2]*dpy[ef1]); */

/*     } */
/*   } else if(mydim==3) { */
/*     for (i=0; i<face_order; i++) { */
/*       myf = jel_face[i]; */
/*       ica = iface_n[myf-1]-1; */
/*       icb = iface_n[myf]-1; */
/*       jcnt=0; */
/*       for(j=ica;j<icb;j++) { */
/* 	ipf[jcnt] = jface_n[j]; */
/* 	xf[jcnt] = xn[ipf[jcnt]-1]; */
/* 	yf[jcnt] = yn[ipf[jcnt]-1]; */
/* 	zf[jcnt] = zn[ipf[jcnt]-1]; */
/* 	jcnt++; */
/*       } */
    
/*       // Get the area and normal vector of the face */
/*       farea = f_area[myf-1]; */
/*       /\* nx = f_norm[(myf-1)*mydim]; *\/ */
/*       /\* ny = f_norm[(myf-1)*mydim+1]; *\/ */
/*       /\* nz = f_norm[(myf-1)*mydim+2]; *\/ */
/*       /\* //  printf("\n\nface %d:\tnx=%f\tny=%f\tnz=%f\n",myf,nx,ny,nz); *\/ */

/*       /\* // Determine proper orientation of basis vectors  Compute n^(\perp)*t.  If + use face_node ordering, if - switch sign *\/ */
/*       /\* tx = (yf[1]-yf[0])*(zf[2]-zf[0]) - (zf[1]-zf[0])*(yf[2]-yf[0]); *\/ */
/*       /\* ty = (zf[1]-zf[0])*(xf[2]-xf[0]) - (xf[1]-xf[0])*(zf[2]-zf[0]); *\/ */
/*       /\* tz = (xf[1]-xf[0])*(yf[2]-yf[0]) - (yf[1]-yf[0])*(xf[2]-xf[0]); *\/ */
/*       /\* mysign = nx*tx + ny*ty + nz*tz; *\/ */
/*       /\* if(mysign<0) { *\/ */
/*       /\* 	nf1 = 0; *\/ */
/*       /\* 	nf2 = 2; *\/ */
/*       /\* 	nf3 = 1; *\/ */
/*       /\* } else { *\/ */
/*       /\* 	nf1 = 0; *\/ */
/*       /\* 	nf2 = 1; *\/ */
/*       /\* 	nf3 = 2; *\/ */
/*       /\* } *\/ */

/*       // Loop through Nodes on element to find corresponding nodes */
/*       for(j=0;j<element_order;j++) { */
/* 	elnd = jel_n[j]; */
/* 	if(ipf[0]==elnd) { */
/* 	  ef1 = j; */
/* 	} */
/* 	if(ipf[1]==elnd) { */
/* 	  ef2 = j; */
/* 	} */
/* 	if(ipf[2]==elnd) { */
/* 	  ef3 = j; */
/* 	} */
/*       } */
/*       //printf("\n\nface %d:\tef1=%d\tef2=%d\tef3=%d\n",myf,ef1,ef2,ef3); */
/*       /\* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 3D */
/*        * phi_fijk = 6*|fijk|*(p(i)(grad(p(j)) x grad(p(k))) - p(j)(grad(p(k)) x grad(p(i))) + p(k)(grad(p(i)) x grad(p(j)))) */
/*        * |fijk| = Area(Face) */
/*        *\/ */
/*       phi[i*mydim+0] = 2*farea*(p[ef1]*(dpy[ef2]*dpz[ef3]-dpz[ef2]*dpy[ef3]) + p[ef2]*(dpy[ef3]*dpz[ef1]-dpz[ef3]*dpy[ef1]) + \ */
/* 				p[ef3]*(dpy[ef1]*dpz[ef2]-dpz[ef1]*dpy[ef2])); */
/*       phi[i*mydim+1] = 2*farea*(p[ef1]*(dpz[ef2]*dpx[ef3]-dpx[ef2]*dpz[ef3]) + p[ef2]*(dpz[ef3]*dpx[ef1]-dpx[ef3]*dpz[ef1]) + \ */
/* 				p[ef3]*(dpz[ef1]*dpx[ef2]-dpx[ef1]*dpz[ef2])); */
/*       phi[i*mydim+2] = 2*farea*(p[ef1]*(dpx[ef2]*dpy[ef3]-dpy[ef2]*dpx[ef3]) + p[ef2]*(dpx[ef3]*dpy[ef1]-dpy[ef3]*dpx[ef1]) + \ */
/* 				p[ef3]*(dpx[ef1]*dpy[ef2]-dpy[ef1]*dpx[ef2])); */
		
/*       // Compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i))) */
/*       dphi[i] = 6*farea*(dpx[ef1]*(dpy[ef2]*dpz[ef3]-dpz[ef2]*dpy[ef3]) + dpy[ef1]*(dpz[ef2]*dpx[ef3]-dpx[ef2]*dpz[ef3]) + \ */
/* 			 dpz[ef1]*(dpx[ef2]*dpy[ef3]-dpy[ef2]*dpx[ef3])); */
/*     } */
/*   } else { */
/*     baddimension(); */
/*   } */
	
/*   //if(ie) free(ie); */
/*   if(p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(xf) free(xf); */
/*   if(yf) free(yf); */
/*   if(ipf) free(ipf); */
/*   if(mydim==3) { */
/*     if(dpz) free(dpz); */
/*     if(zf) free(zf); */
/*   } */
	
/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void bdm1_basis(REAL *phi,REAL *dphix,REAL *dphiy,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_face,CSRinc face_n, \ */
/* 	      REAL *f_area,INT face_order,coordinates cn,INT element_order,INT mydim) */
/* { */
/*   /\*************************************************************************** */
/*    *** Build basis functions *************************************************\/ */
	
/*   /\* Compute the Brezzi-Douglas-Marini (BDM) Elements of order 1 for now (only 2D for now). */
/*    *    INPUT: */
/*    *          x,y,z              Coordinates on physical triangle where to compute basis */
/*    *          jel_n	         Nodes belonging to particular element */
/*    *	      jel_face           Faces belonging to particular element */
/*    *          el_face            Ordering of faces on element */
/*    *	      face_n	         Face to Node Map */
/*    *          f_area             Area of faces (length of edge in 2D) */
/*    *          f_norm             Normal vector of faces */
/*    *          face_order         Number of faces per element */
/*    *	      element_order	 Number of nodes per element */
/*    *	      mydim		 Dimension of Problem */
/*    *          cn                 Coordinates of nodes */
/*    *    OUTPUT: */
/*    *          phi(2*face_order,mydim)                   Basis functions at particular point (2*face_order for */
/*    *                                                    each face, 12 total in 2D) (from reference triangle) */
/*    *          Dphix(2*face_order,mydim),Dphiy           Derivatives of Basis functions at particular point */
/*    *\/ */
	
/*   INT i,j,ica,icb,jcnt; */
/*   REAL a1,a2,a3,a4; */
/*   REAL* p; */
/*   REAL* dpx; */
/*   REAL* dpy; */
/*   REAL* dpz; */
/*   INT* ipf = calloc(mydim,sizeof(INT)); */
/*   INT myf; */
/*   REAL farea; */
/*   INT elnd,ef1,ef2; */
	
/*   /\* Get Linear Basis Functions for particular element *\/ */
/*   p = calloc(element_order,sizeof(REAL)); */
/*   dpx = calloc(element_order,sizeof(REAL)); */
/*   dpy = calloc(element_order,sizeof(REAL)); */
	
/*   if(mydim==2) { */
/*     lin_tri_2D(p,dpx,dpy,x,y,jel_n,cn.x,cn.y,element_order); */
/*   } else if(mydim==3) { */
/*     dpz = calloc(element_order,sizeof(REAL)); */
/*     lin_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_n,cn.x,cn.y,cn.z,element_order); */
/*   } else { */
/*     printf("You have now entered the twilight zone, in bdm1_basis()"); */
/*     exit(0); */
/*   } */
	
/*   // Go through each face and find the corresponding nodes */
/*   if(mydim==2) { */
/*     for (i=0; i<face_order; i++) { */
/*       myf = jel_face[i]; */
/*       ica = face_n.IA[myf-1]-1; */
/*       icb = face_n.IA[myf]-1; */
/*       jcnt=0; */
/*       for(j=ica;j<icb;j++) { */
/* 	ipf[jcnt] = face_n.JA[j]; */
/* 	jcnt++; */
/*       } */
    
/*       // Get the area of the face */
/*       farea = f_area[myf-1]; */

/*       // Loop through Nodes on element to find corresponding nodes and orientation */
/*       for(j=0;j<element_order;j++) { */
/* 	elnd = jel_n[j]; */
/* 	if(ipf[0]==elnd) { */
/* 	  ef1 = j; */
/* 	} */
/* 	if(ipf[1]==elnd) { */
/* 	  ef2 = j; */
/* 	} */
/*       } */
     
/*       /\* Now, with the linear basis functions p, dpx, and dpy, the BDM1 elements are in 2D */
/*        * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i))) */
/*        * psi_fij = alpha*|fij|*curl(p(i)p(j)) */
/*        * |fij| = |eij| */
/*        *\/ */
/*       phi[i*mydim*2] = farea*(p[ef1]*dpy[ef2] - p[ef2]*dpy[ef1]); */
/*       phi[i*mydim*2+1] = farea*(-p[ef1]*dpx[ef2] + p[ef2]*dpx[ef1]); */
/*       phi[i*mydim*2+2] = -6*farea*(p[ef1]*dpy[ef2] + p[ef2]*dpy[ef2]); */
/*       phi[i*mydim*2+3] = 6*farea*(p[ef1]*dpx[ef2] + p[ef2]*dpx[ef2]); */
		
/*       a1 = dpx[ef1]*dpx[ef2]; */
/*       a2 = dpx[ef1]*dpy[ef2]; */
/*       a3 = dpy[ef1]*dpx[ef2]; */
/*       a4 = dpy[ef1]*dpy[ef2]; */

/*       dphix[i*mydim*2] = farea*(a2-a3); */
/*       dphix[i*mydim*2+1] = 0.0; */
/*       dphix[i*mydim*2+2] = -6*farea*(a2+a3); */
/*       dphix[i*mydim*2+3] = 12*farea*a1; */
/*       dphiy[i*mydim*2] = 0.0; */
/*       dphiy[i*mydim*2+1] = farea*(a2-a3); */
/*       dphiy[i*mydim*2+2] = -12*farea*a4; */
/*       dphiy[i*mydim*2+3] = 6*farea*(a2+a3); */

/*     } */
/*   } else { */
/*     baddimension(); // 3D not implemented */
/*   } */
	
/*   //if(ie) free(ie); */
/*   if(p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(ipf) free(ipf); */
/*   if(mydim==3) { */
/*     if(dpz) free(dpz); */
/*   } */
	
/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */



/* /\****************************************************************************************************************************\/ */
/* void ned_basis3D(REAL *phi,REAL *cphi,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_ed,INT *ied_n,INT *jed_n,REAL *xn, \ */
/* 		 REAL *yn,REAL *zn,INT element_order,INT mydim)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Build basis functions *************************************************\/ */
	
/*   /\* Compute the Nedelec Elements */
/*    *    INPUT: */
/*    *          x,y,z            Coordinates on physical triangle where to compute basis */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *			jel_ed			 Edges belonging to particular element */
/*    *			ied_n,jed_n		 Edge to Node Map */
/*    *			element_order	 Number of nodes per element */
/*    *			mydim			 Dimension of Problem */
/*    *          xn,yn,zn         Coordinates of nodes */
/*    *    OUTPUT: */
/*    *          phi(6,3)         Basis functions at particular point (3 for */
/*    *                           each edge, 18 total) (from reference triangle) */
/*    * */
/*    *\/ */
	
/*   INT i,ica,n1,n2,k,ihi,ilo; */
/*   INT mark1 = -1; */
/*   INT mark2 = -1; */
/*   REAL* p; */
/*   REAL* dpx; */
/*   REAL* dpy; */
/*   REAL* dpz; */
/*   INT* ip = calloc(element_order,sizeof(INT)); */
/*   //REAL p[element_order],dpx[element_order],dpy[element_order],dpz[element_order]; */
/*   //INT ip[element_order]; */
	
	
/*   /\* Get Linear Basis Functions for particular element *\/ */
/*   p = calloc(element_order,sizeof(REAL)); */
/*   dpx = calloc(element_order,sizeof(REAL)); */
/*   dpy = calloc(element_order,sizeof(REAL)); */
/*   dpz = calloc(element_order,sizeof(REAL)); */
/*   lin_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_n,xn,yn,zn,element_order); */
	
/*   REAL elen; */
/*   REAL* tau_edge = calloc(mydim,sizeof(REAL)); */
/*   REAL* midpt = calloc(mydim,sizeof(REAL)); */
/*   //REAL tau_edge[mydim],midpt[mydim]; */
	
/*   INT nedge = 3*(mydim-1); */
/*   INT ie[nedge]; */
	
/*   // Get Nodes correspond to linear basis elements */
/*   for (i=0; i<element_order; i++) { */
/*     ip[i] = jel_n[i]; */
/*   } */
	
/* #define max(a,b) ((a) > (b) ? (a) : (b)) */
	
	
/*   /\* Now, with the linear basis functions p, dpx, and dpy, the nedelec elements are  */
/*    * phi_eij = |eij|*(p(i)grad(p(j)) - p(j)grad(p(i))) */
/*    * |eij| = sqrt(|xj-xi|^2)   */
/*    *\/ */
/*   //INT* ie = calloc((mydim-1)*3,sizeof(INT)); */
/*   // Go through each edge and get length and find the corresponding nodes */
/*   for (i=0; i<nedge; i++) { */
/*     ie[i] = jel_ed[i]; */
/*     ica = ied_n[ie[i]-1]; */
/*     n1 = jed_n[ica-1]; */
/*     n2 = jed_n[ica]; */
/*     edge_stats(&elen,tau_edge,midpt,ie[i],xn,yn,zn,ied_n,jed_n,mydim); */
		
/*     // Find out which linear basis elements line up with nodes on this edge */
/*     for (k=0; k<element_order; k++) { */
/*       if (ip[k]==n1) { */
/* 	mark1=k; */
/*       } */
/*       if (ip[k]==n2) { */
/* 	mark2=k; */
/*       } */
/*     } */
/*     // Make sure orientation is correct always go from i->j if nj > ni */
/*     if (max(n1,n2)==n1) { */
/*       ihi = mark1; */
/*       ilo = mark2; */
/*     } else { */
/*       ihi = mark2; */
/*       ilo = mark1; */
/*     } */
		
/*     phi[i*mydim+0] = elen*(p[ilo]*dpx[ihi] - p[ihi]*dpx[ilo]); */
/*     phi[i*mydim+1] = elen*(p[ilo]*dpy[ihi] - p[ihi]*dpy[ilo]); */
/*     phi[i*mydim+2] = elen*(p[ilo]*dpz[ihi] - p[ihi]*dpz[ilo]); */
		
/*     cphi[i*mydim+0] = 2*elen*(dpy[ilo]*dpz[ihi]-dpy[ihi]*dpz[ilo]); */
/*     cphi[i*mydim+1] = 2*elen*(dpx[ihi]*dpz[ilo]-dpx[ilo]*dpz[ihi]); */
/*     cphi[i*mydim+2] = 2*elen*(dpx[ilo]*dpy[ihi]-dpx[ihi]*dpy[ilo]); */
		
/*     /\* Curls not needed in 3D as <curl u, curl v> operator can be found from laplacian matrix *\/ */
/*   } */
	
/*   //if(ie) free(ie); */
/*   if(p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(dpz) free(dpz); */
/*   if(ip) free(ip); */
/*   if(tau_edge) free(tau_edge); */
/*   if(midpt) free(midpt); */
	
/*   return; */
/* } */
/* /\***************************************************************************\/ */

/****************************************************************************************************************************/
// INTERPOLATION ROUTINES 
/****************************************************************************************************************************/

/* /\****************************************************************************************************************************\/ */
/* void H1_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_dof,REAL *xn,REAL *yn,REAL *zn,INT ndof,INT element_order,INT mydim,INT nun)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate a finite-element approximation to any other point in the given element using H1 elements ************************************************* */
/*    *** so far only P0,P1,or P2 ************************************************************************************************************ */
/*    *    INPUT: */
/*    *			u				 Approximation to interpolate */
/*    *          x,y,z            Coordinates where to compute value */
/*    *          jel_dof			 DOF belonging to particular element */
/*    *          xn,yn,zn         Coordinates of nodes */
/*    *			element_order	 Number of DOF per element */
/*    *			ndof			 Total DOF for each unknown */
/*    *			mydim			 Dimension of Problem */
/*    *			nun				 Number of unknowns in u (u1,u2,etc?) 1 is a scalar... */
/*    *    OUTPUT: */
/*    *          val				 Value of approximation at given values */
/*    * */
/*    *\/ */
	
	
/*   INT i,nd,j; */
/*   REAL coef; */
/*   REAL* p = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpx = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpy = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpz = calloc(element_order,sizeof(REAL)); */
	
/*   if(element_order==3) { // P1 elements 2D */
/*     lin_tri_2D(p,dpx,dpy,x,y,jel_dof,xn,yn,element_order); */
/*   } else if(element_order==6) { // P2 elements 2D */
/*     quad_tri_2D(p,dpx,dpy,x,y,jel_dof,xn,yn,element_order); */
/*   } else if(element_order==4) { // P1 Elements 3D */
/*     lin_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_dof,xn,yn,zn,element_order); */
/*   } else if(element_order==10) { // P2 Elements 3D */
/*     quad_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_dof,xn,yn,zn,element_order); */
/*   } */
	
/*   for (i=0; i<nun; i++) { */
/*     coef = 0.0; */
		
/*     for (j=0; j<element_order; j++) {  */
/*       nd = i*ndof + jel_dof[j]-1; */
/*       coef = coef + u[nd]*p[j]; */
/*     } */
/*     val[i] = coef; */
/*   } */
	
/*   if (p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(dpz) free(dpz); */
/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void H1Grad_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_dof,REAL *xn,REAL *yn,REAL *zn,INT ndof,INT element_order,INT mydim,INT nun)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate the gradient of a finite-element approximation to any other point in the given element using H1 elements *********************** */
/*    *** so far only P0,P1,or P2 ************************************************************************************************************ */
/*    *    INPUT: */
/*    *			u				 Approximation to interpolate */
/*    *          x,y,z            Coordinates where to compute value */
/*    *          jel_dof			 DOF belonging to particular element */
/*    *          xn,yn,zn         Coordinates of nodes */
/*    *			element_order	 Number of DOF per element */
/*    *			ndof			 Total DOF for each unknown */
/*    *			mydim			 Dimension of Problem */
/*    *			nun				 Number of unknowns in u (u1,u2,etc?) 1 is a scalar... */
/*    *    OUTPUT: */
/*    *          val				 Value of gradient of approximation at given values (val(mydim,nun)) */
/*    * */
/*    *\/ */
	
	
/*   INT i,nd,j; */
/*   REAL* coef = calloc(mydim,sizeof(REAL)); */
/*   REAL* p = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpx = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpy = calloc(element_order,sizeof(REAL)); */
/*   REAL* dpz = calloc(element_order,sizeof(REAL)); */
	
/*   if(element_order==3) { // P1 elements 2D */
/*     lin_tri_2D(p,dpx,dpy,x,y,jel_dof,xn,yn,element_order); */
/*   } else if(element_order==6) { // P2 elements 2D */
/*     quad_tri_2D(p,dpx,dpy,x,y,jel_dof,xn,yn,element_order); */
/*   } else if(element_order==4) { // P1 Elements 3D */
/*     lin_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_dof,xn,yn,zn,element_order); */
/*   } else if(element_order==10) { // P2 Elements 3D */
/*     quad_tet_3D(p,dpx,dpy,dpz,x,y,z,jel_dof,xn,yn,zn,element_order); */
/*   } */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nun; i++) { */
/*       coef[0] = 0.0; */
/*       coef[1] = 0.0; */
/*       for (j=0; j<element_order; j++) {  */
/* 	nd = i*ndof + jel_dof[j]-1; */
/* 	coef[0] = coef[0] + u[nd]*dpx[j]; */
/* 	coef[1] = coef[1] + u[nd]*dpy[j]; */
/*       } */
/*       val[i] = coef[0]; */
/*       val[nun+i] = coef[1]; */
/*     } */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nun; i++) { */
/*       coef[0] = 0.0; */
/*       coef[1] = 0.0; */
/*       coef[2] = 0.0; */
/*       for (j=0; j<element_order; j++) {  */
/* 	nd = i*ndof + jel_dof[j]-1; */
/* 	coef[0] = coef[0] + u[nd]*dpx[j]; */
/* 	coef[1] = coef[1] + u[nd]*dpy[j]; */
/* 	coef[2] = coef[2] + u[nd]*dpz[j]; */
/*       } */
/*       val[i] = coef[0]; */
/*       val[nun+i] = coef[1]; */
/*       val[2*nun+i] = coef[2]; */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
/*   if (coef) free(coef); */
/*   if (p) free(p); */
/*   if(dpx) free(dpx); */
/*   if(dpy) free(dpy); */
/*   if(dpz) free(dpz); */
	
/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void Ned_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_ed,INT *ied_n,INT *jed_n,REAL *xn,REAL *yn,REAL *zn,INT element_order,INT edge_order,INT mydim)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate a finite-element approximation to any other point in the given element using Nedelec elements ************************************************* */
/*    *** so far only first-order ************************************************************************************************************ */
/*    *    INPUT: */
/*    *			u				 Approximation to interpolate */
/*    *          x,y,z            Coordinates where to compute value */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *			jel_ed			 Edges for particular element */
/*    *			ied_n,jed_n		 Edge to Node map */
/*    *          xn,yn,zn         Coordinates of nodes */
/*    *			element_order	 Number of Nodes per element */
/*    *			edge_order		 Number of Edges per element */
/*    *			mydim			 Dimension of Problem */
/*    *    OUTPUT: */
/*    *          val				 Value of approximation at given values (vector) */
/*    * */
/*    *\/ */
	
	
/*   INT ed,j; */
/*   REAL coef1,coef2,coef3; */
/*   REAL* phi = calloc(edge_order*mydim,sizeof(REAL)); */
/*   REAL* cphi2d= calloc(edge_order,sizeof(REAL)); */
/*   REAL* cphi = calloc(edge_order*mydim,sizeof(REAL)); */
	
/*   if (mydim==2) { */
/*     ned_basis(phi,cphi2d,x,y,jel_n,jel_ed,ied_n,jed_n,xn,yn,element_order,mydim); */
/*   } else if (mydim==3) { */
/*     ned_basis3D(phi,cphi,x,y,z,jel_n,jel_ed,ied_n,jed_n,xn,yn,zn,element_order,mydim); */
/*   } else { */
/*     printf("Dimension not 2 or 3!"); */
/*     exit(0); */
/*   } */
	
/*   coef1 = 0.0; */
/*   coef2 = 0.0; */
/*   coef3 = 0.0; */
		
/*   if (mydim==2) { */
/*     for (j=0; j<edge_order; j++) {  */
/*       ed =jel_ed[j]-1; */
/*       coef1 = coef1 + u[ed]*phi[j*mydim+0]; */
/*       coef2 = coef2 + u[ed]*phi[j*mydim+1]; */
/*     } */
/*     val[0] = coef1; */
/*     val[1] = coef2; */
/*   } else if (mydim==3) { */
/*     for (j=0; j<edge_order; j++) {  */
/*       ed =jel_ed[j]-1; */
/*       coef1 = coef1 + u[ed]*phi[j*mydim+0]; */
/*       coef2 = coef2 + u[ed]*phi[j*mydim+1]; */
/*       coef3 = coef3 + u[ed]*phi[j*mydim+2]; */
/*     } */
/*     val[0] = coef1; */
/*     val[1] = coef2; */
/*     val[2] = coef3; */
/*   } */
	
/*   if (phi) free(phi); */
/*   if(cphi) free(cphi); */
/*   if(cphi2d) free(cphi2d); */

/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void RT_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_f,INT *if_n,INT *jf_n,REAL *xn,REAL *yn,REAL *zn,INT element_order,INT face_order,INT mydim,REAL* f_area,REAL* f_norm)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate a finite-element approximation to any other point in the given element using Raviart Thomas elements *************** */
/*    *** so far only first-order ************************************************************************************************************ */
/*    *    INPUT: */
/*    *	      u				 Approximation to interpolate */
/*    *          x,y,z                      Coordinates where to compute value */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *	      jel_f			 Faces for particular element */
/*    *	      if_n,jf_n	         	 Face to Node map */
/*    *          xn,yn,zn                   Coordinates of nodes */
/*    *	      element_order	         Number of Nodes per element */
/*    *   	      face_order		 Number of Faces per element */
/*    *	      mydim			 Dimension of Problem */
/*    *    OUTPUT: */
/*    *          val				 Value of approximation at given values (vector) */
/*    * */
/*    *\/ */
	
	
/*   INT face,j; */
/*   REAL coef1,coef2,coef3; */
/*   REAL* phi = calloc(face_order*mydim,sizeof(REAL)); */
/*   REAL* dphi= calloc(face_order,sizeof(REAL)); */
  
/*   rt_basis(phi,dphi,x,y,z,jel_n,jel_f,if_n,jf_n,f_area,f_norm,face_order,xn,yn,zn,element_order,mydim); */
    
	
/*   coef1 = 0.0; */
/*   coef2 = 0.0; */
/*   coef3 = 0.0; */
		
/*   if (mydim==2) { */
/*     for (j=0; j<face_order; j++) {  */
/*       face =jel_f[j]-1; */
/*       coef1 = coef1 + u[face]*phi[j*mydim+0]; */
/*       coef2 = coef2 + u[face]*phi[j*mydim+1]; */
/*     } */
/*     val[0] = coef1; */
/*     val[1] = coef2; */
/*   } else if (mydim==3) { */
/*     for (j=0; j<face_order; j++) {  */
/*       face =jel_f[j]-1; */
/*       coef1 = coef1 + u[face]*phi[j*mydim+0]; */
/*       coef2 = coef2 + u[face]*phi[j*mydim+1]; */
/*       coef3 = coef3 + u[face]*phi[j*mydim+2]; */
/*     } */
/*     val[0] = coef1; */
/*     val[1] = coef2; */
/*     val[2] = coef3; */
/*   } */
	
/*   if (phi) free(phi); */
/*   if(dphi) free(dphi); */

/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void NedCurl_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_ed,INT *ied_n,INT *jed_n,REAL *xn,REAL *yn,REAL *zn,INT element_order,INT edge_order,INT mydim)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate the curl of a finite-element approximation to any other point in the given element using Nedelec elements ************************************************* */
/*    *** so far only first-order ************************************************************************************************************ */
/*    *    INPUT: */
/*    *			u				 Approximation to interpolate */
/*    *          x,y,z            Coordinates where to compute value */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *			jel_ed			 Edges for particular element */
/*    *			ied_n,jed_n		 Edge to Node map */
/*    *          xn,yn,zn         Coordinates of nodes */
/*    *			element_order	 Number of Nodes per element */
/*    *			edge_order		 Number of Edges per element */
/*    *			mydim			 Dimension of Problem */
/*    *    OUTPUT: */
/*    *          val				 Value of approximation at given values (vector in 3D scalar in 1D) */
/*    * */
/*    *\/ */
	
	
/*   INT ed,j; */
/*   REAL coef1,coef2,coef3; */
/*   REAL* phi = calloc(edge_order*mydim,sizeof(REAL)); */
/*   REAL* cphi2d= calloc(edge_order,sizeof(REAL)); */
/*   REAL* cphi = calloc(mydim+1,sizeof(REAL)); */
	
/*   if (mydim==2) { */
/*     ned_basis(phi,cphi2d,x,y,jel_n,jel_ed,ied_n,jed_n,xn,yn,element_order,mydim); */
/*   } else if (mydim==3) { */
/*     ned_basis3D(phi,cphi,x,y,z,jel_n,jel_ed,ied_n,jed_n,xn,yn,zn,element_order,mydim); */
/*   } else { */
/*     printf("Dimension not 2 or 3!"); */
/*     exit(0); */
/*   } */
	
/*   coef1 = 0.0; */
/*   coef2 = 0.0; */
/*   coef3 = 0.0; */
	
/*   if (mydim==2) { */
/*     for (j=0; j<edge_order; j++) {  */
/*       ed =jel_ed[j]-1; */
/*       coef1 = coef1 + u[ed]*cphi2d[j]; */
/*     } */
/*     val[0] = coef1; */
/*   } else if (mydim==3) { */
/*     for (j=0; j<edge_order; j++) {  */
/*       ed =jel_ed[j]-1; */
/*       coef1 = coef1 + u[ed]*cphi[j*mydim+0]; */
/*       coef2 = coef2 + u[ed]*cphi[j*mydim+1]; */
/*       coef3 = coef3 + u[ed]*cphi[j*mydim+2]; */
/*     } */
/*     val[0] = coef1; */
/*     val[1] = coef2; */
/*     val[2] = coef3; */
/*   } */
	
/*   if (phi) free(phi); */
/*   if(cphi2d) free(cphi2d); */
/*   if(cphi) free(cphi); */
	
/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* void RTDiv_interpolation(REAL* val,REAL *u,REAL x,REAL y,REAL z,INT *jel_n,INT *jel_f,INT *if_n,INT *jf_n,REAL *xn,REAL *yn,REAL *zn,INT element_order,INT face_order,INT mydim,REAL* f_area,REAL* f_norm)  */
/* { */
	
/*   /\*************************************************************************** */
/*    *** Interpolate the divergence a finite-element approximation to any other point in the given element using Raviart Thomas elements *************** */
/*    *** so far only first-order ************************************************************************************************************ */
/*    *    INPUT: */
/*    *	      u				 Approximation to interpolate */
/*    *          x,y,z                      Coordinates where to compute value */
/*    *          jel_n			 Nodes belonging to particular element */
/*    *	      jel_f			 Faces for particular element */
/*    *	      if_n,jf_n	         	 Face to Node map */
/*    *          xn,yn,zn                   Coordinates of nodes */
/*    *	      element_order	         Number of Nodes per element */
/*    *   	      face_order		 Number of Faces per element */
/*    *	      mydim			 Dimension of Problem */
/*    *    OUTPUT: */
/*    *          val				 Value of approximation at given values (vector) */
/*    * */
/*    *\/ */
	
	
/*   INT face,j; */
/*   REAL coef; */
/*   REAL* phi = calloc(face_order*mydim,sizeof(REAL)); */
/*   REAL* dphi= calloc(face_order,sizeof(REAL)); */
  
/*   rt_basis(phi,dphi,x,y,z,jel_n,jel_f,if_n,jf_n,f_area,f_norm,face_order,xn,yn,zn,element_order,mydim); */
    
/*   coef = 0.0; */
		
/*   for (j=0; j<face_order; j++) {  */
/*     face =jel_f[j]-1; */
/*     coef = coef + u[face]*dphi[j]; */
/*   } */
/*   val[0] = coef; */
    	
/*   if (phi) free(phi); */
/*   if(dphi) free(dphi); */

/*   return; */
/* } */
/* /\****************************************************************************************************************************\/ */

/* /\****************************************************************************************************************************\/ */
/* // Error and Norm ROUTINES */
/* /\****************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void L2error(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT test, \ */
/* 	     REAL param,REAL time)  */
/* { */
	
/*   /\* Computes the L2 Norm of the Error using a high order quadrature or the Mass matrix if given */
/*    * */
/*    * Input:		u				Numerical Solution at DOF */
/*    *	       		xn,yn,zn		Coordinates of Nodes */
/*    *	       		iel_n,jel_n		Element to Node Map */
/*    *	       		nelm			Number of Elements */
/*    *	       		nve			Number of Vertices per element */
/*    *   			mydim			Dimension */
/*    *   			nq1d			Number of Quadrature Points per dimension */
/*    *	       		element_order	        Number of Nodes per Element */
/*    *		       	test			Which solution to compare to */
/*    * */
/*    * Output:		err				L2 Norm of the Error */
/*    * */
/*    *\/ */
	
/*   INT i,j,rowa,iq;				/\* Loop Indices *\/ */
	
/*   REAL el2 = 0;			/\* Error *\/ */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL uh,utrue;			/\* Numerical Solution and True Solution at Quadrature Nodes *\/ */
	
/*   for (i=0; i<nelm; i++) { */
		
/*     // Get Vertices (must be vertices) */
/*     rowa = iel_n[i]-1; */
/*     for (j=0; j<nve; j++) {  // First ones are always vertices  */
/*       ipv[j] = jel_n[rowa+j]; */
/*       ipn[j] = ipv[j]; */
/*     } */
/*     for (j=nve; j<element_order; j++) { */
/*       ipn[j] = jel_n[rowa+j]; */
/*     } */
/*     if (mydim==2) { */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
/*     } else if (mydim==3) { */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
/*     } else { */
/*       printf("Bad Dimension"); */
/*       return; */
/*     } */
		
/*     for (iq=0;iq<nq;iq++) { */
/*       x = xq[iq]; */
/*       y = yq[iq]; */
/*       if (mydim>2) { z = zq[iq]; } */
/*       w = wq[iq]; */
			
/*       if (element_order==1) { */
/* 	uh=u[i]; */
/*       } else { */
/* 	H1_interpolation(&uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
/*       } */

/*       // get true solution at quadrature node */
/*       getknownfunction(&utrue,x,y,z,time,mydim,1,param,test); */
			
/*       el2 = el2 + ( uh - utrue )*( uh - utrue)*w; */
/*     } */
		
/*   }	 */
	
/*   el2 = sqrt ( el2 ); */
 
	
/*   *err = el2; */
	
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(ipn) free(ipn); */
/*   if(ipv) free(ipv); */
	
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void L2normMass(REAL *norm,REAL *u,INT* iM,INT* jM,REAL* M,INT imin,INT imax)  */
/* { */
	
/*   /\* Computes the L2 Norm of a vector using the Mass matrix if given for any type of element... */
/*    * */
/*    * Input:		u		       	Numerical Solution at DOF */
/*    *                    iM,jM,M                 Computes L2 norm much faster by <Mu,u> */
/*    *                    ndof                    Number of DOF */
/*    * */
/*    * Output:		norm		       	L2 Norm */
/*    * */
/*    *\/ */
		
/*   REAL el2 = 0.0;			/\* Error *\/ */

/*   INT i,j,jk,rowa,rowb; */
/*   REAL ui,mij; */

/*   for(i=imin;i<imax;i++) { */
/*     rowa=iM[i]-1; */
/*     rowb=iM[i+1]-1; */
/*     ui=0.0; */
/*     for (jk=rowa ; jk< rowb; jk++) { */
/*       j=jM[jk]-1; */
/*       mij=M[jk]; */
/*       ui+=mij*u[j]; */
/*     } */
/*     el2+=ui*u[i]; */
/*   } */
/*   *norm = sqrt(el2); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void L2normMassassemble(REAL *norm,REAL *u,INT* iel_n,INT* jel_n,INT* iel_dof,INT* jel_dof,REAL* xn,REAL* yn,REAL* zn,REAL* xq,REAL* yq,REAL* zq, \ */
/* 			REAL* wq,INT nelm,INT nq1d,INT mydim,INT dof_order,INT element_order,INT* idof_n,INT* jdof_n,REAL* f_area,REAL* f_norm, \ */
/* 			INT elementtype)  */
/* { */

/*   /\* Computes the L2 Norm of a vector using the Mass matrix if given for any type of element... */
/*    * */
/*    * Input:		u		       	Numerical Solution at DOF */
/*    *                    iM,jM,M                 Computes L2 norm much faster by <Mu,u> */
/*    *                    ndof                    Number of DOF */
/*    * */
/*    * Output:		norm		       	L2 Norm */
/*    * */
/*    *\/ */
		
/*   INT i,j,k,rowa,rowb,ndof,jcntr,elm; */
/*   REAL sum = 0.0; */
/*   REAL* MLoc = calloc(dof_order*dof_order,sizeof(REAL)); */
/*   INT* ipdof = calloc(dof_order,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */

/*   /\* Loop over all Elements *\/ */
/*   for (i=0; i<nelm; i++) { */
		
/*     ndof = iel_dof[i+1]-iel_dof[i]; */
		
/*     // Zero out local matrices */
/*     for (j=0; j<ndof*ndof; j++) { */
/*       MLoc[j]=0.0; */
/*     } */
		
/*     // Find DOF for given Element if not H1 elements */
/*     rowa = iel_dof[i]-1; */
/*     rowb = iel_dof[i+1]-1; */
/*     jcntr = 0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipdof[jcntr] = jel_dof[j]; */
/*       jcntr++; */
/*     } */

/*     //Find Nodes for given Element if not H1 elements */
/*     rowa = iel_n[i]-1; */
/*     rowb = iel_n[i+1]-1; */
/*     jcntr=0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipn[jcntr] = jel_n[j]; */
/*       jcntr++; */
/*     } */
    
		
/*     // Compute Local Stiffness Matrix for given Element */
/*     elm = i+1; */

/*     if(elementtype==0) { // H1 Elements Linears or Quadratics */
/*       H1_massL(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ndof,nq1d,mydim,elm); */
/*     } else if(elementtype==1) { // Nedelec Elements */
/*       if (mydim==2) { */
/* 	ned_massL2D(MLoc,xn,yn,xq,yq,wq,ipn,ipdof,idof_n,jdof_n,ndof,element_order,nq1d,mydim,elm); */
/*       } else if (mydim==3) { */
/* 	ned_massL3D(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ipdof,idof_n,jdof_n,ndof,element_order,nq1d,mydim,elm); */
/*       } else { */
/* 	printf("Your dimension isn't 2 or 3.  Welcome to the Twilight Zone..."); */
/*       } */
/*     } else if(elementtype==2) { // Raviart-Thomas Elements */
/*       rt_massL(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ipdof,idof_n,jdof_n,ndof,f_area,f_norm,element_order,nq1d,mydim,elm); */
/*     } */

/*     for(j=0;j<ndof;j++) { */
/*       for(k=0;k<ndof;k++) { */
/* 	sum+=u[ipdof[j]-1]*MLoc[j*ndof+k]*u[ipdof[k]-1]; */
/*       } */
/*     } */
/*   } */

    
/*   *norm = sqrt(sum); */

/*   if(MLoc) free(MLoc); */
/*   if(ipn) free(ipn); */
/*   if(ipdof) free(ipdof); */

/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void L2innerprodMassassemble(REAL *norm,REAL *u,REAL *v,INT* iel_n,INT* jel_n,INT* iel_dof,INT* jel_dof,REAL* xn,REAL* yn,REAL* zn,REAL* xq,REAL* yq,REAL* zq, \ */
/* 			REAL* wq,INT nelm,INT nq1d,INT mydim,INT dof_order,INT element_order,INT* idof_n,INT* jdof_n,REAL* f_area,REAL* f_norm, \ */
/* 			INT elementtype)  */
/* { */

/*   /\* Computes the L2 inner product of two vectors using the Mass matrix if given for any type of element... */
/*    * */
/*    * Input:		u		       	Numerical Solution at DOF */
/*    *                    iM,jM,M                 Computes L2 norm much faster by <Mu,u> */
/*    *                    ndof                    Number of DOF */
/*    * */
/*    * Output:		norm		       	L2 Norm */
/*    * */
/*    *\/ */
		
/*   INT i,j,k,rowa,rowb,ndof,jcntr,elm; */
/*   REAL sum = 0.0; */
/*   REAL* MLoc = calloc(dof_order*dof_order,sizeof(REAL)); */
/*   INT* ipdof = calloc(dof_order,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */

/*   /\* Loop over all Elements *\/ */
/*   for (i=0; i<nelm; i++) { */
		
/*     ndof = iel_dof[i+1]-iel_dof[i]; */
		
/*     // Zero out local matrices */
/*     for (j=0; j<ndof*ndof; j++) { */
/*       MLoc[j]=0.0; */
/*     } */
		
/*     // Find DOF for given Element if not H1 elements */
/*     rowa = iel_dof[i]-1; */
/*     rowb = iel_dof[i+1]-1; */
/*     jcntr = 0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipdof[jcntr] = jel_dof[j]; */
/*       jcntr++; */
/*     } */

/*     //Find Nodes for given Element if not H1 elements */
/*     rowa = iel_n[i]-1; */
/*     rowb = iel_n[i+1]-1; */
/*     jcntr=0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipn[jcntr] = jel_n[j]; */
/*       jcntr++; */
/*     } */
    
		
/*     // Compute Local Stiffness Matrix for given Element */
/*     elm = i+1; */

/*     if(elementtype==0) { // H1 Elements Linears or Quadratics */
/*       H1_massL(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ndof,nq1d,mydim,elm); */
/*     } else if(elementtype==1) { // Nedelec Elements */
/*       if (mydim==2) { */
/* 	ned_massL2D(MLoc,xn,yn,xq,yq,wq,ipn,ipdof,idof_n,jdof_n,ndof,element_order,nq1d,mydim,elm); */
/*       } else if (mydim==3) { */
/* 	ned_massL3D(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ipdof,idof_n,jdof_n,ndof,element_order,nq1d,mydim,elm); */
/*       } else { */
/* 	printf("Your dimension isn't 2 or 3.  Welcome to the Twilight Zone..."); */
/*       } */
/*     } else if(elementtype==2) { // Raviart-Thomas Elements */
/*       rt_massL(MLoc,xn,yn,zn,xq,yq,zq,wq,ipn,ipdof,idof_n,jdof_n,ndof,f_area,f_norm,element_order,nq1d,mydim,elm); */
/*     } */

/*     for(j=0;j<ndof;j++) { */
/*       for(k=0;k<ndof;k++) { */
/* 	sum+=v[ipdof[j]-1]*MLoc[j*ndof+k]*u[ipdof[k]-1]; */
/*       } */
/*     } */
/*   } */

/*   *norm = sum; */
/*   //\*norm = sqrt(sum); */

/*   if(MLoc) free(MLoc); */
/*   if(ipn) free(ipn); */
/*   if(ipdof) free(ipdof); */

/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void H1seminormLapassemble(REAL *norm,REAL *u,CSRinc el_dof,CSRinc el_n,coordinates cn,INT nelm,INT nq1d,INT mydim,INT dof_order,INT element_order, \ */
/* 			   CSRinc dof_n,REAL* f_area,REAL* f_norm,INT elementtype)  */
/* { */

/*   /\* Computes the H1 semi-Norm of a vector using the laplacian-like matrix if given for any type of element... */
/*    *          Nodal - <grad u, grad u> - |u|_1 */
/*    *          RT - <div u, div u>      - |u|_(H(div)) */
/*    *          Nedelec - <curl u, curl u>  - |u|_(H(curl)) */
/*    * */
/*    * Input:		u		       	Numerical Solution at DOF */
/*    *                    ndof                    Number of DOF */
/*    * */
/*    * Output:		norm		       	L2 Norm */
/*    * */
/*    *\/ */
		
/*   INT i,j,k,rowa,rowb,ndof,jcntr,elm; */
/*   REAL sum = 0.0; */
/*   REAL* MLoc = calloc(dof_order*dof_order,sizeof(REAL)); */
/*   INT* ipdof = calloc(dof_order,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   REAL nq; */
/*   qcoordinates cq; */

/*   /\* Get Quadrature Nodes *\/ */
/*   nq = pow(nq1d,mydim); */
/*   allocateqcoords(&cq,nq1d,nelm,mydim); */
/*   if (mydim==2) { */
/*     get_quadrature(cq.x,cq.y,NULL,cq.w,cn.x,cn.y,NULL,el_n.IA,el_n.JA,nelm,element_order,nq1d,mydim); */
/*   } else { */
/*     get_quadrature(cq.x,cq.y,cq.z,cq.w,cn.x,cn.y,cn.z,el_n.IA,el_n.JA,nelm,element_order,nq1d,mydim); */
/*   } */

/*   /\* Loop over all Elements *\/ */
/*   for (i=0; i<nelm; i++) { */
		
/*     ndof = el_dof.IA[i+1]-el_dof.IA[i]; */
		
/*     // Zero out local matrices */
/*     for (j=0; j<ndof*ndof; j++) { */
/*       MLoc[j]=0.0; */
/*     } */
		
/*     // Find DOF for given Element if not H1 elements */
/*     rowa = el_dof.IA[i]-1; */
/*     rowb = el_dof.IA[i+1]-1; */
/*     jcntr = 0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipdof[jcntr] = el_dof.JA[j]; */
/*       jcntr++; */
/*     } */

/*     //Find Nodes for given Element if not H1 elements */
/*     rowa = el_n.IA[i]-1; */
/*     rowb = el_n.IA[i+1]-1; */
/*     jcntr=0; */
/*     for (j=rowa; j<rowb; j++) { */
/*       ipn[jcntr] = el_n.JA[j]; */
/*       jcntr++; */
/*     } */
    
		
/*     // Compute Local Stiffness Matrix for given Element */
/*     elm = i+1; */

/*     if(elementtype==0) { // H1 Elements Linears or Quadratics */
/*       H1_LapL(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ndof,nq1d,mydim,elm); */
/*     } else if(elementtype==1) { // Nedelec Elements */
/*       if (mydim==2) { */
/* 	ned_curlcurlL2D(MLoc,cn.x,cn.y,cq.x,cq.y,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,element_order,nq1d,mydim,elm); */
/*       } else if (mydim==3) { */
/* 	ned_curlcurlL3D(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,element_order,nq1d,mydim,elm); */
/*       } else { */
/* 	printf("Your dimension isn't 2 or 3.  Welcome to the Twilight Zone..."); */
/*       } */
/*     } else if(elementtype==2) { // Raviart-Thomas Elements */
/*       rt_divdivL(MLoc,cn.x,cn.y,cn.z,cq.x,cq.y,cq.z,cq.w,ipn,ipdof,dof_n.IA,dof_n.JA,ndof,f_area,f_norm,element_order,nq1d,mydim,elm); */
/*     } */

/*     for(j=0;j<ndof;j++) { */
/*       for(k=0;k<ndof;k++) { */
/* 	sum+=u[ipdof[j]-1]*MLoc[j*ndof+k]*u[ipdof[k]-1]; */
/*       } */
/*     } */
/*   } */

/*   if(sum<0.0) { */
/*     printf("Your H1 Semi Norm Squared is negative!  Outputting the square itself\n"); */
/*     *norm = sum; */
/*   } else { */
/*     *norm = sqrt(sum); */
/*   } */

/*   if(MLoc) free(MLoc); */
/*   if(ipn) free(ipn); */
/*   if(ipdof) free(ipdof); */
/*   freeqcoords(cq); */

/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\***************************************************************************\/ */
/* void H1semierror(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT test,REAL param,REAL time)  */
/* { */

/*   /\* Computes the H1 Semi-Norm of the Error using a high order quadrature  */
/*    * */
/*    * Input:		u				Numerical Solution at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				test			Which solution to compare to (This is the gradient of the solution so it's a vector) */
/*    * */
/*    * Output:		err				H1 Norm of the Error */
/*    * */
/*    *\/ */
	
/*   INT i,j,rowa,iq;				/\* Loop Indices *\/ */

/*   REAL el = 0.0;			 */
/*   REAL elx = 0.0;			 */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* ugh = calloc(mydim,sizeof(REAL));			/\* Grad of approximation interpolated at Quadrature Nodes *\/ */
/*   REAL* ugtrue = calloc(mydim,sizeof(REAL));		/\* Grad of true solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d);		 */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	H1Grad_interpolation(ugh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
			
/* 	// get true solution at quadrature node */
/* 	getknownfunction(ugtrue,x,y,z,time,mydim,2,param,test); */
			
/* 	elx = elx + ( ugh[0] - ugtrue[0] )*( ugh[0] - ugtrue[0] )*w; */
/* 	ely = ely + ( ugh[1] - ugtrue[1] )*( ugh[1] - ugtrue[1] )*w; */
/*       } */
		
/*     } */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */

/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(ugh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(ugtrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + ( ugh[0] - ugtrue[0] )*( ugh[0] - ugtrue[0] )*w; */
/* 	ely = ely + ( ugh[1] - ugtrue[1] )*( ugh[1] - ugtrue[1] )*w; */
/* 	elz = elz + ( ugh[2] - ugtrue[2] )*( ugh[2] - ugtrue[2] )*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  YOu are in the twilight zone"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(ipn) free(ipn); */
/*   if(ipv) free(ipv); */
	
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2errorRT(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_f,INT *jel_f,INT *if_n,INT *jf_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT face_order,INT test,REAL param,REAL time,REAL* f_area,REAL* f_norm) */
/* { */
	
/*   /\* Computes the L2 error of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature PoINTs per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipf = calloc(face_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
/*   REAL* utrue = calloc(mydim,sizeof(REAL));		/\* True Solution at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/*       } */
			
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/* 	elz = elz + (uh[2]-utrue[2])*(uh[2]-utrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if(ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipf) free(ipf); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2errorNed(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order,INT test,REAL param,REAL time)   */
/* { */
	
/*   /\* Computes the L2 error of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
/*   REAL* utrue = calloc(mydim,sizeof(REAL));		/\* True Solution at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */

/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/*       } */
			
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// get true solution at quadrature node */
/* 	getknownfunction(utrue,x,y,z,time,mydim,2,param,test); */
				
/* 	elx = elx + (uh[0]-utrue[0])*(uh[0]-utrue[0])*w; */
/* 	ely = ely + (uh[1]-utrue[1])*(uh[1]-utrue[1])*w; */
/* 	elz = elz + (uh[2]-utrue[2])*(uh[2]-utrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */


/* /\*******************************************************************************************************************************************************\/ */
/* void NedsemiCurlerror(REAL *err,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order,INT test,REAL param,REAL time)  */
/* { */
	
/*   /\* Computes the H curl Semi-Norm Error of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature PoINTs per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uch = calloc(mydim,sizeof(REAL));			/\* Curl of approximation interpolated at Quadrature Nodes *\/ */
/*   REAL* uctrue = calloc(mydim,sizeof(REAL));		/\* Curl of true solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uch,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
/* 	getknownfunction(uctrue,x,y,z,time,mydim,1,param,test); */
				
/* 	elx = elx + (uch[0]-uctrue[0])*(uch[0]-uctrue[0])*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uch,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	// True Solution */
/* 	getknownfunction(uctrue,x,y,z,time,mydim,2,param,test); */

				
/* 	elx = elx + (uch[0]-uctrue[0])*(uch[0]-uctrue[0])*w; */
/* 	ely = ely + (uch[1]-uctrue[1])*(uch[1]-uctrue[1])*w; */
/* 	elz = elz + (uch[2]-uctrue[2])*(uch[2]-uctrue[2])*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *err = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2norminterp(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el2 = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL uh;			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   for (i=0; i<nelm; i++) { */
		
/*     // Get Vertices (must be vertices) */
/*     rowa = iel_n[i]-1; */
/*     for (j=0; j<nve; j++) {  // First ones are always vertices  */
/*       ipv[j] = jel_n[rowa+j]; */
/*       ipn[j] = ipv[j]; */
/*     } */
/*     for (j=nve; j<element_order; j++) { */
/*       ipn[j] = jel_n[rowa+j]; */
/*     } */
		
/*     // Get quadrature nodes for given order of accuracy */
/*     if (mydim==2) { */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
/*     } else if (mydim==3) { */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
/*     } else { */
/*       printf("Bad Dimension"); */
/*       return; */
/*     } */
		
/*     for (iq=0;iq<nq;iq++) { */
/*       x = xq[iq]; */
/*       y = yq[iq]; */
/*       if (mydim>2) { z = zq[iq]; } */
/*       w = wq[iq]; */
			
/*       if (element_order==1) { */
/* 	uh=u[i]; */
/*       } else { */
/* 	H1_interpolation(&uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
/*       } */
			
/*       el2 = el2 + uh*uh*w; */
/*     } */
		
/*   }	 */
	
/*   el2 = sqrt ( el2 ); */
	
/*   *norm = el2; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void H1seminorminterp(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order)  */
/* { */
	
/*   /\* Computes the H1 Semi-Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				nelm			Number of Elements */
/*    *				ndof			Number of DOF */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    * */
/*    * Output:		norm				H1 semi-norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Grad of approximation interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	H1Grad_interpolation(uh,u,x,y,z,ipn,xn,yn,zn,0,element_order,mydim,1); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2normNed(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
		
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
			
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
		
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	Ned_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
			
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void L2normRT(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_f,INT *jel_f,INT *if_n,INT *jf_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT face_order,REAL* f_area,REAL* f_norm)  */
/* { */
	
/*   /\* Computes the L2 Norm of an approximation defined on Nedelec elements using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0;	 */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipf = calloc(face_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Numerical Solution interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2) { */
/*     for (i=0; i<nelm; i++) { */
		
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
		
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
		
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
			
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
			
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/*       } */
		
/*     }	 */
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       // Get rest of nodes */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_f[i]-1; */
/*       for (j=0; j<face_order; j++) { */
/* 	ipf[j] = jel_f[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	RT_interpolation(uh,u,x,y,z,ipn,ipf,if_n,jf_n,xn,yn,zn,element_order,face_order,mydim,f_area,f_norm); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
			
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3\n"); */
/*     exit(0); */
/*   } */

	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipf) free(ipf); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */

/* /\*******************************************************************************************************************************************************\/ */
/* void NedsemiCurlnorm(REAL *norm,REAL *u,REAL *xn,REAL *yn,REAL *zn,INT *iel_n,INT *jel_n,INT *iel_ed,INT *jel_ed,INT *ied_n,INT *jed_n,INT nelm,INT nve,INT mydim,INT nq1d,INT element_order,INT edge_order)  */
/* { */
	
/*   /\* Computes the H curl Semi-Norm of an approximation using a high order quadrature  */
/*    * */
/*    * Input:		u				Approximation at DOF */
/*    *				xn,yn,zn		Coordinates of Nodes */
/*    *				iel_n,jel_n		Element to Node Map */
/*    *				iel_ed,jel_ed	Element to Edge Map */
/*    *				ied_n,jed_n		Edge to Node Map */
/*    *				nelm			Number of Elements */
/*    *				nve				Number of vertices per element */
/*    *				mydim			Dimension */
/*    *				nq1d			Number of Quadrature Points per dimension */
/*    *				element_order	Number of Nodes per Element */
/*    *				edge_order		Number of Edges per Element */
/*    * */
/*    * Output:		norm			L2 Norm */
/*    * */
/*    *\/ */
	
/*   INT i,j,iq,rowa;				/\* Loop Indices *\/ */
	
/*   REAL el = 0.0; */
/*   REAL elx = 0.0;			/\* Error *\/ */
/*   REAL ely = 0.0; */
/*   REAL elz = 0.0; */
	
/*   INT nq = (INT) pow(nq1d,mydim); */
/*   REAL* xq = calloc(nq,sizeof(REAL)); */
/*   REAL* yq = calloc(nq,sizeof(REAL)); */
/*   REAL* zq = calloc(nq,sizeof(REAL)); */
/*   REAL* wq = calloc(nq,sizeof(REAL)); */
/*   INT* ipv = calloc(nve,sizeof(INT)); */
/*   INT* ipn = calloc(element_order,sizeof(INT)); */
/*   INT* ipe = calloc(edge_order,sizeof(INT)); */
	
/*   REAL x,y,z,w; */
	
/*   REAL* uh = calloc(mydim,sizeof(REAL));			/\* Curl of approximation interpolated at Quadrature Nodes *\/ */
	
/*   if (mydim==2 ) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_2D_tri(xq,yq,wq,xn,yn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = 0; */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/*       } */
/*     } */
		
/*   } else if (mydim==3) { */
/*     for (i=0; i<nelm; i++) { */
			
/*       // Get Vertices (must be vertices) */
/*       rowa = iel_n[i]-1; */
/*       for (j=0; j<nve; j++) {  // First ones are always vertices  */
/* 	ipv[j] = jel_n[rowa+j]; */
/* 	ipn[j] = ipv[j]; */
/*       } */
/*       for (j=nve; j<element_order; j++) { */
/* 	ipn[j] = jel_n[rowa+j]; */
/*       } */
/*       // Get edges */
/*       rowa = iel_ed[i]-1; */
/*       for (j=0; j<edge_order; j++) { */
/* 	ipe[j] = jel_ed[rowa+j]; */
/*       } */
			
/*       // Get quadrature nodes for given order of accuracy */
/*       quad_3D_tet(xq,yq,zq,wq,xn,yn,zn,ipv,nq1d); */
			
/*       for (iq=0;iq<nq;iq++) { */
/* 	x = xq[iq]; */
/* 	y = yq[iq]; */
/* 	z = zq[iq];  */
/* 	w = wq[iq]; */
				
/* 	NedCurl_interpolation(uh,u,x,y,z,ipn,ipe,ied_n,jed_n,xn,yn,zn,element_order,edge_order,mydim); */
				
/* 	elx = elx + uh[0]*uh[0]*w; */
/* 	ely = ely + uh[1]*uh[1]*w; */
/* 	elz = elz + uh[2]*uh[2]*w; */
/*       } */
/*     } */
/*   } else { */
/*     printf("Dimension not 2 or 3.  What is this the Twighlight Zone??"); */
/*     exit(0); */
/*   } */
	
	
/*   el = sqrt ( elx + ely + elz ); */
	
/*   *norm = el; */
	
/*   if (ipv) free(ipv); */
/*   if(ipn) free(ipn); */
/*   if(ipe) free(ipe); */
/*   if(xq) free(xq); */
/*   if(yq) free(yq); */
/*   if(zq) free(zq); */
/*   if(wq) free(wq); */
/*   return; */
/* } */
/* /\*******************************************************************************************************************************************************\/ */
