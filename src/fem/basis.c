/*! \file src/fem/basis.c
 *
 * \brief Compute the basis functions for triangles or tetrahedra or 1D FEM
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/1/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/14/2016
 *
 * \note Typically this involves DOF defined on either the
 *  vertices, edges, or faces.  In most cases, the basis elements are
 *  defined using the standard Lagrange finite-element basis functions
 *  using barycentric coordinates.  See PX_H1_basis for details on the Lagrange
 *  basis functions.
 *
 */

#include "hazmath.h"

/*******************************************************************************************************/
/*!
 * \fn void PX_H1_basis(REAL *p,REAL *dp,REAL *x,INT *dof,INT porder,trimesh *mesh)
 *
 * \brief Compute Standard Lagrange Finite Element Basis Functions (PX) at a particular point in 1, 2 or 3D
 *        For now, we only assume constants, Linears, or Quadratic Elements (P0 or P1 or P2)
 *
 * \param x       Coordinate on where to compute basis function
 * \param dof       DOF on element
 * \param porder  Order of elements
 * \param mesh    Mesh struct
 *
 * \return p      Basis functions (1 for each DOF on element)
 * \return dp     Derivatives of basis functions (i.e., gradient)
 *
 *  \note For 1D we just compute the basis functions directly
 *        For P1: for element between x_1 and x_2:  x1 ------ x2
 *             phi_1 = (x_2 - x)/(x_2 - x_1)
 *             phi_2 = (x-x_1)/(x_2 - x_1)
 *
 *        For P2: for element between x_1 and x_2: x1 ---- x3 ---- x2
 *             phi_1 = (x_2 - x)/(x_2 - x_1)*(2(x_2 - x)/(x_2 - x_1) - 1)
 *             phi_2 = (x-x_1)/(x_2 - x_1)*(2(x-x_1)/(x_2 - x_1) - 1)
 *             phi_3 = 4(x_2 - x)/(x_2 - x_1)(x-x_1)/(x_2 - x_1)
 *
 *        For 2D, we show an illustration here.
 *
 *        The physical element:
 *
 *        In this picture, we don't mean to suggest that the bottom of
 *        the physical triangle is horizontal.  However, we do assume that
 *        each of the sides is a straight line, and that the intermediate
 *        points are exactly halfway on each side.
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
 *      Reference element T3:
 *
 *       In this picture of the reference element, we really do assume
 *       that one side is vertical, one horizontal, of length 1.
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
 *     Determine the (R,S) coordinates corresponding to (X,Y).
 *
 *     What is happening here is that we are solving the linear system:
 *
 *      ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
 *      ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
 *
 *     by computing the inverse of the coefficient matrix and multiplying
 *     it by the right hand side to get R and S.
 *
 *    The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
 *    for R and S.
 *
 *    For quadratic elements:
 *
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
 *
 */
void PX_H1_basis(REAL *p,REAL *dp,REAL *x,INT *dof,INT porder,trimesh *mesh)
{
  REAL dp1r,dp2r,dp3r,dp4r,dp5r,dp6r,dp7r,dp8r,dp9r,dp10r;
  REAL dp1s,dp2s,dp3s,dp4s,dp5s,dp6s,dp7s,dp8s,dp9s,dp10s;
  REAL dp1t,dp2t,dp3t,dp4t,dp5t,dp6t,dp7t,dp8t,dp9t,dp10t;
  REAL onemrst;
  INT i;

  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT dim = mesh->dim;

  REAL* xp = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* yp = NULL;//(REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* zp = NULL;
  coordinates* cv = mesh->cv;

  // P0 elements are trivial and we just need to return a 1 for each element:
  if(porder==0) {
    p[0] = 1.0;
  } else {
    // The remaining depend on the dimension
    if(dim==1) {
      // Get Physical Coordinates of Vertices
      for (i=0; i<v_per_elm; i++) {
        xp[i] = cv->x[dof[i]-1];
      }
      // Get Barycentric Coordinates
      REAL oneoverh = 1.0/(xp[1]-xp[0]);
      REAL lam1 = (xp[1]-x[0])*oneoverh;
      REAL lam2 = (x[0]-xp[0])*oneoverh;

      // Now Get basis functions
      if(porder==1) {
        p[0] = lam1;
        p[1] = lam2;
        dp[0] = -oneoverh;
        dp[1] = oneoverh;
      } else if(porder==2) {
        p[0] = lam1*(2*lam1-1);
        p[1] = lam2*(2*lam2-1);
        p[2] = 4*lam1*lam2;
        dp[0] = -oneoverh*(2*lam1-1) - 2*lam1*oneoverh;
        dp[1] = oneoverh*(2*lam2-1) + 2*lam2*oneoverh;
        dp[2] = 4*lam1*oneoverh - 4*lam2*oneoverh;
      } else {
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }
    } else if(dim==2) {
      // Get Physical Coordinates of Vertices
      yp = (REAL *) calloc(v_per_elm,sizeof(REAL));
      for (i=0; i<v_per_elm; i++) {
        xp[i] = cv->x[dof[i]-1];
        yp[i] = cv->y[dof[i]-1];
      }

      // Get coordinates on reference triangle
      REAL det = (xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0]);

      REAL r = ((yp[2]-yp[0])*(x[0]-xp[0]) + (xp[0]-xp[2])*(x[1]-yp[0]))/det;

      REAL drdx = (yp[2]-yp[0])/det;
      REAL drdy = (xp[0]-xp[2])/det;

      REAL s = ((yp[0]-yp[1])*(x[0]-xp[0]) + (xp[1]-xp[0])*(x[1]-yp[0]))/det;

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
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }

      /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
       *  to (X,Y) using the chain rule.
       */
      dp[0*dim] = dp1r * drdx + dp1s * dsdx;
      dp[0*dim+1] = dp1r * drdy + dp1s * dsdy;
      dp[1*dim] = dp2r * drdx + dp2s * dsdx;
      dp[1*dim+1] = dp2r * drdy + dp2s * dsdy;
      dp[2*dim] = dp3r * drdx + dp3s * dsdx;
      dp[2*dim+1] = dp3r * drdy + dp3s * dsdy;
      if(porder==2) {
        dp[3*dim] = dp4r * drdx + dp4s * dsdx;
        dp[3*dim+1] = dp4r * drdy + dp4s * dsdy;
        dp[4*dim] = dp5r * drdx + dp5s * dsdx;
        dp[4*dim+1] = dp5r * drdy + dp5s * dsdy;
        dp[5*dim] = dp6r * drdx + dp6s * dsdx;
        dp[5*dim+1] = dp6r * drdy + dp6s * dsdy;
      }
    } else if (dim==3) {
      // Get Nodes and Physical Coordinates
      yp = (REAL *) calloc(v_per_elm,sizeof(REAL));
      zp = (REAL *) calloc(v_per_elm,sizeof(REAL));
      for (i=0; i<v_per_elm; i++) {
        xp[i] = cv->x[dof[i]-1];
        yp[i] = cv->y[dof[i]-1];
        zp[i] = cv->z[dof[i]-1];
      }

      // Get coordinates on reference triangle
      REAL det = (xp[3]-xp[0])*((yp[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(zp[1]-zp[0])) \
          - (xp[2]-xp[0])*((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0])) \
          + (xp[1]-xp[0])*((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0]));

      REAL r = ((xp[3]-xp[0])*((x[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(x[2]-zp[0])) \
          - (xp[2]-xp[0])*((x[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(x[2]-zp[0])) \
          + (x[0]-xp[0])*((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0])))/det;

      REAL drdx = ((yp[2]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[2]-zp[0]))/det;
      REAL drdy = ((xp[3]-xp[0])*(zp[2]-zp[0]) - (xp[2]-xp[0])*(zp[3]-zp[0]))/det;
      REAL drdz = ((xp[2]-xp[0])*(yp[3]-yp[0]) - (xp[3]-xp[0])*(yp[2]-yp[0]))/det;

      REAL s = ((xp[3]-xp[0])*((yp[1]-yp[0])*(x[2]-zp[0])-(x[1]-yp[0])*(zp[1]-zp[0])) \
          - (x[0]-xp[0])*((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0])) \
          + (xp[1]-xp[0])*((x[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(x[2]-zp[0])))/det;

      REAL dsdx = -((yp[1]-yp[0])*(zp[3]-zp[0])-(yp[3]-yp[0])*(zp[1]-zp[0]))/det;
      REAL dsdy = ((xp[1]-xp[0])*(zp[3]-zp[0]) - (xp[3]-xp[0])*(zp[1]-zp[0]))/det;
      REAL dsdz = ((xp[3]-xp[0])*(yp[1]-yp[0]) - (xp[1]-xp[0])*(yp[3]-yp[0]))/det;

      REAL t = ((x[0]-xp[0])*((yp[1]-yp[0])*(zp[2]-zp[0])-(yp[2]-yp[0])*(zp[1]-zp[0])) \
          - (xp[2]-xp[0])*((yp[1]-yp[0])*(x[2]-zp[0])-(x[1]-yp[0])*(zp[1]-zp[0])) \
          + (xp[1]-xp[0])*((yp[2]-yp[0])*(x[2]-zp[0])-(x[1]-yp[0])*(zp[2]-zp[0])))/det;

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
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }

      /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
       *  to (X,Y) using the chain rule.
       */
      dp[0*dim] = dp1r * drdx + dp1s * dsdx + dp1t * dtdx;
      dp[0*dim+1] = dp1r * drdy + dp1s * dsdy + dp1t * dtdy;
      dp[0*dim+2] = dp1r * drdz + dp1s * dsdz + dp1t * dtdz;
      dp[1*dim] = dp2r * drdx + dp2s * dsdx + dp2t * dtdx;
      dp[1*dim+1] = dp2r * drdy + dp2s * dsdy + dp2t * dtdy;
      dp[1*dim+2] = dp2r * drdz + dp2s * dsdz + dp2t * dtdz;
      dp[2*dim] = dp3r * drdx + dp3s * dsdx + dp3t * dtdx;
      dp[2*dim+1] = dp3r * drdy + dp3s * dsdy + dp3t * dtdy;
      dp[2*dim+2] = dp3r * drdz + dp3s * dsdz + dp3t * dtdz;
      dp[3*dim] = dp4r * drdx + dp4s * dsdx + dp4t * dtdx;
      dp[3*dim+1] = dp4r * drdy + dp4s * dsdy + dp4t * dtdy;
      dp[3*dim+2] = dp4r * drdz + dp4s * dsdz + dp4t * dtdz;
      if(porder==2) {
        dp[4*dim] = dp5r * drdx + dp5s * dsdx + dp5t * dtdx;
        dp[4*dim+1] = dp5r * drdy + dp5s * dsdy + dp5t * dtdy;
        dp[4*dim+2] = dp5r * drdz + dp5s * dsdz + dp5t * dtdz;
        dp[5*dim] = dp6r * drdx + dp6s * dsdx + dp6t * dtdx;
        dp[5*dim+1] = dp6r * drdy + dp6s * dsdy + dp6t * dtdy;
        dp[5*dim+2] = dp6r * drdz + dp6s * dsdz + dp6t * dtdz;
        dp[6*dim] = dp7r * drdx + dp7s * dsdx + dp7t * dtdx;
        dp[6*dim+1] = dp7r * drdy + dp7s * dsdy + dp7t * dtdy;
        dp[6*dim+2] = dp7r * drdz + dp7s * dsdz + dp7t * dtdz;
        dp[7*dim] = dp8r * drdx + dp8s * dsdx + dp8t * dtdx;
        dp[7*dim+1] = dp8r * drdy + dp8s * dsdy + dp8t * dtdy;
        dp[7*dim+2] = dp8r * drdz + dp8s * dsdz + dp8t * dtdz;
        dp[8*dim] = dp9r * drdx + dp9s * dsdx + dp9t * dtdx;
        dp[8*dim+1] = dp9r * drdy + dp9s * dsdy + dp9t * dtdy;
        dp[8*dim+2] = dp9r * drdz + dp9s * dsdz + dp9t * dtdz;
        dp[9*dim] = dp10r * drdx + dp10s * dsdx + dp10t * dtdx;
        dp[9*dim+1] = dp10r * drdy + dp10s * dsdy + dp10t * dtdy;
        dp[9*dim+2] = dp10r * drdz + dp10s * dsdz + dp10t * dtdz;
      }
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }
  }

  if(xp) free(xp);
  if(yp) free(yp);
  if(zp) free(zp);

  return;
}
/*******************************************************************************************************/

/*******************************************************************************************************/
/*!
 * \fn void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh *mesh)
 *
 * \brief Compute Standard Quadratic Finite Element Basis Functions (P2) at a particular point
 *        Also compute the 2nd derivatives for some reason...
 *
 * \param x,y,z           Coordinate on where to compute basis function
 * \param dof             DOF for the given element (in this case vertices and their global numbering)
 * \param porder          Order of elements
 * \param mesh            Mesh struct
 *
 * \return p              Basis functions (1 for each DOF on element)
 * \return dpx,dpy        Derivatives of basis functions (i.e., gradient)
 * \return dpxx,dpyy,dpxy 2nd Derivatives of basis functions
 *
 *
 */
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh *mesh) 
{
  REAL dp1r,dp2r,dp3r,dp4r,dp5r,dp6r;
  REAL dp1s,dp2s,dp3s,dp4s,dp5s,dp6s;
  REAL onemrs;
  INT i;
  INT v_per_elm = mesh->v_per_elm;
  REAL* xp = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* yp = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL dp1rr,dp2rr,dp3rr,dp4rr,dp5rr,dp6rr;
  REAL dp1ss,dp2ss,dp3ss,dp4ss,dp5ss,dp6ss;
  REAL dp1rs,dp2rs,dp3rs,dp4rs,dp5rs,dp6rs;

  coordinates* cv = mesh->cv;

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
/*******************************************************************************************************/

/*******************************************************************************************************/
/*!
 * \fn void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
 *
 * \brief Compute Nedelec Finite Element Basis Functions (zeroth order) at a particular point in 2 or 3D
 *
 * \param x         Coordinate on where to compute basis function
 * \param v_on_elm  Vertices on element
 * \param dof       DOF on element
 * \param mesh      Mesh struct
 *
 * \return phi      Basis functions (dim for each edge from reference triangle)
 * \return cphi     Curl of basis functions (1 for each edge in 2D, dim for each edge in 3D)
 *
 */
void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT ed_per_elm = mesh->ed_per_elm;
  INT dim = mesh->dim;

  INT i,k,ica,n1,n2,ihi,ilo;
  INT mark1 = -1;
  INT mark2 = -1;
  REAL* p;
  REAL* dp;

  /* Get Linear Basis Functions for particular element */
  p = (REAL *) calloc(v_per_elm,sizeof(REAL));
  dp = (REAL *) calloc(v_per_elm*dim,sizeof(REAL));
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  REAL elen;

  /* Now, with the linear basis functions p, dpx, and dpy, the nedelec elements are
   * phi_eij = |eij|*(p(i)grad(p(j)) - p(j)grad(p(i)))
   * |eij| = sqrt(|xj-xi|^2)
   */

  // Go through each edge and get length and find the corresponding nodes
  for (i=0; i<ed_per_elm; i++) {
    ica = mesh->ed_v->IA[dof[i]-1];
    n1 = mesh->ed_v->JA[ica-1];
    n2 = mesh->ed_v->JA[ica];
    elen = mesh->ed_len[dof[i]-1];

    // Find out which linear basis elements line up with nodes on this edge
    for (k=0; k<v_per_elm; k++) {
      if (v_on_elm[k]==n1) {
        mark1=k;
      }
      if (v_on_elm[k]==n2) {
        mark2=k;
      }
    }
    // Make sure orientation is correct always go from i->j if nj > ni
    if (MAX(n1,n2)==n1) {
      ihi = mark1;
      ilo = mark2;
    } else {
      ihi = mark2;
      ilo = mark1;
    }

    phi[i*dim+0] = elen*(p[ilo]*dp[ihi*dim] - p[ihi]*dp[ilo*dim]);
    phi[i*dim+1] = elen*(p[ilo]*dp[ihi*dim+1] - p[ihi]*dp[ilo*dim+1]);
    if(dim==3) phi[i*dim+2] = elen*(p[ilo]*dp[ihi*dim+2] - p[ihi]*dp[ilo*dim+2]);

    /* Now compute Curls
     * In 2D curl v = (-dy,dx)*(v1,v2)^T = (dx,dy)(0 1;-1 0)(v1,v2)^T = div (Jv)
     * curl phi_eij = |eij|*(grad(p(i))*(J*grad(p(j)))-grad(p(j))*(J*grad(p(i)))
     * This results from the fact that the p's are linear...
     *
     * In 3D, technically the curls are not needed in 3D as <curl u, curl v> operator can be found from Laplacian matrix.
     * We compute them anyway
     */

    if(dim==2) {
      cphi[i] = 2*elen*(dp[ilo*dim]*dp[ihi*dim+1] - dp[ilo*dim+1]*dp[ihi*dim]);
    } else if(dim==3) {
      cphi[i*dim+0] = 2*elen*(dp[ilo*dim+1]*dp[ihi*dim+2]-dp[ihi*dim+1]*dp[ilo*dim+2]);
      cphi[i*dim+1] = 2*elen*(dp[ihi*dim]*dp[ilo*dim+2]-dp[ilo*dim]*dp[ihi*dim+2]);
      cphi[i*dim+2] = 2*elen*(dp[ilo*dim]*dp[ihi*dim+1]-dp[ihi*dim]*dp[ilo*dim+1]);
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }
  }

  if(p) free(p);
  if(dp) free(dp);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
 *
 * \brief Compute Raviart-Thomas Finite Element Basis Functions (zeroth order) at a particular point in 2 or 3D
 *
 * \param x         Coordinate on where to compute basis function
 * \param v_on_elm  Vertices on element
 * \param dof       DOF on element
 * \param mesh      Mesh struct
 *
 * \return phi      Basis functions (dim for each face from reference triangle)
 * \return dphi     Div of basis functions (1 for each face)
 *
 */
void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
{
  // Flag for erros
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;

  INT i,j,ica,icb,jcnt;
  REAL* p;
  REAL* dp;
  INT* ipf = (INT *) calloc(dim,sizeof(INT));
  INT myf;
  REAL farea;
  INT elnd,ef1,ef2,ef3;

  /* Get Linear Basis Functions for particular element */
  p = (REAL *) calloc(v_per_elm,sizeof(REAL));
  dp = (REAL *) calloc(v_per_elm*dim,sizeof(REAL));
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = mesh->f_v->IA[myf-1]-1;
      icb = mesh->f_v->IA[myf]-1;
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area and normal vector of the face
      farea = mesh->f_area[myf-1];

      // Loop through Nodes on element to find corresponding nodes and get correct orienation
      for(j=0;j<v_per_elm;j++) {
        elnd = v_on_elm[j];
        if(ipf[0]==elnd) {
          ef1 = j;
        }
        if(ipf[1]==elnd) {
          ef2 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 2D
       * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i)))
       * |fij| = |eij|
       */
      phi[i*dim+0] = farea*(p[ef1]*dp[ef2*dim+1] - p[ef2]*dp[ef1*dim+1]);
      phi[i*dim+1] = farea*(-p[ef1]*dp[ef2*dim] + p[ef2]*dp[ef1*dim]);

      // Compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i)))
      dphi[i] = 2*farea*(dp[ef1*dim]*dp[ef2*dim+1] - dp[ef2*dim]*dp[ef1*dim+1]);
    }
  } else if(dim==3) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = mesh->f_v->IA[myf-1]-1;
      icb = mesh->f_v->IA[myf]-1;
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area
      farea = mesh->f_area[myf-1];

      // Loop through Nodes on element to find corresponding nodes for correct orienation
      for(j=0;j<v_per_elm;j++) {
        elnd = v_on_elm[j];
        if(ipf[0]==elnd) {
          ef1 = j;
        }
        if(ipf[1]==elnd) {
          ef2 = j;
        }
        if(ipf[2]==elnd) {
          ef3 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 3D
       * phi_fijk = 6*|fijk|*(p(i)(grad(p(j)) x grad(p(k))) - p(j)(grad(p(k)) x grad(p(i))) + p(k)(grad(p(i)) x grad(p(j))))
       * |fijk| = Area(Face)
       */
      phi[i*dim+0] = 2*farea*(p[ef1]*(dp[ef2*dim+1]*dp[ef3*dim+2]-dp[ef2*dim+2]*dp[ef3*dim+1])
          + p[ef2]*(dp[ef3*dim+1]*dp[ef1*dim+2]-dp[ef3*dim+2]*dp[ef1*dim+1])
          + p[ef3]*(dp[ef1*dim+1]*dp[ef2*dim+2]-dp[ef1*dim+2]*dp[ef2*dim+1]));
      phi[i*dim+1] = 2*farea*(p[ef1]*(dp[ef2*dim+2]*dp[ef3*dim]-dp[ef2*dim]*dp[ef3*dim+2])
          + p[ef2]*(dp[ef3*dim+2]*dp[ef1*dim]-dp[ef3*dim]*dp[ef1*dim+2])
          + p[ef3]*(dp[ef1*dim+2]*dp[ef2*dim]-dp[ef1*dim]*dp[ef2*dim+2]));
      phi[i*dim+2] = 2*farea*(p[ef1]*(dp[ef2*dim]*dp[ef3*dim+1]-dp[ef2*dim+1]*dp[ef3*dim])
          + p[ef2]*(dp[ef3*dim]*dp[ef1*dim+1]-dp[ef3*dim+1]*dp[ef1*dim])
          + p[ef3]*(dp[ef1*dim]*dp[ef2*dim+1]-dp[ef1*dim+1]*dp[ef2*dim]));

      // Compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i)))
      dphi[i] = 6*farea*(dp[ef1*dim]*(dp[ef2*dim+1]*dp[ef3*dim+2]-dp[ef2*dim+2]*dp[ef3*dim+1])
          + dp[ef1*dim+1]*(dp[ef2*dim+2]*dp[ef3*dim]-dp[ef2*dim]*dp[ef3*dim+2])
          + dp[ef1*dim+2]*(dp[ef2*dim]*dp[ef3*dim+1]-dp[ef2*dim+1]*dp[ef3*dim]));
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(p) free(p);
  if(dp) free(dp);
  if(ipf) free(ipf);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void bdm1_basis(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
 *
 * \brief Brezzi-Douglas-Marini (BDM) Elements of order 1.
 *
 * \note ONLY in 2D for now.
 * \note This has NOT been tested.
 *
 * \param x         Coordinate on where to compute basis function
 * \param v_on_elm  Vertices on element
 * \param dof       DOF on element
 * \param mesh      Mesh struct
 *
 * \return phi      Basis functions (2*f_per_elm for each face, 12 total in 2D)
 * \return dphix    Div of basis functions
 *
 */
void bdm1_basis(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;

  INT i,j,ica,icb,jcnt;
  REAL a1,a2,a3,a4;
  REAL* p;
  REAL* dp;
  INT* ipf = (INT *) calloc(dim,sizeof(INT));
  INT myf;
  REAL farea;
  INT elnd,ef1,ef2;

  /* Get Linear Basis Functions for particular element */
  p = (REAL *) calloc(v_per_elm,sizeof(REAL));
  dp = (REAL *) calloc(v_per_elm*dim,sizeof(REAL));
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = mesh->f_v->IA[myf-1]-1;
      icb = mesh->f_v->IA[myf]-1;
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area and normal vector of the face
      farea = mesh->f_area[myf-1];

      // Loop through Nodes on element to find corresponding nodes and get correct orienation
      for(j=0;j<v_per_elm;j++) {
        elnd = v_on_elm[j];
        if(ipf[0]==elnd) {
          ef1 = j;
        }
        if(ipf[1]==elnd) {
          ef2 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the BDM1 elements are in 2D
       * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i)))
       * psi_fij = alpha*|fij|*curl(p(i)p(j))
       * |fij| = |eij|
       */
      phi[i*dim*2] = farea*(p[ef1]*dp[ef2*dim+1] - p[ef2]*dp[ef1*dim+1]);
      phi[i*dim*2+1] = farea*(-p[ef1]*dp[ef2*dim] + p[ef2]*dp[ef1*dim]);
      phi[i*dim*2+2] = -6*farea*(p[ef1]*dp[ef2*dim+1] + p[ef2]*dp[ef2*dim+1]);
      phi[i*dim*2+3] = 6*farea*(p[ef1]*dp[ef2*dim] + p[ef2]*dp[ef2*dim]);

      a1 = dp[ef1*dim]*dp[ef2*dim];
      a2 = dp[ef1*dim]*dp[ef2*dim+1];
      a3 = dp[ef1*dim+1]*dp[ef2*dim];
      a4 = dp[ef1*dim+1]*dp[ef2*dim+1];

      dphix[i*dim*2] = farea*(a2-a3);
      dphix[i*dim*2+1] = 0.0;
      dphix[i*dim*2+2] = -6*farea*(a2+a3);
      dphix[i*dim*2+3] = 12*farea*a1;
      dphiy[i*dim*2] = 0.0;
      dphiy[i*dim*2+1] = farea*(a2-a3);
      dphiy[i*dim*2+2] = -12*farea*a4;
      dphiy[i*dim*2+3] = 6*farea*(a2+a3);
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__); // 3D not implemented
  }

  if(p) free(p);
  if(dp) free(dp);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void bubble_face_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh)
 *
 * \brief Compute Bubble Element Finite Element Basis Functions at a particular point in 2 or 3D
 *
 * \param x         Coordinate on where to compute basis function
 * \param v_on_elm  Vertices on element
 * \param dof       DOF on element
 * \param mesh      Mesh struct
 *
 * \return phi      Basis functions (dim for each face from reference triangle)
 * \return dphi     Tensor from gradient of basis functions
 */
void bubble_face_basis(REAL *phi, REAL *dphi, REAL *x, INT *v_on_elm, INT *dof, trimesh *mesh)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT dof_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;
  INT i,j;

  /* Get Linear Basis Functions for particular element */
  REAL* p = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* dp = (REAL *) calloc(v_per_elm*dim,sizeof(REAL));
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // face to vertex map
  INT* fv = (INT *)calloc(dim,sizeof(INT));
  INT elnd,ef1,ef2,ef3;//face endpoint vertex tracking numbers

  REAL gradp; 

  if(dim==2){
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i]-1,mesh->f_v,fv);
      // Find orientation of face
      for(j=0;j<v_per_elm;j++){
        elnd = v_on_elm[j];
        if(fv[0]==elnd) {
          ef1 = j;
        }
        if(fv[1]==elnd) {
          ef2 = j;
        }
      }
      
      // Multiply basis function by normal vector
      phi[i*dim] = mesh->f_norm[dim*(dof[i]-1)] * 4*p[ef1]*p[ef2];
      phi[i*dim+1] = mesh->f_norm[dim*(dof[i]-1)+1] * 4*p[ef1]*p[ef2];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 4*(p[ef1]*dp[ef2*dim+j] + dp[ef1*dim+j]*p[ef2]);
      
        dphi[i*dim*dim + j*dim + 0] = gradp * mesh->f_norm[dim*(dof[i]-1)+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * mesh->f_norm[dim*(dof[i]-1)+1];
      }

    }
  } else if(dim==3) {
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i]-1,mesh->f_v,fv);
      // Find orientation of face
      for(j=0;j<v_per_elm;j++){
        elnd = v_on_elm[j];
        if(fv[0]==elnd) {
          ef1 = j;
        }
        if(fv[1]==elnd) {
          ef2 = j;
        }
        if(fv[2]==elnd) {
          ef3 = j;
        }
      }
      
      // Multiply basis function by normal vector
      phi[i*dim] = mesh->f_norm[dim*(dof[i]-1)] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+1] = mesh->f_norm[dim*(dof[i]-1)+1] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+2] = mesh->f_norm[dim*(dof[i]-1)+2] * 8*p[ef1]*p[ef2]*p[ef3];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 8*(p[ef1]*p[ef2]*dp[ef3*dim+j] + p[ef1]*dp[ef2*dim+j]*p[ef3] + dp[ef1*dim+j]*p[ef2]*p[ef3]);
      
        dphi[i*dim*dim + j*dim + 0] = gradp * mesh->f_norm[dim*(dof[i]-1)+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * mesh->f_norm[dim*(dof[i]-1)+1];
        dphi[i*dim*dim + j*dim + 2] = gradp * mesh->f_norm[dim*(dof[i]-1)+2];
      }

    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }


  if(p) free(p);
  if(dp) free(dp);
  if(fv) free(fv);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
 * \fn void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh,fespace *FE)
 *
 * \brief Grabs the basis function of a FEM space at a particular point in 2 or 3D
 *
 * \param x         Coordinate on where to compute basis function
 * \param v_on_elm  Vertices on element
 * \param dof       DOF on element
 * \param mesh      Mesh struct
 * \param FE        Fespace struct
 *
 * \return phi      Basis functions
 * \return dphi     Derivatives of basis functions (depends on type)
 *
 */
void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,trimesh *mesh,fespace *FE)
{
  // Flag for erros
  SHORT status;

  // Mesh and FEM Data
  INT FEtype = FE->FEtype;
  INT offset = FE->dof_per_elm/mesh->dim;

  if(FEtype>=0 && FEtype<10) { // PX elements

    PX_H1_basis(phi,dphi,x,dof,FEtype,mesh);

  } else if(FEtype==20) { // Nedelec elements

    ned_basis(phi,dphi,x,v_on_elm,dof,mesh);

  } else if(FEtype==30) { // Raviart-Thomas elements

    rt_basis(phi,dphi,x,v_on_elm,dof,mesh);

  } else if(FEtype==60) { // Vector element

    PX_H1_basis(phi, dphi, x, dof, 1, mesh);
    PX_H1_basis(phi+offset, dphi+FE->dof_per_elm, x, dof, 1, mesh);

  } else if(FEtype==61) { // Bubble element

    bubble_face_basis(phi,dphi,x,v_on_elm,dof,mesh);

  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
/****************************************************************************************************************************/
