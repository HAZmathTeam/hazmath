/*! \file src/fem/quadrature.c
 *
 * \brief Computes quadrature nodes and weights for each element or for entire domain
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/11/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/14/2016
 * \note Updated on 11/3/2018 for 0-1 fix.
 *
 */

#include "hazmath.h"

/*!
 * \fn struct qcoords *alloc_quadrature(INT nq1d,INT nsimplex,INT dim)
 *
 * \brief Allocates memory and properties of quadrature struct qcoords
 *
 * \param nq1d       Number of quadrature nodes on an simplex in 1D direction
 * \param nsimplex   Number of simplicies to get quadrature on
 * \param dim        Dimension of problem
 *
 * \return Q     Quadrature struct
 *
 * \note Will replace older routines based on old struct
 *
 */
struct quadrature *alloc_quadrature(INT nq1d,INT nsimplex,INT dim)
{

  struct quadrature *Q = malloc(sizeof(struct quadrature));
  assert(Q != NULL);

  INT nq = (INT) pow(nq1d,dim);
  INT nq_tot = nq*nsimplex;
  Q->dim = dim;
  Q->nq1d = nq1d;
  Q->nq_simplex = nq;
  Q->nq = nq_tot;
  Q->x = (REAL *) calloc(nq_tot*dim,sizeof(REAL));
  Q->w = (REAL *) calloc(nq_tot,sizeof(REAL));

  return Q;
}
/******************************************************************************/

/*!
 * \fn struct qcoords *alloc_quadrature_bdry(INT nq1d,INT nregion,INT dim,INT ed_or_f)
 *
 * \brief Allocates memory and properties of quadrature coordinates struct.
 *        Assumes we are allocated on a boundary, so dimension is 1 or 2 less
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param nregion    Number of "elements" (faces or edges) to get quadrature on
 * \param dim     Dimension of problem
 * \param ed_or_f Whether we are computing quadrature on faces or edges
 *
 * \return Q     Quadrature struct
 *
 * \note Will replace older routines based on old struct
 *
 */
struct quadrature *alloc_quadrature_bdry(INT nq1d,INT nregion,INT dim,INT ed_or_f)
{
  // Flag for errors
  SHORT status;

  struct quadrature *Q = malloc(sizeof(struct quadrature));
  assert(Q != NULL);

  INT nq = 0;

  switch (ed_or_f)
  {
    case 1:
    // Compute quadrature on edges
    nq = nq1d;
    break;

    case 2:
    nq = (INT) pow(nq1d,dim-1);
    break;

    default:
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  INT nq_tot = nq*nregion;
  Q->dim = dim;
  Q->nq1d = nq1d;
  Q->nq_simplex = nq;
  Q->nq = nq_tot;
  Q->x = (REAL *) calloc(nq_tot*dim,sizeof(REAL));
  Q->w = (REAL *) calloc(nq_tot,sizeof(REAL));

  return Q;
}
/******************************************************************************/

/*!
 * \fn void free_quadrature(quadrature* Q)
 *
 * \brief Frees memory of arrays of quadrature struct
 *
 * \return Q  truct for quadratures to be freed
 *
 */
void free_quadrature(quadrature* Q)
{
  if (Q==NULL) return;

  if(Q->x) {
    free(Q->x);
    Q->x = NULL;
  }

  if(Q->w) {
    free(Q->w);
    Q->w = NULL;
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void quad_refelm(quadrature *cqelm,INT nq1d,INT dim)
 *
 * \brief Computes quadrature weights and nodes for the reference element using nq1d^(dim)
 *        quadrature nodes on simplex
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param dim     Dimension of problem
 * \param elm     Index of current element
 *
 * \return cq_ref Quadrature struct on reference element
 *
 */
void quad_refelm(quadrature *cqelm,INT nq1d,INT dim)
{
  // Flag for errors
  SHORT status;

  /* Loop indices */
  INT q;

  /* Total Number of Quadrature Nodes */
  INT nq = pow(nq1d,dim);

  /* Gaussian points for reference element */
  REAL* gp = (REAL *) calloc(dim*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  /* Points on Reference Triangle */
  REAL r[dim];

  // Get coordinates of vertices for given element
  if(dim==1) {

    // Get Quad Nodes and Weights
    quad1d(gp,gw,nq1d);

    // Map to Real Interval
    // x = 0.5*(x2-x1)*r + 0.5*(x1+x2)
    // w = 0.5*(x2-x1)*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      cqelm->x[q] = 0.5*(1+r[0]);
      cqelm->w[q] = 0.5*gw[q];
    }
  } else if(dim==2) {

    // Get Quad Nodes and Weights
    triquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s) + x2*r + x3*s
    // y = y1*(1-r-s) + y2*r + y3*s
    // w = 2*Element Area*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      r[1] = gp[nq+q];
      cqelm->x[q*dim+0] = r[0];
      cqelm->x[q*dim+1] = r[1];
      cqelm->w[q] = gw[q];
    }
  } else if(dim==3) {

    // Get Quad Nodes and Weights
    tetquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s-t) + x2*r + x3*s + x4*t
    // y = y1*(1-r-s-t) + y2*r + y3*s + y4*t
    // z = z1*(1-r-s-t) + z2*r + z3*s + z4*t3*s
    // w = 6*Element Vol*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      r[1] = gp[nq+q];
      r[2] = gp[2*nq+q];
      cqelm->x[q*dim+0] = r[0];
      cqelm->x[q*dim+1] = r[1];
      cqelm->x[q*dim+2] = r[2];
      cqelm->w[q] = gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);

  return;
}
/*****************************************************************************/

/*!
 * \fn void quad_elm_local(quadrature *cqelm,simplex_local_data *loc_data,INT nq1d)
 *
 * \brief Computes quadrature weights and nodes for SINGLE element using nq1d^(dim)
 *        quadrature nodes on simplex
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param loc_data Mesh data on current simplex (element)
 *
 * \return cq_elm Quadrature struct on element
 *
 */
void quad_elm_local(quadrature *cqelm,simplex_local_data *loc_data,INT nq1d)
{
  // Flag for errors
  SHORT status;

  /* Loop indices */
  INT q,j;

  /* Dimension */
  INT dim = loc_data->dim;

  /* Total Number of Quadrature Nodes */
  INT nq = pow(nq1d,dim);

  /* Coordinates of vertices of element */
  REAL* xv = loc_data->xv;

  /* Gaussian points for reference element */
  REAL* gp = (REAL *) calloc(dim*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  /* Points on Reference Triangle */
  REAL r[dim];

  REAL e_vol = loc_data->vol; /* Area/Volume of Element */
  REAL voldim=0.0;

  // only 1D, 2D, and 3D quadrature is implemented
  if(dim==1) {

    // Scale the quadrature weight appropriately
    voldim = 0.5*e_vol;

    // Get Quad Nodes and Weights on ref element
    quad1d(gp,gw,nq1d);

    // Map to Real Interval
    // x = 0.5*(x2-x1)*r + 0.5*(x1+x2)
    // w = 0.5*(x2-x1)*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      cqelm->x[q] = 0.5*(xv[0]*(1-r[0]) + xv[1]*(1+r[0]));
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==2) {
    voldim = 2.0*e_vol;

    // Get Quad Nodes and Weights
    triquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s) + x2*r + x3*s
    // y = y1*(1-r-s) + y2*r + y3*s
    // w = 2*Element Area*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      r[1] = gp[nq+q];
      for(j=0;j<dim;j++) cqelm->x[q*dim+j] = xv[0*dim+j]*(1-r[0]-r[1]) + xv[1*dim+j]*r[0]+ xv[2*dim+j]*r[1];
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==3) {

    voldim = 6.0*e_vol;

    // Get Quad Nodes and Weights
    tetquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s-t) + x2*r + x3*s + x4*t
    // y = y1*(1-r-s-t) + y2*r + y3*s + y4*t
    // z = z1*(1-r-s-t) + z2*r + z3*s + z4*t3*s
    // w = 6*Element Vol*wref
    for (q=0; q<nq; q++) {
      r[0] = gp[q];
      r[1] = gp[nq+q];
      r[2] = gp[2*nq+q];
      for(j=0;j<dim;j++) cqelm->x[q*dim+j] = xv[0*dim+j]*(1-r[0]-r[1]-r[2]) + xv[1*dim+j]*r[0]+ xv[2*dim+j]*r[1] + xv[3*dim+j]*r[2];
      cqelm->w[q] = voldim*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);

  return;
}
/******************************************************************************/

/*!
* \fn void quad_edge_local(quadrature *cqedge,simplex_local_data *loc_data,INT nq1d,INT edge)
*
* \brief Computes quadrature weights and nodes for SINGLE Edge using nq1d
*          quadrature nodes on a line/surface.  Can be used to compute integrals on
*          1D boundaries (curves).
*
* \param nq1d    Number of quadrature nodes on an edge in 1D direction
* \param loc_data    Mesh data on element that contains edge
* \param edge     Index of current edge
*
* \return cqedge Quadrature struct on edge
*
* \note An edge integral is always a 1D integral.
*
*/
void quad_edge_local(quadrature *cqedge,simplex_local_data *loc_data,INT nq1d,INT edge)
{

  INT q,i,j,k; /* Loop indices */
  INT dim = loc_data->dim;

  // Test for Simpson's Rule
  INT nqdum=1;
  if(nq1d==-1) {
    nqdum = -1;
    nq1d = 3;
  }

  // Get total quadrature points
  INT nq = nq1d;

  /* Coordinates of vertices of edge/face */
  INT v_per_ed = 2;
  INT v_per_elm = dim+1;
  INT* v_on_ed = loc_data->v_on_ed + edge*v_per_ed;
  INT* v_on_elm = loc_data->local_v;
  REAL* xved = (REAL *) calloc(v_per_ed*dim,sizeof(REAL));
  for(i=0;i<v_per_ed;i++) {
    for(j=0;j<v_per_elm;j++) {
      if(v_on_ed[i]==v_on_elm[j]) {
        for(k=0;i<dim;k++) {
          xved[i*dim+k] = loc_data->xv[j*dim+k];
        }
      }
    }
  }

  // Gaussian points for reference element
  // (edges: [-1,1]
  REAL* gp = (REAL *) calloc(nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r;      	/* Points on Reference Edge */

  REAL w = 0.5*loc_data->ed_len[edge]; /* Jacobian = 1/2 |e| */

  // Get Quad Nodes and Weights
  if(nqdum==-1) { // Simpsons Rule
    gp[0] = -1.0;
    gp[1] = 0.0;
    gp[2] = 1.0;
    gw[0] = 1.0/3.0;
    gw[1] = 4.0/3.0;
    gw[2] = gw[0];
  } else {
    quad1d(gp,gw,nq1d);
  }

  // Map to Real Edge
  // Edges: x = 0.5(x1*(1-r) + x2*(1+r))
  //        y = 0.5(y1*(1-r) + y2*(1+r))
  //        z = 0.5(z1*(1-r) + z2*(1+r))
  //        w = 0.5*Edge Length*wref
  for (q=0; q<nq1d; q++) {
    r = gp[q];
    for(k=0;k<dim;k++) cqedge->x[q*dim+k] = 0.5*xved[0*dim+k]*(1-r) + 0.5*xved[1*dim+k]*(1+r);
    cqedge->w[q] = w*gw[q];
  }

  if(gp) free(gp);
  if(gw) free(gw);
  if(xved) free(xved);

  return;
}
/******************************************************************************/

/*!
* \fn void quad_face(quadrature *cqface,simplex_local_data *loc_data,INT nq1d,INT face)
*
* \brief Computes quadrature weights and nodes for SINGLE Face using nq1d^(dim-1)
*          quadrature nodes on a line/surface.  Can be used to compute integrals on
*          1D/2D boundaries (curves/surfaces).
*
* \param nq1d    Number of quadrature nodes on a face in 1D direction
* \param loc_data    Mesh data on element that contains face
* \param face     Index of current face
*
* \return cq_bdry Quadrature struct on face
*
* \note A face integral is a 1D integral in 2D and a 2D integral in 3D.
*
*/
void quad_face_local(quadrature *cqface,simplex_local_data *loc_data,INT nq1d,INT face)
{
  // Flag for errors
  SHORT status;

  INT q,i,j,k; /* Loop indices */
  INT dim = loc_data->dim;
  // Face is an edge in 2D
  INT nq = (INT) pow(nq1d,dim-1);

  /* Coordinates of vertices of face */
  INT v_per_face = dim;
  INT v_per_elm = dim+1;
  INT* v_on_face = loc_data->v_on_f + face*v_per_face;
  INT* v_on_elm = loc_data->local_v;
  REAL* xvf = (REAL *) calloc(v_per_face*dim,sizeof(REAL));
  for(i=0;i<v_per_face;i++) {
    for(j=0;j<v_per_elm;j++) {
      if(v_on_face[i]==v_on_elm[j]) {
        for(k=0;i<dim;k++) {
          xvf[i*dim+k] = loc_data->xv[j*dim+k];
        }
      }
    }
  }

  // Gaussian points for reference element
  // (edges: [-1,1] faces: tri[(0,0),(1,0),(0,1)])
  REAL* gp = (REAL *) calloc((dim-1)*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r,s;      	/* Points on Reference Face */

  REAL w = 0.0;
  if(dim==2) { // Faces are Edges in 2D
    w = 0.5*loc_data->f_area[face]; /* Jacobian = 1/2 |e| */
  } else if(dim==3) {
    w = 2*loc_data->f_area[face]; /* Jacobian = 2*|f| */
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get Quad Nodes and Weights
  if(dim==2) { // face is an edge
    quad1d(gp,gw,nq1d);
  } else if(dim==3) { // face is a face
    triquad_(gp,gw,nq1d);
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Map to Real Face
  if(dim==2) {
    // Edges: x = 0.5(x1*(1-r) + x2*(1+r))
    //        y = 0.5(y1*(1-r) + y2*(1+r))
    //        z = 0.5(z1*(1-r) + z2*(1+r))
    //        w = 0.5*Edge Length*wref
    for (q=0; q<nq1d; q++) {
      r = gp[q];
      for(k=0;k<dim;k++) cqface->x[q*dim+k] = 0.5*xvf[0*dim+k]*(1-r) + 0.5*xvf[1*dim+k]*(1-r);
      cqface->w[q] = w*gw[q];
    }
  } else if(dim==3) {
    // Faces: x = x1*(1-r-s) + x2*r + x3*s
    //        y = y1*(1-r-s) + y2*r + y3*s
    //        z = z1*(1-r-s) + z2*r + z3*s
    //        w = 2*Element Area*wref
    for (q=0; q<nq1d*nq1d; q++) {
      r = gp[q];
      s = gp[nq1d*nq1d+q];
      for(k=0;k<dim;k++) cqface->x[q*dim+k] = xvf[0*dim+j]*(1-r-s) + xvf[1*dim+j]*r + xvf[2*dim+j]*s;
      cqface->w[q] = w*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  if(xvf) free(xvf);

  return;
}
/******************************************************************************/

/*********************************************************************/
/**** Integration Routines for a given function ********/
/*********************************************************************/

/*!
* \fn REAL integrate_elm_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *elm_data,REAL time)
*
* \brief Integrate a given scalar function over an element
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param cq     Precomputed quadrature points and weights on the element
* \param loc_data  Local Mesh Information on element
* \param time   If needed for function
*
* \return integral Integral of scalar function over element
*
* \note If cq is given, we will just use these precomputed values
*       Otherwise, we will use those allocated by elm_data
*/
REAL integrate_elm_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *elm_data,REAL time)
{

  // Loop indices
  INT quad;
  INT dim = elm_data->dim;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = NULL;//(REAL *) calloc(dim,sizeof(REAL));

  INT nq_per_elm;

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    nq_per_elm = cq->nq_simplex;
    for (quad=0;quad<nq_per_elm;quad++) {
      qx = cq->x + quad*dim;
      w = cq->w[quad];
      (*expr)(uval,qx,dim,time,&(elm_data->flag));
      integral += w*uval[comp];
    }
  } else { // use quadrature saved on the element
    quadrature *cqelm = elm_data->quad_local;
    nq_per_elm = cqelm->nq_simplex;
    for (quad=0;nq_per_elm;quad++) {
      qx = cqelm->x + quad*dim;
      w = cqelm->w[quad];
      (*expr)(uval,qx,dim,time,&(elm_data->flag));
      integral += w*uval[comp];
    }
  }

  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_face_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *face_data,REAL time)
*
* \brief Integrate a given scalar function over a face
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param cq     Precomputed quadrature points and weights on the face
* \param loc_data  Local Mesh Information on face
* \param time   If needed for function
*
* \return integral Integral of scalar function over element
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will use the one stored in face_data
*/
REAL integrate_face_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *face_data,REAL time)
{

  // Loop indices
  INT quad;
  INT dim = face_data->dim;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = NULL;

  INT nq_per_face;

  // Quadrature on face
  if(cq) { // assuming quadrature is given
    nq_per_face = cq->nq_simplex;
    for (quad=0;quad<nq_per_face;quad++) {
      qx = cq->x + quad*dim;
      w = cq->w[quad];
      (*expr)(uval,qx,time,dim,&(face_data->flag));
      integral += w*uval[comp];
    }
  } else { // assemble quadrature again
    quadrature *cqface = face_data->quad_local;
    nq_per_face = cqface->nq_simplex;
    for (quad=0;nq_per_face;quad++) {
      qx = cqface->x + quad*dim;
      w = cqface->w[quad];
      (*expr)(uval,qx,time,dim,&(face_data->flag));
      integral += w*uval[comp];
    }
  }

  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_edge_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *edge_data,REAL time)
*
* \brief Integrate a given scalar function over an edge
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param cq     Precomputed quadrature points and weights on the edge
* \param loc_data  Local Mesh Information on edge
* \param time   If needed for function
*
* \return integral Integral of scalar function over edge
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will use the one stored in edge_data
*/
REAL integrate_edge_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *edge_data,REAL time)
{

  // Loop indices
  INT quad;
  INT dim = edge_data->dim;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = NULL;

  INT nq_per_edge;

  // Quadrature on face
  if(cq) { // assuming quadrature is given
    nq_per_edge = cq->nq_simplex;
    for (quad=0;quad<nq_per_edge;quad++) {
      qx = cq->x + quad*dim;
      w = cq->w[quad];
      (*expr)(uval,qx,time,dim,&(edge_data->flag));
      integral += w*uval[comp];
    }
  } else { // assemble quadrature again
    quadrature *cqedge = edge_data->quad_local;
    nq_per_edge = cqedge->nq_simplex;
    for (quad=0;nq_per_edge;quad++) {
      qx = cqedge->x + quad*dim;
      w = cqedge->w[quad];
      (*expr)(uval,qx,time,dim,&(edge_data->flag));
      integral += w*uval[comp];
    }
  }

  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_edge_vector_tangent_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *edge_data,REAL time)
*
* \brief Integrate the tangential component of a given vector function along an edge
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate (this is the start of the vector)
* \param cq     Precomputed quadrature points and weights on the edge
* \param loc_data  Local Mesh Information on edge
* \param time   If needed for function
*
* \return integral Integral of vector function along edge
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will use the one stored in edge_data
*
* \note This is a tangential line integral in any dimension
*
*/
REAL integrate_edge_vector_tangent_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *edge_data,REAL time)
{

  // Loop indices
  INT quad,j;
  INT dim = edge_data->dim;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = NULL;

  INT nq_per_edge;

  // Quadrature on face
  if(cq) { // assuming quadrature is given
    nq_per_edge = cq->nq_simplex;
    for (quad=0;quad<nq_per_edge;quad++) {
      qx = cq->x + quad*dim;
      w = cq->w[quad];
      (*expr)(uval,qx,time,dim,&(edge_data->flag));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*edge_data->ed_tau[j];
    }
  } else { // assemble quadrature again
    quadrature *cqedge = edge_data->quad_local;
    nq_per_edge = cqedge->nq_simplex;
    for (quad=0;nq_per_edge;quad++) {
      qx = cqedge->x + quad*dim;
      w = cqedge->w[quad];
      (*expr)(uval,qx,time,dim,&(edge_data->flag));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*edge_data->ed_tau[j];
      integral += w*uval[comp];
    }
  }

  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_face_vector_normal_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *face_data,REAL time)
*
* \brief Integrate the normal component of a given vector function across a face
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param cq     Precomputed quadrature points and weights on the face
* \param loc_data  Local Mesh Information on face
* \param time   If needed for function
*
* \return integral Integral of vector function across the face
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will use the one stored in face_data
*/
REAL integrate_face_vector_normal_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,quadrature *cq,simplex_local_data *face_data,REAL time)
{

  // Loop indices
  INT quad,j;
  INT dim = face_data->dim;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = NULL;

  INT nq_per_face;

  // Quadrature on face
  if(cq) { // assuming quadrature is given
    nq_per_face = cq->nq_simplex;
    for (quad=0;quad<nq_per_face;quad++) {
      qx = cq->x + quad*dim;
      w = cq->w[quad];
      (*expr)(uval,qx,time,dim,&(face_data->flag));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*face_data->f_norm[dim+j];
    }
  } else { // assemble quadrature again
    quadrature *cqface = face_data->quad_local;
    nq_per_face = cqface->nq_simplex;
    for (quad=0;nq_per_face;quad++) {
      qx = cqface->x + quad*dim;
      w = cqface->w[quad];
      (*expr)(uval,qx,time,dim,&(face_data->flag));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*face_data->f_norm[dim+j];
    }
  }

  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_domain_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,INT nq1d,quadrature *cq,simplex_local_data *elm_data,mesh_struct *mesh,REAL time)
*
* \brief Integrate a given scalar function over the entire mesh
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights (if given - can be NULL)
* \param elm_data Local Mesh Data
* \param mesh   Mesh Information
* \param time   If needed for function
*
* \return integral Integral of scalar function over domain
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*/
REAL integrate_domain_local(void (*expr)(REAL *,REAL *,INT,REAL,void *),INT nun,INT comp,INT nq1d,quadrature *cq,simplex_local_data *elm_data,mesh_struct *mesh,REAL time)
{
  // Loop indices
  INT elm;

  // Integral to return
  REAL integral = 0.0;
  // Loop over all elements and call integrate_elm
  for(elm=0;elm<mesh->nelm;elm++) {
    get_elmlocaldata(elm_data,mesh,elm);
    integral += integrate_elm_local(expr,nun,comp,cq,elm_data,time);
  }

  return integral;
}
/******************************************************************************/

/*********************************************************************/
/**** Tables of Quadrature Nodes and Weights 1D and 2D and 3D ********/
/*********************************************************************/

/************************************************************************************/
/*!
* \fn void quad1d(REAL *gaussp, REAL *gaussc, INT ng1d)
*
* \brief 1D Quadrature on Reference Element [-1,1]
*
* \note Up to 5 Gaussian points on the reference domain
*
* \param ng1d            Number of Gaussian points in 1 direction
*
* \return gaussp         x coordinates of the Gaussian points
* \return gaussc         Weights of the Gaussian points
*
*/
void quad1d(REAL *gaussp, REAL *gaussc, INT ng1d)
{
  // Check for errors
  if(ng1d<1) {
    ng1d = 1;
    printf("HEY!  You better have at least 1 quadrature point...\n");
    printf("Forcing you to use 1 quadrature point.\n\n");
  }
  if(ng1d>5) {
    ng1d = 5;
    printf("Getting greedy with the quadrature aren't we??\n  Knocking you down to 5 quadrature points...\n");
  }

  switch (ng1d) {
    case 1:
    gaussp[0] = 0.0;
    gaussc[0] = 2.0;
    break;
    case 2:
    gaussp[0] =     -0.57735026918962576450915;
    gaussp[1] =      0.57735026918962576450915;
    gaussc[0] =       1.0000000000000000000000;
    gaussc[1] =       1.0000000000000000000000;
    break;

  default:// case 3:
    gaussp[  0]=     -0.77459666924148337703585;
    gaussp[  1]=                            0.0;
    gaussp[  2]=      0.77459666924148337703585;
    gaussc[  0]=      0.55555555555555555555556;
    gaussc[  1]=      0.88888888888888888888889;
    gaussc[  2]=      0.55555555555555555555556;
    break;
  case 4:
    gaussp[  0]=     -0.86113631159405257522395;
    gaussp[  1]=     -0.33998104358485626480267;
    gaussp[  2]=      0.33998104358485626480267;
    gaussp[  3]=      0.86113631159405257522395;
    gaussc[  0]=      0.34785484513745385315771;
    gaussc[  1]=      0.65214515486254614656693;
    gaussc[  2]=      0.65214515486254614289181;
    gaussc[  3]=      0.34785484513745385738355;
    break;
  case 5:
    gaussp[  0]=     -0.90617984593866399279763;
    gaussp[  1]=     -0.53846931010568309103631;
    gaussp[  2]=                            0.0;
    gaussp[  3]=      0.53846931010568309103631;
    gaussp[  4]=      0.90617984593866399279763;
    gaussc[  0]=      0.23692688505618908751426;
    gaussc[  1]=      0.47862867049936646804129;
    gaussc[  2]=      0.56888888888888888903082;
    gaussc[  3]=      0.47862867049936646790566;
    gaussc[  4]=      0.23692688505618908750796;
    break;
  }
  return;
}
/************************************************************************************/

/************************************************************************************/
/*!
* \fn void triquad_(REAL *gp, REAL *gc, INT ng1d)
*
* \brief 2D Quadrature on Reference Triangle ( Vertices (0,0),(1,0),(0,1) )
*
* \note Up to 7x7 Gaussian points on a triangle.
*
* \param ng1d        Number of Gaussian points in 1 direction
*
* \return gp         (x,y) coordinates of the Gaussian points
* \return gc         Weights of the Gaussian points
*
*/
void triquad_(REAL *gp, REAL *gc, INT ng1d)
{
  // Check for errors
  if(ng1d<1) {
    ng1d = 1;
    printf("HEY!  You better have at least 1 quadrature point...\n");
    printf("Forcing you to use 1 quadrature point.\n\n");
  }
  if(ng1d>7) {
    ng1d = 7;
    printf("Getting greedy with the quadrature aren't we??\n  Knocking you down to 7x7 quadrature points...\n");
  }

  INT ngauss = ng1d*ng1d;
  switch (ng1d) {
  case 1:
    gp[0]=   0.3333333333333333333333333E+00;
    gp[ngauss+0]=   0.3333333333333333333333333E+00;
    gc[0]=     0.5000000000000000000000000E+00;
    break;
  default:
    gp[0]=   0.7503111022260811817747560E-01;
    gp[ngauss+0]=   0.6449489742783178098197284E+00;
    gc[0]=   0.9097930912801141530281550E-01;
    gp[1]=   0.2800199154990740720027960E+00;
    gp[ngauss+1]=   0.6449489742783178098197284E+00;
    gc[1]=   0.9097930912801141530281550E-01;
    gp[2]=   0.1785587282636164231170351E+00;
    gp[ngauss+2]=   0.1550510257216821901802716E+00;
    gc[2]=   0.1590206908719885846971845E+00;
    gp[3]=   0.6663902460147013867026933E+00;
    gp[ngauss+3]=   0.1550510257216821901802716E+00;
    gc[3]=   0.1590206908719885846971845E+00;
    break;
    case 3:
    gp[0]=   0.2393113228708061896738275E-01;
    gp[ngauss+0]=   0.7876594617608470560252419E+00;
    gc[0]=   0.1939638330595947848163369E-01;
    gp[1]=   0.1061702691195764719873791E+00;
    gp[ngauss+1]=   0.7876594617608470560252419E+00;
    gc[1]=   0.3103421328953516557061390E-01;
    gp[2]=   0.1884094059520723250073754E+00;
    gp[ngauss+2]=   0.7876594617608470560252419E+00;
    gc[2]=   0.1939638330595947848163369E-01;
    gp[3]=   0.6655406783916450412865057E-01;
    gp[ngauss+3]=   0.4094668644407347108649263E+00;
    gc[3]=   0.6367808509988506852608905E-01;
    gp[4]=   0.2952665677796326445675369E+00;
    gp[ngauss+4]=   0.4094668644407347108649263E+00;
    gc[4]=   0.1018849361598161096417425E+00;
    gp[5]=   0.5239790677201007850064232E+00;
    gp[ngauss+5]=   0.4094668644407347108649263E+00;
    gc[5]=   0.6367808509988506852608905E-01;
    gp[6]=   0.1027176548096262680160926E+00;
    gp[ngauss+6]=   0.8858795951270394739554614E-01;
    gc[6]=   0.5581442048304434188116615E-01;
    gp[7]=   0.4557060202436480263022269E+00;
    gp[ngauss+7]=   0.8858795951270394739554614E-01;
    gc[7]=   0.8930307277287094700986584E-01;
    gp[8]=   0.8086943856776697845883612E+00;
    gp[ngauss+8]=   0.8858795951270394739554614E-01;
    gc[8]=   0.5581442048304434188116615E-01;
    break;
    case 4:
    gp[0]=   0.9703785126946112175961724E-02;
    gp[ngauss+0]=   0.8602401356562194478479129E+00;
    gc[0]=   0.5423225910525254453528332E-02;
    gp[1]=   0.4612207990645204861943693E-01;
    gp[ngauss+1]=   0.8602401356562194478479129E+00;
    gc[1]=   0.1016725956447878663340924E-01;
    gp[2]=   0.9363778443732850353265015E-01;
    gp[ngauss+2]=   0.8602401356562194478479129E+00;
    gc[2]=   0.1016725956447878663340924E-01;
    gp[3]=   0.1300560792168344399761254E+00;
    gp[ngauss+3]=   0.8602401356562194478479129E+00;
    gc[3]=   0.5423225910525254453528332E-02;
    gp[4]=   0.2891208422438901271682135E-01;
    gp[ngauss+4]=   0.5835904323689168200566977E+00;
    gc[4]=   0.2258404928236993135842135E-01;
    gp[5]=   0.1374191041345743684268067E+00;
    gp[ngauss+5]=   0.5835904323689168200566977E+00;
    gc[5]=   0.4233972452174628905480675E-01;
    gp[6]=   0.2789904634965088115164956E+00;
    gp[ngauss+6]=   0.5835904323689168200566977E+00;
    gc[6]=   0.4233972452174628905480675E-01;
    gp[7]=   0.3874974834066941672264810E+00;
    gp[ngauss+7]=   0.5835904323689168200566977E+00;
    gc[7]=   0.2258404928236993135842135E-01;
    gp[8]=   0.5021012321136977210504381E-01;
    gp[ngauss+8]=   0.2768430136381238276800460E+00;
    gc[8]=   0.3538806789808594616389450E-01;
    gp[9]=   0.2386486597314429209541046E+00;
    gp[ngauss+9]=   0.2768430136381238276800460E+00;
    gc[9]=   0.6634421610704973423180774E-01;
    gp[10]=   0.4845083266304332513658494E+00;
    gp[ngauss+10]=   0.2768430136381238276800460E+00;
    gc[10]=   0.6634421610704973423180774E-01;
    gp[11]=   0.6729468631505064002149102E+00;
    gp[ngauss+11]=   0.2768430136381238276800460E+00;
    gc[11]=   0.3538806789808594616389450E-01;
    gp[12]=   0.6546699455501446386445480E-01;
    gp[ngauss+12]=   0.5710419611451768219312119E-01;
    gc[12]=   0.2356836819338233236742181E-01;
    gp[13]=   0.3111645522443570344411343E+00;
    gp[ngauss+13]=   0.5710419611451768219312119E-01;
    gc[13]=   0.4418508852236172573671028E-01;
    gp[14]=   0.6317312516411252833657445E+00;
    gp[ngauss+14]=   0.5710419611451768219312119E-01;
    gc[14]=   0.4418508852236172573671028E-01;
    gp[15]=   0.8774288093304678539424240E+00;
    gp[ngauss+15]=   0.5710419611451768219312119E-01;
    gc[15]=   0.2356836819338233236742181E-01;
    break;
  case 5:
    gp[0]=   0.4622288465046428525209780E-02;
    gp[ngauss+0]=   0.9014649142011735738765011E+00;
    gc[0]=   0.1865552166877838436183754E-02;
    gp[1]=   0.2273848306376403459813202E-01;
    gp[ngauss+1]=   0.9014649142011735738765011E+00;
    gc[1]=   0.3768701695327620376776386E-02;
    gp[2]=   0.4926754289941321306174945E-01;
    gp[ngauss+2]=   0.9014649142011735738765011E+00;
    gc[2]=   0.4479406797281358559372037E-02;
    gp[3]=   0.7579660273506239152536687E-01;
    gp[ngauss+3]=   0.9014649142011735738765011E+00;
    gc[3]=   0.3768701695327620376776386E-02;
    gp[4]=   0.9391279733377999759828912E-01;
    gp[ngauss+4]=   0.9014649142011735738765011E+00;
    gc[4]=   0.1865552166877838436183754E-02;
    gp[5]=   0.1428579439557138533782080E-01;
    gp[ngauss+5]=   0.6954642733536360945146148E+00;
    gc[5]=   0.8755499182163831736919008E-02;
    gp[6]=   0.7027629200828172118339338E-01;
    gp[ngauss+6]=   0.6954642733536360945146148E+00;
    gc[6]=   0.1768745211048346587741866E-01;
    gp[7]=   0.1522678633231819527426926E+00;
    gp[ngauss+7]=   0.6954642733536360945146148E+00;
    gc[7]=   0.2102296748732207512195789E-01;
    gp[8]=   0.2342594346380821843019918E+00;
    gp[ngauss+8]=   0.6954642733536360945146148E+00;
    gc[8]=   0.1768745211048346587741866E-01;
    gp[9]=   0.2902499322507925201475644E+00;
    gp[ngauss+9]=   0.6954642733536360945146148E+00;
    gc[9]=   0.8755499182163831736919008E-02;
    gp[10]=   0.2636464494447091747928250E-01;
    gp[ngauss+10]=   0.4379748102473861440050125E+00;
    gc[10]=   0.1734150643136570012831098E-01;
    gp[11]=   0.1296959367822541214837999E+00;
    gp[ngauss+11]=   0.4379748102473861440050125E+00;
    gc[11]=   0.3503250450337172031698701E-01;
    gp[12]=   0.2810125948763069279974937E+00;
    gp[ngauss+12]=   0.4379748102473861440050125E+00;
    gc[12]=   0.4163896521519496780744188E-01;
    gp[13]=   0.4323292529703597345111876E+00;
    gp[ngauss+13]=   0.4379748102473861440050125E+00;
    gc[13]=   0.3503250450337172031698701E-01;
    gp[14]=   0.5356605448081429385157050E+00;
    gp[ngauss+14]=   0.4379748102473861440050125E+00;
    gc[14]=   0.1734150643136570012831098E-01;
    gp[15]=   0.3762125234511119174727219E-01;
    gp[ngauss+15]=   0.1980134178736081725357921E+00;
    gc[15]=   0.1980408313204735378039807E-01;
    gp[16]=   0.1850707102673894331855753E+00;
    gp[ngauss+16]=   0.1980134178736081725357921E+00;
    gc[16]=   0.4000728738616042409551254E-01;
    gp[17]=   0.4009932910631959137321039E+00;
    gp[ngauss+17]=   0.1980134178736081725357921E+00;
    gc[17]=   0.4755189705795400973985440E-01;
    gp[18]=   0.6169158718590023942786325E+00;
    gp[ngauss+18]=   0.1980134178736081725357921E+00;
    gc[18]=   0.4000728738616042409551254E-01;
    gp[19]=   0.7643653297812806357169357E+00;
    gp[ngauss+19]=   0.1980134178736081725357921E+00;
    gc[19]=   0.1980408313204735378039807E-01;
    gp[20]=   0.4504259356980372309546899E-01;
    gp[ngauss+20]=   0.3980985705146874234080669E-01;
    gc[20]=   0.1146508035159254779675419E-01;
    gp[21]=   0.2215786095523792017723069E+00;
    gp[ngauss+21]=   0.3980985705146874234080669E-01;
    gc[21]=   0.2316122192949838634362829E-01;
    gp[22]=   0.4800950714742656288295967E+00;
    gp[ngauss+22]=   0.3980985705146874234080669E-01;
    gc[22]=   0.2752898566446981099359601E-01;
    gp[23]=   0.7386115333961520558868864E+00;
    gp[ngauss+23]=   0.3980985705146874234080669E-01;
    gc[23]=   0.2316122192949838634362829E-01;
    gp[24]=   0.9151475493787275345637243E+00;
    gp[ngauss+24]=   0.3980985705146874234080669E-01;
    gc[24]=   0.1146508035159254779675419E-01;
    break;
  case 6:
    gp[0]=   0.2466697152670243054005080E-02;
    gp[ngauss+0]=   0.9269456713197411148518740E+00;
    gc[0]=   0.7485425612363183140950521E-03;
    gp[1]=   0.1237506041744003817266408E-01;
    gp[ngauss+1]=   0.9269456713197411148518740E+00;
    gc[1]=   0.1576221754023588582963184E-02;
    gp[2]=   0.2781108211536058069826929E-01;
    gp[ngauss+2]=   0.9269456713197411148518740E+00;
    gc[2]=   0.2044386591544858980950415E-02;
    gp[3]=   0.4524324656489830444985675E-01;
    gp[ngauss+3]=   0.9269456713197411148518740E+00;
    gc[3]=   0.2044386591544858980950415E-02;
    gp[4]=   0.6067926826281884697546196E-01;
    gp[ngauss+4]=   0.9269456713197411148518740E+00;
    gc[4]=   0.1576221754023588582963184E-02;
    gp[5]=   0.7058763152758864209412095E-01;
    gp[ngauss+5]=   0.9269456713197411148518740E+00;
    gc[5]=   0.7485425612363183140950521E-03;
    gp[6]=   0.7791874701286432033793818E-02;
    gp[ngauss+6]=   0.7692338620300545009168834E+00;
    gc[6]=   0.3765298212691672929234314E-02;
    gp[7]=   0.3909070073282424404541472E-01;
    gp[ngauss+7]=   0.7692338620300545009168834E+00;
    gc[7]=   0.7928667333796484710025645E-02;
    gp[8]=   0.8785045497599719116592015E-01;
    gp[ngauss+8]=   0.7692338620300545009168834E+00;
    gc[8]=   0.1028361722876633011482835E-01;
    gp[9]=   0.1429156829939483079171965E+00;
    gp[ngauss+9]=   0.7692338620300545009168834E+00;
    gc[9]=   0.1028361722876633011482835E-01;
    gp[10]=   0.1916754372371212550377019E+00;
    gp[ngauss+10]=   0.7692338620300545009168834E+00;
    gc[10]=   0.7928667333796484710025645E-02;
    gp[11]=   0.2229742632686590670493228E+00;
    gp[ngauss+11]=   0.7692338620300545009168834E+00;
    gc[11]=   0.3765298212691672929234314E-02;
    gp[12]=   0.1490156336667116035714823E-01;
    gp[ngauss+12]=   0.5586715187715501320813933E+00;
    gc[12]=   0.8451535796943121648933723E-02;
    gp[13]=   0.7475897346264909767772818E-01;
    gp[ngauss+13]=   0.5586715187715501320813933E+00;
    gc[13]=   0.1779657599702627725499296E-01;
    gp[14]=   0.1680095191211918575326299E+00;
    gp[ngauss+14]=   0.5586715187715501320813933E+00;
    gc[14]=   0.2308246365135823315636558E-01;
    gp[15]=   0.2733189621072580103859768E+00;
    gp[ngauss+15]=   0.5586715187715501320813933E+00;
    gc[15]=   0.2308246365135823315636558E-01;
    gp[16]=   0.3665695077658007702408785E+00;
    gp[ngauss+16]=   0.5586715187715501320813933E+00;
    gc[16]=   0.1779657599702627725499296E-01;
    gp[17]=   0.4264269178617787075614584E+00;
    gp[ngauss+17]=   0.5586715187715501320813933E+00;
    gc[17]=   0.8451535796943121648933723E-02;
    gp[18]=   0.2238687297803063445050099E-01;
    gp[ngauss+18]=   0.3369846902811542990970530E+00;
    gc[18]=   0.1206060640426510907696060E-01;
    gp[19]=   0.1123116817809536957220250E+00;
    gp[ngauss+19]=   0.3369846902811542990970530E+00;
    gc[19]=   0.2539627158904765582035844E-01;
    gp[20]=   0.2524035680765180133752919E+00;
    gp[ngauss+20]=   0.3369846902811542990970530E+00;
    gc[20]=   0.3293939890078669916221938E-01;
    gp[21]=   0.4106117416423276875276552E+00;
    gp[ngauss+21]=   0.3369846902811542990970530E+00;
    gc[21]=   0.3293939890078669916221938E-01;
    gp[22]=   0.5507036279378920051809220E+00;
    gp[ngauss+22]=   0.3369846902811542990970530E+00;
    gc[22]=   0.2539627158904765582035844E-01;
    gp[23]=   0.6406284367408150664524460E+00;
    gp[ngauss+23]=   0.3369846902811542990970530E+00;
    gc[23]=   0.1206060640426510907696060E-01;
    gp[24]=   0.2876533301255912843698113E-01;
    gp[ngauss+24]=   0.1480785996684842918499769E+00;
    gc[24]=   0.1161087476699751443083611E-01;
    gp[25]=   0.1443114869504166464557392E+00;
    gp[ngauss+25]=   0.1480785996684842918499769E+00;
    gc[25]=   0.2444926225805781423674754E-01;
    gp[26]=   0.3243183045887760364106504E+00;
    gp[ngauss+26]=   0.1480785996684842918499769E+00;
    gc[26]=   0.3171111159070397975276155E-01;
    gp[27]=   0.5276030957427396717393727E+00;
    gp[ngauss+27]=   0.1480785996684842918499769E+00;
    gc[27]=   0.3171111159070397975276155E-01;
    gp[28]=   0.7076099133810990616942839E+00;
    gp[ngauss+28]=   0.1480785996684842918499769E+00;
    gc[28]=   0.2444926225805781423674754E-01;
    gp[29]=   0.8231560673189565797130420E+00;
    gp[ngauss+29]=   0.1480785996684842918499769E+00;
    gc[29]=   0.1161087476699751443083611E-01;
    gp[30]=   0.3277536661445989520154516E-01;
    gp[ngauss+30]=   0.2931642715978489197205028E-01;
    gc[30]=   0.6194265352658849860014235E-02;
    gp[31]=   0.1644292415948274481657064E+00;
    gp[ngauss+31]=   0.2931642715978489197205028E-01;
    gc[31]=   0.1304339433008283128737061E-01;
    gp[32]=   0.3695299243723766991833510E+00;
    gp[ngauss+32]=   0.2931642715978489197205028E-01;
    gc[32]=   0.1691750568001266068034231E-01;
    gp[33]=   0.6011536484678384088445987E+00;
    gp[ngauss+33]=   0.2931642715978489197205028E-01;
    gc[33]=   0.1691750568001266068034231E-01;
    gp[34]=   0.8062543312453876598622433E+00;
    gp[ngauss+34]=   0.2931642715978489197205028E-01;
    gc[34]=   0.1304339433008283128737061E-01;
    gp[35]=   0.9379082062257552128264046E+00;
    gp[ngauss+35]=   0.2931642715978489197205028E-01;
    gc[35]=   0.6194265352658849860014235E-02;
    break;
  case 7:
    gp[0]=   0.1431659581332948445688999E-02;
    gp[ngauss+0]=   0.9437374394630778535343478E+00;
    gc[0]=   0.3375907567113747844459523E-03;
    gp[1]=   0.7271058658560282492949164E-02;
    gp[ngauss+1]=   0.9437374394630778535343478E+00;
    gc[1]=   0.7292426106515660112115446E-03;
    gp[2]=   0.1671433656946750295425480E-01;
    gp[ngauss+2]=   0.9437374394630778535343478E+00;
    gc[2]=   0.9955000916249671892041120E-03;
    gp[3]=   0.2813128026846107323282610E-01;
    gp[ngauss+3]=   0.9437374394630778535343478E+00;
    gc[3]=   0.1089695284831588119968157E-02;
    gp[4]=   0.3954822396745464351139739E-01;
    gp[ngauss+4]=   0.9437374394630778535343478E+00;
    gc[4]=   0.9955000916249671892041120E-03;
    gp[5]=   0.4899150187836186397270303E-01;
    gp[ngauss+5]=   0.9437374394630778535343478E+00;
    gc[5]=   0.7292426106515660112115446E-03;
    gp[6]=   0.5483090095558919801996319E-01;
    gp[ngauss+6]=   0.9437374394630778535343478E+00;
    gc[6]=   0.3375907567113747844459523E-03;
    gp[7]=   0.4586412541637882763079511E-02;
    gp[ngauss+7]=   0.8197593082631076350124201E+00;
    gc[7]=   0.1774485071438049608442658E-02;
    gp[8]=   0.2329329894998979644828858E-01;
    gp[ngauss+8]=   0.8197593082631076350124201E+00;
    gc[8]=   0.3833132573484684075609282E-02;
    gp[9]=   0.5354544045728325221282914E-01;
    gp[ngauss+9]=   0.8197593082631076350124201E+00;
    gc[9]=   0.5232667115687632726379401E-02;
    gp[10]=   0.9012034586844618249378997E-01;
    gp[ngauss+10]=   0.8197593082631076350124201E+00;
    gc[10]=   0.5727787200652742623468366E-02;
    gp[11]=   0.1266952512796091127747508E+00;
    gp[ngauss+11]=   0.8197593082631076350124201E+00;
    gc[11]=   0.5232667115687632726379401E-02;
    gp[12]=   0.1569473927869025685392914E+00;
    gp[ngauss+12]=   0.8197593082631076350124201E+00;
    gc[12]=   0.3833132573484684075609282E-02;
    gp[13]=   0.1756542791952544822245004E+00;
    gp[ngauss+13]=   0.8197593082631076350124201E+00;
    gc[13]=   0.1774485071438049608442658E-02;
    gp[14]=   0.8972904006716703697492974E-02;
    gp[ngauss+14]=   0.6473752828868303626260922E+00;
    gc[14]=   0.4297910087982423247056434E-02;
    gp[15]=   0.4557124628029494113490578E-01;
    gp[ngauss+15]=   0.6473752828868303626260922E+00;
    gc[15]=   0.9284078756888546352613518E-02;
    gp[16]=   0.1047568427084817262927626E+00;
    gp[ngauss+16]=   0.6473752828868303626260922E+00;
    gc[16]=   0.1267383600209279954896470E-01;
    gp[17]=   0.1763123585565848186869539E+00;
    gp[ngauss+17]=   0.6473752828868303626260922E+00;
    gc[17]=   0.1387304677156393168637868E-01;
    gp[18]=   0.2478678744046879110811452E+00;
    gp[ngauss+18]=   0.6473752828868303626260922E+00;
    gc[18]=   0.1267383600209279954896470E-01;
    gp[19]=   0.3070534708328746962390020E+00;
    gp[ngauss+19]=   0.6473752828868303626260922E+00;
    gc[19]=   0.9284078756888546352613518E-02;
    gp[20]=   0.3436518131064529336764148E+00;
    gp[ngauss+20]=   0.6473752828868303626260922E+00;
    gc[20]=   0.4297910087982423247056434E-02;
    gp[21]=   0.1392289515659608599518407E-01;
    gp[ngauss+21]=   0.4528463736694446169985514E+00;
    gc[21]=   0.6935542753734072742362389E-02;
    gp[22]=   0.7071107454632530338118497E-01;
    gp[ngauss+22]=   0.4528463736694446169985514E+00;
    gc[22]=   0.1498172921938941357190799E-01;
    gp[23]=   0.1625469900128696646167676E+00;
    gp[ngauss+23]=   0.4528463736694446169985514E+00;
    gc[23]=   0.2045178462250981418417614E-01;
    gp[24]=   0.2735768131652776915007243E+00;
    gp[ngauss+24]=   0.4528463736694446169985514E+00;
    gc[24]=   0.2238695250460706899401921E-01;
    gp[25]=   0.3846066363176857183846809E+00;
    gp[ngauss+25]=   0.4528463736694446169985514E+00;
    gc[25]=   0.2045178462250981418417614E-01;
    gp[26]=   0.4764425517842300796202636E+00;
    gp[ngauss+26]=   0.4528463736694446169985514E+00;
    gc[26]=   0.1498172921938941357190799E-01;
    gp[27]=   0.5332307311739592970062645E+00;
    gp[ngauss+27]=   0.4528463736694446169985514E+00;
    gc[27]=   0.6935542753734072742362389E-02;
    gp[28]=   0.1868274434884273534596982E-01;
    gp[ngauss+28]=   0.2657898227845894684767894E+00;
    gc[28]=   0.8247603013529574038759331E-02;
    gp[29]=   0.9488521701286283095347225E-01;
    gp[ngauss+29]=   0.2657898227845894684767894E+00;
    gc[29]=   0.1781596040067579543578056E-01;
    gp[30]=   0.2181172683502983220175448E+00;
    gp[ngauss+30]=   0.2657898227845894684767894E+00;
    gc[30]=   0.2432083637489711574352285E-01;
    gp[31]=   0.3671050886077052657616053E+00;
    gp[ngauss+31]=   0.2657898227845894684767894E+00;
    gc[31]=   0.2662209772138335648260900E-01;
    gp[32]=   0.5160929088651122095056658E+00;
    gp[ngauss+32]=   0.2657898227845894684767894E+00;
    gc[32]=   0.2432083637489711574352285E-01;
    gp[33]=   0.6393249602025477005697384E+00;
    gp[ngauss+33]=   0.2657898227845894684767894E+00;
    gc[33]=   0.1781596040067579543578056E-01;
    gp[34]=   0.7155274328665677961772408E+00;
    gp[ngauss+34]=   0.2657898227845894684767894E+00;
    gc[34]=   0.8247603013529574038759331E-02;
    gp[35]=   0.2252791561566364109969225E-01;
    gp[ngauss+35]=   0.1146790531609042319096402E+00;
    gc[35]=   0.7154643779096141969509067E-02;
    gp[36]=   0.1144139277467613129097319E+00;
    gp[ngauss+36]=   0.1146790531609042319096402E+00;
    gc[36]=   0.1545501766273406746103162E-01;
    gp[37]=   0.2630088665758011781230588E+00;
    gp[ngauss+37]=   0.1146790531609042319096402E+00;
    gc[37]=   0.2109787781815243944540462E-01;
    gp[38]=   0.4426604734195478840451799E+00;
    gp[ngauss+38]=   0.1146790531609042319096402E+00;
    gc[38]=   0.2309417967090930466923013E-01;
    gp[39]=   0.6223120802632945899673009E+00;
    gp[ngauss+39]=   0.1146790531609042319096402E+00;
    gc[39]=   0.2109787781815243944540462E-01;
    gp[40]=   0.7709070190923344551806278E+00;
    gp[ngauss+40]=   0.1146790531609042319096402E+00;
    gc[40]=   0.1545501766273406746103162E-01;
    gp[41]=   0.8627930312234321269906675E+00;
    gp[ngauss+41]=   0.1146790531609042319096402E+00;
    gc[41]=   0.7154643779096141969509067E-02;
    gp[42]=   0.2487403237606075687067163E-01;
    gp[ngauss+42]=   0.2247938643871249810882550E-01;
    gc[42]=   0.3623466079725786927077026E-02;
    gp[43]=   0.1263292970196692449335864E+00;
    gp[ngauss+43]=   0.2247938643871249810882550E-01;
    gc[43]=   0.7827186648495094067212433E-02;
    gp[44]=   0.2903993060879903088904502E+00;
    gp[ngauss+44]=   0.2247938643871249810882550E-01;
    gc[44]=   0.1068501060131496739994062E-01;
    gp[45]=   0.4887603067806437509455873E+00;
    gp[ngauss+45]=   0.2247938643871249810882550E-01;
    gc[45]=   0.1169603676441935436310196E-01;
    gp[46]=   0.6871213074732971930007243E+00;
    gp[ngauss+46]=   0.2247938643871249810882550E-01;
    gc[46]=   0.1068501060131496739994062E-01;
    gp[47]=   0.8511913165416182569575881E+00;
    gp[ngauss+47]=   0.2247938643871249810882550E-01;
    gc[47]=   0.7827186648495094067212433E-02;
    gp[48]=   0.9526465811852267450205029E+00;
    gp[ngauss+48]=   0.2247938643871249810882550E-01;
    gc[48]=   0.3623466079725786927077026E-02;
    break;
  }
  return;
}
/************************************************************************************/

/************************************************************************************/
/*!
* \fn void tetquad_(REAL *gp, REAL *gc, INT ng1d)
*
* \brief 3D Quadrature on Reference Tetrahedron ( Vertices (0,0,0),(0,1,0),(1,0,0),(0,0,1) )
*
* \note Up to 5x5x5 Gaussian points on a tetrahedron.
*
* \param ng1d        Number of Gaussian points in 1 direction
*
* \return gp         (x,y) coordinates of the Gaussian points
* \return gc         Weights of the Gaussian points
*
*/
void tetquad_(REAL *gp, REAL *gc, INT ng1d)
{
  // Check for errors
  if(ng1d<1) {
    ng1d = 1;
    printf("HEY!  You better have at least 1 quadrature point...\n");
    printf("Forcing you to use 1 quadrature point.\n\n");
  }
  if(ng1d>5) {
    ng1d = 5;
    printf("Getting greedy with the quadrature aren't we??\n  Knocking you down to 5x5x5 quadrature points...\n");
  }

  INT ngauss=ng1d*ng1d*ng1d;
  switch (ng1d) {
  case 1:
    gp[0]=    0.2500000000000000000000000E+00;
    gp[ngauss+  0]=   0.2500000000000000000000000E+00;
    gp[2*ngauss+  0]=   0.2500000000000000000000000E+00;
    gc[  0]=   0.1666666666666666666666667E+00;
    break;
  default:
    gp[  0]=   0.3420279323676641430060446E-01;
    gp[ngauss+  0]=   0.2939988006316228658907916E+00;
    gp[2*ngauss+  0]=   0.5441518440112252887999262E+00;
    gc[  0]=   0.9169429921479743922682354E-02;
    gp[  1]=   0.1276465621203854310086777E+00;
    gp[ngauss+  1]=   0.2939988006316228658907916E+00;
    gp[2*ngauss+  1]=   0.5441518440112252887999262E+00;
    gc[  1]=   0.9169429921479743922682354E-02;
    gp[  2]=   0.8139566701467025507670959E-01;
    gp[ngauss+  2]=   0.7067972415939690306926744E-01;
    gp[2*ngauss+  2]=   0.5441518440112252887999262E+00;
    gc[  2]=   0.1602704059847661372315674E-01;
    gp[  3]=   0.3037727648147075530540967E+00;
    gp[ngauss+  3]=   0.7067972415939690306926744E-01;
    gp[2*ngauss+  3]=   0.5441518440112252887999262E+00;
    gc[  3]=   0.1602704059847661372315674E-01;
    gp[  4]=   0.6583868706004440993602967E-01;
    gp[ngauss+  4]=   0.5659331650728008805355130E+00;
    gp[2*ngauss+  4]=   0.1225148226554413778667404E+00;
    gc[  4]=   0.2115700645452406117825615E-01;
    gp[  5]=   0.2457133252117133316617169E+00;
    gp[ngauss+  5]=   0.5659331650728008805355130E+00;
    gp[2*ngauss+  5]=   0.1225148226554413778667404E+00;
    gc[  5]=   0.2115700645452406117825615E-01;
    gp[  6]=   0.1566826373368183090793373E+00;
    gp[ngauss+  6]=   0.1360549768028460171710947E+00;
    gp[2*ngauss+  6]=   0.1225148226554413778667404E+00;
    gc[  6]=   0.3697985635885291450923809E-01;
    gp[  7]=   0.5847475632048942958828276E+00;
    gp[ngauss+  7]=   0.1360549768028460171710947E+00;
    gp[2*ngauss+  7]=   0.1225148226554413778667404E+00;
    gc[  7]=   0.3697985635885291450923809E-01;
    break;
  case      3:
    gp[  0]=   0.7059631139554788090931258E-02;
    gp[ngauss+  0]=   0.2323578005798646935907283E+00;
    gp[2*ngauss+  0]=   0.7050022098884983831223985E+00;
    gc[  0]=   0.5809353158373849777015838E-03;
    gp[  1]=   0.3131999476581846164343659E-01;
    gp[ngauss+  1]=   0.2323578005798646935907283E+00;
    gp[2*ngauss+  1]=   0.7050022098884983831223985E+00;
    gc[  1]=   0.9294965053398159643225341E-03;
    gp[  2]=   0.5558035839208213519594193E-01;
    gp[ngauss+  2]=   0.2323578005798646935907283E+00;
    gp[2*ngauss+  2]=   0.7050022098884983831223985E+00;
    gc[  2]=   0.5809353158373849777015838E-03;
    gp[  3]=   0.1963330293548449033824621E-01;
    gp[ngauss+  3]=   0.1207918201339025431243858E+00;
    gp[2*ngauss+  3]=   0.7050022098884983831223985E+00;
    gc[  3]=   0.1907203414981785439750756E-02;
    gp[  4]=   0.8710298498879953687660784E-01;
    gp[ngauss+  4]=   0.1207918201339025431243858E+00;
    gp[2*ngauss+  4]=   0.7050022098884983831223985E+00;
    gc[  4]=   0.3051525463970856703601209E-02;
    gp[  5]=   0.1545726670421145834149695E+00;
    gp[ngauss+  5]=   0.1207918201339025431243858E+00;
    gp[2*ngauss+  5]=   0.7050022098884983831223985E+00;
    gc[  5]=   0.1907203414981785439750756E-02;
    gp[  6]=   0.3030148117427580438384368E-01;
    gp[ngauss+  6]=   0.2613325228673484212751636E-01;
    gp[2*ngauss+  6]=   0.7050022098884983831223985E+00;
    gc[  6]=   0.1671681131483704306363840E-02;
    gp[  7]=   0.1344322689123833873750426E+00;
    gp[ngauss+  7]=   0.2613325228673484212751636E-01;
    gp[2*ngauss+  7]=   0.7050022098884983831223985E+00;
    gc[  7]=   0.2674689810373926890182144E-02;
    gp[  8]=   0.2385630566504909703662415E+00;
    gp[ngauss+  8]=   0.2613325228673484212751636E-01;
    gp[2*ngauss+  8]=   0.7050022098884983831223985E+00;
    gc[  8]=   0.1671681131483704306363840E-02;
    gp[ 9]=   0.1562693925790164701335452E-01;
    gp[ngauss+ 9]=   0.5143386621740919113570116E+00;
    gp[2*ngauss+ 9]=   0.3470037660383518847217635E+00;
    gc[ 9]=   0.2836648695630920164099731E-02;
    gp[ 10]=   0.6932878589377810196061241E-01;
    gp[ngauss+ 10]=   0.5143386621740919113570116E+00;
    gp[2*ngauss+ 10]=   0.3470037660383518847217635E+00;
    gc[ 10]=   0.4538637913009472262559569E-02;
    gp[ 11]=   0.1230306325296545569078703E+00;
    gp[ngauss+ 11]=   0.5143386621740919113570116E+00;
    gp[2*ngauss+ 11]=   0.3470037660383518847217635E+00;
    gc[ 11]=   0.2836648695630920164099731E-02;
    gp[ 12]=   0.4345955565380246496495420E-01;
    gp[ngauss+ 12]=   0.2673803204118844564054627E+00;
    gp[2*ngauss+ 12]=   0.3470037660383518847217635E+00;
    gc[ 12]=   0.9312682379470454264458097E-02;
    gp[ 13]=   0.1928079567748818294363869E+00;
    gp[ngauss+ 13]=   0.2673803204118844564054627E+00;
    gp[2*ngauss+ 13]=   0.3470037660383518847217635E+00;
    gc[ 13]=   0.1490029180715272682313296E-01;
    gp[ 14]=   0.3421563578959611939078196E+00;
    gp[ngauss+ 14]=   0.2673803204118844564054627E+00;
    gp[2*ngauss+ 14]=   0.3470037660383518847217635E+00;
    gc[ 14]=   0.9312682379470454264458097E-02;
    gp[ 15]=   0.6707424175205852430583444E-01;
    gp[ngauss+ 15]=   0.5784760393614263759525826E-01;
    gp[2*ngauss+ 15]=   0.3470037660383518847217635E+00;
    gc[ 15]=   0.8162650766546684183041063E-02;
    gp[ 16]=   0.2975743150127527388414891E+00;
    gp[ngauss+ 16]=   0.5784760393614263759525826E-01;
    gp[2*ngauss+ 16]=   0.3470037660383518847217635E+00;
    gc[ 16]=   0.1306024122647469469286570E-01;
    gp[ 17]=   0.5280743882734469533771438E+00;
    gp[ngauss+ 17]=   0.5784760393614263759525826E-01;
    gp[2*ngauss+ 17]=   0.3470037660383518847217635E+00;
    gc[ 17]=   0.8162650766546684183041063E-02;
    gp[ 18]=   0.2218430264081972545955689E-01;
    gp[ngauss+ 18]=   0.7301650280476316250995886E+00;
    gp[2*ngauss+ 18]=   0.7299402407314973215583798E-01;
    gc[ 18]=   0.3047877090518187685409914E-02;
    gp[ 19]=   0.9842047393960932137228673E-01;
    gp[ngauss+ 19]=   0.7301650280476316250995886E+00;
    gp[2*ngauss+ 19]=   0.7299402407314973215583798E-01;
    gc[ 19]=   0.4876603344829100296655863E-02;
    gp[ 20]=   0.1746566452383989172850166E+00;
    gp[ngauss+ 20]=   0.7301650280476316250995886E+00;
    gp[2*ngauss+ 20]=   0.7299402407314973215583798E-01;
    gc[ 20]=   0.3047877090518187685409914E-02;
    gp[ 21]=   0.6169601860914648993801939E-01;
    gp[ngauss+ 21]=   0.3795782302805905833418882E+00;
    gp[2*ngauss+ 21]=   0.7299402407314973215583798E-01;
    gc[ 21]=   0.1000614257217611647115416E-01;
    gp[ 22]=   0.2737138728231298422511369E+00;
    gp[ngauss+ 22]=   0.3795782302805905833418882E+00;
    gp[2*ngauss+ 22]=   0.7299402407314973215583798E-01;
    gc[ 22]=   0.1600982811548178635384666E-01;
    gp[ 23]=   0.4857317270371131945642544E+00;
    gp[ngauss+ 23]=   0.3795782302805905833418882E+00;
    gp[2*ngauss+ 23]=   0.7299402407314973215583798E-01;
    gc[ 23]=   0.1000614257217611647115416E-01;
    gp[ 24]=   0.9521987984171492384049553E-01;
    gp[ngauss+ 24]=   0.8212156786344242164387440E-01;
    gp[2*ngauss+ 24]=   0.7299402407314973215583798E-01;
    gc[ 24]=   0.8770474929651058804317146E-02;
    gp[ 25]=   0.4224422040317039231001438E+00;
    gp[ngauss+ 25]=   0.8212156786344242164387440E-01;
    gp[2*ngauss+ 25]=   0.7299402407314973215583798E-01;
    gc[ 25]=   0.1403275988744169408690743E-01;
    gp[ 26]=   0.7496645282216929223597921E+00;
    gp[ngauss+ 26]=   0.8212156786344242164387440E-01;
    gp[2*ngauss+ 26]=   0.7299402407314973215583798E-01;
    gc[ 26]=   0.8770474929651058804317146E-02;
    break;
  case 4:
    gp[0]=   0.1981013974700432744909435E-02;
    gp[ngauss +  0]=   0.1756168039625049658342796E+00;
    gp[2*ngauss +  0]=   0.7958514178967728633033780E+00;
    gc[0]=   0.5614254026695104258530883E-04;
    gp[1]=   0.9415757216553928623888957E-02;
    gp[ngauss+  1]=   0.1756168039625049658342796E+00;
    gp[2*ngauss+ 1]=   0.7958514178967728633033780E+00;
    gc[  1]=   0.1052539187783914959760637E-03;
    gp[2]=   0.1911602092416824223845346E-01;
    gp[ngauss+  2]=   0.1756168039625049658342796E+00;
    gp[2*ngauss+  2]=   0.7958514178967728633033780E+00;
    gc[  2]=   0.1052539187783914959760637E-03;
    gp[3]=   0.2655076416602173811743298E-01;
    gp[ngauss+  3]=   0.1756168039625049658342796E+00;
    gp[2*ngauss+  3]=   0.7958514178967728633033780E+00;
    gc[  3]=   0.5614254026695104258530883E-04;
    gp[4]=   0.5902361000058098432934300E-02;
    gp[ngauss+  4]=   0.1191391592971236390275109E+00;
    gp[2*ngauss+  4]=   0.7958514178967728633033780E+00;
    gc[  4]=   0.2337955152791078423282521E-03;
    gp[5]=   0.2805391526296907513510520E-01;
    gp[ngauss+  5]=   0.1191391592971236390275109E+00;
    gp[2*ngauss+  5]=   0.7958514178967728633033780E+00;
    gc[  5]=   0.4383110215343271014270785E-03;
    gp[6]=   0.5695550754313442253400591E-01;
    gp[ngauss+  6]=   0.1191391592971236390275109E+00;
    gp[2*ngauss+  6]=   0.7958514178967728633033780E+00;
    gc[  6]=   0.4383110215343271014270785E-03;
    gp[7]=   0.7910706180604539923617681E-01;
    gp[ngauss+  7]=   0.1191391592971236390275109E+00;
    gp[2*ngauss+  7]=   0.7958514178967728633033780E+00;
    gc[  7]=   0.2337955152791078423282521E-03;
    gp[8]=   0.1025032546082947250520215E-01;
    gp[ngauss+  8]=   0.5651710869940735217362115E-01;
    gp[2*ngauss+  8]=   0.7958514178967728633033780E+00;
    gc[  8]=   0.3663457985554326676012287E-03;
    gp[ 9]=   0.4871978550500959094728183E-01;
    gp[ngauss+ 9]=   0.5651710869940735217362115E-01;
    gp[2*ngauss+ 9]=   0.7958514178967728633033780E+00;
    gc[ 9]=   0.6868112975047707265125985E-03;
    gp[ 10]=   0.9891168789881019357571906E-01;
    gp[ngauss+ 10]=   0.5651710869940735217362115E-01;
    gp[2*ngauss+ 10]=   0.7958514178967728633033780E+00;
    gc[ 10]=   0.6868112975047707265125985E-03;
    gp[ 11]=   0.1373811479429903120177987E+00;
    gp[ngauss+ 11]=   0.5651710869940735217362115E-01;
    gp[2*ngauss+ 11]=   0.7958514178967728633033780E+00;
    gc[ 11]=   0.3663457985554326676012287E-03;
    gp[ 12]=   0.1336499411296589418346618E-01;
    gp[ngauss+ 12]=   0.1165774066892339709191637E-01;
    gp[2*ngauss+ 12]=   0.7958514178967728633033780E+00;
    gc[ 12]=   0.2439854216206053939189319E-03;
    gp[ 13]=   0.6352380214147103185255417E-01;
    gp[ngauss+ 13]=   0.1165774066892339709191637E-01;
    gp[2*ngauss+ 13]=   0.7958514178967728633033780E+00;
    gc[ 13]=   0.4574146739399300507499519E-03;
    gp[ 14]=   0.1289670392928327077521515E+00;
    gp[ngauss+ 14]=   0.1165774066892339709191637E-01;
    gp[2*ngauss+ 14]=   0.7958514178967728633033780E+00;
    gc[ 14]=   0.4574146739399300507499519E-03;
    gp[ 15]=   0.1791258473213378454212395E+00;
    gp[ngauss+ 15]=   0.1165774066892339709191637E-01;
    gp[2*ngauss+ 15]=   0.7958514178967728633033780E+00;
    gc[ 15]=   0.2439854216206053939189319E-03;
    gp[ 16]=   0.4686469274784633447665181E-02;
    gp[ngauss+ 16]=   0.4154553003749570180402003E+00;
    gp[2*ngauss+ 16]=   0.5170472951043675023405734E+00;
    gc[ 16]=   0.3722170752562633273986158E-03;
    gp[ 17]=   0.2227478324623351755095983E-01;
    gp[ngauss+ 17]=   0.4154553003749570180402003E+00;
    gp[2*ngauss+ 17]=   0.5170472951043675023405734E+00;
    gc[ 17]=   0.6978185458062600471342006E-03;
    gp[ 18]=   0.4522262127444196206826646E-01;
    gp[ngauss+ 18]=   0.4154553003749570180402003E+00;
    gp[2*ngauss+ 18]=   0.5170472951043675023405734E+00;
    gc[ 18]=   0.6978185458062600471342006E-03;
    gp[ 19]=   0.6281093524589084617156111E-01;
    gp[ngauss+ 19]=   0.4154553003749570180402003E+00;
    gp[2*ngauss+ 19]=   0.5170472951043675023405734E+00;
    gc[ 19]=   0.3722170752562633273986158E-03;
    gp[ 20]=   0.1396316928033901864590850E-01;
    gp[ngauss+ 20]=   0.2818465778637800603501812E+00;
    gp[2*ngauss+ 20]=   0.5170472951043675023405734E+00;
    gc[ 20]=   0.1550031090353912216078243E-02;
    gp[ 21]=   0.6636692804612728658298820E-01;
    gp[ngauss+ 21]=   0.2818465778637800603501812E+00;
    gp[2*ngauss+ 21]=   0.5170472951043675023405734E+00;
    gc[ 21]=   0.2905939875758179218047843E-02;
    gp[ 22]=   0.1347391989857251507262572E+00;
    gp[ngauss+ 22]=   0.2818465778637800603501812E+00;
    gp[2*ngauss+ 22]=   0.5170472951043675023405734E+00;
    gc[ 22]=   0.2905939875758179218047843E-02;
    gp[ 23]=   0.1871429577515134186633369E+00;
    gp[ngauss+ 23]=   0.2818465778637800603501812E+00;
    gp[2*ngauss+ 23]=   0.5170472951043675023405734E+00;
    gc[ 23]=   0.1550031090353912216078243E-02;
    gp[ 24]=   0.2424911481807401304158134E-01;
    gp[ngauss+ 24]=   0.1337020822679903798291838E+00;
    gp[2*ngauss+ 24]=   0.5170472951043675023405734E+00;
    gc[ 24]=   0.2428820659384971875777962E-02;
    gp[ 25]=   0.1152560157370177676747899E+00;
    gp[ngauss+ 25]=   0.1337020822679903798291838E+00;
    gp[2*ngauss+ 25]=   0.5170472951043675023405734E+00;
    gc[ 25]=   0.4553461442867277241397578E-02;
    gp[ 26]=   0.2339946068906243501554529E+00;
    gp[ngauss+ 26]=   0.1337020822679903798291838E+00;
    gp[2*ngauss+ 26]=   0.5170472951043675023405734E+00;
    gc[ 26]=   0.4553461442867277241397578E-02;
    gp[ 27]=   0.3250015078095681047886615E+00;
    gp[ngauss+ 27]=   0.1337020822679903798291838E+00;
    gp[2*ngauss+ 27]=   0.5170472951043675023405734E+00;
    gc[ 27]=   0.2428820659384971875777962E-02;
    gp[ 28]=   0.3161746210173187993001242E-01;
    gp[ngauss+ 28]=   0.2757862597439698206385973E-01;
    gp[2*ngauss+ 28]=   0.5170472951043675023405734E+00;
    gc[ 28]=   0.1617588723434511855713981E-02;
    gp[ 29]=   0.1502777621740505836344576E+00;
    gp[ngauss+ 29]=   0.2757862597439698206385973E-01;
    gp[2*ngauss+ 29]=   0.5170472951043675023405734E+00;
    gc[ 29]=   0.3032594380369393047795663E-02;
    gp[ 30]=   0.3050963167471849319611093E+00;
    gp[ngauss+ 30]=   0.2757862597439698206385973E-01;
    gp[2*ngauss+ 30]=   0.5170472951043675023405734E+00;
    gc[ 30]=   0.3032594380369393047795663E-02;
    gp[ 31]=   0.4237566168195036356655545E+00;
    gp[ngauss+ 31]=   0.2757862597439698206385973E-01;
    gp[2*ngauss+ 31]=   0.5170472951043675023405734E+00;
    gc[ 31]=   0.1617588723434511855713981E-02;
    gp[ 32]=   0.7388454838611978023539041E-02;
    gp[ngauss+ 32]=   0.6549862048169314047901757E+00;
    gp[2*ngauss+ 32]=   0.2386007375518623050589814E+00;
    gc[ 32]=   0.7780094259316945643908710E-03;
    gp[ 33]=   0.3511731762334666143239029E-01;
    gp[ngauss+ 33]=   0.6549862048169314047901757E+00;
    gp[2*ngauss+ 33]=   0.2386007375518623050589814E+00;
    gc[ 33]=   0.1458582752694612457632238E-02;
    gp[ 34]=   0.7129574000785962871845260E-01;
    gp[ngauss+ 34]=   0.6549862048169314047901757E+00;
    gp[2*ngauss+ 34]=   0.2386007375518623050589814E+00;
    gc[ 34]=   0.1458582752694612457632238E-02;
    gp[ 35]=   0.9902460279259431212730384E-01;
    gp[ngauss+ 35]=   0.6549862048169314047901757E+00;
    gp[2*ngauss+ 35]=   0.2386007375518623050589814E+00;
    gc[ 35]=   0.7780094259316945643908710E-03;
    gp[ 36]=   0.2201363960428823146375467E-01;
    gp[ngauss+ 36]=   0.4443453247774830496819952E+00;
    gp[2*ngauss+ 36]=   0.2386007375518623050589814E+00;
    gc[ 36]=   0.3239880378814602487009997E-02;
    gp[ 37]=   0.1046308045343487533720147E+00;
    gp[ngauss+ 37]=   0.4443453247774830496819952E+00;
    gp[2*ngauss+ 37]=   0.2386007375518623050589814E+00;
    gc[ 37]=   0.6074005640321836238018089E-02;
    gp[ 38]=   0.2124231331363058918870087E+00;
    gp[ngauss+ 38]=   0.4443453247774830496819952E+00;
    gp[2*ngauss+ 38]=   0.2386007375518623050589814E+00;
    gc[ 38]=   0.6074005640321836238018089E-02;
    gp[ 39]=   0.2950402980663664137952687E+00;
    gp[ngauss+ 39]=   0.4443453247774830496819952E+00;
    gp[2*ngauss+ 39]=   0.2386007375518623050589814E+00;
    gc[ 39]=   0.3239880378814602487009997E-02;
    gp[ 40]=   0.3822995078056706336853632E-01;
    gp[ngauss+ 40]=   0.2107880663979872074525160E+00;
    gp[2*ngauss+ 40]=   0.2386007375518623050589814E+00;
    gc[ 40]=   0.5076729393991831949508787E-02;
    gp[ 41]=   0.1817069135037572184823919E+00;
    gp[ngauss+ 41]=   0.2107880663979872074525160E+00;
    gp[2*ngauss+ 41]=   0.2386007375518623050589814E+00;
    gc[ 41]=   0.9517660952894889436594081E-02;
    gp[ 42]=   0.3689042825463932690061107E+00;
    gp[ngauss+ 42]=   0.2107880663979872074525160E+00;
    gp[2*ngauss+ 42]=   0.2386007375518623050589814E+00;
    gc[ 42]=   0.9517660952894889436594081E-02;
    gp[ 43]=   0.5123812452695834241199663E+00;
    gp[ngauss+ 43]=   0.2107880663979872074525160E+00;
    gp[2*ngauss+ 43]=   0.2386007375518623050589814E+00;
    gc[ 43]=   0.5076729393991831949508787E-02;
    gp[ 44]=   0.4984652136888425922032195E-01;
    gp[ngauss+ 44]=   0.4347909280428757352601284E-01;
    gp[2*ngauss+ 44]=   0.2386007375518623050589814E+00;
    gc[ 44]=   0.3381089578564921959659668E-02;
    gp[ 45]=   0.2369204605788584548781286E+00;
    gp[ngauss+ 45]=   0.4347909280428757352601284E-01;
    gp[2*ngauss+ 45]=   0.2386007375518623050589814E+00;
    gc[ 45]=   0.6338739326589163169268273E-02;
    gp[ 46]=   0.4809997090649916665368772E+00;
    gp[ngauss+ 46]=   0.4347909280428757352601284E-01;
    gp[2*ngauss+ 46]=   0.2386007375518623050589814E+00;
    gc[ 46]=   0.6338739326589163169268273E-02;
    gp[ 47]=   0.6680736482749658621946838E+00;
    gp[ngauss+ 47]=   0.4347909280428757352601284E-01;
    gp[2*ngauss+ 47]=   0.2386007375518623050589814E+00;
    gc[ 47]=   0.3381089578564921959659668E-02;
    gp[ 48]=   0.9233146216573625006194481E-02;
    gp[ngauss+ 48]=   0.8185180164205332861703353E+00;
    gp[2*ngauss+ 48]=   0.4850054944699732929706726E-01;
    gc[ 48]=   0.6013729287201758834679817E-03;
    gp[ 49]=   0.4388513368935080907940955E-01;
    gp[ngauss+ 49]=   0.8185180164205332861703353E+00;
    gp[2*ngauss+ 49]=   0.4850054944699732929706726E-01;
    gc[ 49]=   0.1127431304213664877060578E-02;
    gp[ 50]=   0.8909630044311857545318785E-01;
    gp[ngauss+ 50]=   0.8185180164205332861703353E+00;
    gp[2*ngauss+ 50]=   0.4850054944699732929706726E-01;
    gc[ 50]=   0.1127431304213664877060578E-02;
    gp[ 51]=   0.1237482879158957595264029E+00;
    gp[ngauss+ 51]=   0.8185180164205332861703353E+00;
    gp[2*ngauss+ 51]=   0.4850054944699732929706726E-01;
    gc[ 51]=   0.6013729287201758834679817E-03;
    gp[ 52]=   0.2750983225384828197777377E-01;
    gp[ngauss+ 52]=   0.5552859757470136190763871E+00;
    gp[2*ngauss+ 52]=   0.4850054944699732929706726E-01;
    gc[ 52]=   0.2504309443009021240723960E-02;
    gp[ 53]=   0.1307542020795333691342280E+00;
    gp[ngauss+ 53]=   0.5552859757470136190763871E+00;
    gp[2*ngauss+ 53]=   0.4850054944699732929706726E-01;
    gc[ 53]=   0.4694984969634420460775905E-02;
    gp[ 54]=   0.2654592727264556824923177E+00;
    gp[ngauss+ 54]=   0.5552859757470136190763871E+00;
    gp[2*ngauss+ 54]=   0.4850054944699732929706726E-01;
    gc[ 54]=   0.4694984969634420460775905E-02;
    gp[ 55]=   0.3687036425521407696487719E+00;
    gp[ngauss+ 55]=   0.5552859757470136190763871E+00;
    gp[2*ngauss+ 55]=   0.4850054944699732929706726E-01;
    gc[ 55]=   0.2504309443009021240723960E-02;
    gp[ 56]=   0.4777490464781690413678532E-01;
    gp[ngauss+ 56]=   0.2634159753661122469767895E+00;
    gp[2*ngauss+ 56]=   0.4850054944699732929706726E-01;
    gc[ 56]=   0.3924126780763078895076854E-02;
    gp[ 57]=   0.2270740686096784331853873E+00;
    gp[ngauss+ 57]=   0.2634159753661122469767895E+00;
    gp[2*ngauss+ 57]=   0.4850054944699732929706726E-01;
    gc[ 57]=   0.7356805009082974006098323E-02;
    gp[ 58]=   0.4610094065772119905407560E+00;
    gp[ngauss+ 58]=   0.2634159753661122469767895E+00;
    gp[2*ngauss+ 58]=   0.4850054944699732929706726E-01;
    gc[ 58]=   0.7356805009082974006098323E-02;
    gp[ 59]=   0.6403085705390735195893580E+00;
    gp[ngauss+ 59]=   0.2634159753661122469767895E+00;
    gp[2*ngauss+ 59]=   0.4850054944699732929706726E-01;
    gc[ 59]=   0.3924126780763078895076854E-02;
    gp[ 60]=   0.6229180934845267994089097E-01;
    gp[ngauss+ 60]=   0.5433461122723448458170192E-01;
    gp[2*ngauss+ 60]=   0.4850054944699732929706726E-01;
    gc[ 60]=   0.2613459007507404913181355E-02;
    gp[ 61]=   0.2960729004920768122935820E+00;
    gp[ngauss+ 61]=   0.5433461122723448458170192E-01;
    gp[2*ngauss+ 61]=   0.4850054944699732929706726E-01;
    gc[ 61]=   0.4899614459888755644422873E-02;
    gp[ 62]=   0.6010919388336913738276488E+00;
    gp[ngauss+ 62]=   0.5433461122723448458170192E-01;
    gp[2*ngauss+ 62]=   0.4850054944699732929706726E-01;
    gc[ 62]=   0.4899614459888755644422873E-02;
    gp[ 63]=   0.8348730299773155061803398E+00;
    gp[ngauss+ 63]=   0.5433461122723448458170192E-01;
    gp[2*ngauss+ 63]=   0.4850054944699732929706726E-01;
    gc[ 63]=   0.2613459007507404913181355E-02;
    break;
  case 5:
    gp[0]=   0.6884703934122676876048399E-03;
    gp[ngauss+  0]=   0.1342694011463441143815669E+00;
    gp[2*ngauss+  0]=   0.8510542129470164181162242E+00;
    gc[  0]=   0.7674555521798018092669516E-05;
    gp[1]=   0.3386801256323271594835794E-02;
    gp[ngauss+  1]=   0.1342694011463441143815669E+00;
    gp[2*ngauss+  1]=   0.8510542129470164181162242E+00;
    gc[  1]=   0.1550378001720072350948701E-04;
    gp[2]=   0.7338192953319733751104449E-02;
    gp[ngauss+  2]=   0.1342694011463441143815669E+00;
    gp[2*ngauss+  2]=   0.8510542129470164181162242E+00;
    gc[  2]=   0.1842749657758906164643500E-04;
    gp[3]=   0.1128958465031619590737310E-01;
    gp[ngauss+  3]=   0.1342694011463441143815669E+00;
    gp[2*ngauss+  3]=   0.8510542129470164181162242E+00;
    gc[  3]=   0.1550378001720072350948701E-04;
    gp[4]=   0.1398791551322719981460406E-01;
    gp[ngauss+  4]=   0.1342694011463441143815669E+00;
    gp[2*ngauss+  4]=   0.8510542129470164181162242E+00;
    gc[  4]=   0.7674555521798018092669516E-05;
    gp[5]=   0.2127808889925481860954521E-02;
    gp[ngauss+  5]=   0.1035864735618886456835587E+00;
    gp[2*ngauss+  5]=   0.8510542129470164181162242E+00;
    gc[  5]=   0.3601859320129832317508089E-04;
    gp[6]=   0.1046735762433882115148403E-01;
    gp[ngauss+  6]=   0.1035864735618886456835587E+00;
    gp[2*ngauss+  6]=   0.8510542129470164181162242E+00;
    gc[  6]=   0.7276308627071362096789688E-04;
    gp[7]=   0.2267965674554746810010855E-01;
    gp[ngauss+  7]=   0.1035864735618886456835587E+00;
    gp[2*ngauss+  7]=   0.8510542129470164181162242E+00;
    gc[  7]=   0.8648481349327657670994506E-04;
    gp[  8]=   0.3489195586675611504873307E-01;
    gp[ngauss+  8]=   0.1035864735618886456835587E+00;
    gp[2*ngauss+  8]=   0.8510542129470164181162242E+00;
    gc[  8]=   0.7276308627071362096789688E-04;
    gp[ 9]=   0.4323150460116945433926258E-01;
    gp[ngauss+ 9]=   0.1035864735618886456835587E+00;
    gp[2*ngauss+ 9]=   0.8510542129470164181162242E+00;
    gc[ 9]=   0.3601859320129832317508089E-04;
    gp[ 10]=   0.3926902791626685426815838E-02;
    gp[ngauss+ 10]=   0.6523450282167806813349548E-01;
    gp[2*ngauss+ 10]=   0.8510542129470164181162242E+00;
    gc[ 10]=   0.7133992621705575404868601E-04;
    gp[ 11]=   0.1931766338160684304489376E-01;
    gp[ngauss+ 11]=   0.6523450282167806813349548E-01;
    gp[2*ngauss+ 11]=   0.8510542129470164181162242E+00;
    gc[ 11]=   0.1441175999536500782735814E-03;
    gp[ 12]=   0.4185564211565275687514017E-01;
    gp[ngauss+ 12]=   0.6523450282167806813349548E-01;
    gp[2*ngauss+ 12]=   0.8510542129470164181162242E+00;
    gc[ 12]=   0.1712954245332319660698772E-03;
    gp[ 13]=   0.6439362084969867070538657E-01;
    gp[ngauss+ 13]=   0.6523450282167806813349548E-01;
    gp[2*ngauss+ 13]=   0.8510542129470164181162242E+00;
    gc[ 13]=   0.1441175999536500782735814E-03;
    gp[ 14]=   0.7978438143967882832346449E-01;
    gp[ngauss+ 14]=   0.6523450282167806813349548E-01;
    gp[2*ngauss+ 14]=   0.8510542129470164181162242E+00;
    gc[ 14]=   0.7133992621705575404868601E-04;
    gp[ 15]=   0.5603527040461490761502753E-02;
    gp[ngauss+ 15]=   0.2949326437223589592796798E-01;
    gp[2*ngauss+ 15]=   0.8510542129470164181162242E+00;
    gc[ 15]=   0.8147053631288434281442928E-04;
    gp[ 16]=   0.2756550260123100869305021E-01;
    gp[ngauss+ 16]=   0.2949326437223589592796798E-01;
    gp[2*ngauss+ 16]=   0.8510542129470164181162242E+00;
    gc[ 16]=   0.1645829871568117858991650E-03;
    gp[ 17]=   0.5972626134037384297790392E-01;
    gp[ngauss+ 17]=   0.2949326437223589592796798E-01;
    gp[2*ngauss+ 17]=   0.8510542129470164181162242E+00;
    gc[ 17]=   0.1956201925721807731387582E-03;
    gp[ 18]=   0.9188702007951667726275762E-01;
    gp[ngauss+ 18]=   0.2949326437223589592796798E-01;
    gp[2*ngauss+ 18]=   0.8510542129470164181162242E+00;
    gc[ 18]=   0.1645829871568117858991650E-03;
    gp[ 19]=   0.1138489956402861951943051E+00;
    gp[ngauss+ 19]=   0.2949326437223589592796798E-01;
    gp[2*ngauss+ 19]=   0.8510542129470164181162242E+00;
    gc[ 19]=   0.8147053631288434281442928E-04;
    gp[ 20]=   0.6708904550162072916648035E-02;
    gp[ngauss+ 20]=   0.5929510490997780154719579E-02;
    gp[2*ngauss+ 20]=   0.8510542129470164181162242E+00;
    gc[ 20]=   0.4716533650593665851425186E-04;
    gp[ 21]=   0.3300320039388486633331973E-01;
    gp[ngauss+ 21]=   0.5929510490997780154719579E-02;
    gp[2*ngauss+ 21]=   0.8510542129470164181162242E+00;
    gc[ 21]=   0.9528121850813989619054645E-04;
    gp[ 22]=   0.7150813828099290086452812E-01;
    gp[ngauss+ 22]=   0.5929510490997780154719579E-02;
    gp[2*ngauss+ 22]=   0.8510542129470164181162242E+00;
    gc[ 22]=   0.1132494350422471987715167E-03;
    gp[ 23]=   0.1100130761681009353957365E+00;
    gp[ngauss+ 23]=   0.5929510490997780154719579E-02;
    gp[2*ngauss+ 23]=   0.8510542129470164181162242E+00;
    gc[ 23]=   0.9528121850813989619054645E-04;
    gp[ 24]=   0.1363073720118237288124082E+00;
    gp[ngauss+ 24]=   0.5929510490997780154719579E-02;
    gp[2*ngauss+ 24]=   0.8510542129470164181162242E+00;
    gc[ 24]=   0.4716533650593665851425186E-04;
    gp[ 25]=   0.1690216171511836227042953E-02;
    gp[ngauss+ 25]=   0.3296355447210387441801802E+00;
    gp[2*ngauss+ 25]=   0.6343334726308867723471639E+00;
    gc[ 25]=   0.5980139538929241281078968E-04;
    gp[ 26]=   0.8314702139567988954602070E-02;
    gp[ngauss+ 26]=   0.3296355447210387441801802E+00;
    gp[2*ngauss+ 26]=   0.6343334726308867723471639E+00;
    gc[ 26]=   0.1208079967893718785024756E-03;
    gp[ 27]=   0.1801549132403724173632794E-01;
    gp[ngauss+ 27]=   0.3296355447210387441801802E+00;
    gp[2*ngauss+ 27]=   0.6343334726308867723471639E+00;
    gc[ 27]=   0.1435900757693728527211053E-03;
    gp[ 28]=   0.2771628050850649451805380E-01;
    gp[ngauss+ 28]=   0.3296355447210387441801802E+00;
    gp[2*ngauss+ 28]=   0.6343334726308867723471639E+00;
    gc[ 28]=   0.1208079967893718785024756E-03;
    gp[ 29]=   0.3434076647656264724561292E-01;
    gp[ngauss+ 29]=   0.3296355447210387441801802E+00;
    gp[2*ngauss+ 29]=   0.6343334726308867723471639E+00;
    gc[ 29]=   0.5980139538929241281078968E-04;
    gp[ 30]=   0.5223836827337728335993803E-02;
    gp[ngauss+ 30]=   0.2543080057465078161577979E+00;
    gp[2*ngauss+ 30]=   0.6343334726308867723471639E+00;
    gc[ 30]=   0.2806627859136634173965294E-03;
    gp[ 31]=   0.2569768765504614119336387E-01;
    gp[ngauss+ 31]=   0.2543080057465078161577979E+00;
    gp[2*ngauss+ 31]=   0.6343334726308867723471639E+00;
    gc[ 31]=   0.5669819026601681133017675E-03;
    gp[ 32]=   0.5567926081130270574751912E-01;
    gp[ngauss+ 32]=   0.2543080057465078161577979E+00;
    gp[2*ngauss+ 32]=   0.6343334726308867723471639E+00;
    gc[ 32]=   0.6739038517854064343439875E-03;
    gp[ 33]=   0.8566083396755927030167437E-01;
    gp[ngauss+ 33]=   0.2543080057465078161577979E+00;
    gp[2*ngauss+ 33]=   0.6343334726308867723471639E+00;
    gc[ 33]=   0.5669819026601681133017675E-03;
    gp[ 34]=   0.1061346847952676831590444E+00;
    gp[ngauss+ 34]=   0.2543080057465078161577979E+00;
    gp[2*ngauss+ 34]=   0.6343334726308867723471639E+00;
    gc[ 34]=   0.2806627859136634173965294E-03;
    gp[ 35]=   0.9640668162164327438527513E-02;
    gp[ngauss+ 35]=   0.1601527279383079979472859E+00;
    gp[2*ngauss+ 35]=   0.6343334726308867723471639E+00;
    gc[ 35]=   0.5558924060985351094218042E-03;
    gp[ 36]=   0.4742546281705090567357321E-01;
    gp[ngauss+ 36]=   0.1601527279383079979472859E+00;
    gp[2*ngauss+ 36]=   0.6343334726308867723471639E+00;
    gc[ 36]=   0.1122987976685449209851463E-02;
    gp[ 37]=   0.1027568997154026148527751E+00;
    gp[ngauss+ 37]=   0.1601527279383079979472859E+00;
    gp[2*ngauss+ 37]=   0.6343334726308867723471639E+00;
    gc[ 37]=   0.1334762043455590017807318E-02;
    gp[ 38]=   0.1580883366137543240319770E+00;
    gp[ngauss+ 38]=   0.1601527279383079979472859E+00;
    gp[2*ngauss+ 38]=   0.6343334726308867723471639E+00;
    gc[ 38]=   0.1122987976685449209851463E-02;
    gp[ 39]=   0.1958731312686409022670227E+00;
    gp[ngauss+ 39]=   0.1601527279383079979472859E+00;
    gp[2*ngauss+ 39]=   0.6343334726308867723471639E+00;
    gc[ 39]=   0.5558924060985351094218042E-03;
    gp[ 40]=   0.1375683270031391679650191E-01;
    gp[ngauss+ 40]=   0.7240687888633139719987719E-01;
    gp[2*ngauss+ 40]=   0.6343334726308867723471639E+00;
    gc[ 40]=   0.6348317815652551040684828E-03;
    gp[ 41]=   0.6767416394121158260026687E-01;
    gp[ngauss+ 41]=   0.7240687888633139719987719E-01;
    gp[2*ngauss+ 41]=   0.6343334726308867723471639E+00;
    gc[ 41]=   0.1282457630459549363510960E-02;
    gp[ 42]=   0.1466298242413909152264795E+00;
    gp[ngauss+ 42]=   0.7240687888633139719987719E-01;
    gp[2*ngauss+ 42]=   0.6343334726308867723471639E+00;
    gc[ 42]=   0.1524304625709161315839557E-02;
    gp[ 43]=   0.2255854845415702478526920E+00;
    gp[ngauss+ 43]=   0.7240687888633139719987719E-01;
    gp[2*ngauss+ 43]=   0.6343334726308867723471639E+00;
    gc[ 43]=   0.1282457630459549363510960E-02;
    gp[ 44]=   0.2795028157824679136564570E+00;
    gp[ngauss+ 44]=   0.7240687888633139719987719E-01;
    gp[2*ngauss+ 44]=   0.6343334726308867723471639E+00;
    gc[ 44]=   0.6348317815652551040684828E-03;
    gp[ 45]=   0.1647056877436847659039497E-01;
    gp[ngauss+ 45]=   0.1455713218307138008948628E-01;
    gp[2*ngauss+ 45]=   0.6343334726308867723471639E+00;
    gc[ 45]=   0.3675200380073265705508936E-03;
    gp[ 46]=   0.8102388069429512304981796E-01;
    gp[ngauss+ 46]=   0.1455713218307138008948628E-01;
    gp[2*ngauss+ 46]=   0.6343334726308867723471639E+00;
    gc[ 46]=   0.7424468824279099317121437E-03;
    gp[ 47]=   0.1755546975930209237816749E+00;
    gp[ngauss+ 47]=   0.1455713218307138008948628E-01;
    gp[2*ngauss+ 47]=   0.6343334726308867723471639E+00;
    gc[ 47]=   0.8824581727683867844078806E-03;
    gp[ 48]=   0.2700855144917467245135319E+00;
    gp[ngauss+ 48]=   0.1455713218307138008948628E-01;
    gp[2*ngauss+ 48]=   0.6343334726308867723471639E+00;
    gc[ 48]=   0.7424468824279099317121437E-03;
    gp[ 49]=   0.3346388264116733709729549E+00;
    gp[ngauss+ 49]=   0.1455713218307138008948628E-01;
    gp[2*ngauss+ 49]=   0.6343334726308867723471639E+00;
    gc[ 49]=   0.3675200380073265705508936E-03;
    gp[ 50]=   0.2820121115434851485096519E-02;
    gp[ngauss+ 50]=   0.5499960157369496423867172E+00;
    gp[2*ngauss+ 50]=   0.3898863870655193282408954E+00;
    gc[ 50]=   0.1664075540527897608455068E-03;
    gp[ 51]=   0.1387305805468257439256853E-01;
    gp[ngauss+ 51]=   0.5499960157369496423867172E+00;
    gp[2*ngauss+ 51]=   0.3898863870655193282408954E+00;
    gc[ 51]=   0.3361687988193032947249193E-03;
    gp[ 52]=   0.3005879859876551468619370E-01;
    gp[ngauss+ 52]=   0.5499960157369496423867172E+00;
    gp[2*ngauss+ 52]=   0.3898863870655193282408954E+00;
    gc[ 52]=   0.3995638084945832988385996E-03;
    gp[ 53]=   0.4624453914284845497981886E-01;
    gp[ngauss+ 53]=   0.5499960157369496423867172E+00;
    gp[2*ngauss+ 53]=   0.3898863870655193282408954E+00;
    gc[ 53]=   0.3361687988193032947249193E-03;
    gp[ 54]=   0.5729747608209617788729088E-01;
    gp[ngauss+ 54]=   0.5499960157369496423867172E+00;
    gp[2*ngauss+ 54]=   0.3898863870655193282408954E+00;
    gc[ 54]=   0.1664075540527897608455068E-03;
    gp[ 55]=   0.8715957632321213455687685E-02;
    gp[ngauss+ 55]=   0.4243122204826401923058484E+00;
    gp[2*ngauss+ 55]=   0.3898863870655193282408954E+00;
    gc[ 55]=   0.7809919386245131847165196E-03;
    gp[ 56]=   0.4287652242081133138941515E-01;
    gp[ngauss+ 56]=   0.4243122204826401923058484E+00;
    gp[2*ngauss+ 56]=   0.3898863870655193282408954E+00;
    gc[ 56]=   0.1577723579854277460560149E-02;
    gp[ 57]=   0.9290069622592023972662808E-01;
    gp[ngauss+ 57]=   0.4243122204826401923058484E+00;
    gp[2*ngauss+ 57]=   0.3898863870655193282408954E+00;
    gc[ 57]=   0.1875252089225373929164060E-02;
    gp[ 58]=   0.1429248700310291480638410E+00;
    gp[ngauss+ 58]=   0.4243122204826401923058484E+00;
    gp[2*ngauss+ 58]=   0.3898863870655193282408954E+00;
    gc[ 58]=   0.1577723579854277460560149E-02;
    gp[ 59]=   0.1770854348195192659975685E+00;
    gp[ngauss+ 59]=   0.4243122204826401923058484E+00;
    gp[2*ngauss+ 59]=   0.3898863870655193282408954E+00;
    gc[ 59]=   0.7809919386245131847165196E-03;
    gp[ 60]=   0.1608542878080594201063844E-01;
    gp[ngauss+ 60]=   0.2672143938543263687711941E+00;
    gp[2*ngauss+ 60]=   0.3898863870655193282408954E+00;
    gc[ 60]=   0.1546865169503060314679319E-02;
    gp[ 61]=   0.7912925657314306568907516E-01;
    gp[ngauss+ 61]=   0.2672143938543263687711941E+00;
    gp[2*ngauss+ 61]=   0.3898863870655193282408954E+00;
    gc[ 61]=   0.3124905049696835174549313E-02;
    gp[ 62]=   0.1714496095400771514939552E+00;
    gp[ngauss+ 62]=   0.2672143938543263687711941E+00;
    gp[2*ngauss+ 62]=   0.3898863870655193282408954E+00;
    gc[ 62]=   0.3714202410295569084604130E-02;
    gp[ 63]=   0.2637699625070112372988353E+00;
    gp[ngauss+ 63]=   0.2672143938543263687711941E+00;
    gp[2*ngauss+ 63]=   0.3898863870655193282408954E+00;
    gc[ 63]=   0.3124905049696835174549313E-02;
    gp[ 64]=   0.3268137902993483609772720E+00;
    gp[ngauss+ 64]=   0.2672143938543263687711941E+00;
    gp[2*ngauss+ 64]=   0.3898863870655193282408954E+00;
    gc[ 64]=   0.1546865169503060314679319E-02;
    gp[ 65]=   0.2295323819139559290243112E-01;
    gp[ngauss+ 65]=   0.1208106817883721533703770E+00;
    gp[2*ngauss+ 65]=   0.3898863870655193282408954E+00;
    gc[ 65]=   0.1766527408224395001754439E-02;
    gp[ 66]=   0.1129141596895874545450887E+00;
    gp[ngauss+ 66]=   0.1208106817883721533703770E+00;
    gp[2*ngauss+ 66]=   0.3898863870655193282408954E+00;
    gc[ 66]=   0.3568656484883993820498823E-02;
    gp[ 67]=   0.2446514655730542591943638E+00;
    gp[ngauss+ 67]=   0.1208106817883721533703770E+00;
    gp[2*ngauss+ 67]=   0.3898863870655193282408954E+00;
    gc[ 67]=   0.4241636883961948877319364E-02;
    gp[ 68]=   0.3763887714565210638436389E+00;
    gp[ngauss+ 68]=   0.1208106817883721533703770E+00;
    gp[2*ngauss+ 68]=   0.3898863870655193282408954E+00;
    gc[ 68]=   0.3568656484883993820498823E-02;
    gp[ 69]=   0.4663496929547129254862965E+00;
    gp[ngauss+ 69]=   0.1208106817883721533703770E+00;
    gp[2*ngauss+ 69]=   0.3898863870655193282408954E+00;
    gc[ 69]=   0.1766527408224395001754439E-02;
    gp[ 70]=   0.2748109949881235672570692E-01;
    gp[ngauss+ 70]=   0.2428853571607680625473734E-01;
    gp[2*ngauss+ 70]=   0.3898863870655193282408954E+00;
    gc[ 70]=   0.1022687015780539027115484E-02;
    gp[ 71]=   0.1351881260230007058889347E+00;
    gp[ngauss+ 71]=   0.2428853571607680625473734E-01;
    gp[2*ngauss+ 71]=   0.3898863870655193282408954E+00;
    gc[ 71]=   0.2065984730200281882763274E-02;
    gp[ 72]=   0.2929125386092019327521836E+00;
    gp[ngauss+ 72]=   0.2428853571607680625473734E-01;
    gp[2*ngauss+ 72]=   0.3898863870655193282408954E+00;
    gc[ 72]=   0.2455589959537547058785597E-02;
    gp[ 73]=   0.4506369511954031596154326E+00;
    gp[ngauss+ 73]=   0.2428853571607680625473734E-01;
    gp[2*ngauss+ 73]=   0.3898863870655193282408954E+00;
    gc[ 73]=   0.2065984730200281882763274E-02;
    gp[ 74]=   0.5583439777195915087786603E+00;
    gp[ngauss+ 74]=   0.2428853571607680625473734E-01;
    gp[2*ngauss+ 74]=   0.3898863870655193282408954E+00;
    gc[ 74]=   0.1022687015780539027115484E-02;
    gp[ 75]=   0.3820412379430865050367239E-02;
    gp[ngauss+ 75]=   0.7450784917211248190869680E+00;
    gp[2*ngauss+ 75]=   0.1734803207716957231045924E+00;
    gc[ 75]=   0.2354307468301136510627038E-03;
    gp[ 76]=   0.1879380372800047934136409E-01;
    gp[ngauss+ 76]=   0.7450784917211248190869680E+00;
    gp[2*ngauss+ 76]=   0.1734803207716957231045924E+00;
    gc[ 76]=   0.4756062416607821967306727E-03;
    gp[ 77]=   0.4072059375358972890421977E-01;
    gp[ngauss+ 77]=   0.7450784917211248190869680E+00;
    gp[2*ngauss+ 77]=   0.1734803207716957231045924E+00;
    gc[ 77]=   0.5652964877443147112891440E-03;
    gp[ 78]=   0.6264738377917897846707544E-01;
    gp[ngauss+ 78]=   0.7450784917211248190869680E+00;
    gp[2*ngauss+ 78]=   0.1734803207716957231045924E+00;
    gc[ 78]=   0.4756062416607821967306727E-03;
    gp[ 79]=   0.7762077512774859275807230E-01;
    gp[ngauss+ 79]=   0.7450784917211248190869680E+00;
    gp[2*ngauss+ 79]=   0.1734803207716957231045924E+00;
    gc[ 79]=   0.2354307468301136510627038E-03;
    gp[ 80]=   0.1180749020134916839035834E-01;
    gp[ngauss+ 80]=   0.5748149081269930263556251E+00;
    gp[2*ngauss+ 80]=   0.1734803207716957231045924E+00;
    gc[ 80]=   0.1104934907704599594434288E-02;
    gp[ 81]=   0.5808473832803965156390707E-01;
    gp[ngauss+ 81]=   0.5748149081269930263556251E+00;
    gp[2*ngauss+ 81]=   0.1734803207716957231045924E+00;
    gc[ 81]=   0.2232138094997411949593085E-02;
    gp[ 82]=   0.1258523855506556252698913E+00;
    gp[ngauss+ 82]=   0.5748149081269930263556251E+00;
    gp[2*ngauss+ 82]=   0.1734803207716957231045924E+00;
    gc[ 82]=   0.2653076672955636507052818E-02;
    gp[ 83]=   0.1936200327732715989758754E+00;
    gp[ngauss+ 83]=   0.5748149081269930263556251E+00;
    gp[2*ngauss+ 83]=   0.1734803207716957231045924E+00;
    gc[ 83]=   0.2232138094997411949593085E-02;
    gp[ 84]=   0.2398972808999620821494242E+00;
    gp[ngauss+ 84]=   0.5748149081269930263556251E+00;
    gp[2*ngauss+ 84]=   0.1734803207716957231045924E+00;
    gc[ 84]=   0.1104934907704599594434288E-02;
    gp[ 85]=   0.2179089788247223673946365E-01;
    gp[ngauss+ 85]=   0.3619947996757470286840036E+00;
    gp[2*ngauss+ 85]=   0.1734803207716957231045924E+00;
    gc[ 85]=   0.2188480109418990024542297E-02;
    gp[ 86]=   0.1071962440664831064485623E+00;
    gp[ngauss+ 86]=   0.3619947996757470286840036E+00;
    gp[2*ngauss+ 86]=   0.1734803207716957231045924E+00;
    gc[ 86]=   0.4421065701079485357792387E-02;
    gp[ 87]=   0.2322624397762786241057020E+00;
    gp[ngauss+ 87]=   0.3619947996757470286840036E+00;
    gp[2*ngauss+ 87]=   0.1734803207716957231045924E+00;
    gc[ 87]=   0.5254794184744129496514531E-02;
    gp[ 88]=   0.3573286354860741417628417E+00;
    gp[ngauss+ 88]=   0.3619947996757470286840036E+00;
    gp[2*ngauss+ 88]=   0.1734803207716957231045924E+00;
    gc[ 88]=   0.4421065701079485357792387E-02;
    gp[ 89]=   0.4427339816700850114719404E+00;
    gp[ngauss+ 89]=   0.3619947996757470286840036E+00;
    gp[2*ngauss+ 89]=   0.1734803207716957231045924E+00;
    gc[ 89]=   0.2188480109418990024542297E-02;
    gp[ 90]=   0.3109470542044839223481313E-01;
    gp[ngauss+ 90]=   0.1636619866237947995192818E+00;
    gp[2*ngauss+ 90]=   0.1734803207716957231045924E+00;
    gc[ 90]=   0.2499254732643923770293104E-02;
    gp[ 91]=   0.1529645840847571531666495E+00;
    gp[ngauss+ 91]=   0.1636619866237947995192818E+00;
    gp[2*ngauss+ 91]=   0.1734803207716957231045924E+00;
    gc[ 91]=   0.5048878136564868823849832E-02;
    gp[ 92]=   0.3314288463022547386880629E+00;
    gp[ngauss+ 92]=   0.1636619866237947995192818E+00;
    gp[2*ngauss+ 92]=   0.1734803207716957231045924E+00;
    gc[ 92]=   0.6001000045085251254970781E-02;
    gp[ 93]=   0.5098931085197523242094763E+00;
    gp[ngauss+ 93]=   0.1636619866237947995192818E+00;
    gp[2*ngauss+ 93]=   0.1734803207716957231045924E+00;
    gc[ 93]=   0.5048878136564868823849832E-02;
    gp[ 94]=   0.6317629871840610851413126E+00;
    gp[ngauss+ 94]=   0.1636619866237947995192818E+00;
    gp[2*ngauss+ 94]=   0.1734803207716957231045924E+00;
    gc[ 94]=   0.2499254732643923770293104E-02;
    gp[ 95]=   0.3722858998892505406031969E-01;
    gp[ngauss+ 95]=   0.3290363028030459202550237E-01;
    gp[2*ngauss+ 95]=   0.1734803207716957231045924E+00;
    gc[ 95]=   0.1446881238470051746865003E-02;
    gp[ 96]=   0.1831390812910861357644692E+00;
    gp[ngauss+ 96]=   0.3290363028030459202550237E-01;
    gp[2*ngauss+ 96]=   0.1734803207716957231045924E+00;
    gc[ 96]=   0.2922922163836161298883409E-02;
    gp[ 97]=   0.3968080244739998424349526E+00;
    gp[ngauss+ 97]=   0.3290363028030459202550237E-01;
    gp[2*ngauss+ 97]=   0.1734803207716957231045924E+00;
    gc[ 97]=   0.3474129413013635216495632E-02;
    gp[ 98]=   0.6104769676569135491054361E+00;
    gp[ngauss+ 98]=   0.3290363028030459202550237E-01;
    gp[2*ngauss+ 98]=   0.1734803207716957231045924E+00;
    gc[ 98]=   0.2922922163836161298883409E-02;
    gp[99]=   0.7563874589590746308095855E+00;
    gp[ngauss+99]=   0.3290363028030459202550237E-01;
    gp[2*ngauss+99]=   0.1734803207716957231045924E+00;
    gc[99]=   0.1446881238470051746865003E-02;
    gp[100]=   0.4462454629928929415083641E-02;
    gp[ngauss+100]=   0.8702932130946322704376958E+00;
    gp[2*ngauss+100]=   0.3457893991821509152445743E-01;
    gc[100]=   0.1525364704986189692495816E-03;
    gp[101]=   0.2195221042407078662784791E-01;
    gp[ngauss+101]=   0.8702932130946322704376958E+00;
    gp[2*ngauss+101]=   0.3457893991821509152445743E-01;
    gc[101]=   0.3081470811558820321245740E-03;
    gp[102]=   0.4756392349357631901892338E-01;
    gp[ngauss+102]=   0.8702932130946322704376958E+00;
    gp[2*ngauss+102]=   0.3457893991821509152445743E-01;
    gc[102]=   0.3662577305079262619620616E-03;
    gp[103]=   0.7317563656308185140999884E-01;
    gp[ngauss+103]=   0.8702932130946322704376958E+00;
    gp[2*ngauss+103]=   0.3457893991821509152445743E-01;
    gc[103]=   0.3081470811558820321245740E-03;
    gp[104]=   0.9066539235722370862276311E-01;
    gp[ngauss+104]=   0.8702932130946322704376958E+00;
    gp[2*ngauss+104]=   0.3457893991821509152445743E-01;
    gc[104]=   0.1525364704986189692495816E-03;
    gp[105]=   0.1379180676948294852564964E-01;
    gp[ngauss+105]=   0.6714158560300755951647965E+00;
    gp[2*ngauss+105]=   0.3457893991821509152445743E-01;
    gc[105]=   0.7158915019438693925839182E-03;
    gp[106]=   0.6784621232925240815339392E-01;
    gp[ngauss+106]=   0.6714158560300755951647965E+00;
    gp[2*ngauss+106]=   0.3457893991821509152445743E-01;
    gc[106]=   0.1446210706378584148049987E-02;
    gp[107]=   0.1470026020258546566553730E+00;
    gp[ngauss+107]=   0.6714158560300755951647965E+00;
    gp[2*ngauss+107]=   0.3457893991821509152445743E-01;
    gc[107]=   0.1718938401647664926715154E-02;
    gp[108]=   0.2261589917224569051573522E+00;
    gp[ngauss+108]=   0.6714158560300755951647965E+00;
    gp[2*ngauss+108]=   0.3457893991821509152445743E-01;
    gc[108]=   0.1446210706378584148049987E-02;
    gp[109]=   0.2802133972822263647850965E+00;
    gp[ngauss+109]=   0.6714158560300755951647965E+00;
    gp[2*ngauss+109]=   0.3457893991821509152445743E-01;
    gc[109]=   0.7158915019438693925839182E-03;
    gp[110]=   0.2545298347097098436579517E-01;
    gp[ngauss+110]=   0.4228301055981501231453074E+00;
    gp[2*ngauss+110]=   0.3457893991821509152445743E-01;
    gc[110]=   0.1417924532550925506744889E-02;
    gp[111]=   0.1252111887766239334716453E+00;
    gp[ngauss+111]=   0.4228301055981501231453074E+00;
    gp[2*ngauss+111]=   0.3457893991821509152445743E-01;
    gc[111]=   0.2864425173708486951862259E-02;
    gp[112]=   0.2712954772418173926651176E+00;
    gp[ngauss+112]=   0.4228301055981501231453074E+00;
    gp[2*ngauss+112]=   0.3457893991821509152445743E-01;
    gc[112]=   0.3404601008703135370818105E-02;
    gp[113]=   0.4173797657070108518585898E+00;
    gp[ngauss+113]=   0.4228301055981501231453074E+00;
    gp[2*ngauss+113]=   0.3457893991821509152445743E-01;
    gc[113]=   0.2864425173708486951862259E-02;
    gp[114]=   0.5171379710126638009644400E+00;
    gp[ngauss+114]=   0.4228301055981501231453074E+00;
    gp[2*ngauss+114]=   0.3457893991821509152445743E-01;
    gc[114]=   0.1417924532550925506744889E-02;
    gp[115]=   0.3632034932062158323429498E-01;
    gp[ngauss+115]=   0.1911663237939562572118897E+00;
    gp[2*ngauss+115]=   0.3457893991821509152445743E-01;
    gc[115]=   0.1619276585269326374535570E-02;
    gp[116]=   0.1786711612964319811195394E+00;
    gp[ngauss+116]=   0.1911663237939562572118897E+00;
    gp[2*ngauss+116]=   0.3457893991821509152445743E-01;
    gc[116]=   0.3271187222988250904745399E-02;
    gp[117]=   0.3871273681439143256318264E+00;
    gp[ngauss+117]=   0.1911663237939562572118897E+00;
    gp[2*ngauss+117]=   0.3457893991821509152445743E-01;
    gc[117]=   0.3888070605322794358683006E-02;
    gp[118]=   0.5955835749913966701441134E+00;
    gp[ngauss+118]=   0.1911663237939562572118897E+00;
    gp[2*ngauss+118]=   0.3457893991821509152445743E-01;
    gc[118]=   0.3271187222988250904745399E-02;
    gp[119]=   0.7379343869672070680293578E+00;
    gp[ngauss+119]=   0.1911663237939562572118897E+00;
    gp[2*ngauss+119]=   0.3457893991821509152445743E-01;
    gc[119]=   0.1619276585269326374535570E-02;
    gp[120]=   0.4348506843299289873538162E-01;
    gp[ngauss+120]=   0.3843327439633327330290728E-01;
    gp[2*ngauss+120]=   0.3457893991821509152445743E-01;
    gc[120]=   0.9374398217669952625390973E-03;
    gp[121]=   0.2139166561255058407993538E+00;
    gp[ngauss+121]=   0.3843327439633327330290728E-01;
    gp[2*ngauss+121]=   0.3457893991821509152445743E-01;
    gc[121]=   0.1893772314860302438326724E-02;
    gp[122]=   0.4634938928427258175863176E+00;
    gp[ngauss+122]=   0.3843327439633327330290728E-01;
    gp[2*ngauss+122]=   0.3457893991821509152445743E-01;
    gc[122]=   0.2250901574461454072738044E-02;
    gp[123]=   0.7130711295599457943732815E+00;
    gp[ngauss+123]=   0.3843327439633327330290728E-01;
    gp[2*ngauss+123]=   0.3457893991821509152445743E-01;
    gc[123]=   0.1893772314860302438326724E-02;
    gp[124]=   0.8835027172524587364372537E+00;
    gp[ngauss+124]=   0.3843327439633327330290728E-01;
    gp[2*ngauss+124]=   0.3457893991821509152445743E-01;
    gc[124]=   0.9374398217669952625390973E-03;
    break;
  }
  return;
}
/******************************************************************************/


/*****************************************************************************/
/******* SOON TO BE DEPRECATED ROUTINES **************************************/
/*****************************************************************************/

/*!
 * \fn struct qcoordinates *allocateqcoords(INT nq1d,INT nelm,INT mydim)
 *
 * \brief Allocates memory and properties of quadrature coordinates struct.
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param nelm    Number of elements to get quadrature on
 * \param mydim   Dimension of problem
 *
 * \return A      Quadrature struct
 *
 * \note TODO: Will be deprecated soon
 *
 */
struct qcoordinates *allocateqcoords(INT nq1d,INT nelm,INT mydim)
{
  // Flag for errors
  SHORT status;

  struct qcoordinates *A = malloc(sizeof(struct qcoordinates));
  assert(A != NULL);

  INT nq = 0;
  switch (mydim)
  {
    case 1:
    nq = nq1d;
    A->x = (REAL *) calloc(nq*nelm,sizeof(REAL));
    A->y = NULL;
    A->z = NULL;
    break;
    case 2:
    nq = nq1d*nq1d;
    A->x = (REAL *) calloc(nq*nelm,sizeof(REAL));
    A->y = (REAL *) calloc(nq*nelm,sizeof(REAL));
    A->z = NULL;
    break;
    case 3:
    nq = nq1d*nq1d*nq1d;
    A->x = (REAL *) calloc(nq*nelm,sizeof(REAL));
    A->y = (REAL *) calloc(nq*nelm,sizeof(REAL));
    A->z = (REAL *) calloc(nq*nelm,sizeof(REAL));
    break;
    case 4:
    nq = nq1d;
    A->z = (REAL *) calloc(nq*nelm,sizeof(REAL));
    default:
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  A->w = (REAL *) calloc(nq*nelm,sizeof(REAL));
  A->n = nq*nelm;
  A->nq_per_elm = nq;
  A->nq1d = nq1d;

  return A;
}
/******************************************************************************/

/*!
* \fn struct qcoordinates *allocateqcoords_bdry(INT nq1d,INT nregion,INT dim,INT ed_or_f)
*
* \brief Allocates memory and properties of quadrature coordinates struct.
*        Assumes we are allocated on a boundary, so dimension is 1 or 2 less
*
* \param nq1d    Number of quadrature nodes on an element in 1D direction
* \param nregion    Number of "elements" (faces or edges) to get quadrature on
* \param dim     Dimension of problem
* \param ed_or_f Whether we are computing quadrature on faces or edges
*
* \return A      Quadrature struct
*
* \note TODO: will be deprecated soon
*
*/
struct qcoordinates *allocateqcoords_bdry(INT nq1d,INT nregion,INT dim,INT ed_or_f)
{
  // Flag for errors
  SHORT status;

  struct qcoordinates *A = malloc(sizeof(struct qcoordinates));
  assert(A != NULL);

  INT nq = 0;
  switch (ed_or_f)
  {
    case 1:
    // Compute quadrature on edges
    nq = nq1d;
    break;
    case 2:
    if(dim==2) { // face is edge in 2D
      nq = nq1d;
    } else {
      nq = nq1d*nq1d;
    }
    break;
    default:
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  A->x = (REAL *) calloc(nq*nregion,sizeof(REAL));
  A->y = (REAL *) calloc(nq*nregion,sizeof(REAL));
  if(dim==2) {
    A->z = NULL;
  } else if(dim==3) {
    A->z = (REAL *) calloc(nq*nregion,sizeof(REAL));
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  A->w = (REAL *) calloc(nq*nregion,sizeof(REAL));
  A->n = nq*nregion;
  A->nq_per_elm = nq;
  A->nq1d = nq1d;

  return A;
}
/******************************************************************************/

/*!
 * \fn void free_qcoords(qcoordinates* A)
 *
 * \brief Frees memory of arrays of quadrature struct
 *
 * \return A       Struct for quadratures to be freed
 *
 * \note TODO: will be deprecated soon
 */
void free_qcoords(qcoordinates* A)
{
  if (A==NULL) return;

  if(A->x) {
    free(A->x);
    A->x = NULL;
  }

  if(A->y) {
    free(A->y);
    A->y = NULL;
  }

  if(A->z) {
    free(A->z);
    A->z = NULL;
  }

  if(A->w) {
    free(A->w);
    A->w = NULL;
  }

  return;
}
/******************************************************************************/

/*!
 * \fn qcoordinates* get_quadrature(mesh_struct *mesh,INT nq1d)
 *
 * \brief Computes quadrature weights and nodes for entire domain using nq1d^(dim)
 *        quadrature nodes per element
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param mesh    Mesh struct
 *
 * \return cq_all      Quadrature struct
 *
 * \note TODO: This will be deprecated soon
 *
 */
qcoordinates* get_quadrature(mesh_struct *mesh,INT nq1d)
{
  INT i,j;

  INT dim = mesh->dim;
  INT nelm = mesh->nelm;
  INT nq = (INT) pow(nq1d,dim);
  qcoordinates *cq_all = allocateqcoords(nq1d,nelm,dim);
  qcoordinates *cqelm = allocateqcoords(nq1d,1,dim);

  for (i=0; i<nelm; i++) {
    quad_elm(cqelm,mesh,nq1d,i);
    for (j=0; j<nq; j++) {
      cq_all->x[i*nq+j] = cqelm->x[j];
      cq_all->w[i*nq+j] = cqelm->w[j];
      if(dim>1) cq_all->y[i*nq+j] = cqelm->y[j];
      if(dim>2) cq_all->z[i*nq+j] = cqelm->z[j];
    }
  }

  if(cqelm) {
    free_qcoords(cqelm);
    free(cqelm);
    cqelm = NULL;
  }

  return cq_all;
}
/******************************************************************************/

/*!
 * \fn void quad_elm(qcoordinates *cqelm,mesh_struct *mesh,INT nq1d,INT elm)
 *
 * \brief Computes quadrature weights and nodes for SINGLE element using nq1d^(dim)
 *        quadrature nodes on simplex
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param mesh    Mesh struct
 * \param elm     Index of current element
 *
 * \return cq_elm Quadrature struct on element
 *
 * \note TODO Will be deprecated soon
 *
 */
void quad_elm(qcoordinates *cqelm,mesh_struct *mesh,INT nq1d,INT elm)
{
  // Flag for errors
  SHORT status;

  /* Loop indices */
  INT q,j;

  /* Dimension */
  INT dim = mesh->dim;

  /* Total Number of Quadrature Nodes */
  INT nq = pow(nq1d,dim);

  /* Coordinates of vertices of element */
  INT v_per_elm = mesh->v_per_elm;
  coordinates *cvelm = allocatecoords(v_per_elm,dim);

  /* Gaussian points for reference element */
  REAL* gp = (REAL *) calloc(dim*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  /* Points on Reference Triangle */
  REAL r,s,t;

  REAL e_vol = mesh->el_vol[elm]; /* Area/Volume of Element */
  REAL voldim=0.0;

  // Get vertices of element
  INT* thiselm_v = (INT *) calloc(v_per_elm,sizeof(INT));
  iCSRmat* el_v = mesh->el_v;
  get_incidence_row(elm,el_v,thiselm_v);

  // Get coordinates of vertices for given element
  if(dim==1) {
    for (j=0; j<v_per_elm; j++) {
      cvelm->x[j] = mesh->cv->x[thiselm_v[j]];
    }
    voldim = 0.5*e_vol;

    // Get Quad Nodes and Weights
    quad1d(gp,gw,nq1d);

    // Map to Real Interval
    // x = 0.5*(x2-x1)*r + 0.5*(x1+x2)
    // w = 0.5*(x2-x1)*wref
    for (q=0; q<nq; q++) {
      r = gp[q];
      cqelm->x[q] = 0.5*(cvelm->x[0]*(1-r) + cvelm->x[1]*(1+r));
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==2) {
    for (j=0; j<v_per_elm; j++) {
      cvelm->x[j] = mesh->cv->x[thiselm_v[j]];
      cvelm->y[j] = mesh->cv->y[thiselm_v[j]];
    }
    voldim = 2.0*e_vol;

    // Get Quad Nodes and Weights
    triquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s) + x2*r + x3*s
    // y = y1*(1-r-s) + y2*r + y3*s
    // w = 2*Element Area*wref
    for (q=0; q<nq; q++) {
      r = gp[q];
      s = gp[nq+q];
      cqelm->x[q] = cvelm->x[0]*(1-r-s) + cvelm->x[1]*r + cvelm->x[2]*s;
      cqelm->y[q] = cvelm->y[0]*(1-r-s) + cvelm->y[1]*r + cvelm->y[2]*s;
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==3) {
    for (j=0; j<v_per_elm; j++) {
      cvelm->x[j] = mesh->cv->x[thiselm_v[j]];
      cvelm->y[j] = mesh->cv->y[thiselm_v[j]];
      cvelm->z[j] = mesh->cv->z[thiselm_v[j]];
    }
    voldim = 6.0*e_vol;

    // Get Quad Nodes and Weights
    tetquad_(gp,gw,nq1d);

    // Map to Real Triangle
    // x = x1*(1-r-s-t) + x2*r + x3*s + x4*t
    // y = y1*(1-r-s-t) + y2*r + y3*s + y4*t
    // z = z1*(1-r-s-t) + z2*r + z3*s + z4*t3*s
    // w = 6*Element Vol*wref
    for (q=0; q<nq; q++) {
      r = gp[q];
      s = gp[nq+q];
      t = gp[2*nq+q];
      cqelm->x[q] = cvelm->x[0]*(1-r-s-t) + cvelm->x[1]*r + cvelm->x[2]*s + cvelm->x[3]*t;
      cqelm->y[q] = cvelm->y[0]*(1-r-s-t) + cvelm->y[1]*r + cvelm->y[2]*s + cvelm->y[3]*t;
      cqelm->z[q] = cvelm->z[0]*(1-r-s-t) + cvelm->z[1]*r + cvelm->z[2]*s + cvelm->z[3]*t;
      cqelm->w[q] = voldim*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  if(cvelm) {
    free_coords(cvelm);
    free(cvelm);
    cvelm=NULL;
  }
  if(thiselm_v) free(thiselm_v);

  return;
}
/******************************************************************************/

/*!
 * \fn qcoordinates* get_quadrature_boundary(mesh_struct *mesh,INT nq1d,INT ed_or_f)
 *
 * \brief Computes quadrature weights and nodes for all faces (surface integral)
 *        or edges (line integral) in entire domain using nq1d quadrature nodes
 *        per 1D direction.
 *
 * \param nq1d    Number of quadrature nodes on an element in 1D direction
 * \param mesh    Mesh struct
 * \param ed_or_f Whether we do an edge integral (1) or face integral (2)
 *
 * \return cq_all Quadrature struct on boundary
 *
 * \note This will be deprecated
 */
qcoordinates* get_quadrature_boundary(mesh_struct *mesh,INT nq1d,INT ed_or_f)
{
  INT i,j;

  INT dim = mesh->dim;
  INT nedge = mesh->nedge;
  INT nface = mesh->nface;
  INT ndof=0;
  INT nq=0;

  if(ed_or_f==1) { // 1D integral -> get quadrature on edges
    ndof=nedge;
    nq = nq1d;
  } else if(ed_or_f==2) { // 2D integral -> get quadrature on faces
    ndof=nface;
    if(dim==2) { // face is an edge
      nq = nq1d;
    } else {
      nq = nq1d*nq1d;
    }
  }

  // Allocate quadrature points
  qcoordinates *cq_all = allocateqcoords_bdry(nq1d,ndof,dim,ed_or_f);
  qcoordinates *cq_surf = allocateqcoords_bdry(nq1d,1,dim,ed_or_f);

  if(ed_or_f==1) { // Edges
    for (i=0; i<ndof; i++) {
      quad_edge(cq_surf,mesh,nq1d,i);
      for (j=0; j<nq; j++) {
        cq_all->x[i*nq+j] = cq_surf->x[j];
        cq_all->y[i*nq+j] = cq_surf->y[j];
        cq_all->w[i*nq+j] = cq_surf->w[j];
        if(dim==3) {
          cq_all->z[i*nq+j] = cq_surf->z[j];
        }
      }
    }
  } else { // Faces
    for (i=0; i<ndof; i++) {
      quad_face(cq_surf,mesh,nq1d,i);
      for (j=0; j<nq; j++) {
        cq_all->x[i*nq+j] = cq_surf->x[j];
        cq_all->y[i*nq+j] = cq_surf->y[j];
        cq_all->w[i*nq+j] = cq_surf->w[j];
        if(dim==3) {
          cq_all->z[i*nq+j] = cq_surf->z[j];
        }
      }
    }
  }

  if(cq_surf) {
    free_qcoords(cq_surf);
    free(cq_surf);
    cq_surf = NULL;
  }

  return cq_all;
}
/******************************************************************************/

/*!
* \fn void quad_edge(qcoordinates *cqbdry,mesh_struct *mesh,INT nq1d,INT dof)
*
* \brief Computes quadrature weights and nodes for SINGLE Edge using nq1d
*          quadrature nodes on a line/surface.  Can be used to compute integrals on
*          1D boundaries (curves).
*
* \param nq1d    Number of quadrature nodes on an edge in 1D direction
* \param mesh    Mesh struct
* \param dof     Index of current edge
*
* \return cq_bdry Quadrature struct on edge
*
* \note An edge integral is always a 1D integral.
*
* \note This will be deprecated soon
*
*/
void quad_edge(qcoordinates *cqbdry,mesh_struct *mesh,INT nq1d,INT dof)
{
  // Flag for errors
  SHORT status;

  INT q,j; /* Loop indices */
  INT dim = mesh->dim;

  // Test for Simpson's Rule
  INT nqdum=1;
  if(nq1d==-1) {
    nqdum = -1;
    nq1d = 3;
  }

  // Get total quadrature points
  INT nq = nq1d;

  /* Coordinates of vertices of edge/face */
  coordinates *cvdof = allocatecoords(2,dim);

  // Add in for Simpson's Rule if using edges


  // Gaussian points for reference element
  // (edges: [-1,1]
  REAL* gp = (REAL *) calloc(nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r;      	/* Points on Reference Edge */

  REAL w = 0.5*mesh->ed_len[dof]; /* Jacobian = 1/2 |e| */

  // Get coordinates of vertices for given edge
  INT* thisdof_v = (INT *) calloc(2,sizeof(INT));
  iCSRmat* dof_v = NULL;
  dof_v = mesh->ed_v;
  get_incidence_row(dof,dof_v,thisdof_v);
  if(dim==2) {
    for (j=0; j<2; j++) {
      cvdof->x[j] = mesh->cv->x[thisdof_v[j]];
      cvdof->y[j] = mesh->cv->y[thisdof_v[j]];
    }
  } else if(dim==3) {
    for (j=0; j<2; j++) {
      cvdof->x[j] = mesh->cv->x[thisdof_v[j]];
      cvdof->y[j] = mesh->cv->y[thisdof_v[j]];
      cvdof->z[j] = mesh->cv->z[thisdof_v[j]];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get Quad Nodes and Weights
  if(nqdum==-1) { // Simpsons Rule
    gp[0] = -1.0;
    gp[1] = 0.0;
    gp[2] = 1.0;
    gw[0] = 1.0/3.0;
    gw[1] = 4.0/3.0;
    gw[2] = gw[0];
  } else {
    quad1d(gp,gw,nq1d);
  }

  // Map to Real Edge
  // Edges: x = 0.5(x1*(1-r) + x2*(1+r))
  //        y = 0.5(y1*(1-r) + y2*(1+r))
  //        z = 0.5(z1*(1-r) + z2*(1+r))
  //        w = 0.5*Edge Length*wref
  if(dim==2) {
    for (q=0; q<nq1d; q++) {
      r = gp[q];
      cqbdry->x[q] = 0.5*cvdof->x[0]*(1-r) + 0.5*cvdof->x[1]*(1+r);
      cqbdry->y[q] = 0.5*cvdof->y[0]*(1-r) + 0.5*cvdof->y[1]*(1+r);
      cqbdry->w[q] = w*gw[q];
    }
  } else if(dim==3) {
    for (q=0; q<nq1d; q++) {
      r = gp[q];
      cqbdry->x[q] = 0.5*cvdof->x[0]*(1-r) + 0.5*cvdof->x[1]*(1+r);
      cqbdry->y[q] = 0.5*cvdof->y[0]*(1-r) + 0.5*cvdof->y[1]*(1+r);
      cqbdry->z[q] = 0.5*cvdof->z[0]*(1-r) + 0.5*cvdof->z[1]*(1+r);
      cqbdry->w[q] = w*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  if(cvdof) {
    free_coords(cvdof);
    free(cvdof);
    cvdof=NULL;
  }
  if(thisdof_v) free(thisdof_v);

  return;
}
/******************************************************************************/

/*!
* \fn void quad_face(qcoordinates *cqbdry,mesh_struct *mesh,INT nq1d,INT dof)
*
* \brief Computes quadrature weights and nodes for SINGLE Face using nq1d^(dim-1)
*          quadrature nodes on a line/surface.  Can be used to compute integrals on
*          1D/2D boundaries (curves/surfaces).
*
* \param nq1d    Number of quadrature nodes on a face in 1D direction
* \param mesh    Mesh struct
* \param dof     Index of current face
*
* \return cq_bdry Quadrature struct on face
*
* \note A face integral is a 1D integral in 2D and a 2D integral in 3D.
*
* \note This will be deprecated soon.
*
*/
void quad_face(qcoordinates *cqbdry,mesh_struct *mesh,INT nq1d,INT dof)
{
  // Flag for errors
  SHORT status;

  INT q,j; /* Loop indices */
  INT dim = mesh->dim;
  INT nq = 0;

  if(dim==2) { // face is edge in 2D
    nq = nq1d;
  } else if(dim==3) {
    nq = nq1d*nq1d;
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  /* Coordinates of vertices of face */
  coordinates *cvdof = allocatecoords(dim,dim);

  // Gaussian points for reference element
  // (edges: [-1,1] faces: tri[(0,0),(1,0),(0,1)])
  REAL* gp = (REAL *) calloc((dim-1)*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r,s;      	/* Points on Reference Face */

  REAL w = 0.0;
  if(dim==2) { // Faces are Edges in 2D
    w = 0.5*mesh->f_area[dof]; /* Jacobian = 1/2 |e| */
  } else if(dim==3) {
    w = 2*mesh->f_area[dof]; /* Jacobian = 2*|f| */
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get coordinates of vertices for given edge/face
  INT* thisdof_v = (INT *) calloc(dim,sizeof(INT));
  iCSRmat* dof_v = NULL;
  dof_v = mesh->f_v;
  get_incidence_row(dof,dof_v,thisdof_v);
  if(dim==2) {
    for (j=0; j<dim; j++) {
      cvdof->x[j] = mesh->cv->x[thisdof_v[j]];
      cvdof->y[j] = mesh->cv->y[thisdof_v[j]];
    }
  } else if(dim==3){
    for (j=0; j<dim; j++) {
      cvdof->x[j] = mesh->cv->x[thisdof_v[j]];
      cvdof->y[j] = mesh->cv->y[thisdof_v[j]];
      cvdof->z[j] = mesh->cv->z[thisdof_v[j]];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get Quad Nodes and Weights
  if(dim==2) { // face is an edge
    quad1d(gp,gw,nq1d);
  } else if(dim==3) { // face is a face
    triquad_(gp,gw,nq1d);
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Map to Real Face
  if(dim==2) {
    // Edges: x = 0.5(x1*(1-r) + x2*(1+r))
    //        y = 0.5(y1*(1-r) + y2*(1+r))
    //        z = 0.5(z1*(1-r) + z2*(1+r))
    //        w = 0.5*Edge Length*wref
    for (q=0; q<nq1d; q++) {
      r = gp[q];
      cqbdry->x[q] = 0.5*cvdof->x[0]*(1-r) + 0.5*cvdof->x[1]*(1+r);
      cqbdry->y[q] = 0.5*cvdof->y[0]*(1-r) + 0.5*cvdof->y[1]*(1+r);
      cqbdry->w[q] = w*gw[q];
    }
  } else if(dim==3) {
    // Faces: x = x1*(1-r-s) + x2*r + x3*s
    //        y = y1*(1-r-s) + y2*r + y3*s
    //        z = z1*(1-r-s) + z2*r + z3*s
    //        w = 2*Element Area*wref
    for (q=0; q<nq1d*nq1d; q++) {
      r = gp[q];
      s = gp[nq1d*nq1d+q];
      cqbdry->x[q] = cvdof->x[0]*(1-r-s) + cvdof->x[1]*r + cvdof->x[2]*s;
      cqbdry->y[q] = cvdof->y[0]*(1-r-s) + cvdof->y[1]*r + cvdof->y[2]*s;
      cqbdry->z[q] = cvdof->z[0]*(1-r-s) + cvdof->z[1]*r + cvdof->z[2]*s;
      cqbdry->w[q] = w*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  if(cvdof) {
    free_coords(cvdof);
    free(cvdof);
    cvdof=NULL;
  }
  if(thisdof_v) free(thisdof_v);

  return;
}
/******************************************************************************/

/*!
 * \fn void dump_qcoords(qcoordinates *q)
 *
 * \brief Dump the quadrature data to file for plotting purposes
 *
 * \param q           Quadrature struct
 *
 * \return qcoord.dat File with quadrature in format: qcoord(nq,dim+1)
 *
 */
void dump_qcoords(qcoordinates *q)
{
  // Loop indices
  INT i;

  FILE* fid = fopen("output/qcoord.dat","w");

  if (q->z) {
    for (i=0; i<q->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\t%25.16e\t%25.16e\n",q->x[i],q->y[i],q->z[i],q->w[i]);
    }
  } else if (q->y) {
    for (i=0; i<q->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\t%25.16e\n",q->x[i],q->y[i],q->w[i]);
    }
  } else {
    for (i=0; i<q->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\n",q->x[i],q->w[i]);
    }
  }

  fclose(fid);

  return;
}
/******************************************************************************/

/*!
* \fn REAL integrate_elm(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT elm)
*
* \brief Integrate a given scalar function over an element
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
* \param elm    Element to integrate over (assumes counting at 0)
*
* \return integral Integral of scalar function over element
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*/
REAL integrate_elm(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT elm)
{

  // Loop indices
  INT quad;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim>1) qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*expr)(uval,qx,time,&(mesh->el_flag[elm]));
      integral += w*uval[comp];
    }
  } else { // assemble quadrature again
    qcoordinates *cqelm = allocateqcoords(nq1d,1,mesh->dim);
    quad_elm(cqelm,mesh,nq1d,elm);
    for (quad=0;quad<cqelm->nq_per_elm;quad++) {
      qx[0] = cqelm->x[quad];
      qx[1] = cqelm->y[quad];
      if(mesh->dim==3) qx[2] = cqelm->z[quad];
      w = cqelm->w[quad];
      (*expr)(uval,qx,time,&(mesh->el_flag[elm]));
      integral += w*uval[comp];
    }
    free_qcoords(cqelm);
    free(cqelm);
  }

  if(qx) free(qx);
  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_domain(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time)
*
* \brief Integrate a given scalar function over the entire mesh
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
*
* \return integral Integral of scalar function over domain
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*/
REAL integrate_domain(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time)
{
  // Loop indices
  INT elm;

  // Integral to return
  REAL integral = 0.0;
  // Loop over all elements and call integrate_elm
  for(elm=0;elm<mesh->nelm;elm++) {
    integral += integrate_elm(expr,nun,comp,nq1d,cq,mesh,time,elm);
  }

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_face(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT face)
*
* \brief Integrate a given scalar function over a face
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
* \param face   Face to integrate over (assumes counting at 0)
*
* \return integral Integral of scalar function over face
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*
* \note In 3D, this is an area (2D) integral.  In 2D, this is a line (1D) integral.
*
*/
REAL integrate_face(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT face)
{
  // Loop indices
  INT quad;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];
      (*expr)(uval,qx,time,&(mesh->f_flag[face]));
      integral += w*uval[comp];
    }
  } else { // assemble quadrature again
    qcoordinates *cqface = allocateqcoords_bdry(nq1d,1,mesh->dim,2);
    quad_face(cqface,mesh,nq1d,face);
    for (quad=0;quad<cqface->nq_per_elm;quad++) {
      qx[0] = cqface->x[quad];
      qx[1] = cqface->y[quad];
      if(mesh->dim==3) qx[2] = cqface->z[quad];
      w = cqface->w[quad];
      (*expr)(uval,qx,time,&(mesh->f_flag[face]));
      integral += w*uval[comp];
    }
    free_qcoords(cqface);
    free(cqface);
  }

  if(qx) free(qx);
  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_edge(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT edge)
*
* \brief Integrate a given scalar function over an edge
*
* \param expr   Function to be integrated
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
* \param edge   Edge to integrate over (assumes counting at 0)
*
* \return integral Integral of scalar function over edge
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*
* \note This is a line integral in any dimension
*
*/
REAL integrate_edge(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT edge)
{
  // Loop indices
  INT quad;

  // Function at quadrature node
  REAL* uval = (REAL *) calloc(nun,sizeof(REAL));;

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(mesh->dim,sizeof(REAL));

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[edge*cq->nq_per_elm+quad];
      qx[1] = cq->y[edge*cq->nq_per_elm+quad];
      if(mesh->dim==3) qx[2] = cq->z[edge*cq->nq_per_elm+quad];
      w = cq->w[edge*cq->nq_per_elm+quad];
      (*expr)(uval,qx,time,&(mesh->ed_flag[edge]));
      integral += w*uval[comp];
    }
  } else { // assemble quadrature again
    qcoordinates *cqedge = allocateqcoords_bdry(nq1d,1,mesh->dim,1);
    quad_edge(cqedge,mesh,nq1d,edge);
    for (quad=0;quad<cqedge->nq_per_elm;quad++) {
      qx[0] = cqedge->x[quad];
      qx[1] = cqedge->y[quad];
      if(mesh->dim==3) qx[2] = cqedge->z[quad];
      w = cqedge->w[quad];
      (*expr)(uval,qx,time,&(mesh->ed_flag[edge]));
      integral += w*uval[comp];
    }
    free_qcoords(cqedge);
    free(cqedge);
  }

  if(qx) free(qx);
  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_edge_vector_tangent(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT edge)
*
* \brief Integrate the tangential component of a given vector function along an edge
*
* \param expr   Function to be integrated (vector size dim)
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate (this is the start of the vector)
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
* \param edge   Edge to integrate over (assumes counting at 0)
*
* \return integral Integral of vector function along the edge
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*
* \note This is a tangential line integral in any dimension
*
*/
REAL integrate_edge_vector_tangent(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT edge)
{
  // Loop indices
  INT quad;
  INT j;
  INT dim = mesh->dim;

  // Function at quadrature node
  REAL* uval= (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[edge*cq->nq_per_elm+quad];
      qx[1] = cq->y[edge*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[edge*cq->nq_per_elm+quad];
      w = cq->w[edge*cq->nq_per_elm+quad];
      (*expr)(uval,qx,time,&(mesh->ed_flag[edge]));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*mesh->ed_tau[edge*dim+j];
    }
  } else { // assemble quadrature again
    qcoordinates *cqedge = allocateqcoords_bdry(nq1d,1,dim,1);
    quad_edge(cqedge,mesh,nq1d,edge);
    for (quad=0;quad<cqedge->nq_per_elm;quad++) {
      qx[0] = cqedge->x[quad];
      qx[1] = cqedge->y[quad];
      if(dim==3) qx[2] = cqedge->z[quad];
      w = cqedge->w[quad];
      (*expr)(uval,qx,time,&(mesh->ed_flag[edge]));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*mesh->ed_tau[edge*dim+j];
    }
    free_qcoords(cqedge);
    free(cqedge);
  }

  if(qx) free(qx);
  if(uval) free(uval);

  return integral;
}
/******************************************************************************/

/*!
* \fn REAL integrate_face_vector_normal(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT face)
*
* \brief Integrate the normal component of a given vector function across a face
*
* \param expr   Function to be integrated (vector size dim)
* \param nun    Number of unknowns (all scalar components) in expr (could be from block function)
* \param comp   Index of unknown we want to integrate (start of vector)
* \param nq1d   Number of quadrature points per direction (2*nq1d-1 is order of quadrature)
* \param cq     Precomputed quadrature points and weights
* \param mesh   Mesh Information
* \param time   If needed for function
* \param face   face to integrate over (assumes counting at 0)
*
* \return integral Integral of vector function across the face
*
* \note If cq is given, we will just use these precomputed values instead of reallocating the quadrature.
*       Otherwise, we will allocate a new set of quadrature based on nq1d
*
*/
REAL integrate_face_vector_normal(void (*expr)(REAL *,REAL *,REAL,void *),INT nun,INT comp,INT nq1d,qcoordinates *cq,mesh_struct *mesh,REAL time,INT face)
{
  // Loop indices
  INT quad;
  INT j;
  INT dim = mesh->dim;

  // Function at quadrature node
  REAL* uval= (REAL *) calloc(nun,sizeof(REAL));

  // Integral value
  REAL integral = 0.0;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Quadrature on elm
  if(cq) { // assuming quadrature is given
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];
      (*expr)(uval,qx,time,&(mesh->f_flag[face]));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*mesh->f_norm[face*dim+j];
    }
  } else { // assemble quadrature again
    qcoordinates *cqface = allocateqcoords_bdry(nq1d,1,dim,2);
    quad_face(cqface,mesh,nq1d,face);
    for (quad=0;quad<cqface->nq_per_elm;quad++) {
      qx[0] = cqface->x[quad];
      qx[1] = cqface->y[quad];
      if(dim==3) qx[2] = cqface->z[quad];
      w = cqface->w[quad];
      (*expr)(uval,qx,time,&(mesh->f_flag[face]));
      for(j=0; j<dim; j++) integral += w*uval[comp+j]*mesh->f_norm[face*dim+j];
    }
    free_qcoords(cqface);
    free(cqface);
  }

  if(qx) free(qx);
  if(uval) free(uval);

  return integral;
}
