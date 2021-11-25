/*! \file src/fem/basis.c
*
* \brief Compute the basis functions for triangles or tetrahedra or 1D FEM
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/1/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
* \note Typically this involves DOF defined on either the
*  vertices, edges, or faces.  In most cases, the basis elements are
*  defined using the standard Lagrange finite-element basis functions
*  using barycentric coordinates.  See PX_basis for details on the Lagrange
*  basis functions.
*
* \note Updated on 9/26/2018 for 0-1 fix.
* \note This version will assume only local data is provided, though versions
*       that take in global mesh data will be kept for posterity.
*
*/

#include "hazmath.h"

/*!
* \fn void P1_basis_ref(REAL *lam,REAL *dlam,REAL *x,INT dim)
*
* \brief Compute Standard Lagrange Finite Element Basis Functions (P1) at a particular point
*        in the reference element.  This will be used to construct all element types
*        on the physical elements later.
*
* \param x       Coordinate on where to compute basis function (on ref element)
* \param dim     Dimension of problem
*
* \return lam      Basis functions (1 for each vertex on element)
* \return dlam     Derivatives of basis functions (i.e., gradient)
*
*  \note Functions are defined as follows:
*
*        lam0 = 1 - x_1 - x_2 - ... - x_dim
*        lami = x_i      i = 1,...,dim
*        dlam0 = (-1,-1,...,-1)^T
*        dlami = e_i
*
*/
void P1_basis_ref(REAL *lam,REAL *dlam,REAL *x,INT dim)
{

  // Loop Counters
  INT idim,jdim;

  // Allocate here
  lam = (REAL *) calloc(dim+1,sizeof(REAL));
  dlam = (REAL *) calloc((dim+1)*dim,sizeof(REAL));
  lam[0] = 1.0;
  for(idim=0;idim<dim;idim++) {
    lam[0] -= x[idim];
    dlam[idim] = -1.0;
  }
  for(idim=1;idim<dim+1;idim++) {
    lam[idim] = x[idim-1];
    for(jdim=0;jdim<dim;jdim++) {
      if(idim-1==jdim) {
        dlam[idim*dim+jdim] = 1.0;
      } else {
        dlam[idim*dim+jdim] = 0.0;
      }
    }
  }

  return;
}
/*****************************************************************************/

/*!
* \fn void P1_basis_physical(REAL *lam,REAL *dlam,REAL *x,simplex_local_data* elm_data)
*
* \brief Compute Standard Lagrange Finite Element Basis Functions (P1) at a
*        particular point in the physical element given the local element data.
* \note  This function will not be used often, as we usually just compute the
*        bases on quadrature points, which are fixed ahead of time.  But in case
*        you need to create a basis function at a random point in an element
*        (not associated with quadrature), you have this option.
*
* \param x         Coordinate on where to compute basis function (on ref element)
* \param elm_data  Local mesh data on element
*
* \return lam      Basis functions (1 for each vertex on element)
* \return dlam     Derivatives of basis functions (i.e., gradient)
*
*  \note Functions are computed by mapping from reference simplex directly:
*        B*xr = x - x0 => xr = B^{-1} (x-xv0)
*        B is included in elm_data as ref_map
*        2D Example
*
*      ( xv1-xv0  xv2-xv0 ) * ( xr ) = ( x - xv0 )
*      ( yv1-yv0  yv2-yv0 )   ( yr )   ( y - yv0 )
*
*      lam = 1-xr-yr; xr; yr
*      dlam = rows of Binv map provided in gradlams of elm_data
*
*/
void P1_basis_physical(REAL *lam,REAL *dlam,REAL *x,simplex_local_data* elm_data)
{

  // Loop counters
  INT i,j;
  // Temporarily grab stuff from elm_data
  INT dim = elm_data->dim;
  REAL* binv = elm_data->gradlams;
  REAL* xv = elm_data->xv;

  // Get coordinate of x on ref_elm using inverse reference map
  // Note that we ignore first rwo of B^(-1), since this is grad(lam0) and not part
  // of the original B^(-1)
  REAL* xr = (REAL *) calloc(dim,sizeof(REAL));
  for(i=0;i<dim;i++) {
    xr[i] = 0.0;
    for(j=0;j<dim;j++) {
      xr[i] += binv[(i+1)*dim+j]*(x[j]-xv[0*dim+j]);
    }
  }

  // Get lambdas -> lam(x) = lamr(xr)
  lam[0] = 1.0;
  for(i=0;i<dim;i++) {
    lam[0] -= xr[i];
    lam[i+1] = xr[i];
  }

  // Get gradients
  // dlam[0] = sum of rows of B^{-1}
  // Then rest is just B^{-1}
  REAL dlam0val;
  for(i=0;i<dim;i++) {
    dlam0val = 0.0;
    for(j=0;j<dim;j++) {
      dlam0val += binv[j*dim+i];
      dlam[(i+1)*dim+j] = binv[i*dim+j];
    }
    dlam[i] = dlam0val;
  }

  if(xr) free(xr);

  return;
}
/*****************************************************************************/

/*!
* \fn void compute_refelm_mapping(REAL* ref_map,REAL* lamgrads,REAL *xv,INT dim)
*
* \brief Compute the mappings to move from reference element to physical elements
*
* \param xv      Coordinates of vertices on element
* \param dim     Dimension of problem
*
* \return ref_map  B matrix that defines the mapping: B*xr = x - x0
* \return lamgrads Special matrix that contains the gradients of the P1 basis functions as rows
*
*  \note Assuming physical coordinates xv, let x be coordinate on physical element
*        and xr be coordinate on reference element, then,
*        we define the mapping as B*xr = x - x0 (where x0 is the corner point of simplex)
*        or
*      (    |        |    ) * ( |  ) = (   |   )
*      ( xv1-xv0  xv2-xv0 )   ( xr )   ( x-xv0 )
*      (    |        |    )   ( |  )   (   |   )
*        Note that lam(x) = hat(lam)(B^(-1)(x)), so after chain rule
*
*        B^{-1}^T grad(hat(lam)) = grad(lam)
*
*        But hat(lam0) = 1 - x_0 - x_1 - ... - x_d
*            hat(lami) = x_i i = 1,..,d
*        so
*        B^{-1}^T grad(hat(lam0)) = (-(1,1,...1)*B^{-1})^T = sumofrows(B^{-1}) = grad(lam0)
*        B^{-1}^T hat(lami) = B^{-1}^T ei = grad(lami) = (B^{-1})_rowi
*
*        lamgrads stores the following then,
*        grad(lam)(x) = ( --- grad(lam0) --- ) =  ( sumrows(B^{-1}) )
*                     = ( --- grad(lam1) --- ) =  (   B^{-1)_row1   )
*                     = (         |          ) =  (      |          )
*                     = ( --- grad(lamd) --- ) =  (   B^{-1}_rowd   )
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
*    |        2
*    |       / \
*    |      /   \
*    Y    e20    e12
*    |    /       \
*    |   /         \
*    |  0----e01-----1
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
*    1  2
*    |  |\
*    |  | \
*    S e20 e12
*    |  |   \
*    |  |    \
*    0  0-e01-1
*    |
*    +--0--R--1-------->
*
*     Determine the (R,S) coordinates corresponding to (X,Y).
*
*     What is happening here is that we are solving the linear system:
*
*      ( X1-X0  X2-X0 ) * ( R ) = ( X - X0 )
*      ( Y1-Y0  Y2-Y0 )   ( S )   ( Y - Y0 )
*
*     by computing the inverse of the coefficient matrix and multiplying
*     it by the right hand side to get R and S.
*
*    The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
*    for R and S.
*
*/
void compute_refelm_mapping(REAL* ref_map,REAL* lamgrads,REAL *xv,INT dim)
{

  // Loop counters
  INT i,j;

  // Get B
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      ref_map[i*dim+j] = xv[(j+1)*dim+i]-xv[0];
    }
  }

  // Invert B
  REAL* binv = (REAL *) calloc(dim*dim,sizeof(REAL));
  void* wrk = malloc(dim*sizeof(INT)+dim*(dim+1)*sizeof(REAL));
  ddense_inv(binv,dim,ref_map,wrk);

  // Place all grad(lam) in one matrix as rows
  // First row is sum of rows of B-1
  // Then rest is just B^{-1}
  REAL dlam0val;
  for(i=0;i<dim;i++) {
    dlam0val = 0.0;
    for(j=0;j<dim;j++) {
      dlam0val += binv[j*dim+i];
      lamgrads[(i+1)*dim+j] = binv[i*dim+j];
    }
    lamgrads[i] = dlam0val;
  }

  if(binv) free(binv);

  return;
}
/******************************************************************************/

/*!
* \fn void PX_basis(REAL *p,REAL *dp,INT porder,INT dim,REAL* lam,REAL* dlam)
*
* \brief Compute Standard Lagrange Finite Element Basis Functions (PX) at a particular point in 1, 2 or 3D
*        For now, we only assume constants, Linears, or Quadratic Elements (P0 or P1 or P2)
*
* \param porder   Order of elements
* \param dim      Dimension of problem
* \param lam,dlam P1 basis functions at quadrature point
*
* \return p      Basis functions (1 for each DOF on element)
* \return dp     Derivatives of basis functions (i.e., gradient)
*
* \note For P0 elements, we just return p=1
*       For P1 elements, we copy the given input Functions:
*          lam0 = 1 - x_1 - x_2 - ... - x_dim
*          lam_i = x_i   i=1,...dim
*       For P2 elements we get
*          lam_i = 2*lam_i(lam_i - 1/2)  i = 0,...,dim
*          lam_k = 4*lam_i*lam_j         k = dim+1,...,dim+1 + (dim*(dim+1))/2, i = 0,...dim-1, j=i+1,...dim
*
*
*    For quadratic elements:
*
*    |
*    1  2
*    |  |\
*    |  | \
*    S  4  5
*    |  |   \
*    |  |    \
*    0  0--3--1
*    |
*    +--0--R--1-------->
*
*/
void PX_basis(REAL *p,REAL *dp,INT porder,INT dim,REAL* lam, REAL* dlam)
{

  // Loop Counters
  INT idim,jdim,kdim,cntr;

  // Flag for Errors
  short status;

  switch(porder) {
    case 0: // P0 elements are trivial and we just need to return a 1 for each element:
    p[0] = 1.0;
    break;

    case 1: // P1 elements - just copy from input
    for(idim=0;idim<dim+1;idim++) {
      p[idim] = lam[idim];
      for(jdim=0;jdim<dim;jdim++) {
        dp[idim*dim+jdim] = dlam[idim*dim+jdim];
      }
    }
    break;

    case 2: // P2 elements - add edges dof
    cntr=dim+1;
    for(idim=0;idim<dim+1;idim++) {
      p[idim] = 2*lam[idim]*(lam[idim]-0.5);
      for(jdim=0;jdim<dim;jdim++) {
        dp[idim*dim+jdim] = 2*(dlam[idim*dim+jdim]*(lam[idim]-0.5) + lam[idim]*dlam[idim*dim+jdim]);
      }
      for(jdim=idim;jdim<dim;jdim++) {
        p[cntr] = 4*lam[idim]*lam[jdim+1];
        for(kdim=0;kdim<dim;kdim++) {
          dp[cntr*dim+kdim] = 4*(dlam[idim*dim+kdim]*lam[jdim+1] + lam[idim]*dlam[(jdim+1)*dim+kdim]);
        }
        cntr++;
      }
    }
    break;

    default:
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
    break;
  }
  return;
}
/******************************************************************************/

/*!
* \fn void ned0_basis(REAL *phi,REAL *cphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_ed,REAL* ed_len)
*
* \brief Compute Nedelec Finite Element Basis Functions (zeroth order) at a particular point in 2 or 3D
*
* \param lam,dlam  P1 bases at quadrature point
* \param dim       Dimension of problem
* \param v_on_elm  Vertices on element
* \param v_on_ed   Vertices on each edge (DoF) on the elements (dim+1 edges with 2 vertices each -> (dim+1 x 2) matrix)
* \param ed_len    Lenght of each edge on element
*
* \return phi      Basis functions (dim for each edge from reference triangle)
* \return cphi     Curl of basis functions (1 for each edge in 2D, dim for each edge in 3D)
*
*/
void ned0_basis(REAL *phi,REAL *cphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_ed,REAL* ed_len)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = dim+1;
  INT ed_per_elm = dim*(dim+1)/2;
  INT v_per_ed = 2;

  INT i,k,ve1,ve2,ihi,ilo;
  INT mark1 = -1;
  INT mark2 = -1;

  REAL elen;

  /* Now, with the linear basis functions p, dpx, and dpy, the nedelec elements are
  * phi_eij = |eij|*(p(i)grad(p(j)) - p(j)grad(p(i)))
  * |eij| = sqrt(|xj-xi|^2)
  */

  // Go through each edge and get length and find the corresponding nodes
  for (i=0; i<ed_per_elm; i++) {
    ve1 = v_on_ed[i*v_per_ed+0];
    ve2 = v_on_ed[i*v_per_ed+1];
    elen = ed_len[i];

    // Find out which linear basis elements line up with nodes on this edge
    for (k=0; k<v_per_elm; k++) {
      if (v_on_elm[k]==ve1) {
        mark1=k;
      }
      if (v_on_elm[k]==ve2) {
        mark2=k;
      }
    }
    // Make sure orientation is correct always go from i->j if nj > ni
    if (MAX(ve1,ve2)==ve1) {
      ihi = mark1;
      ilo = mark2;
    } else {
      ihi = mark2;
      ilo = mark1;
    }

    phi[i*dim+0] = elen*(lam[ilo]*dlam[ihi*dim] - lam[ihi]*dlam[ilo*dim]);
    phi[i*dim+1] = elen*(lam[ilo]*dlam[ihi*dim+1] - lam[ihi]*dlam[ilo*dim+1]);
    if(dim==3) phi[i*dim+2] = elen*(lam[ilo]*dlam[ihi*dim+2] - lam[ihi]*dlam[ilo*dim+2]);

    /* Now compute Curls
    * In 2D curl v = (-dy,dx)*(v1,v2)^T = (dx,dy)(0 1;-1 0)(v1,v2)^T = div (Jv)
    * curl phi_eij = |eij|*(grad(p(i))*(J*grad(p(j)))-grad(p(j))*(J*grad(p(i)))
    * This results from the fact that the p's are linear...
    *
    * In 3D, technically the curls are not needed in 3D as <curl u, curl v> operator can be found from Laplacian matrix.
    * We compute them anyway
    */

    if(dim==2) {
      cphi[i] = 2*elen*(dlam[ilo*dim]*dlam[ihi*dim+1] - dlam[ilo*dim+1]*dlam[ihi*dim]);
    } else if(dim==3) {
      cphi[i*dim+0] = 2*elen*(dlam[ilo*dim+1]*dlam[ihi*dim+2]-dlam[ihi*dim+1]*dlam[ilo*dim+2]);
      cphi[i*dim+1] = 2*elen*(dlam[ihi*dim]*dlam[ilo*dim+2]-dlam[ilo*dim]*dlam[ihi*dim+2]);
      cphi[i*dim+2] = 2*elen*(dlam[ilo*dim]*dlam[ihi*dim+1]-dlam[ihi*dim]*dlam[ilo*dim+1]);
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }
  }

  return;
}
/******************************************************************************/

/*!
* \fn void rt0_basis_local(REAL *phi,REAL *dphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* f_area)
*
* \brief Compute Raviart-Thomas Finite Element Basis Functions (zeroth order) at a particular point in 2 or 3D
*
* \param lam,dlam  P1 bases at quadrature point
* \param dim       Dimension of problem
* \param v_on_elm  Vertices on element
* \param v_on_face   Vertices on each face (DoF) ordered of element (dim+1 faces per elm with dim vertices per face -> (dim+1)xdim matrix)
* \param f_area    Area of each face on element
*
* \return phi      Basis functions (dim for each face from reference triangle)
* \return dphi     Div of basis functions (1 for each face)
*
*/
void rt0_basis(REAL *phi,REAL *dphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* f_area)
{
  // Flag for erros
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = dim+1;
  INT f_per_elm = dim+1;
  INT v_per_f = dim;

  INT i,j;
  INT vel,ef1,ef2,ef3;

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {

      // Loop through vertices on element to find corresponding verties on face and get correct orientation
      for(j=0;j<v_per_elm;j++) {
        vel = v_on_elm[j];
        if(v_on_face[i*v_per_f+0]==vel) {
          ef1 = j;
        }
        if(v_on_face[i*v_per_f+1]==vel) {
          ef2 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 2D
      * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i)))
      * |fij| = |eij|
      */
      phi[i*dim+0] = f_area[i]*(lam[ef1]*dlam[ef2*dim+1] - lam[ef2]*dlam[ef1*dim+1]);
      phi[i*dim+1] = f_area[i]*(-lam[ef1]*dlam[ef2*dim] + lam[ef2]*dlam[ef1*dim]);

      // Compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i)))
      dphi[i] = 2*f_area[i]*(dlam[ef1*dim]*dlam[ef2*dim+1] - dlam[ef2*dim]*dlam[ef1*dim+1]);
    }
  } else if(dim==3) {
    for (i=0; i<f_per_elm; i++) {

      // Loop through Nodes on element to find corresponding nodes for correct orienation
      for(j=0;j<v_per_elm;j++) {
        vel = v_on_elm[j];
        if(v_on_face[i*v_per_f+0]==vel) {
          ef1 = j;
        }
        if(v_on_face[i*v_per_f+1]==vel) {
          ef2 = j;
        }
        if(v_on_face[i*v_per_f+2]==vel) {
          ef3 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the RT elements are in 3D
      * phi_fijk = 2*|fijk|*(p(i)(grad(p(j)) x grad(p(k))) + p(j)(grad(p(k)) x grad(p(i))) + p(k)(grad(p(i)) x grad(p(j))))
      * |fijk| = Area(Face)
      */
      phi[i*dim+0] = 2*f_area[i]*(lam[ef1]*(dlam[ef2*dim+1]*dlam[ef3*dim+2]-dlam[ef2*dim+2]*dlam[ef3*dim+1])
      + lam[ef2]*(dlam[ef3*dim+1]*dlam[ef1*dim+2]-dlam[ef3*dim+2]*dlam[ef1*dim+1])
      + lam[ef3]*(dlam[ef1*dim+1]*dlam[ef2*dim+2]-dlam[ef1*dim+2]*dlam[ef2*dim+1]));
      phi[i*dim+1] = 2*f_area[i]*(lam[ef1]*(dlam[ef2*dim+2]*dlam[ef3*dim]-dlam[ef2*dim]*dlam[ef3*dim+2])
      + lam[ef2]*(dlam[ef3*dim+2]*dlam[ef1*dim]-dlam[ef3*dim]*dlam[ef1*dim+2])
      + lam[ef3]*(dlam[ef1*dim+2]*dlam[ef2*dim]-dlam[ef1*dim]*dlam[ef2*dim+2]));
      phi[i*dim+2] = 2*f_area[i]*(lam[ef1]*(dlam[ef2*dim]*dlam[ef3*dim+1]-dlam[ef2*dim+1]*dlam[ef3*dim])
      + lam[ef2]*(dlam[ef3*dim]*dlam[ef1*dim+1]-dlam[ef3*dim+1]*dlam[ef1*dim])
      + lam[ef3]*(dlam[ef1*dim]*dlam[ef2*dim+1]-dlam[ef1*dim+1]*dlam[ef2*dim]));

      // Compute divs div(phi_fijk) = 6*|fijk| grad(pi) dot (grad(pj) cross grad(pk)) (
      dphi[i] = 6*f_area[i]*(dlam[ef1*dim]*(dlam[ef2*dim+1]*dlam[ef3*dim+2]-dlam[ef2*dim+2]*dlam[ef3*dim+1])
      + dlam[ef1*dim+1]*(dlam[ef2*dim+2]*dlam[ef3*dim]-dlam[ef2*dim]*dlam[ef3*dim+2])
      + dlam[ef1*dim+2]*(dlam[ef2*dim]*dlam[ef3*dim+1]-dlam[ef2*dim+1]*dlam[ef3*dim]));
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  return;
}
/******************************************************************************/

/*!
* \fn void bdm1_basis(REAL *phi,REAL *dphix,REAL *dphiy,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* farea)
* \brief Brezzi-Douglas-Marini (BDM) Elements of order 1.
*
* \note ONLY in 2D for now.
* \note This has NOT been tested.
*
* \param lam,dlam  P1 bases at quadrature point
* \param dim       Dimension of problem
* \param v_on_elm  Vertices on element
* \param v_on_face   Vertices on each face (DoF) on the elm
* \param farea    Area of each face on element
*
* \return phi      Basis functions (2*f_per_elm for each face, 12 total in 2D)
* \return dphix    Div of basis functions
*
*/
void bdm1_basis(REAL *phi,REAL *dphix,REAL *dphiy,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* farea)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = dim+1;
  INT f_per_elm = dim+1;
  INT v_per_f = dim;

  INT i,j;
  REAL a1,a2,a3,a4;
  INT vel,ef1,ef2;

  // Go through each face and find the corresponding nodes
  if(dim==2) {

    for (i=0; i<f_per_elm; i++) {
      // Loop through vertices on element to find corresponding nodes and get correct orienation
      for(j=0;j<v_per_elm;j++) {
        vel = v_on_elm[j];
        if(v_on_face[i*v_per_f+0]==vel) {
          ef1 = j;
        }
        if(v_on_face[i*v_per_f+1]==vel) {
          ef2 = j;
        }
      }

      /* Now, with the linear basis functions p, dpx, and dpy, the BDM1 elements are in 2D
      * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i)))
      * psi_fij = alpha*|fij|*curl(p(i)p(j))
      * |fij| = |eij|
      */
      phi[i*dim*2] = farea[i]*(lam[ef1]*dlam[ef2*dim+1] - lam[ef2]*dlam[ef1*dim+1]);
      phi[i*dim*2+1] = farea[i]*(-lam[ef1]*dlam[ef2*dim] + lam[ef2]*dlam[ef1*dim]);
      phi[i*dim*2+2] = -6*farea[i]*(lam[ef1]*dlam[ef2*dim+1] + lam[ef2]*dlam[ef2*dim+1]);
      phi[i*dim*2+3] = 6*farea[i]*(lam[ef1]*dlam[ef2*dim] + lam[ef2]*dlam[ef2*dim]);

      a1 = dlam[ef1*dim]*dlam[ef2*dim];
      a2 = dlam[ef1*dim]*dlam[ef2*dim+1];
      a3 = dlam[ef1*dim+1]*dlam[ef2*dim];
      a4 = dlam[ef1*dim+1]*dlam[ef2*dim+1];

      dphix[i*dim*2] = farea[i]*(a2-a3);
      dphix[i*dim*2+1] = 0.0;
      dphix[i*dim*2+2] = -6*farea[i]*(a2+a3);
      dphix[i*dim*2+3] = 12*farea[i]*a1;
      dphiy[i*dim*2] = 0.0;
      dphiy[i*dim*2+1] = farea[i]*(a2-a3);
      dphiy[i*dim*2+2] = -12*farea[i]*a4;
      dphiy[i*dim*2+3] = 6*farea[i]*(a2+a3);
    }

  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__); // 3D not implemented
  }

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
* \fn void face_bubble_basis(REAL *phi, REAL *dphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* f_area,REAL* f_norm)
*
* \brief Compute Bubble Element Finite Element Basis Functions at a particular point in 2 or 3D
*
* \param lam,dlam  P1 bases at quadrature point
* \param dim       Dimension of problem
* \param v_on_elm  Vertices on element
* \param v_on_face   Vertices on each face (DoF) on the element
* \param f_area    Area of each face on element
* \param f_norm    Normal vector of each face on element (f_per_elm*dim)
*
* \return phi      Basis functions (dim for each face from reference triangle)
* \return dphi     Tensor from gradient of basis functions
*/
void face_bubble_basis(REAL *phi, REAL *dphi,REAL* lam,REAL* dlam,INT dim,INT* v_on_elm,INT* v_on_face,REAL* f_area,REAL* f_norm)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = dim+1;
  INT dof_per_elm = dim+1;
  INT v_per_f = dim;
  INT i,j;

  // face to vertex map
  INT vel,ef1,ef2,ef3;//face endpoint vertex tracking numbers

  REAL gradp;

  if(dim==2){
    for (i=0;i<dof_per_elm;i++) {
      for(j=0;j<v_per_elm;j++){
        vel = v_on_elm[j];
        if(v_on_face[i*v_per_f+0]==vel) {
          ef1 = j;
        }
        if(v_on_face[i*v_per_f+1]==vel) {
          ef2 = j;
        }
      }

      // Multiply basis function by normal vector
      /* phi[i*dim] =   ABS(mesh->f_norm[dim*(dof[i])]  )* 4*p[ef1]*p[ef2]; */
      /* phi[i*dim+1] = ABS(mesh->f_norm[dim*(dof[i])+1]) * 4*p[ef1]*p[ef2]; */
      phi[i*dim] =   f_norm[i*dim+0]*4.0*lam[ef1]*lam[ef2];
      phi[i*dim+1] = f_norm[i*dim+1]*4.0*lam[ef1]*lam[ef2];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 4.0*(lam[ef1]*dlam[ef2*dim+j] + dlam[ef1*dim+j]*lam[ef2]);

        /* dphi[i*dim*dim + j*dim + 0] = gradp * ABS(mesh->f_norm[dim*(dof[i])+0]); */
        /* dphi[i*dim*dim + j*dim + 1] = gradp * ABS(mesh->f_norm[dim*(dof[i])+1]); */
        dphi[i*dim*dim + j*dim + 0] = gradp * f_norm[i*dim+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * f_norm[i*dim+1];
      }

    }
  } else if(dim==3) {
    for (i=0;i<dof_per_elm;i++) {
      for(j=0;j<v_per_elm;j++){
        vel = v_on_elm[j];
        if(v_on_face[i*v_per_f+0]==vel) {
          ef1 = j;
        }
        if(v_on_face[i*v_per_f+1]==vel) {
          ef2 = j;
        }
        if(v_on_face[i*v_per_f+2]==vel) {
          ef3 = j;
        }
      }

      // Multiply basis function by normal vector
      phi[i*dim+0] = f_norm[i*dim+0]*8.0*lam[ef1]*lam[ef2]*lam[ef3];
      phi[i*dim+1] = f_norm[i*dim+1]*8.0*lam[ef1]*lam[ef2]*lam[ef3];
      phi[i*dim+2] = f_norm[i*dim+2]*8.0*lam[ef1]*lam[ef2]*lam[ef3];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 8.0*(lam[ef1]*lam[ef2]*dlam[ef3*dim+j] + lam[ef1]*dlam[ef2*dim+j]*lam[ef3] + dlam[ef1*dim+j]*lam[ef2]*lam[ef3]);

        dphi[i*dim*dim + j*dim + 0] = gradp * f_norm[i*dim+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * f_norm[i*dim+1];
        dphi[i*dim*dim + j*dim + 2] = gradp * f_norm[i*dim+2];
      }
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  return;
}
/******************************************************************************/

/*!
* \fn void get_FEM_basis_on_elm(REAL *phi,REAL *dphi,simplex_local_data *simplex_data,fe_local_data *fe_data,REAL *lam, REAL* dlam,INT space_index)
*
* \brief Grabs the basis function of a given FEM space at a particular point
*        on a given element, given the local mesh/fem data on that simplex for
*        that space.  We assume the P1 basis functions are predefined at the x
*        point of interest and input.
*
* \note This is the "main" routine for grabbing basis functions.  When, new
*       spaces are added, this routine should be updated.  Other routines Will
*       call this one for instance if you want a basis function at a random point
*       versus a quadrature point.
*
* \param simplex_data    Local mesh data
* \param fe_data         Local FE data
* \param lam, dlam       P1 Basis functions defined at x on element
* \param space_index     Which FE space in fe_block we're considering
*
* \return phi      Basis functions
* \return dphi     Derivatives of basis functions (depends on type)
*
*/
void get_FEM_basis_on_elm(REAL *phi,REAL *dphi,simplex_local_data *simplex_data,fe_local_data *fe_data,REAL *lam, REAL* dlam,INT space_index)
{
  // Flag for erros
  SHORT status;
  INT i;

  // Mesh and FEM Data
  INT dim = simplex_data->dim;

  INT fe_type = fe_data->fe_types[space_index];
  INT dof_per_elm = fe_data->n_dof_per_space[space_index];

  INT dim_offset = dof_per_elm/dim;

  if(fe_type>=0 && fe_type<10) { // PX elements - only P0, P1, and P2 implemented

    PX_basis(phi,dphi,fe_type,dim,lam,dlam);

  } else if(fe_type==20) { // Nedelec elements

    ned0_basis(phi,dphi,lam,dlam,dim,simplex_data->local_v,simplex_data->v_on_ed,simplex_data->ed_len);

  } else if(fe_type==30) { // Raviart-Thomas elements

    rt0_basis(phi,dphi,lam,dlam,dim,simplex_data->local_v,simplex_data->v_on_f,simplex_data->f_area);

  } else if(fe_type==60) { // Vector element

    for(i=0;i<dim;i++) PX_basis(phi+i*dim_offset,dphi+i*dof_per_elm,fe_type,dim,lam,dlam);

  } else if(fe_type==61) { // Bubble element

    face_bubble_basis(phi,dphi,lam,dlam,dim,simplex_data->local_v,simplex_data->v_on_f,simplex_data->f_area,simplex_data->f_norm);

  } else if(fe_type==99) { // Constraint Single DOF Space

    phi[0] = 1.0;

  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
/******************************************************************************/

/*!
* \fn void get_FEM_basis_at_quadpt(simplex_local_data *simplex_data,fe_local_data *fe_data,INT space_index,INT quadpt)
*
* \brief Grabs the basis function of a given FEM space at a particular quadrature point
*        on a given element, given the local fem data on that simplex and the predefined
*        quadrature point saved in simplex_data.  Note we will store the result`
*        inside fe_data directly
*
* \param simplex_data    Local mesh data
* \param fe_data         Local FE data
* \param space_index     Which FE space in fe_block we're considering
* \param quadpt          Index of quadrature point basis is to be evaluated at
*
* \note ordering of quad point is determined by get_quadrature routine
*
* \return fe_data->phi[space_index]      Basis functions
* \return fe_data->dphi[space_index]     Derivatives of basis functions (depends on type)
*
*/
void get_FEM_basis_at_quadpt(simplex_local_data *simplex_data,fe_local_data *fe_data,INT space_index,INT quadpt)
{

  INT dim = simplex_data->dim;
  // P1 basis functions at quadpt
  REAL* lam = simplex_data->lams + quadpt*(dim+1);
  REAL* dlam = simplex_data->gradlams + quadpt*((dim+1)*dim);

  // Call general function at any x
  get_FEM_basis_on_elm(fe_data->phi[space_index],fe_data->dphi[space_index],simplex_data,fe_data,lam,dlam,space_index);

  return;
}
/******************************************************************************/

/*!
* \fn void get_FEM_basis_at_x(REAL *phi,REAL *dphi,simplex_local_data *simplex_data,fe_local_data *fe_data,INT space_index,REAL *x)
*
* \brief Grabs the basis function of a given FEM space at a particular point x
*        on a given element, given the local fem data on that simplex. Note we will store the result`
*        inside fe_data directly, overriding what was there before.
*
* \param simplex_data    Local mesh data
* \param fe_data         Local FE data
* \param space_index     Which FE space in fe_block we're considering
* \param x               Physical coordinate to compute basis function on
*
* \note This routine might not be used to often, unless you need to Interpolate
*       to a random point in the element.  This "rebuilds" the P1 basis functions
*       on the physical element, though it does use the preset reference element
*       mapping provided inside simplex_data
*
* \return fe_data->phi[space_index]      Basis functions
* \return fe_data->dphi[space_index]     Derivatives of basis functions (depends on type)
*
*/
void get_FEM_basis_at_x(simplex_local_data *simplex_data,fe_local_data *fe_data,INT space_index,REAL *x)
{

  INT dim = simplex_data->dim;
  // P1 basis functions at x
  REAL *lam = (REAL *) calloc(dim+1,sizeof(REAL));
  REAL *dlam = (REAL *) calloc((dim+1)*dim,sizeof(REAL));
  P1_basis_physical(lam,dlam,x,simplex_data);

  // Call general function at any x
  get_FEM_basis_on_elm(fe_data->phi[space_index],fe_data->dphi[space_index],simplex_data,fe_data,lam,dlam,space_index);

  if(lam) free(lam);
  if(dlam) free(dlam);
  
  return;
}
/******************************************************************************/


/********************** OLD STUFF **************************************/
/*             Uses global mesh data and is inefficient                */

/*******************************************************************************************************/
/*!
* \fn void PX_H1_basis(REAL *p,REAL *dp,REAL *x,INT *dof,INT porder,mesh_struct *mesh)
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
*        For P1: for element between x_0 and x_1:  x0 ------ x1
*             phi_0 = (x_1 - x)/(x_1 - x_0)
*             phi_1 = (x-x_0)/(x_1 - x_0)
*
*        For P2: for element between x_0 and x_1: x0 ---- x2 ---- x1
*             phi_0 = (x_1 - x)/(x_1 - x_0)*(2(x_1 - x)/(x_1 - x_0) - 1)
*             phi_1 = (x-x_0)/(x_1 - x_0)*(2(x-x_0)/(x_1 - x_0) - 1)
*             phi_2 = 4(x_1 - x)/(x_1 - x_0)(x-x_0)/(x_1 - x_0)
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
*    |        2
*    |       / \
*    |      /   \
*    Y    e20    e12
*    |    /       \
*    |   /         \
*    |  0----e01-----1
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
*    1  2
*    |  |\
*    |  | \
*    S e20 e12
*    |  |   \
*    |  |    \
*    0  0-e01-1
*    |
*    +--0--R--1-------->
*
*     Determine the (R,S) coordinates corresponding to (X,Y).
*
*     What is happening here is that we are solving the linear system:
*
*      ( X1-X0  X2-X0 ) * ( R ) = ( X - X0 )
*      ( Y1-Y0  Y2-Y0 )   ( S )   ( Y - Y0 )
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
*    1  2
*    |  |\
*    |  | \
*    S  4  5
*    |  |   \
*    |  |    \
*    0  0--3--1
*    |
*    +--0--R--1-------->
*
*/
void PX_H1_basis(REAL *p,REAL *dp,REAL *x,INT *dof,INT porder,mesh_struct *mesh)
{
  REAL dp0r,dp1r,dp2r,dp3r,dp4r,dp5r,dp6r,dp7r,dp8r,dp9r;
  REAL dp0s,dp1s,dp2s,dp3s,dp4s,dp5s,dp6s,dp7s,dp8s,dp9s;
  REAL dp0t,dp1t,dp2t,dp3t,dp4t,dp5t,dp6t,dp7t,dp8t,dp9t;
  REAL onemrst;

  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT dim = mesh->dim;

  // Store coordinates of vertices of element (2 in 1D, 3 in 2D, 4 in 3D)
  REAL xv0,xv1,xv2,xv3;
  REAL yv0,yv1,yv2,yv3;
  REAL zv0,zv1,zv2,zv3;
  coordinates* cv = mesh->cv;

  // P0 elements are trivial and we just need to return a 1 for each element:
  if(porder==0) {
    p[0] = 1.0;
  } else {
    // The remaining depend on the dimension
    if(dim==1) {

      // Get Physical Coordinates of Vertices
      xv0 = cv->x[dof[0]];
      xv1 = cv->x[dof[1]];

      // Get Barycentric Coordinates
      REAL oneoverh = 1.0/(xv1-xv0);
      REAL lam0 = (xv1-x[0])*oneoverh;
      REAL lam1 = (x[0]-xv0)*oneoverh;

      // Now Get basis functions
      if(porder==1) {
        p[0] = lam0;
        p[1] = lam1;
        dp[0] = -oneoverh;
        dp[1] = oneoverh;
      } else if(porder==2) {
        p[0] = lam0*(2*lam0-1);
        p[1] = lam1*(2*lam1-1);
        p[2] = 4*lam0*lam1;
        dp[0] = -oneoverh*(2*lam0-1) - 2*lam0*oneoverh;
        dp[1] = oneoverh*(2*lam1-1) + 2*lam1*oneoverh;
        dp[2] = 4*lam0*oneoverh - 4*lam1*oneoverh;
      } else {
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }
    } else if(dim==2) {

      // Get Physical Coordinates of Vertices
      xv0 = cv->x[dof[0]];
      xv1 = cv->x[dof[1]];
      xv2 = cv->x[dof[2]];
      yv0 = cv->y[dof[0]];
      yv1 = cv->y[dof[1]];
      yv2 = cv->y[dof[2]];

      // Get coordinates on reference triangle
      REAL det = (xv1-xv0)*(yv2-yv0) - (xv2-xv0)*(yv1-yv0);

      REAL r = ((yv2-yv0)*(x[0]-xv0) + (xv0-xv2)*(x[1]-yv0))/det;

      REAL drdx = (yv2-yv0)/det;
      REAL drdy = (xv0-xv2)/det;

      REAL s = ((yv0-yv1)*(x[0]-xv0) + (xv1-xv0)*(x[1]-yv0))/det;

      REAL dsdx = (yv0-yv1)/det;
      REAL dsdy = (xv1-xv0)/det;

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
        dp0r = -1;
        dp1r = 1;
        dp2r = 0;
        dp0s = -1;
        dp1s = 0;
        dp2s = 1;
      } else if(porder==2) {
        p[0] = 2*onemrst*(onemrst-0.5);
        p[1] = 2*r*(r-0.5);
        p[2] = 2*s*(s-0.5);
        p[3] = 4*r*onemrst;
        p[4] = 4*s*onemrst;
        p[5] = 4*r*s;
        dp0r = 4*r+4*s-3;
        dp1r = 4*r - 1;
        dp2r = 0;
        dp3r = 4-8*r-4*s;
        dp4r = -4*s;
        dp5r = 4*s;
        dp0s = dp0r;
        dp1s = 0;
        dp2s = 4*s-1;
        dp3s = -4*r;
        dp4s = 4-4*r-8*s;
        dp5s = 4*r;
      } else {
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }

      /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
      *  to (X,Y) using the chain rule.
      */
      dp[0*dim] = dp0r * drdx + dp0s * dsdx;
      dp[0*dim+1] = dp0r * drdy + dp0s * dsdy;
      dp[1*dim] = dp1r * drdx + dp1s * dsdx;
      dp[1*dim+1] = dp1r * drdy + dp1s * dsdy;
      dp[2*dim] = dp2r * drdx + dp2s * dsdx;
      dp[2*dim+1] = dp2r * drdy + dp2s * dsdy;
      if(porder==2) {
        dp[3*dim] = dp3r * drdx + dp3s * dsdx;
        dp[3*dim+1] = dp3r * drdy + dp3s * dsdy;
        dp[4*dim] = dp4r * drdx + dp4s * dsdx;
        dp[4*dim+1] = dp4r * drdy + dp4s * dsdy;
        dp[5*dim] = dp5r * drdx + dp5s * dsdx;
        dp[5*dim+1] = dp5r * drdy + dp5s * dsdy;
      }
    } else if (dim==3) {

      // Get Nodes and Physical Coordinates
      xv0 = cv->x[dof[0]];
      xv1 = cv->x[dof[1]];
      xv2 = cv->x[dof[2]];
      xv3 = cv->x[dof[3]];
      yv0 = cv->y[dof[0]];
      yv1 = cv->y[dof[1]];
      yv2 = cv->y[dof[2]];
      yv3 = cv->y[dof[3]];
      zv0 = cv->z[dof[0]];
      zv1 = cv->z[dof[1]];
      zv2 = cv->z[dof[2]];
      zv3 = cv->z[dof[3]];

      // Get coordinates on reference triangle
      REAL det = (xv3-xv0)*((yv1-yv0)*(zv2-zv0)-(yv2-yv0)*(zv1-zv0)) \
      - (xv2-xv0)*((yv1-yv0)*(zv3-zv0)-(yv3-yv0)*(zv1-zv0)) \
      + (xv1-xv0)*((yv2-yv0)*(zv3-zv0)-(yv3-yv0)*(zv2-zv0));

      REAL r = ((xv3-xv0)*((x[1]-yv0)*(zv2-zv0)-(yv2-yv0)*(x[2]-zv0)) \
      - (xv2-xv0)*((x[1]-yv0)*(zv3-zv0)-(yv3-yv0)*(x[2]-zv0)) \
      + (x[0]-xv0)*((yv2-yv0)*(zv3-zv0)-(yv3-yv0)*(zv2-zv0)))/det;

      REAL drdx = ((yv2-yv0)*(zv3-zv0)-(yv3-yv0)*(zv2-zv0))/det;
      REAL drdy = ((xv3-xv0)*(zv2-zv0) - (xv2-xv0)*(zv3-zv0))/det;
      REAL drdz = ((xv2-xv0)*(yv3-yv0) - (xv3-xv0)*(yv2-yv0))/det;

      REAL s = ((xv3-xv0)*((yv1-yv0)*(x[2]-zv0)-(x[1]-yv0)*(zv1-zv0)) \
      - (x[0]-xv0)*((yv1-yv0)*(zv3-zv0)-(yv3-yv0)*(zv1-zv0)) \
      + (xv1-xv0)*((x[1]-yv0)*(zv3-zv0)-(yv3-yv0)*(x[2]-zv0)))/det;

      REAL dsdx = -((yv1-yv0)*(zv3-zv0)-(yv3-yv0)*(zv1-zv0))/det;
      REAL dsdy = ((xv1-xv0)*(zv3-zv0) - (xv3-xv0)*(zv1-zv0))/det;
      REAL dsdz = ((xv3-xv0)*(yv1-yv0) - (xv1-xv0)*(yv3-yv0))/det;

      REAL t = ((x[0]-xv0)*((yv1-yv0)*(zv2-zv0)-(yv2-yv0)*(zv1-zv0)) \
      - (xv2-xv0)*((yv1-yv0)*(x[2]-zv0)-(x[1]-yv0)*(zv1-zv0)) \
      + (xv1-xv0)*((yv2-yv0)*(x[2]-zv0)-(x[1]-yv0)*(zv2-zv0)))/det;

      REAL dtdx = ((yv1-yv0)*(zv2-zv0)-(yv2-yv0)*(zv1-zv0))/det;
      REAL dtdy = ((xv2-xv0)*(zv1-zv0) - (xv1-xv0)*(zv2-zv0))/det;
      REAL dtdz = ((xv1-xv0)*(yv2-yv0) - (xv2-xv0)*(yv1-yv0))/det;

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
        dp0r = -1;
        dp1r = 1;
        dp2r = 0;
        dp3r = 0;
        dp0s = -1;
        dp1s = 0;
        dp2s = 1;
        dp3s = 0;
        dp0t = -1;
        dp1t = 0;
        dp2t = 0;
        dp3t = 1;
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

        dp0r = 4*r+4*s+4*t-3;
        dp1r = 4*r-1;
        dp2r = 0;
        dp3r = 0;
        dp4r = 4*onemrst - 4*r;
        dp5r = -4*s;
        dp6r = -4*t;
        dp7r = 4*s;
        dp8r = 4*t;
        dp9r = 0;
        dp0s = 4*r+4*s+4*t-3;
        dp1s = 0;
        dp2s = 4*s-1;
        dp3s = 0;
        dp4s = -4*r;
        dp5s = 4*onemrst - 4*s;
        dp6s = -4*t;
        dp7s = 4*r;
        dp8s = 0;
        dp9s = 4*t;
        dp0t = 4*r+4*s+4*t-3;
        dp1t = 0;
        dp2t = 0;
        dp3t = 4*t-1;
        dp4t = -4*r;
        dp5t = -4*s;
        dp6t = 4*onemrst - 4*t;
        dp7t = 0;
        dp8t = 4*r;
        dp9t = 4*s;
      } else {
        status = ERROR_FE_TYPE;
        check_error(status, __FUNCTION__);
      }

      /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
      *  to (X,Y) using the chain rule.
      */
      dp[0*dim] = dp0r * drdx + dp0s * dsdx + dp0t * dtdx;
      dp[0*dim+1] = dp0r * drdy + dp0s * dsdy + dp0t * dtdy;
      dp[0*dim+2] = dp0r * drdz + dp0s * dsdz + dp0t * dtdz;
      dp[1*dim] = dp1r * drdx + dp1s * dsdx + dp1t * dtdx;
      dp[1*dim+1] = dp1r * drdy + dp1s * dsdy + dp1t * dtdy;
      dp[1*dim+2] = dp1r * drdz + dp1s * dsdz + dp1t * dtdz;
      dp[2*dim] = dp2r * drdx + dp2s * dsdx + dp2t * dtdx;
      dp[2*dim+1] = dp2r * drdy + dp2s * dsdy + dp2t * dtdy;
      dp[2*dim+2] = dp2r * drdz + dp2s * dsdz + dp2t * dtdz;
      dp[3*dim] = dp3r * drdx + dp3s * dsdx + dp3t * dtdx;
      dp[3*dim+1] = dp3r * drdy + dp3s * dsdy + dp3t * dtdy;
      dp[3*dim+2] = dp3r * drdz + dp3s * dsdz + dp3t * dtdz;
      if(porder==2) {
        dp[4*dim] = dp4r * drdx + dp4s * dsdx + dp4t * dtdx;
        dp[4*dim+1] = dp4r * drdy + dp4s * dsdy + dp4t * dtdy;
        dp[4*dim+2] = dp4r * drdz + dp4s * dsdz + dp4t * dtdz;
        dp[5*dim] = dp5r * drdx + dp5s * dsdx + dp5t * dtdx;
        dp[5*dim+1] = dp5r * drdy + dp5s * dsdy + dp5t * dtdy;
        dp[5*dim+2] = dp5r * drdz + dp5s * dsdz + dp5t * dtdz;
        dp[6*dim] = dp6r * drdx + dp6s * dsdx + dp6t * dtdx;
        dp[6*dim+1] = dp6r * drdy + dp6s * dsdy + dp6t * dtdy;
        dp[6*dim+2] = dp6r * drdz + dp6s * dsdz + dp6t * dtdz;
        dp[7*dim] = dp7r * drdx + dp7s * dsdx + dp7t * dtdx;
        dp[7*dim+1] = dp7r * drdy + dp7s * dsdy + dp7t * dtdy;
        dp[7*dim+2] = dp7r * drdz + dp7s * dsdz + dp7t * dtdz;
        dp[8*dim] = dp8r * drdx + dp8s * dsdx + dp8t * dtdx;
        dp[8*dim+1] = dp8r * drdy + dp8s * dsdy + dp8t * dtdy;
        dp[8*dim+2] = dp8r * drdz + dp8s * dsdz + dp8t * dtdz;
        dp[9*dim] = dp9r * drdx + dp9s * dsdx + dp9t * dtdx;
        dp[9*dim+1] = dp9r * drdy + dp9s * dsdy + dp9t * dtdy;
        dp[9*dim+2] = dp9r * drdz + dp9s * dsdz + dp9t * dtdz;
      }
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }
  }

  return;
}
/*******************************************************************************************************/

/*******************************************************************************************************/
/*!
* \fn void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,mesh_struct *mesh)
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
*/
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,mesh_struct *mesh)
{
  REAL dp0r,dp1r,dp2r,dp3r,dp4r,dp5r;
  REAL dp0s,dp1s,dp2s,dp3s,dp4s,dp5s;
  REAL onemrs;
  REAL xv0,xv1,xv2,yv0,yv1,yv2;
  REAL dp0rr,dp1rr,dp2rr,dp3rr,dp4rr,dp5rr;
  REAL dp0ss,dp1ss,dp2ss,dp3ss,dp4ss,dp5ss;
  REAL dp0rs,dp1rs,dp2rs,dp3rs,dp4rs,dp5rs;

  coordinates* cv = mesh->cv;

  // Get Nodes and Physical Coordinates of vertices only
  xv0 = cv->x[dof[0]];
  xv1 = cv->x[dof[1]];
  xv2 = cv->x[dof[2]];
  yv0 = cv->y[dof[0]];
  yv1 = cv->y[dof[1]];
  yv2 = cv->y[dof[2]];

  // Get coordinates on reference triangle
  REAL det = (xv1-xv0)*(yv2-yv0) - (xv2-xv0)*(yv1-yv0);

  REAL r = ((yv2-yv0)*(x-xv0) + (xv0-xv2)*(y-yv0))/det;

  REAL drdx = (yv2-yv0)/det;
  REAL drdy = (xv0-xv2)/det;

  REAL s = ((yv0-yv1)*(x-xv0) + (xv1-xv0)*(y-yv0))/det;

  REAL dsdx = (yv0-yv1)/det;
  REAL dsdy = (xv1-xv0)/det;

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
  dp0r = 4*r+4*s-3;
  dp0rr = 4.0;
  dp0rs = 4.0;
  dp1r = 4*r - 1;
  dp1rr = 4.0;
  dp1rs = 0.0;
  dp2r = 0.0;
  dp2rr = 0.0;
  dp2rs = 0.0;
  dp3r = 4.0-8*r-4*s;
  dp3rr = -8.0;
  dp3rs = -4.0;
  dp4r = -4.0*s;
  dp4rr = 0.0;
  dp4rs = -4.0;
  dp5r = 4.0*s;
  dp5rr = 0.0;
  dp5rs = 4.0;
  dp0s = dp0r;
  dp0ss = 4.0;
  dp1s = 0.0;
  dp1ss = 0.0;
  dp2s = 4.0*s-1.0;
  dp2ss = 4.0;
  dp3s = -4.0*r;
  dp3ss = 0;
  dp4s = 4.0-4*r-8*s;
  dp4ss = -8.0;
  dp5s = 4.0*r;
  dp5ss = 0.0;


  /*  We need to convert the derivative information from (R(X,Y),S(X,Y))
  *  to (X,Y) using the chain rule.
  */

  dpx[0] = dp0r * drdx + dp0s * dsdx;
  dpy[0] = dp0r * drdy + dp0s * dsdy;
  dpx[1] = dp1r * drdx + dp1s * dsdx;
  dpy[1] = dp1r * drdy + dp1s * dsdy;
  dpx[2] = dp2r * drdx + dp2s * dsdx;
  dpy[2] = dp2r * drdy + dp2s * dsdy;
  dpx[3] = dp3r * drdx + dp3s * dsdx;
  dpy[3] = dp3r * drdy + dp3s * dsdy;
  dpx[4] = dp4r * drdx + dp4s * dsdx;
  dpy[4] = dp4r * drdy + dp4s * dsdy;
  dpx[5] = dp5r * drdx + dp5s * dsdx;
  dpy[5] = dp5r * drdy + dp5s * dsdy;
  dpxx[0] = dp0rr*drdx*drdx + 2*dp0rs*drdx*dsdx + dp0ss*dsdx*dsdx;
  dpxx[1] = dp1rr*drdx*drdx + 2*dp1rs*drdx*dsdx + dp1ss*dsdx*dsdx;
  dpxx[2] = dp2rr*drdx*drdx + 2*dp2rs*drdx*dsdx + dp2ss*dsdx*dsdx;
  dpxx[3] = dp3rr*drdx*drdx + 2*dp3rs*drdx*dsdx + dp3ss*dsdx*dsdx;
  dpxx[4] = dp4rr*drdx*drdx + 2*dp4rs*drdx*dsdx + dp4ss*dsdx*dsdx;
  dpxx[5] = dp5rr*drdx*drdx + 2*dp5rs*drdx*dsdx + dp5ss*dsdx*dsdx;

  dpyy[0] = dp0rr*drdy*drdy + 2*dp0rs*drdy*dsdy + dp0ss*dsdy*dsdy;
  dpyy[1] = dp1rr*drdy*drdy + 2*dp1rs*drdy*dsdy + dp1ss*dsdy*dsdy;
  dpyy[2] = dp2rr*drdy*drdy + 2*dp2rs*drdy*dsdy + dp2ss*dsdy*dsdy;
  dpyy[3] = dp3rr*drdy*drdy + 2*dp3rs*drdy*dsdy + dp3ss*dsdy*dsdy;
  dpyy[4] = dp4rr*drdy*drdy + 2*dp4rs*drdy*dsdy + dp4ss*dsdy*dsdy;
  dpyy[5] = dp5rr*drdy*drdy + 2*dp5rs*drdy*dsdy + dp5ss*dsdy*dsdy;

  dpxy[0] = dp0rr*drdy*drdx + dp0rs*drdy*dsdx + dp0rs*drdx*dsdy + dp0ss*dsdx*dsdy;
  dpxy[1] = dp1rr*drdy*drdx + dp1rs*drdy*dsdx + dp1rs*drdx*dsdy + dp1ss*dsdx*dsdy;
  dpxy[2] = dp2rr*drdy*drdx + dp2rs*drdy*dsdx + dp2rs*drdx*dsdy + dp2ss*dsdx*dsdy;
  dpxy[3] = dp3rr*drdy*drdx + dp3rs*drdy*dsdx + dp3rs*drdx*dsdy + dp3ss*dsdx*dsdy;
  dpxy[4] = dp4rr*drdy*drdx + dp4rs*drdy*dsdx + dp4rs*drdx*dsdy + dp4ss*dsdx*dsdy;
  dpxy[5] = dp5rr*drdy*drdx + dp5rs*drdy*dsdx + dp5rs*drdx*dsdy + dp5ss*dsdx*dsdy;

  return;
}
/*******************************************************************************************************/

/*******************************************************************************************************/
/*!
* \fn void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
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
void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
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

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in mesh to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = mesh->dwork;
  REAL* dp = mesh->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  REAL elen;

  /* Now, with the linear basis functions p, dpx, and dpy, the nedelec elements are
  * phi_eij = |eij|*(p(i)grad(p(j)) - p(j)grad(p(i)))
  * |eij| = sqrt(|xj-xi|^2)
  */

  // Go through each edge and get length and find the corresponding nodes
  for (i=0; i<ed_per_elm; i++) {
    ica = mesh->ed_v->IA[dof[i]];
    n1 = mesh->ed_v->JA[ica];
    n2 = mesh->ed_v->JA[ica+1];
    elen = mesh->ed_len[dof[i]];

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

  // Clean up working array for next person to use.
  array_set(v_per_elm*(dim+1),mesh->dwork,0.0);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
* \fn void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
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
void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
{
  // Flag for erros
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;

  INT i,j,ica,icb,jcnt;
  INT ipf[dim];
  INT myf;
  REAL farea;
  INT elnd,ef1,ef2,ef3;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in mesh to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = mesh->dwork;
  REAL* dp = mesh->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = mesh->f_v->IA[myf];
      icb = mesh->f_v->IA[myf+1];
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area and normal vector of the face
      farea = mesh->f_area[myf];

      // Loop through Nodes on element to find corresponding nodes and get correct orientation
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
      ica = mesh->f_v->IA[myf];
      icb = mesh->f_v->IA[myf+1];
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area
      farea = mesh->f_area[myf];

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
      * phi_fijk = 2*|fijk|*(p(i)(grad(p(j)) x grad(p(k))) + p(j)(grad(p(k)) x grad(p(i))) + p(k)(grad(p(i)) x grad(p(j))))
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

      // Compute divs div(phi_fijk) = 6*|fijk| grad(pi) dot (grad(pj) cross grad(pk)) (
      dphi[i] = 6*farea*(dp[ef1*dim]*(dp[ef2*dim+1]*dp[ef3*dim+2]-dp[ef2*dim+2]*dp[ef3*dim+1])
      + dp[ef1*dim+1]*(dp[ef2*dim+2]*dp[ef3*dim]-dp[ef2*dim]*dp[ef3*dim+2])
      + dp[ef1*dim+2]*(dp[ef2*dim]*dp[ef3*dim+1]-dp[ef2*dim+1]*dp[ef3*dim]));
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Clean up working array for next person to use.
  array_set(v_per_elm*(dim+1),mesh->dwork,0.0);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
* \fn void bdm1_basis_global(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
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
void bdm1_basis_global(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT f_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;

  INT i,j,ica,icb,jcnt;
  REAL a1,a2,a3,a4;
  INT ipf[dim];
  INT myf;
  REAL farea;
  INT elnd,ef1,ef2;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in mesh to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = mesh->dwork;
  REAL* dp = mesh->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = mesh->f_v->IA[myf];
      icb = mesh->f_v->IA[myf+1];
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = mesh->f_v->JA[j];
        jcnt++;
      }

      // Get the area and normal vector of the face
      farea = mesh->f_area[myf];

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

  // Clean up working array for next person to use.
  array_set(v_per_elm*(dim+1),mesh->dwork,0.0);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
* \fn void bubble_face_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh)
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
void bubble_face_basis(REAL *phi, REAL *dphi, REAL *x, INT *v_on_elm, INT *dof, mesh_struct *mesh)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  INT v_per_elm = mesh->v_per_elm;
  INT dof_per_elm = mesh->f_per_elm;
  INT dim = mesh->dim;
  INT i,j;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in mesh to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = mesh->dwork;
  REAL* dp = mesh->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,mesh);

  // face to vertex map
  INT fv[dim];// = (INT *)calloc(dim,sizeof(INT));
  INT elnd,ef1,ef2,ef3;//face endpoint vertex tracking numbers

  REAL gradp;

  if(dim==2){
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i],mesh->f_v,fv);
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
      /* phi[i*dim] =   ABS(mesh->f_norm[dim*(dof[i])]  )* 4*p[ef1]*p[ef2]; */
      /* phi[i*dim+1] = ABS(mesh->f_norm[dim*(dof[i])+1]) * 4*p[ef1]*p[ef2]; */
      phi[i*dim] =   mesh->f_norm[dim*(dof[i])]* 4*p[ef1]*p[ef2];
      phi[i*dim+1] = mesh->f_norm[dim*(dof[i])+1]* 4*p[ef1]*p[ef2];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 4*(p[ef1]*dp[ef2*dim+j] + dp[ef1*dim+j]*p[ef2]);

        /* dphi[i*dim*dim + j*dim + 0] = gradp * ABS(mesh->f_norm[dim*(dof[i])+0]); */
        /* dphi[i*dim*dim + j*dim + 1] = gradp * ABS(mesh->f_norm[dim*(dof[i])+1]); */
        dphi[i*dim*dim + j*dim + 0] = gradp * mesh->f_norm[dim*(dof[i])+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * mesh->f_norm[dim*(dof[i])+1];
      }

    }
  } else if(dim==3) {
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i],mesh->f_v,fv);
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
      phi[i*dim] = mesh->f_norm[dim*(dof[i])] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+1] = mesh->f_norm[dim*(dof[i])+1] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+2] = mesh->f_norm[dim*(dof[i])+2] * 8*p[ef1]*p[ef2]*p[ef3];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 8*(p[ef1]*p[ef2]*dp[ef3*dim+j] + p[ef1]*dp[ef2*dim+j]*p[ef3] + dp[ef1*dim+j]*p[ef2]*p[ef3]);

        dphi[i*dim*dim + j*dim + 0] = gradp * mesh->f_norm[dim*(dof[i])+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * mesh->f_norm[dim*(dof[i])+1];
        dphi[i*dim*dim + j*dim + 2] = gradp * mesh->f_norm[dim*(dof[i])+2];
      }

    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Clean up working array for next person to use.
  array_set(v_per_elm*(dim+1),mesh->dwork,0.0);

  return;
}
/****************************************************************************************************************************/

/****************************************************************************************************************************/
/*!
* \fn void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh,fespace *FE)
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
void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,mesh_struct *mesh,fespace *FE)
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
    if (mesh->dim>1){
      PX_H1_basis(phi+offset, dphi+FE->dof_per_elm, x, dof, 1, mesh);
    }
    if (mesh->dim>2){
      PX_H1_basis(phi+offset*2, dphi+FE->dof_per_elm*2, x, dof, 1, mesh);
    }

  } else if(FEtype==61) { // Bubble element

    bubble_face_basis(phi,dphi,x,v_on_elm,dof,mesh);

  } else if(FEtype==99) { // Constraint Single DOF Space

    phi[0] = 1.0;

  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
/****************************************************************************************************************************/
