/*! \file src/fem/zcommented_out_basis.c
 *
 * \brief Commented-out basis and assembly functions preserved for reference.
 *        These are not compiled (wrapped in #if 0).
 *
 */

#include "hazmath.h"

#if 0 /* quad_tri_2D_2der — dead code, never called */
/*******************************************************************************************************/
/*!
* \fn void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,scomplex *sc)
*
* \brief Compute Standard Quadratic Finite Element Basis Functions (P2) at a particular point
*        Also compute the 2nd derivatives for some reason...
*
* \param x,y,z           Coordinate on where to compute basis function
* \param dof             DOF for the given element (in this case vertices and their global numbering)
* \param porder          Order of elements
* \param sc              Simplicial complex
*
* \return p              Basis functions (1 for each DOF on element)
* \return dpx,dpy        Derivatives of basis functions (i.e., gradient)
* \return dpxx,dpyy,dpxy 2nd Derivatives of basis functions
*
*/
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,scomplex *sc)
{
  REAL dp0r,dp1r,dp2r,dp3r,dp4r,dp5r;
  REAL dp0s,dp1s,dp2s,dp3s,dp4s,dp5s;
  REAL onemrs;
  REAL xv0,xv1,xv2,yv0,yv1,yv2;
  REAL dp0rr,dp1rr,dp2rr,dp3rr,dp4rr,dp5rr;
  REAL dp0ss,dp1ss,dp2ss,dp3ss,dp4ss,dp5ss;
  REAL dp0rs,dp1rs,dp2rs,dp3rs,dp4rs,dp5rs;

  INT dim = sc->dim;

  // Get Nodes and Physical Coordinates of vertices only
  xv0 = sc->x[dof[0]*dim];
  xv1 = sc->x[dof[1]*dim];
  xv2 = sc->x[dof[2]*dim];
  yv0 = sc->x[dof[0]*dim+1];
  yv1 = sc->x[dof[1]*dim+1];
  yv2 = sc->x[dof[2]*dim+1];

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
#endif /* quad_tri_2D_2der */

#if 0 /* bdm1_basis_global — dead code, never called */
/****************************************************************************************************************************/
/*!
* \fn void bdm1_basis_global(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
*
* \brief Brezzi-Douglas-Marini (BDM) Elements of order 1.
*
* \note ONLY in 2D for now.
* \note This has NOT been tested.
*
* \param x         Coordinate on where to compute basis function
* \param v_on_elm  Vertices on element
* \param dof       DOF on element
* \param sc        Simplicial complex
*
* \return phi      Basis functions (2*f_per_elm for each face, 12 total in 2D)
* \return dphix    Div of basis functions
*
*/
void bdm1_basis_global(REAL *phi,REAL *dphix,REAL *dphiy,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT v_per_elm = (dim + 1);
  INT f_per_elm = (dim + 1);

  INT i,j,ica,icb,jcnt;
  REAL a1,a2,a3,a4;
  INT ipf[dim];
  INT myf;
  REAL farea;
  INT elnd,ef1,ef2;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in fem to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = fem->dwork;
  REAL* dp = fem->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,sc);

  // Go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      myf = dof[i];
      ica = fem->f_v->IA[myf];
      icb = fem->f_v->IA[myf+1];
      jcnt=0;
      for(j=ica;j<icb;j++) {
        ipf[jcnt] = fem->f_v->JA[j];
        jcnt++;
      }

      // Get the area and normal vector of the face
      farea = fem->f_area[myf];

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
  array_set(v_per_elm*(dim+1),fem->dwork,0.0);

  return;
}
#endif /* bdm1_basis_global */

#if 0 /* FEM_RHS_Local_face — dead code, never called */
/******************************************************************************************************/
/*!
* \fn void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,scomplex *sc,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the local assembly of a RHS for any "boundary" bilinear form using various element types
*        (eg. P1, P2, Nedelec, and Raviart-Thomas).
*        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
*        a(u,v)_i, where i denotes a set of faces (or edges) with in a boundary region marked with flag
*
*
*        For this problem we compute local RHS of:
*
*        Lu = f  ---->   a(u,v)_bdry = <f,v>_bdry
*
*        which gives Ax = b,
*
*        A_ij = a( phi_j, phi_i)_bdry
*
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.

* \note Assumes different type of integral for different Element type:
*       Scalar -> <f,v>_bdry
*       Vector -> <f,n*v>_bdry
*
* \note This reallocates the quadrature on each face for every call to this function.
*       This is not optimal, and the quadrature should be set outside if you are
*       building your own local assembly.
*
* \param old_sol                 FE approximation of previous solution if needed
* \param FE                      FE Space
* \param mesh                    Mesh Data
* \param dof_on_f                DOF on the given face
* \param dof_on_elm              DOF on given element
* \param v_on_elm                Vertices on given element
* \param dof_per_f               # of DOF per face
* \param face                    Given Face
* \param elm                     Given Element
* \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
* \param time                    Physical Time if time dependent
*
* \return bLoc                   Local RHS vector
*
*/
void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,scomplex *sc,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Mesh and FE data
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;

  // Loop Indices
  INT j,quad,test,doft;

  // Quadrature Weights and Nodes
  qcoordinates *cq_face = allocateqcoords_bdry(cq->nq1d,1,dim,2);
  quad_face(cq_face,sc,cq->nq1d,face);
  REAL qx[dim+1];
  REAL w;

  // Get normal vector components on face if needed
  REAL nx=0.0,ny=0.0,nz=0.0;
  if(FE->FEtype>=20) {
    nx = fem->f_norm[face*dim];
    ny = fem->f_norm[face*dim+1];
    if(dim==3) nz = fem->f_norm[face*dim+2];
  }

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Value of f
  REAL rhs_val;

  if(FE->scal_or_vec==0) { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq_face->n;quad++) {
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time,&(fem->f_flag[face]));

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,sc,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_face;test++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<FE->dof_per_elm;j++) {
          if(dof_on_f[test]==dof_on_elm[j]) {
            doft = j;
          }
        }
        kij = rhs_val*FE->phi[doft];
        bLoc[test] += w*kij;
      }
    }
  } else { // Vector Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq_face->n;quad++) {
    //for (quad=0;quad<1;quad++) {
    //  qx[0] = fem->f_mid[face*dim];
    //  qx[1] = fem->f_mid[face*dim+1];
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
//        qx[2] = fem->f_mid[face*dim+2];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time,&(fem->f_flag[face]));

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,sc,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_face;test++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<FE->dof_per_elm;j++) {
          if(dof_on_f[test]==dof_on_elm[j]) {
            doft = j;
          }
        }
        kij = rhs_val*(nx*FE->phi[doft*dim] + ny*FE->phi[doft*dim+1]);
        if(dim==3) kij +=rhs_val*nz*FE->phi[doft*dim+2];
        bLoc[test] += w*kij;
      }
    }
  }

  return;
}
#endif /* FEM_RHS_Local_face */

#if 0 /* ned_basis — replaced by ned0_basis called via get_FEM_basis */
/*******************************************************************************************************/
/*!
* \fn void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
*
* \brief Compute Nedelec Finite Element Basis Functions (zeroth order) at a particular point in 2 or 3D
*        Old scomplex-based version. Now replaced by ned0_basis (local data version).
*/
void ned_basis(REAL *phi,REAL *cphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
{
  SHORT status;
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT v_per_elm = (dim + 1);
  INT ed_per_elm = dim*(dim+1)/2;
  INT i,k,ica,n1,n2,ihi,ilo;
  INT mark1 = -1, mark2 = -1;
  REAL* p = fem->dwork;
  REAL* dp = fem->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,sc);
  REAL elen;
  for (i=0; i<ed_per_elm; i++) {
    ica = fem->ed_v->IA[dof[i]];
    n1 = fem->ed_v->JA[ica];
    n2 = fem->ed_v->JA[ica+1];
    elen = fem->ed_len[dof[i]];
    for (k=0; k<v_per_elm; k++) {
      if (v_on_elm[k]==n1) mark1=k;
      if (v_on_elm[k]==n2) mark2=k;
    }
    if (MAX(n1,n2)==n1) { ihi = mark1; ilo = mark2; }
    else { ihi = mark2; ilo = mark1; }
    for(k=0;k<dim;k++)
      phi[i*dim+k] = elen*(p[ilo]*dp[ihi*dim+k] - p[ihi]*dp[ilo*dim+k]);
    if(dim==2) {
      cphi[i] = 2*elen*(dp[ilo*dim]*dp[ihi*dim+1] - dp[ilo*dim+1]*dp[ihi*dim]);
    } else if(dim==3) {
      cphi[i*dim+0] = 2*elen*(dp[ilo*dim+1]*dp[ihi*dim+2]-dp[ihi*dim+1]*dp[ilo*dim+2]);
      cphi[i*dim+1] = 2*elen*(dp[ihi*dim]*dp[ilo*dim+2]-dp[ilo*dim]*dp[ihi*dim+2]);
      cphi[i*dim+2] = 2*elen*(dp[ilo*dim]*dp[ihi*dim+1]-dp[ihi*dim]*dp[ilo*dim+1]);
    } else { status = ERROR_DIM; check_error(status, __FUNCTION__); }
  }
  array_set(v_per_elm*(dim+1),fem->dwork,0.0);
  return;
}
#endif /* ned_basis */

#if 0 /* rt_basis — replaced by rt0_basis called via get_FEM_basis */
/****************************************************************************************************************************/
/*!
* \fn void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
*
* \brief Compute Raviart-Thomas Finite Element Basis Functions (zeroth order).
*        Old scomplex-based version. Now replaced by rt0_basis (local data version).
*/
void rt_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT v_per_elm = (dim + 1);
  INT f_per_elm = (dim + 1);
  INT v_per_f = dim;
  INT i,j,k,c,ica,icb,jcnt;
  INT ipf[dim], ef[dim];
  INT myf;
  REAL farea;
  REAL* p = fem->dwork;
  REAL* dp = fem->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,sc);
  REAL fac = 1.0;
  for (i=2; i<dim; i++) fac *= i;
  for (i=0; i<f_per_elm; i++) {
    myf = dof[i];
    ica = fem->f_v->IA[myf];
    icb = fem->f_v->IA[myf+1];
    jcnt=0;
    for(j=ica;j<icb;j++) { ipf[jcnt] = fem->f_v->JA[j]; jcnt++; }
    farea = fem->f_area[myf];
    for (k=0; k<v_per_f; k++)
      for (j=0; j<v_per_elm; j++)
        if (ipf[k] == v_on_elm[j]) { ef[k] = j; break; }
    INT opp = -1;
    for (j=0; j<v_per_elm; j++) {
      INT on_face = 0;
      for (k=0; k<v_per_f; k++) if (ef[k] == j) { on_face = 1; break; }
      if (!on_face) { opp = j; break; }
    }
    REAL M[dim*dim];
    for (c=0; c<dim; c++)
      for (k=0; k<dim; k++)
        M[c*dim+k] = dp[ef[k]*dim+c];
    REAL sigma = fac * haz_det(dim, M);
    INT vopp = v_on_elm[opp];
    for (c=0; c<dim; c++)
      phi[i*dim+c] = sigma * farea * (x[c] - sc->x[vopp*dim+c]);
    dphi[i] = sigma * dim * farea;
  }
  array_set(v_per_elm*(dim+1),fem->dwork,0.0);
  return;
}
#endif /* rt_basis */

#if 0 /* FE_Interpolation, FE_DerivativeInterpolation, P2_2ndDerivativeInterpolation
       — replaced by fe_interpolation_to_x, fe_dinterp_to_x, fe_2ndderiv_interp_to_x */
/*!
 * \fn void FE_Interpolation(REAL* val,REAL *u,REAL* x,INT *dof_on_elm,INT
 * *v_on_elm,fespace *FE,scomplex *sc)
 *
 * \brief Interpolate a finite-element approximation to any other point in the
 * given element using the given type of elements.
 *
 * \param u 	      Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param val         Pointer to value of approximation at given values
 *
 * \note This will be deprecated soon.
 *
 */
void FE_Interpolation(REAL* val, REAL* u, REAL* x, INT* dof_on_elm,
                      INT* v_on_elm, fespace* FE, scomplex* sc) {
  sc_fem* fem = sc->fem;
  INT dim = sc->dim;
  INT i, j, dof;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  // INT FEtype = FE->FEtype;
  INT scal_or_vec = FE->scal_or_vec;

  REAL coef[dim];

  get_FEM_basis(FE->phi, FE->dphi, x, v_on_elm, dof_on_elm, sc, FE);

  if (scal_or_vec == 0) {  // Scalar Element
    coef[0] = 0.0;
    for (j = 0; j < dof_per_elm; j++) {
      dof = dof_on_elm[j];
      coef[0] += u[dof] * FE->phi[j];
    }
    val[0] = coef[0];
  } else {  // Vector Element
    for (i = 0; i < dim; i++) {
      coef[i] = 0.0;
      for (j = 0; j < dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[i] += u[dof] * FE->phi[j * dim + i];
      }
      val[i] = coef[i];
    }
  }

  return;
}
/******************************************************************************/

/*!
 * \fn void FE_DerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT
 * *dof_on_elm,INT *v_on_elm,fespace *FE,scomplex *sc)
 *
 * \brief Interpolate the "derivative" of a finite-element approximation to any
 * other point in the given element using the given type of elements. Note that
 * for Lagrange Elements this means the Gradient, grad u, for Nedelec it means
 * the Curl, curl u, and for RT it is the Divergence, div u.
 *
 * \param u 	      Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param val         Pointer to value of approximation at given values
 *
 * \note This will be deprecated soon.
 *
 */
void FE_DerivativeInterpolation(REAL* val, REAL* u, REAL* x, INT* dof_on_elm,
                                INT* v_on_elm, fespace* FE, scomplex* sc) {
  sc_fem* fem = sc->fem;
  INT dim = sc->dim;
  INT dof, j, k, i;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT FEtype = FE->FEtype;

  // Basis Functions and its derivatives if necessary
  // REAL coef[dim];
  REAL coef[dim * dim];

  get_FEM_basis(FE->phi, FE->dphi, x, v_on_elm, dof_on_elm, sc, FE);

  if (FEtype == 0 ||
      FEtype == 99) {  // Don't compute derivatives of P0 elements (set to 0)
    for (j = 0; j < dim; j++) {
      val[j] = 0.0;
    }
  } else if ((FEtype < 20 && FEtype != 0) || FEtype == 103) {  // Scalar Element
    for (j = 0; j < dim; j++) {
      coef[j] = 0.0;
      for (k = 0; k < dof_per_elm; k++) {
        dof = dof_on_elm[k];
        coef[j] += u[dof] * FE->dphi[k * dim + j];
      }
      val[j] = coef[j];
    }
  } else if (FEtype == 20) {  // Nedelec
    for (j = 0; j < dim; j++) coef[j] = 0.0;
    if (dim == 2) {  // Curl is scalar
      for (j = 0; j < dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[0] += u[dof] * FE->dphi[j];
      }
      val[0] = coef[0];
    } else if (dim == 3) {  // Curl is vector
      for (j = 0; j < dof_per_elm; j++) {
        dof = dof_on_elm[j];
        coef[0] += u[dof] * FE->dphi[j * dim + 0];
        coef[1] += u[dof] * FE->dphi[j * dim + 1];
        coef[2] += u[dof] * FE->dphi[j * dim + 2];
      }
      val[0] = coef[0];
      val[1] = coef[1];
      val[2] = coef[2];
    }
  } else if (FEtype == 30) {  // Raviart-Thomas (div is scalar)
    coef[0] = 0.0;

    for (j = 0; j < dof_per_elm; j++) {
      dof = dof_on_elm[j];
      coef[0] += u[dof] * FE->dphi[j];
    }
    val[0] = coef[0];
  } else if (FEtype >= 60) {  // Vector Lagrange Functions (i.e. bubbles) ->
                              // gradients of vectors -> tensor
    for (j = 0; j < dim; j++) {
      for (i = 0; i < dim; i++) {
        coef[j * dim + i] = 0.0;
        for (k = 0; k < dof_per_elm; k++) {
          dof = dof_on_elm[k];
          coef[j * dim + i] += u[dof] * FE->dphi[k * dim * dim + j * dim + i];
        }
        val[j * dim + i] = coef[j * dim + i];
      }
    }

  } else {
    check_error(ERROR_FE_TYPE, __FUNCTION__);
  }

  return;
}
/******************************************************************************/

/* blockFE_Interpolation removed — use blockfe_interpolation_to_x/quadpt instead */
/* blockFE_DerivativeInterpolation removed — use blockfe_dinterp_to_x/quadpt instead */
/******************************************************************************/

/******************************************************************************/
/* Some special interpolation routines needed for some error estimators */
/******************************************************************************/
/*!
 * \fn void P2_2ndDerivativeInterpolation(REAL* val,REAL *u,REAL* x,INT
 * *dof_on_elm,INT *v_on_elm,fespace *FE,scomplex *sc)
 *
 * \brief Interpolate the "second derivative" of a P2 finite-element
 * approximation to any other point in the given element
 *
 * \param u 	        Approximation to interpolate
 * \param x           Coordinates where to compute value
 * \param dof_on_elm  DOF belonging to particular element
 * \param v_on_elm    Vertices belonging to particular element
 * \param FE          FE Space
 * \param mesh        Mesh Data
 * \param val         Pointer 2nd derivatives (uxx, uxy, uyy, uxz, uyz, uzz) of
 * approximation at given values
 *
 */
void P2_2ndDerivativeInterpolation(REAL* val, REAL* u, REAL* x, INT* dof_on_elm,
                                   fespace* FE, scomplex* sc) {
  sc_fem* fem = sc->fem;
  INT dim = sc->dim;
  INT dof, j, k;

  // Get FE and Mesh data
  INT dof_per_elm = FE->dof_per_elm;
  INT twoders = 3 * (dim - 1);

  REAL* ddp = (REAL*)calloc(3 * (dim - 1) * dof_per_elm, sizeof(REAL));

  // Basis Functions and its derivatives
  P2_basis_2der(FE->phi, FE->dphi, ddp, x, dof_on_elm, sc);

  for (j = 0; j < twoders; j++) {
    val[j] = 0.0;
    for (k = 0; k < dof_per_elm; k++) {
      dof = dof_on_elm[k];
      val[j] += u[dof] * ddp[k * twoders + j];
    }
  }

  if (ddp) free(ddp);

  return;
}
#endif /* FE_Interpolation etc */

#if 0 /* bubble_face_basis, MINI_basis, get_FEM_basis — scomplex-based, replaced by local-data interface */
/* bubble_face_basis, MINI_basis: scomplex-based versions still needed by get_FEM_basis
 * which is called from FE_Interpolation/FE_DerivativeInterpolation in interp.c */
/****************************************************************************************************************************/
/*!
* \fn void bubble_face_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc)
*
* \brief Compute Bubble Element Finite Element Basis Functions at a particular point in 2 or 3D
*
* \param x         Coordinate on where to compute basis function
* \param v_on_elm  Vertices on element
* \param dof       DOF on element
* \param sc        Simplicial complex
*
* \return phi      Basis functions (dim for each face from reference triangle)
* \return dphi     Tensor from gradient of basis functions
*/
void bubble_face_basis(REAL *phi, REAL *dphi, REAL *x, INT *v_on_elm, INT *dof, scomplex *sc)
{
  // Flag for errors
  SHORT status;

  // Get Mesh Data
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT v_per_elm = (dim + 1);
  INT dof_per_elm = (dim + 1);
  INT i,j;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in fem to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* p = fem->dwork;
  REAL* dp = fem->dwork + v_per_elm;
  PX_H1_basis(p,dp,x,v_on_elm,1,sc);

  // face to vertex map
  INT fv[dim];// = (INT *)calloc(dim,sizeof(INT));
  INT elnd,ef1,ef2,ef3;//face endpoint vertex tracking numbers

  REAL gradp;

  if(dim==2){
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i],fem->f_v,fv);
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
      /* phi[i*dim] =   ABS(fem->f_norm[dim*(dof[i])]  )* 4*p[ef1]*p[ef2]; */
      /* phi[i*dim+1] = ABS(fem->f_norm[dim*(dof[i])+1]) * 4*p[ef1]*p[ef2]; */
      phi[i*dim] =   fem->f_norm[dim*(dof[i])]* 4*p[ef1]*p[ef2];
      phi[i*dim+1] = fem->f_norm[dim*(dof[i])+1]* 4*p[ef1]*p[ef2];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 4*(p[ef1]*dp[ef2*dim+j] + dp[ef1*dim+j]*p[ef2]);

        /* dphi[i*dim*dim + j*dim + 0] = gradp * ABS(fem->f_norm[dim*(dof[i])+0]); */
        /* dphi[i*dim*dim + j*dim + 1] = gradp * ABS(fem->f_norm[dim*(dof[i])+1]); */
        dphi[i*dim*dim + j*dim + 0] = gradp * fem->f_norm[dim*(dof[i])+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * fem->f_norm[dim*(dof[i])+1];
      }

    }
  } else if(dim==3) {
    for (i=0;i<dof_per_elm;i++) {
      get_incidence_row(dof[i],fem->f_v,fv);
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
      phi[i*dim] = fem->f_norm[dim*(dof[i])] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+1] = fem->f_norm[dim*(dof[i])+1] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+2] = fem->f_norm[dim*(dof[i])+2] * 8*p[ef1]*p[ef2]*p[ef3];

      // Gradient
      for(j=0;j<dim;j++) {
        gradp = 8*(p[ef1]*p[ef2]*dp[ef3*dim+j] + p[ef1]*dp[ef2*dim+j]*p[ef3] + dp[ef1*dim+j]*p[ef2]*p[ef3]);

        dphi[i*dim*dim + j*dim + 0] = gradp * fem->f_norm[dim*(dof[i])+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * fem->f_norm[dim*(dof[i])+1];
        dphi[i*dim*dim + j*dim + 2] = gradp * fem->f_norm[dim*(dof[i])+2];
      }

    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Clean up working array for next person to use.
  array_set(v_per_elm*(dim+1),fem->dwork,0.0);

  return;
}
/****************************************************************************************************************************/

/*******************************************************************************************************/
/*!
* \fn void MINI_basis(REAL *p,REAL *dp,REAL *x,INT *v_on_elm,scomplex *sc)
*
* \brief Compute a MINI element (P1 + cubic bubble) at a particular point in 1, 2 or 3D
*
* \param x       Coordinate on where to compute basis function
* \param v_on_elm     vertices on element
* \param sc      Simplicial complex
*
* \return p      Basis functions (1 for each DOF on element)
* \return dp     Derivatives of basis functions (i.e., gradient)
*
*  \note We call the P1 basis routine, then augment with cubic bubble:
*         b = lam1*lam2*lam3 (in 2D)
*        where lam_i are the P1 basis functions.
*
*     DoF of the bubble is eta(f) = c(int_T (f) - f(v1)*int_T lam1 - f(v2)*int_T lam2 - f(v3)*int_T lam3)
*      c = (int_T lam1*lam2*lam3)^{-1}  (again in 2D)
*        = d(2d+1)!/(d! * |V|)  where |V| is volume of element (area in 2D)
*
*/
void MINI_basis(REAL *p,REAL *dp,REAL *x,INT *v_on_elm,scomplex *sc)
{

  // Get Mesh Data
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT v_per_elm = (dim + 1);

  INT i,j,k;
  REAL temp;

  /* Get Linear Basis Functions for particular element */
  // We use the working double arrays in fem to store the values
  // This way we do not need to reallocate each time this is called.
  REAL* lam = fem->dwork;
  REAL* dlam = fem->dwork + v_per_elm;
  PX_H1_basis(lam,dlam,x,v_on_elm,1,sc);

  for(i=0;i<v_per_elm;i++) {
    p[i] = lam[i];
    for(j=0;j<dim;j++) {
      dp[i*dim+j] = dlam[i*dim+j];
    }
  }
  // bubbles
  p[v_per_elm] = 1.0;
  for(i=0;i<v_per_elm;i++) {
    p[v_per_elm] = p[v_per_elm]*lam[i];
  }

  // Gradient of bubble
  for(k=0;k<dim;k++) {
    dp[v_per_elm*dim+k] = 0.0;
    for(i=0;i<v_per_elm;i++) {
      temp = 1.0;
      for(j=0;j<v_per_elm;j++) {
        if(j==i) {
          temp = temp*dlam[i*dim+k];
        } else {
          temp = temp*lam[j];
        }
      }
      dp[v_per_elm*dim+k] += temp;
    }
  }

  return;
}
/*******************************************************************************************************/


/****************************************************************************************************************************/
/*!
* \fn void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc,fespace *FE)
*
* \brief Grabs the basis function of a FEM space at a particular point in 2 or 3D
*
* \param x         Coordinate on where to compute basis function
* \param v_on_elm  Vertices on element
* \param dof       DOF on element
* \param sc        Simplicial complex
* \param FE        Fespace struct
*
* \return phi      Basis functions
* \return dphi     Derivatives of basis functions (depends on type)
*
*/
void get_FEM_basis(REAL *phi,REAL *dphi,REAL *x,INT *v_on_elm,INT *dof,scomplex *sc,fespace *FE)
{
  // Flag for erros
  SHORT status;

  // Mesh and FEM Data
  sc_fem *fem = sc->fem;
  INT FEtype = FE->FEtype;
  INT dim = sc->dim;
  INT v_per_elm = dim+1;
  INT offset = FE->dof_per_elm/dim;

  if(FEtype>=0 && FEtype<10) { // PX elements

    PX_H1_basis(phi,dphi,x,dof,FEtype,sc);

  } else if(FEtype==20) { // Nedelec elements
    // Compute P1 basis (lam/dlam) then call local ned0_basis
    REAL lam[v_per_elm], dlam[v_per_elm*dim];
    PX_H1_basis(lam,dlam,x,v_on_elm,1,sc);
    // Build local v_on_ed and ed_len from sc->fem for edges on this element
    INT ed_per_elm = dim*(dim+1)/2;
    INT v_per_ed = 2;
    INT v_on_ed[ed_per_elm*v_per_ed];
    REAL ed_len[ed_per_elm];
    for(INT i=0; i<ed_per_elm; i++) {
      INT eid = dof[i];
      INT ica = fem->ed_v->IA[eid];
      v_on_ed[i*v_per_ed+0] = fem->ed_v->JA[ica];
      v_on_ed[i*v_per_ed+1] = fem->ed_v->JA[ica+1];
      ed_len[i] = fem->ed_len[eid];
    }
    ned0_basis(phi,dphi,lam,dlam,dim,v_on_elm,v_on_ed,ed_len);

  } else if(FEtype==30) { // Raviart-Thomas elements
    // Compute P1 basis (lam/dlam) then call local rt0_basis
    REAL lam[v_per_elm], dlam[v_per_elm*dim];
    PX_H1_basis(lam,dlam,x,v_on_elm,1,sc);
    // Build local v_on_f, f_area, xv from sc->fem for faces on this element
    INT f_per_elm = dim+1;
    INT v_per_f = dim;
    INT v_on_f[f_per_elm*v_per_f];
    REAL f_area[f_per_elm];
    REAL xv[v_per_elm*dim];
    for(INT i=0; i<f_per_elm; i++) {
      INT fid = dof[i];
      INT ica = fem->f_v->IA[fid];
      for(INT j=0; j<v_per_f; j++)
        v_on_f[i*v_per_f+j] = fem->f_v->JA[ica+j];
      f_area[i] = fem->f_area[fid];
    }
    for(INT i=0; i<v_per_elm; i++)
      for(INT j=0; j<dim; j++)
        xv[i*dim+j] = sc->x[v_on_elm[i]*dim+j];
    rt0_basis(phi,dphi,lam,dlam,dim,v_on_elm,v_on_f,f_area,xv);

  } else if(FEtype==60) { // Vector element

    PX_H1_basis(phi, dphi, x, dof, 1, sc);
    if (dim>1){
      PX_H1_basis(phi+offset, dphi+FE->dof_per_elm, x, dof, 1, sc);
    }
    if (dim>2){
      PX_H1_basis(phi+offset*2, dphi+FE->dof_per_elm*2, x, dof, 1, sc);
    }

  } else if(FEtype==61) { // Bubble element

    bubble_face_basis(phi,dphi,x,v_on_elm,dof,sc);

  } else if(FEtype==99) { // Constraint Single DOF Space

    phi[0] = 1.0;

  } else if(FEtype==103) { // MINI element

    MINI_basis(phi,dphi,x,v_on_elm,sc);

  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
#endif /* bubble_face_basis, MINI_basis, get_FEM_basis */
