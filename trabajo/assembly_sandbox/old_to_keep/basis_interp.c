#include "hazmath.h"
#include "biot_headers.h"
void zpx_h1_basis(REAL *p,REAL *dp,REAL *x,INT  dim, REAL *xel,INT  porder)
{
  REAL dp0r,dp1r,dp2r,dp3r,dp4r,dp5r,dp6r,dp7r,dp8r,dp9r;
  REAL dp0s,dp1s,dp2s,dp3s,dp4s,dp5s,dp6s,dp7s,dp8s,dp9s;
  REAL dp0t,dp1t,dp2t,dp3t,dp4t,dp5t,dp6t,dp7t,dp8t,dp9t;
  REAL onemrst;

  // flag for errors
  short status;

  // get mesh data
  // store coordinates of vertices of element (2 in 1d, 3 in 2d, 4 in 3d)
  REAL xv0,xv1,xv2,xv3;
  REAL yv0,yv1,yv2,yv3;
  REAL zv0,zv1,zv2,zv3;

  // p0 elements are trivial and we just need to return a 1 for each element:
  if(porder==0) {
    p[0] = 1.0;
  } else {
    // the remaining depend on the dimension
    if(dim==1) {

      // get physical coordinates of vertices
      xv0 = xel[0];
      xv1 = xel[1];

      // get barycentric coordinates
      REAL oneoverh = 1.0/(xv1-xv0);
      REAL lam0 = (xv1-x[0])*oneoverh;
      REAL lam1 = (x[0]-xv0)*oneoverh;

      // now get basis functions
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

      // get physical coordinates of vertices
      xv0 = xel[0*dim+0];
      xv1 = xel[1*dim+0];
      xv2 = xel[2*dim+0];
      yv0 = xel[0*dim+1];
      yv1 = xel[1*dim+1];
      yv2 = xel[2*dim+1];

      // get coordinates on reference triangle
      REAL det = (xv1-xv0)*(yv2-yv0) - (xv2-xv0)*(yv1-yv0);

      REAL r = ((yv2-yv0)*(x[0]-xv0) + (xv0-xv2)*(x[1]-yv0))/det;

      REAL drdx = (yv2-yv0)/det;
      REAL drdy = (xv0-xv2)/det;

      REAL s = ((yv0-yv1)*(x[0]-xv0) + (xv1-xv0)*(x[1]-yv0))/det;

      REAL dsdx = (yv0-yv1)/det;
      REAL dsdy = (xv1-xv0)/det;

      /*  get the basis functions for linear elements on each node.
      *  the basis functions can now be evaluated in terms of the
      *  reference coordinates r and s.  it's also easy to determine
      *  the values of the derivatives with respect to r and s.
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

      /*  we need to convert the derivative information from (r(x,y),s(x,y))
      *  to (x,y) using the chain rule.
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

      // get nodes and physical coordinates
      xv0 = xel[0*dim+0];
      xv1 = xel[1*dim+0];
      xv2 = xel[2*dim+0];
      xv3 = xel[3*dim+0];
      yv0 = xel[0*dim+1];
      yv1 = xel[1*dim+1];
      yv2 = xel[2*dim+1];
      yv3 = xel[3*dim+1];
      zv0 = xel[0*dim+2];
      zv1 = xel[1*dim+2];
      zv2 = xel[2*dim+2];
      zv3 = xel[3*dim+2];

      // get coordinates on reference triangle
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

      /*  get the basis functions for linear elements on each node.
      *  the basis functions can now be evaluated in terms of the
      *  reference coordinates r and s.  it's also easy to determine
      *  the values of the derivatives with respect to r and s.
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

      /*  we need to convert the derivative information from (r(x,y),s(x,y))
      *  to (x,y) using the chain rule.
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
/*************************************************************/
void zrt_basis(REAL *phi,REAL *dphi,REAL *x,simplex_data *splex,	\
	       REAL *wrk)
{
  // flag for erros
  short status;

  // get mesh data
  INT  v_per_elm = splex->nv;
  INT  f_per_elm = splex->nf;
  INT  dim = splex->dim;
  INT  *fv=splex->fv;
  INT  *v_on_elm=splex->v;
  INT  i,j,ica,icb,jcnt;
  REAL farea;
  INT  elnd,ef1,ef2,ef3;

  /* get linear basis functions for particular element */
  // we use the working double arrays in mesh to store the values
  // this way we do not need to REALlocate each time this is called.
  REAL* p = wrk;
  REAL* dp = p + splex->nv;
  REAL *wrkend=dp+splex->nv*dim;
  zpx_h1_basis(p,dp,x,splex->dim,splex->x,1);

  // go through each face and find the corresponding nodes
  if(dim==2) {
    for (i=0; i<f_per_elm; i++) {
      farea = splex->f_area[i];
      // loop through nodes on element to find corresponding nodes and get correct orientation
      for(j=0;j<v_per_elm;j++) {
        elnd = v_on_elm[j];
        if(fv[i*dim+0]==elnd) {
          ef1 = j;
        }
        if(fv[i*dim+1]==elnd) {
          ef2 = j;
        }
      }
      /* now, with the linear basis functions p, dpx, and dpy, the rt elements are in 2d
       * phi_fij = |fij|*(p(i)curl(p(j)) - p(j)curl(p(i)))
       * |fij| = |eij|
       */
      phi[i*dim+0] = farea*(p[ef1]*dp[ef2*dim+1] - p[ef2]*dp[ef1*dim+1]);
      phi[i*dim+1] = farea*(-p[ef1]*dp[ef2*dim] + p[ef2]*dp[ef1*dim]);
      
      // compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i)))
      dphi[i] = 2*farea*(dp[ef1*dim]*dp[ef2*dim+1] - dp[ef2*dim]*dp[ef1*dim+1]);
    }
  } else if(dim==3) {
    for (i=0; i<f_per_elm; i++) {
      // get the area
      farea = splex->f_area[i];
      // loop through nodes on element to find corresponding nodes for correct orienation
      for(j=0;j<v_per_elm;j++) {
        elnd = v_on_elm[j];
        if(fv[i*dim+0]==elnd) {
          ef1 = j;
        }
        if(fv[i*dim+1]==elnd) {
          ef2 = j;
        }
        if(fv[i*dim+2]==elnd) {
          ef3 = j;
        }
      }
      /* now, with the linear basis functions p, dpx, and dpy, the rt elements are in 3d
      * phi_fijk = 6*|fijk|*(p(i)(grad(p(j)) x grad(p(k))) - p(j)(grad(p(k)) x grad(p(i))) + p(k)(grad(p(i)) x grad(p(j))))
      * |fijk| = area(face)
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

      // compute divs div(phi_fij) = 2*|fij|(dx(p(i))*dy(p(j)) - dx(p(j))*dy(p(i)))
      dphi[i] = 6*farea*(dp[ef1*dim]*(dp[ef2*dim+1]*dp[ef3*dim+2]-dp[ef2*dim+2]*dp[ef3*dim+1])
      + dp[ef1*dim+1]*(dp[ef2*dim+2]*dp[ef3*dim]-dp[ef2*dim]*dp[ef3*dim+2])
      + dp[ef1*dim+2]*(dp[ef2*dim]*dp[ef3*dim+1]-dp[ef2*dim+1]*dp[ef3*dim]));
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  return;
}
/*****************************************************************************/
void zbubble_face_basis(REAL *phi, REAL *dphi, REAL *x,
			simplex_data *splex, REAL *wrk)
{
  // flag for errors
  short status;

  // get mesh data
  INT  dim = splex->dim;
  INT  v_per_elm = splex->nv;
  INT  dof_per_elm = splex->nf;
  INT  *v_on_elm=splex->v;
  REAL  *f_norm=splex->f_norm;
  INT  i,j;

  /* get linear basis functions for particular element */
  // we use the working double arrays in mesh to store the values
  // this way we do not need to REALlocate each time this is called.
  REAL* p = wrk;
  REAL* dp = p + splex->nv;
  REAL *wrkend=dp+splex->nv*dim;
  zpx_h1_basis(p,dp,x,splex->dim,splex->x,1);
  // face to vertex map
  INT  *fv=splex->fv;// = (INT  *)calloc(dim,sizeof(INT ));
  INT  elnd,ef1,ef2,ef3;//face endpoINT  vertex tracking numbers
  REAL gradp;
  if(dim==2){
    for (i=0;i<dof_per_elm;i++) {
      // find orientation of face
      for(j=0;j<v_per_elm;j++){
        elnd = v_on_elm[j];
        if(fv[i*dim+0]==elnd) {
          ef1 = j;
        }
        if(fv[i*dim+1]==elnd) {
          ef2 = j;
        }
      }
      // multiply basis function by normal vector
      //      print_full_mat(1,dim,(f_norm+dim*i),"iiifnorm");
      phi[i*dim] =   f_norm[dim*i]* 4*p[ef1]*p[ef2];
      phi[i*dim+1] = f_norm[dim*i+1]* 4*p[ef1]*p[ef2];
      // gradient
      for(j=0;j<dim;j++) {
        gradp = 4*(p[ef1]*dp[ef2*dim+j] + dp[ef1*dim+j]*p[ef2]);
	//
        dphi[i*dim*dim + j*dim + 0] = gradp * f_norm[dim*i+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * f_norm[dim*i+1];
      }
    }
  } else if(dim==3) {
    for (i=0;i<dof_per_elm;i++) {
      // find orientation of face
      for(j=0;j<v_per_elm;j++){
        elnd = v_on_elm[j];
        if(fv[i*dim+0]==elnd) {
          ef1 = j;
        }
        if(fv[i*dim+1]==elnd) {
          ef2 = j;
        }
        if(fv[i*dim+2]==elnd) {
          ef3 = j;
        }
      }
      // multiply basis function by normal vector
      phi[i*dim+0] = f_norm[dim*i+0] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+1] = f_norm[dim*i+1] * 8*p[ef1]*p[ef2]*p[ef3];
      phi[i*dim+2] = f_norm[dim*i+2] * 8*p[ef1]*p[ef2]*p[ef3];
      // gradient
      for(j=0;j<dim;j++) {
        gradp = 8*(p[ef1]*p[ef2]*dp[ef3*dim+j] + p[ef1]*dp[ef2*dim+j]*p[ef3] + dp[ef1*dim+j]*p[ef2]*p[ef3]);

        dphi[i*dim*dim + j*dim + 0] = gradp * f_norm[dim*i+0];
        dphi[i*dim*dim + j*dim + 1] = gradp * f_norm[dim*i+1];
        dphi[i*dim*dim + j*dim + 2] = gradp * f_norm[dim*i+2];
      }
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  return;
}
/***************************************************************/
INT zget_fem_basis(REAL *phi,REAL *dphi,REAL *x,	\
		   INT  fetype,simplex_data *splex,REAL *wrk)
// returns the num of dofs for the space of fetype;
{
  // flag for erros
  short status;

  // mesh and fem data
  INT  dim=splex->dim,dof_per_elm=-1;
  dof_per_elm=calc_ndofs(fetype,dim,					\
			  splex->nv,splex->ne,splex->nf,splex->ns);
  if(fetype>=0 && fetype<10) { // px elements
    zpx_h1_basis(phi,dphi,x,splex->dim,splex->x,fetype);
  } else if(fetype==30) { // raviart-thomas elements
    zrt_basis(phi,dphi,x,splex,wrk);
  } else if(fetype==60) { // vector element
    zpx_h1_basis(phi, dphi, x, splex->dim,splex->x, 1);
    if (dim>1){
      zpx_h1_basis(phi + splex->nv, dphi + splex->nv, x,	\
		   dim, splex->x,1);
    }
    if (dim>2){
      zpx_h1_basis(phi + splex->nv*2, dphi + 2*splex->nv, x,	\
		   dim,splex->x,1);
    }
  } else if(fetype==61) { // bubble element
    zbubble_face_basis(phi,dphi,x,splex,wrk);
    /* fprintf(stdout,"\nx=[%.15f,%.15f];  ",x[0],x[1]); */
    /* fprintf(stdout,"phi1(x) = [%.15e,%.15e];  ",phi[0*dim+0],phi[0*dim+1]); */
    /* fprintf(stdout,"phi2(x) = [%.15e,%.15e];  ",phi[1*dim+0],phi[1*dim+1]); */
    /* fprintf(stdout,"phi3(x) = [%.15e,%.15e];\n",phi[2*dim+0],phi[2*dim+1]); */
    /* fprintf(stdout,"dphi1(x) = [%.15e,%.15e\n            %.15e,%.15e];\n  ",dphi[0*dim*dim+0*dim+0],	\ */
    /* 	    dphi[0*dim*dim+1*dim+0],					\ */
    /* 	    dphi[0*dim*dim+0*dim+1],					\ */
    /* 	    dphi[0*dim*dim+1*dim+1]); */

    /* fprintf(stdout,"\n"); */
    //    fprintf(stdout,"grad  = [%e,%e]\n",val_sol[0],val_sol[1]);
  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }
  return dof_per_elm;
}
/*********************interpolations:***************************/
/***************************************************************/
void zfe_interp(REAL* val,REAL *u,REAL* x,	\
		INT fetype,simplex_data *splex, REAL *wrk)
{
  INT  i,j;
  // get fe and mesh data
  INT  dim = splex->dim,dof_per_elm=-1;
  REAL *phi = wrk;
  REAL *dphi = phi + dim*splex->nf; // this is overkill: always assume vector
  REAL *wrkend  = dphi + dim*dim*splex->nf;
  dof_per_elm=zget_fem_basis(phi,dphi,x,fetype,splex,wrkend);
  if(fetype<20) { // scalar element
    val[0] = 0.0;
    for(j=0; j<dof_per_elm; j++)
      val[0] += u[j]*phi[j];
  } else { // vector element
    for(i=0;i<dim;i++) {
      val[i] = 0.0;
      for(j=0; j<dof_per_elm; j++){
	val[i] += u[j]*phi[j*dim+i];
      }
    }
  }
  return;
}
/****************************************************************************************************************************/
void zfe_dinterp(REAL* val,REAL *u,REAL *x,
		 INT fetype, simplex_data *splex, REAL *wrk)
{
  INT  j,k,i;

  // get fe and mesh data
  INT  dim = splex->dim,dof_per_elm=-1;
  REAL *phi = wrk;
  REAL *dphi = phi + dim*splex->nf; // this is overkill: always assume vector
  REAL *wrkend = dphi + dim*dim*splex->nf; // 
  // basis functions and its derivatives if necessary
  dof_per_elm=zget_fem_basis(phi,dphi,x,fetype,splex,wrkend);
  if(fetype==0) { // don't compute derivatives of p0 elements (set to 0)
    for(j=0;j<dim;j++) {
      val[j] = 0.0;
    }
  } else if(fetype<20 && fetype!=0) { // scalar element
    for(j=0;j<dim;j++) {
      val[j] = 0.0;
      for(k=0; k<dof_per_elm; k++)
        val[j] += u[k]*dphi[k*dim+j];
    }
  } else if (fetype==30) { // raviart-thomas (div is scalar)
    val[0] = 0.0;
    for (j=0; j<dof_per_elm; j++) {
      val[0] += u[j]*dphi[j];
      //      fprintf(stdout,"\nRTRT_val=%.15e=(%.5e,%.5e,%.5e) * (%e,%e,%e)",val[0],u[0],u[1],u[2],dphi[0],dphi[1],dphi[2]);
    }
  } else if (fetype>=60) { // vector lagrange functions (i.e. bubbles) -> gradients of vectors -> tensor
    for(j=0;j<dim;j++) {
      for(i=0;i<dim;i++) {
        val[j*dim + i] = 0.0;
        for(k=0; k<dof_per_elm; k++){
	  /* if(fabs(u[k])>1e-10) */
	  /*   //	  if(fabs(dphi[k*dim*dim + j*dim + i])>1e-10) */
	  /*   fprintf(stdout,"\nu(%d)=%.7e:%.7e",k,u[k],dphi[k*dim*dim + j*dim + i]); */
          val[j*dim + i] += u[k]*dphi[k*dim*dim + i*dim + j];
	}
	//	if(fetype==61){
	//	  fprintf(stdout,"\nG(%d,%d)=%.8e (%.15e,%.15e,%.15e)",i,j,val[j*dim + i],u[0],u[1],u[2]);
	//}
      }
    }
  } else {
    check_error(ERROR_FE_TYPE,__FUNCTION__);
  }

  return;
}
/********************************************************************/
void zblockfe_interp(REAL* val,REAL *u,REAL* x,simplex_data *splex)
{
  INT  k;
  INT  dim = splex->dim;

  INT *fetypes=splex->fetypes;
  INT *ndofs=splex->ndofs;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<splex->nspaces;k++) {
    //    fprintf(stdout,"\nfetypefetypefetype=%i",fetypes[k]);
    zfe_interp(val_sol,u_comp,x,fetypes[k],splex,splex->wrk);
    if(fetypes[k]<20) { // scalar
      val_sol++;
    } else { // vector
      val_sol += dim;
    }
    u_comp += ndofs[k];
  }
  return;
}
/**************************************************************************/
void zblockfe_dinterp(REAL* val,REAL *u,REAL* x,
		      simplex_data *splex)
{
  INT  k;
  INT  dim = splex->dim;
  INT *fetypes=splex->fetypes, *ndofs=splex->ndofs;
  REAL* val_sol = val;
  REAL* u_comp = u;

  for(k=0;k<splex->nspaces;k++) {
    zfe_dinterp(val_sol,u_comp,x,fetypes[k],splex,splex->wrk);
    if(fetypes[k]<20) { // scalar
      /* if(fetypes[k]>0){ */
      /* 	fprintf(stdout,"\nfetypefetypefetype=%i(%d)",fetypes[k],ndofs[k]); */
      /* 	fprintf(stdout,"\nx=[%f,%f];  ",x[0],x[1]); */
      /* 	fprintf(stdout,"uDOF= [%e,%e,%e];  ",u_comp[0],u_comp[1],u_comp[2]); */
      /* 	fprintf(stdout,"grad  = [%e,%e]\n",val_sol[0],val_sol[1]); */
      /* } */
      val_sol += dim;
    } else if(fetypes[k]==30) { // div is scalar
      val_sol++;
    } else if(fetypes[k]>=60) { // vectorlagrange (i.e., bubbles) gradient is matrix.
      val_sol+=dim*dim;
    } else {
      check_error(ERROR_FE_TYPE,__FUNCTION__);
    }
    u_comp += ndofs[k];
  }
  return;
}
/**********************************************************************/
/*************************************************************************/
void zbubble_hessian2d(REAL *depsp,REAL *dtrepsp,		\
		       REAL *u, REAL *x, simplex_data *splex)
{
  // computes div(eps(uh)) and div(div(u)I) for uh\in FE space. 
  // 3d also probably works but no time for it. 
  // Flag for errors
  SHORT status;
  // Get Mesh Data
  INT v_per_elm = splex->nv;
  INT dim = splex->dim;
  INT i,j,k,l,ij,ijk,ijkl;
  REAL *wrk=(REAL *)splex->wrk;
  REAL *p = wrk;
  REAL *dp = p + v_per_elm;
  REAL *zdtrepsp=dp+dim*v_per_elm;
  REAL *wrkend=zdtrepsp+dim;
  zpx_h1_basis(p,dp,x,splex->dim,splex->x,1);
  // face to vertex map
  INT *fv=splex->fv;
  REAL *f_norm=splex->f_norm;
  INT *fetypes=splex->fetypes;
  INT *ndofs=splex->ndofs;
  INT elnd,ef1,ef2,ef3;//face endpoint vertex tracking numbers
  REAL gradp,trhessp,hessp_jk;
  if(dim==2){
    memset(depsp,0,dim*sizeof(REAL));
    memset(dtrepsp,0,dim*sizeof(REAL));
    for (i=0;i<ndofs[0];i++) {
      //      fprintf(stdout,"\n%d:%d:bubble_dof(%d)=%e",fv[0],fv[1],i,u[i]);
      //      print_full_mat(1,dim,(f_norm+i*dim),"XXXfnorm");
      // Find orientation of face
      for(j=0;j<v_per_elm;j++){
        elnd = splex->v[j];
        if(fv[i*dim+0]==elnd) ef1 = j;
        if(fv[i*dim+1]==elnd) ef2 = j;
      }
      // Gradient
      for(j=0;j<dim;j++) {
	ij=i*dim + j;
	zdtrepsp[j]=0.;
	trhessp=0.;
	for(k=0;k<dim;k++){
	  // 0.5*hessian
	  hessp_jk=4.*(dp[ef2*dim+k]*dp[ef1*dim+j]+dp[ef2*dim+j]*dp[ef1*dim+k]);
	  if(j==k) trhessp+=hessp_jk; // this is tr(H(phi))
	  //zdtrepsp[ij]+=hessp_jk*mesh->f_norm[dof[i]+k];// this is H(phi)*n
	  zdtrepsp[j]+=(hessp_jk*f_norm[i*dim+k]);// this is H(phi)*n
	  //	  fprintf(stdout,"\n%d:%d:XXXH(%d,%d)=%e",ef1,ef2,j,k,hessp_jk);
	}
	//	zdepsp[ij]=0.5*dtrepsp[ij]+trhessp*mesh->f_norm[dof[i]+j];
	depsp[j]+=u[i]*(0.5*zdtrepsp[j] + trhessp*f_norm[i*dim+j]);
	dtrepsp[j]+=u[i]*zdtrepsp[j];// this is H(phi)*n
      }
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  return;
}
/********************************************************************/
void zquad_elm(qcoordinates *cqelm,simplex_data *splex,INT nq1d)
{
  // Flag for errors
  SHORT status;

  /* Loop indices */
  INT q,j;

  /* Dimension */
  INT dim = splex->dim;

  /* Total Number of Quadrature Nodes */
  INT nq=nq1d;
  for(j=1;j<dim;j++) nq*=nq1d;
  //  INT nq = pow(nq1d,dim);

  /* Coordinates of vertices of element */
  INT v_per_elm = splex->nv;
  //  coordinates *cvelm = allocatecoords(v_per_elm,dim);
  REAL *cvelm=splex->x;
  //  print_full_mat(v_per_elm,dim,cvelm,"zcvelm");
  /* Gaussian points for reference element */
  REAL* gp = (REAL *) calloc(dim*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  /* Points on Reference Triangle */
  REAL r,s,t;

  REAL e_vol = splex->vol; /* Area/Volume of Element */
  REAL voldim=0.0;

  // Get vertices of element
  //  INT* thiselm_v = splex->v;
  /* INT* thiselm_v = (INT *) calloc(v_per_elm,sizeof(INT)); */
  /* iCSRmat* el_v = mesh->el_v; */
  /* get_incidence_row(elm,el_v,thiselm_v); */

  // Get coordinates of vertices for given element
  if(dim==1) {
    /* for (j=0; j<v_per_elm; j++) { */
    /*   //      cvelm->x[j] = mesh->cv->x[thiselm_v[j]]; */
    /*   //      cvelm[0*dim+j] should be splex->x[0*dim+j];  */
    /* } */
    voldim = 0.5*e_vol;

    // Get Quad Nodes and Weights
    quad1d(gp,gw,nq1d);

    // Map to Real Interval
    // x = 0.5*(x2-x1)*r + 0.5*(x1+x2)
    // w = 0.5*(x2-x1)*wref
    for (q=0; q<nq; q++) {
      r = gp[q];
      cqelm->x[q] = 0.5*(cvelm[0*dim+0]*(1-r) + cvelm[1*dim+0]*(1+r));
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==2) {
    /* for (j=0; j<v_per_elm; j++) { */
    /*   cvelm->x[j] = mesh->cv->x[thiselm_v[j]]; */
    /*   cvelm->y[j] = mesh->cv->y[thiselm_v[j]]; */
    /* } */
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
      //      cqelm->x[q] = cvelm->x[0]*(1-r-s) + cvelm->x[1]*r + cvelm->x[2]*s;
      //      cqelm->y[q] = cvelm->y[0]*(1-r-s) + cvelm->y[1]*r + cvelm->y[2]*s;
      cqelm->x[q] =
	cvelm[0*dim+0]*(1-r-s) +		\
	cvelm[1*dim+0]*r +			\
	cvelm[2*dim+0]*s;
      cqelm->y[q] = \
	cvelm[0*dim+1]*(1-r-s) +	     \
	cvelm[1*dim+1]*r +		     \
	cvelm[2*dim+1]*s;
      cqelm->w[q] = voldim*gw[q];
    }
  } else if(dim==3) {
    /* for (j=0; j<v_per_elm; j++) { */
    /*   cvelm->x[j] = mesh->cv->x[thiselm_v[j]]; */
    /*   cvelm->y[j] = mesh->cv->y[thiselm_v[j]]; */
    /*   cvelm->z[j] = mesh->cv->z[thiselm_v[j]]; */
    /* } */
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
      cqelm->x[q] = 
	cvelm[0*dim+0]*(1-r-s-t) +		\
	cvelm[1*dim+0]*r +			\
	cvelm[2*dim+0]*s +			\
	cvelm[3*dim+0]*t;
      cqelm->y[q] =		   \
	cvelm[0*dim+1]*(1-r-s-t) + \
	cvelm[1*dim+1]*r +	   \
	cvelm[2*dim+1]*s +	   \
	cvelm[3*dim+1]*t;
      cqelm->z[q] =				\
	cvelm[0*dim+2]*(1-r-s-t) +		\
	cvelm[1*dim+2]*r +			\
	cvelm[2*dim+2]*s +			\
	cvelm[3*dim+2]*t;
      cqelm->w[q] = voldim*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  /* if(cvelm) { */
  /*   free_coords(cvelm); */
  /*   free(cvelm); */
  /*   cvelm=NULL; */
  /* } */

  return;
}
/*******************************************************************************/
void zquad_face(qcoordinates *cqbdry,INT nq1d, INT dim, REAL *xf, REAL farea)
{
  // Flag for errors
  SHORT status;

  INT q,j; /* Loop indices */
  INT nq = 0;

  if(dim==2) { // face is edge in 2D
    nq = nq1d;
  } else if(dim==3) {
    nq = nq1d*nq1d;
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  /* Coordinates of vertices of face are in xf ordered by rows (every
     vertex corresponds to a row)*/

  // Gaussian points for reference element
  // (edges: [-1,1] faces: tri[(0,0),(1,0),(0,1)])
  REAL* gp = (REAL *) calloc((dim-1)*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r,s;      	/* Points on Reference Face */

  REAL w = 0.0;
  if(dim==2) { // Faces are Edges in 2D
    w = 0.5*farea; /* Jacobian = 1/2 |e| */
  } else if(dim==3) {
    w = 2.*farea; /* Jacobian = 2*|f| */
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get coordinates of vertices for given edge/face
  //  get_incidence_row(dof,dof_v,thisdof_v);
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
      /* cqbdry->x[q] = 0.5*cvdof->x[0]*(1-r) + 0.5*cvdof->x[1]*(1+r); */
      /* cqbdry->y[q] = 0.5*cvdof->y[0]*(1-r) + 0.5*cvdof->y[1]*(1+r); */
      cqbdry->x[q] =				\
	0.5*(xf[0*dim+0]*(1-r) +		\
	     xf[1*dim+0]*(1+r));
      cqbdry->y[q] =				\
	0.5*(xf[0*dim+1]*(1-r) +			\
	     xf[1*dim+1]*(1+r));
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
      /* cqbdry->x[q] = cvdof->x[0]*(1-r-s) + cvdof->x[1]*r + cvdof->x[2]*s; */
      /* cqbdry->y[q] = cvdof->y[0]*(1-r-s) + cvdof->y[1]*r + cvdof->y[2]*s; */
      /* cqbdry->z[q] = cvdof->z[0]*(1-r-s) + cvdof->z[1]*r + cvdof->z[2]*s; */
      cqbdry->x[q] =	      \
	xf[0*dim+0]*(1-r-s) + \
	xf[1*dim+0]*r +       \
	xf[2*dim+0]*s;
      cqbdry->y[q] =	      \
	xf[0*dim+1]*(1-r-s) + \
	xf[1*dim+1]*r +       \
	xf[2*dim+1]*s;
      cqbdry->z[q] =	      \
	xf[0*dim+2]*(1-r-s) + \
	xf[1*dim+2]*r +       \
	xf[2*dim+2]*s;
      cqbdry->w[q] = w*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  return;
}
/*********************************************************************************/
