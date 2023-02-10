/*! \file examples/basic_elliptic/eafe_data.h
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the Data parameters and coefficients
*        for the basic_elliptic examples.  This includes exact solutions,
*        RHS functions, and boundary conditions.
*
* \note We include different examples for 2D and 3D problems, and for FE types
*       P1, P2, RT0, Nedelec0
*/

// Exact Solution (if we have one)
void exactsol(REAL *val,REAL* x,REAL time,	\
	      const INT dim,void *param)
{
  switch(dim){
  case 1:
    *val = -x[0]*x[0];
    break;
  case 2:
    *val = 4*x[0]*x[0]-5e0*x[1]*x[1];
    break;
  case 3:
    *val = 3*x[0]*x[0]+x[1]*x[1]-5e0*x[2]*x[2];
    break;
  default:
    *val = 0e0;
    break;
  }
}
// PDE Coefficients
/*====================================================================*/
void diffusion(REAL *val,REAL *x, const REAL time,	\
	       const INT dim, void *params)
{
  // diffusion coefficient. 
  val[0]=1e0;
  return;
}
/*====================================================================*/
void advection(REAL *val,REAL *x, const REAL time,	\
	       const INT dim, void *param)
{
  // advection coefficient (vector valued);
  memset(val,0,dim*sizeof(REAL));
  //  val[0]=1e0;
  return;
}
/*EOF*/
void reaction(REAL *val,REAL* x,REAL time,	\
	      const INT dim, void *param)
{
  // reaction coefficieent
  val[0] = 0e0;
}
//
void rhs_pde(REAL *val,REAL* x,REAL time,	\
	     const INT dim, void *param)
{
  REAL u;
  // right hand side for the PDE
  exactsol(val,x,time,dim,param);
  val[0]=2e0; // + 2e0*val[0];
}
void rhs_natural_bc(REAL *val,REAL* x,REAL time,	\
		    const INT dim, void *param)
{
  // right hand side for the natural boundary condition
  // (flux.n=rhs_natural_bc)
  val[0] = 0e0;
}

void rhs_essential_bc(REAL *val,REAL* x,REAL time,	\
		      const INT dim,void *param)
{
  // dirichlet boundary conditions: set up through exactsol().
  exactsol(val,x,time,dim,param);
}
