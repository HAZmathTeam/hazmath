/*! \file examples/reaction_diffusion/rd_data.h
*
*  Created by Casey Cavanaugh on 6/5/19.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the Data parameters and coefficients
*        for the reaction_diffusion example.  This includes exact solutions,
*        RHS functions, and boundary conditions, for different types of
*        test problems.
*
*/

/*! \fn void diffusion_coeff(REAL *val,REAL* x,REAL time,void *dim)
*
* \brief Diffusion coefficient: a scalar array of dim^2 x 1, A(x)
*        Example of array ordering in 3D
*
*       | a_11  a_12  a_13 |
*  A =  | a_21  a_22  a_23 |  => val = [a_11, a_12, a_13, a_21, a_22, a_23, a_31, a_32, a_33]
*       | a_31  a_32  a_33 |
*
* \param x    coordinates of where to evaluate function
* \param time time at which to evaluate function
* \param dim  void parameter which can be anything, usually dimension
*
* \return val array of coefficients of size dim^2
*
* \note Uncomment the example you want to try or input one yourself
*/
void diffusion_coeff(REAL *val,REAL* x,REAL time,void *dim) {

  // EXAMPLE 1: 2D POISSON, A=I
  val[0] = 1.0;
  val[1] = 0.0;
  val[2] = 0.0;
  val[3] = 1.0;

  // EXAMPLE 2: 2D ISOTROPIC DIFFUSION
  /* val[0] = x[0] + x[1] + 1;
  val[1] = x[0]*x[1];
  val[2] = x[0]*x[1];
  val[3] = exp(x[0]+x[1]);
  */

  // EXAMPLE 3: 3D ANISOTROPIC DIFFUSION
  /*val[0] = 100;
  val[1] = 0;
  val[2] = 0;
  val[3] = 0;
  val[4] = 1;
  val[5] = 0;
  val[6] = 0;
  val[7] = 0;
  val[8] = 1;
  */
}

/*! \fn void reaction_coeff(REAL *val,REAL* x,REAL time,void *dim)
*
* \brief Scalar coefficient
*
* \param x    coordinates of where to evaluate function
* \param time time at which to evaluate function
* \param dim  void parameter which can be anything, usually dimension
*
* \return val coefficient of size 1
*
* \note Uncomment the example you want to try or input one yourself
*/
void reaction_coeff(REAL *val,REAL* x,REAL time,void *dim) {

  // EXAMPLE 1: 2D POISSON
  *val = 0.0;

  // EXAMPLE 2: 2D ISOTROPIC DIFFUSION
  //*val = exp(x[0]+x[1]);

  // EXAMPLE 3: 3D ANISOTROPIC DIFFUSION
  //*val = 1.0;
}

/*! \fn void pde_coeff(REAL *val,REAL* x,REAL time,void *dim)
*
* \brief Combines all PDE coefficients into one array.
*        Do not change anything here.
*
* \param x    coordinates of where to evaluate function
* \param time time at which to evaluate function
* \param dim  void parameter which can be anything, usually dimension
*
* \return val coefficient of size 1
*
*/
void pde_coeff(REAL *val,REAL* x,REAL time,void *dim){

  INT d = *((INT *) dim);
  diffusion_coeff(val, x, time, &d);
  reaction_coeff(val+(d)*(d),x,time, &d);
}

// Exact Solutions (if you have one)
// Here we include several examples for different test problems or dimensions.
void exactsol_1D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0]);
}
void exactsol_2D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void exactsol_3D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}

// Derivative of Exact Solution (if you have one)
void D_exactsol_1D(REAL *val,REAL* x,REAL time,void *param) {
  *val = M_PI*cos(M_PI*x[0]);
}
void D_exactsol_2D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_exactsol_3D(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}

// Right-hand side functions.  Again we include several examples.
void rhs_1D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_1D(&myu,x,time,param);

  // 1D POISSON
  *val = (mya*M_PI*M_PI + myc)*myu;		// A=I

}
void rhs_2D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  exactsol_2D(&myu,x,time,param);

  // EXAMPLE 1: 2D POISSON
  *val = (2*M_PI*M_PI + myc)*myu; // A=I

  // EXAMPLE 2: ISOTROPIC DIFFUSION
  /*  REAL ss = sin(M_PI*x[0])*sin(M_PI*x[1]);
  REAL cc = cos(M_PI*x[0])*cos(M_PI*x[1]);
  REAL sc = sin(M_PI*x[0])*cos(M_PI*x[1]);
  REAL cs = cos(M_PI*x[0])*sin(M_PI*x[1]);
  REAL myexp = exp(x[0]+x[1]);

  *val = M_PI*M_PI*(x[0] + x[1] + 1 + myexp)*ss - 2*M_PI*M_PI*x[0]*x[1]*cc - M_PI*(myexp + x[1])*sc - M_PI*(1+x[0])*cs + myexp*ss; */

}
void rhs_3D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myc=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  exactsol_3D(&myu,x,time,param);

  // 3D POISSON
  *val = (3*M_PI*M_PI + myc)*myu; // A=I

  // EXAMPLE 3: 3D ANISOTROPIC DIFFUSION
  //*val = (102*M_PI*M_PI + myc)*myu;

}
// Boundary Conditions.  Again we include several examples.
void bc_1D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_1D(&myu,x,time,param);
  *val= myu;
}
void bc_2D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_2D(&myu,x,time,param);
  *val= myu;
}
void bc_3D(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_3D(&myu,x,time,param);
  *val= myu;
}
