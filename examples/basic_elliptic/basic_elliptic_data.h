/*! \file examples/basic_elliptic/basic_elliptic_data.h
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

// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time,void *param) {
  // a(x)
  *val = 1.0;
}
void reaction_coeff(REAL *val,REAL* x,REAL time,void *param) {
  // c(x)
  *val = 1.0;
}

// Exact Solutions (if you have one)
// We have different ones for different dimensions and different D's
void exactsol_1D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 1D - grad grad
  *val = sin(M_PI*x[0]);
}
void exactsol_2D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void exactsol_3D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad grad
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}
void exactsol_2D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = -sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void exactsol_3D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - curl curl
  val[0] = cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = -sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}
void exactsol_2D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1]);
}
void exactsol_3D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad div
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
  val[1] = cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
  val[2] = cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
}

// Derivatives of Exact Solutions (if you have one)
// We have different ones for different dimensions and different D's
void D_exactsol_1D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 1D - grad grad
  *val = M_PI*cos(M_PI*x[0]);
}
void D_exactsol_2D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad grad
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_exactsol_3D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad grad
  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
}
void D_exactsol_2D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - curl curl
  *val = -2*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_exactsol_3D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - curl curl
  val[0] = -2*M_PI*sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
  val[1] = 2*M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]);
  val[2] = 0;
}
void D_exactsol_2D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad div
  *val = 2*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
}
void D_exactsol_3D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad div
  *val = 3*M_PI*cos(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]);
}

// Right-hand Sides
// We have different ones for different dimensions and different D's
void rhs_1D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 1D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_1D_PX(&myu,x,time,param);
  *val = (mya*M_PI*M_PI + myc)*myu;
}
void rhs_2D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_2D_PX(&myu,x,time,param);
  *val = (mya*2*M_PI*M_PI + myc)*myu;
}
void rhs_3D_PX(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad grad
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu=-666.6;
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_3D_PX(&myu,x,time,param);
  *val = (mya*3*M_PI*M_PI + myc)*myu;
}
void rhs_2D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_2D_Ned(myu,x,time,param);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
void rhs_3D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - curl curl
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_3D_Ned(myu,x,time,param);
  val[0] = (2*mya*M_PI*M_PI + myc)*myu[0];
  val[1] = (2*mya*M_PI*M_PI + myc)*myu[1];
  val[2] = (4*mya*M_PI*M_PI + myc)*myu[2];
}
void rhs_2D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 2D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[2];
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_2D_RT(myu,x,time,param);
  val[0] = (mya*2.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*2.0*M_PI*M_PI + myc)*myu[1];
}
void rhs_3D_RT(REAL *val,REAL* x,REAL time,void *param) {
  // 3D - grad div
  REAL myc=-666.6;
  REAL mya=-666.6;
  REAL myu[3];
  reaction_coeff(&myc,x,time,param);
  diffusion_coeff(&mya,x,time,param);
  exactsol_3D_RT(myu,x,time,param);
  val[0] = (mya*3.0*M_PI*M_PI + myc)*myu[0];
  val[1] = (mya*3.0*M_PI*M_PI + myc)*myu[1];
  val[2] = (mya*3.0*M_PI*M_PI + myc)*myu[2];
}

// Boundary Conditions
// We have different ones for different dimensions and different D's
void bc_1D_PX(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_1D_PX(&myu,x,time,param);
  *val= myu;
}
void bc_2D_PX(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_2D_PX(&myu,x,time,param);
  *val= myu;
}
void bc_3D_PX(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu;
  exactsol_3D_PX(&myu,x,time,param);
  *val= myu;
}
void bc_2D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[2];
  exactsol_2D_Ned(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
}
void bc_3D_Ned(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[3];
  exactsol_3D_Ned(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
void bc_2D_RT(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[2];
  exactsol_2D_RT(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
}
void bc_3D_RT(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[3];
  exactsol_3D_RT(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
