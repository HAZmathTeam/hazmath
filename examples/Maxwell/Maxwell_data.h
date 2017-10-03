/*! \file Maxwell_data.h
 *
 *  Created by Adler, Hu, Zikatanov on 09/29/2017.
 *  Copyright 2016_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Maxwell example.  This includes exact solutions,
 *        RHS functions, coefficients, and boundary conditions.
 *
 */

// PDE Coefficients
void permitivity(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void permeability(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void oneovermu(REAL *val,REAL* x,REAL time,void *param) {
  REAL mu = -6.66;
  permeability(&mu,x,time,param);
  *val = 1.0/mu;
}
void gamm(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.05;
}
void oneplusgamm(REAL *val,REAL* x,REAL time,void *param) {
  REAL gam = -666.6;
  gamm(&gam,x,time,param);
  *val = 1.0 + gam;
}

// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time,void *param) {
  // E = exp(r*(|x|+t)) * (r^2/|x|^2-r/|x|^3) * (0,z,-y)^T
  // B = exp(r*|x|+t)) * [ (1/|x|^3) * (r^2-3r/|x|+3/|x|^2) * (z^2+y^2,-xy,-xz)^T +
  //                                                   ( (2r/|x|^2-2/|x|^3),0,0 )^T ]
  // p = 0 
  // r = 0.5*(1 - \sqrt(1 + 4/gamma))
  REAL gam = -666.6;
  gamm(&gam,x,time,param);
  REAL r = 0.5*(1-sqrt(1+4/gam));
  REAL normx = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  REAL myexp = exp(r*(normx + time));

  // ordering (E1,E2,E3,B1,B2,B3,p)
  val[0] = 0.0;
  val[1] = (myexp/(normx*normx))*(r*r - r/normx)*x[2];
  val[2] = -(myexp/(normx*normx))*(r*r - r/normx)*x[1];
  val[3] = myexp*((1/(normx*normx*normx))*(r*r-3*r/normx+3/(normx*normx))*(x[2]*x[2]+x[1]*x[1]) + (2*r/(normx*normx)) - 2/(normx*normx*normx));
  val[4] = myexp*((1/(normx*normx*normx))*(r*r-3*r/normx+3/(normx*normx))*(-x[0]*x[1]));
  val[5] = myexp*((1/(normx*normx*normx))*(r*r-3*r/normx+3/(normx*normx))*(-x[0]*x[2]));
  val[6] = 0.0;
}
void Etrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
void Btrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[3];
  val[1] = myu[4];
  val[2] = myu[5];
}
void ptrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  *val = myu[6];
}

// Right-hand Side
void current_density(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
}
void zero_vec(REAL *val,REAL* x,REAL time,void *param) {
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
}
void one_scal(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void zero_scal(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.0;
}

// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
  val[3] = myu[3];
  val[4] = myu[4];
  val[5] = myu[5];
  val[6] = myu[6];
}
void bc_test(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.0;
}
