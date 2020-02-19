/*! \file Maxwell_data.h
*
*  Created by Adler, Hu, Zikatanov on 09/29/2017.
*  Edited by Casey Cavanaugh 5/20/19
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


// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time,void *param) {

  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  //ordering (E1,E2,E3,B1,B2,B3,p)

  //Real test problem
  val[0] = -cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[2] = 0.0;
  val[3] = -sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[4] = -cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[5] = 2*cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp;
  val[6] = 0.0;

  /* val[0] = myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2]);
  val[1] = myexp*(2*x[2]*x[2] - 2*x[0]*x[1]);
  val[2] = myexp*(x[1]*x[2] - 3*x[1]*x[1]);
  val[3] = myexp*(-6*x[1] - 3*x[2]);
  val[4] = myexp*2*x[1];
  val[5] = myexp*(x[0] - 2*x[1] - 2*x[2]);
  val[6] = 0.0;
  */

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
//NEED TO INPUT -j as RHS
void current_density(REAL *val,REAL* x,REAL time,void *param) {
  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  val[0] = myexp*(3*M_PI + 1/M_PI)*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = myexp*(3*M_PI + 1/M_PI)*-sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = 0.0;

  /* val[0] = -myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2] - 2);
  val[1] = -myexp*(2*x[2]*x[2] - 2*x[0]*x[1] - 4);
  val[2] = -myexp*(x[1]*x[2] - 3*x[1]*x[1] + 6);
  */
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



