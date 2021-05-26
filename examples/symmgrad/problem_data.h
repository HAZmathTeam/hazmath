/*! \file examples/stokes/stokes_data.h
 *
 *  Created by Peter Ohm on 1/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Stokes example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 * \note We include different examples for a 2D and 3D problem.
 *
 */

// Exact Solutions
void exact_sol3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//-sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//-sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  return;
}
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  return;
}

// Gradients of Exact Solution
void Dexact_sol3D(REAL *val, REAL *x, REAL time,void *param) {

  val[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//-M_PI*cos(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = M_PI*cos(M_PI*x[1])*sin(M_PI*x[0])*sin(M_PI*x[2]);//-M_PI*sin(M_PI*x[0])*cos(M_PI*(x[1]-x[2]));
  val[2] = M_PI*cos(M_PI*x[2])*sin(M_PI*x[0])*sin(M_PI*x[1]);//M_PI*sin(M_PI*x[0])*cos(M_PI*(x[1]-x[2]));

  val[3] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//M_PI*sin(M_PI*x[1])*cos(M_PI*(x[0]-x[2]));
  val[4] = M_PI*cos(M_PI*x[1])*sin(M_PI*x[0])*sin(M_PI*x[2]);//M_PI*cos(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[5] = M_PI*cos(M_PI*x[2])*sin(M_PI*x[0])*sin(M_PI*x[1]);//-M_PI*sin(M_PI*x[1])*cos(M_PI*(x[0]-x[2]));

  val[6] = M_PI*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);//-M_PI*sin(M_PI*x[2])*cos(M_PI*(x[0]-x[1]));
  val[7] = M_PI*cos(M_PI*x[1])*sin(M_PI*x[0])*sin(M_PI*x[2]);//M_PI*sin(M_PI*x[2])*cos(M_PI*(x[0]-x[1]));
  val[8] = M_PI*cos(M_PI*x[2])*sin(M_PI*x[0])*sin(M_PI*x[1]);//-M_PI*cos(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));

  return;
}
void Dexact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  return;
}

// RHS
void source3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  double pi2 = pi*pi;
  REAL sinx = sin(pi*x[0]);
  REAL cosx = cos(pi*x[0]);
  REAL siny = sin(pi*x[1]);
  REAL cosy = cos(pi*x[1]);
  REAL sinz = sin(pi*x[2]);
  REAL cosz = cos(pi*x[2]);
  val[0] = pi2*(4.0*sinx*siny*sinz - cosx*sin(pi*(x[1]+x[2])));//-3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2]));
  val[1] = pi2*(4.0*sinx*siny*sinz - cosy*sin(pi*(x[0]+x[2])));//3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = pi2*(4.0*sinx*siny*sinz - cosz*sin(pi*(x[0]+x[1])));//-3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  return;
}
void source2D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  val[0] = 2*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]);
  val[1] = -2*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  val[2] = 0.0;
  return;
}

// Boundary Conditions
void bc2D(REAL *val, REAL *x, REAL time,void *param) {

  exact_sol2D(val,x,time,param);
  return;
}
void bc3D(REAL *val, REAL *x, REAL time,void *param) {

  exact_sol3D(val,x,time,param);
  return;
}
