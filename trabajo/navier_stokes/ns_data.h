/*! \file examples/navier_stokes/ns_data.h
 *
 *  Created by James Adler on 05/20/2020
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Navier-Stokes example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 * \note We include different examples for a 2D and 3D problem.
 *
 */

// Constants
// Reynolds number
void get_reynolds_number(REAL *val){
  *val = 10.0;
  return;
}

// Initial Guess
void initial_guess3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];
  return;
}
void initial_guess2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = 0.0;//sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = 0.0;//-cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];
  return;
}

// Exact Solutions
void exact_sol3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];
  return;
}
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];
  return;
}

// Gradients of Exact Solution
void Dexact_sol3D(REAL *val, REAL *x, REAL time,void *param) {

  val[0] = -M_PI*cos(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = -M_PI*sin(M_PI*x[0])*cos(M_PI*(x[1]-x[2]));
  val[2] = M_PI*sin(M_PI*x[0])*cos(M_PI*(x[1]-x[2]));

  val[3] = M_PI*sin(M_PI*x[1])*cos(M_PI*(x[0]-x[2]));
  val[4] = M_PI*cos(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[5] = -M_PI*sin(M_PI*x[1])*cos(M_PI*(x[0]-x[2]));

  val[6] = -M_PI*sin(M_PI*x[2])*cos(M_PI*(x[0]-x[1]));
  val[7] = M_PI*sin(M_PI*x[2])*cos(M_PI*(x[0]-x[1]));
  val[8] = -M_PI*cos(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));

  val[9] = -1.0;
  val[10] = 0.0;
  val[11] = 0.0;
  return;
}
void Dexact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[4] = -1.0;
  val[5] = 0.0;
  return;
}

// RHS
void source3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL Re = 0.0;
  get_reynolds_number(&Re);
  val[0] = -3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0 +
            Re*pi*sin(pi*x[0])*(cos(pi*x[0]) - 0.5*cos(pi*(x[0]-2*x[1])) - 0.5*cos(pi*(x[0]-2*x[2])));
  val[1] = 3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2])) + Re*pi*sin(pi*x[1])*(
            cos(pi*x[1]) - 0.5*cos(pi*(2*x[0]-x[1])) - 0.5*cos(pi*(x[1]-2*x[2])));
  val[2] = -3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1])) + Re*pi*sin(pi*x[2])*(cos(pi*x[2]) - cos(pi*(x[0]-x[1]))*cos(pi*(x[0]+x[1]-x[2])));
  val[3] = 0.0;
  return;
}
void source2D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL Re = 0.0;
  get_reynolds_number(&Re);
  val[0] = 2*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]) + Re*pi*cos(pi*x[0])*sin(pi*x[0]) - 1.0;
  val[1] = -2*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]) + Re*pi*cos(pi*x[1])*sin(pi*x[1]);
  val[2] = 0.0;
  return;
}

// Boundary Conditions
void bc2D(REAL *val, REAL *x, REAL time,void *param) {

  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
  return;
}
void bc3D(REAL *val, REAL *x, REAL time,void *param) {

  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
  val[3] = 0.0;
  return;
}
