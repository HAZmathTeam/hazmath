/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * Basic tilt-twist setup similar to that used for 
 * COMBINING DEFLATION AND NESTED ITERATION FOR COMPUTING MULTIPLE SOLUTIONS OF NONLINEAR
 * VARIATIONAL PROBLEMS
 * 
 * Periodic boundary conditions in x-axis
 * n(x, 0) = (cos -pi/4, 0, sin(-pi/4))
 * n(x, 1) = (cos pi/4, 0, sin(pi/4))
 * 
 * Initial Guess
 * n(x,y) = (cos(pi/40), sin(pi/40), 0)
 * 
 * NOTE: For K1=1.0, K2=3.0, K3=1.2, the expected Free Energy is 3.59293
*/

// Frank Constants
void get_frank_constants(REAL *val) {
  val[0] = 1.0;
  val[1] = 3.0;
  val[2] = 1.2;
}

void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL theta = M_PI*(1.0/40.0);
  REAL phi = M_PI*(1.0/4.0);

  REAL approx_zero = 0.0000000000001;
  REAL approx_one = 0.9999999999999;

  if(x[1] <= approx_zero) {
    val[0] = cos(-phi);
    val[1] = 0.0;
    val[2] = sin(-phi);
    val[3] = 0.0;
  }
  else if (x[1] >= approx_one) {
    val[0] = cos(phi);
    val[1] = 0.0;
    val[2] = sin(phi);
    val[3] = 0.0;
  }
  else {
    val[0] = cos(theta);
    val[1] = sin(theta);
    val[2] = 0.05;
    val[3] = 0.0;
  }
  return;
}

// Boundary Conditions for the update during each newton step
void bc(REAL *val, REAL *x, REAL time,void *param) {
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
  val[3] = 0.0;
 return;
}

INT characterize_boundary_conditions() {
  // return 0 for Dirichlet and 1 for Periodic
  return 1;
}
