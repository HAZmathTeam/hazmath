/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * 2D nano-patterned solutions as used in
 * ENERGY-MINIMIZATION FINITE-ELEMENT APPROACH FOR THE FRANK–OSEEN MODEL
 * OF NEMATIC LIQUID CRYSTALS∗
 * 
 * Periodic boundary condtions
 * 
 * The initial guess is
 * n1(x, y) = sqrt(2)/2
 * n2(x, y) = sqrt(2)/2
 * 
 * Boundary conditions
 * Nano patterning. See below
*/

#include <math.h>

// Frank Constants
void get_frank_constants(REAL *val) {
  val[0] = 1.0;
  val[1] = 0.62903;
  val[2] = 1.32258;
}

void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL r = 0.25;
  REAL s = 0.95;
  
  REAL X_m = (-s * sin(2*M_PI*(x[0] + r)))/(-s * cos(2*M_PI*(x[0] + r)) - 1);
  REAL X_p = (-s * sin(2*M_PI*(x[0] + r)))/(-s * cos(2*M_PI*(x[0] + r)) + 1);

  REAL approx_zero = 0.0000000000001;
  REAL approx_one = 0.9999999999999;

  if(x[1] <= approx_zero || x[1] >= approx_one) {
    val[0] = 0.0;
    val[1] = cos(r*(M_PI + 2*atan(X_m) - 2*atan(X_p)));
    val[2] = sin(r*(M_PI + 2*atan(X_m) - 2*atan(X_p)));
    val[3] = 0.0;
  }
  else {
    val[0] = 0.0;
    val[1] = cos(r*(M_PI + 2*atan(X_m) - 2*atan(X_p)));
    val[2] = sin(r*(M_PI + 2*atan(X_m) - 2*atan(X_p)));
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
