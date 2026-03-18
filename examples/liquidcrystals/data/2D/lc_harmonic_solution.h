/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * 2D harmonic mapping solutions as used in
 * Error Estimators and Marking Strategies for Electrically Coupled Liquid Crystal Systems
 * 
 * Dirichlet boundary conditions throughout the domain
 * 
 * True Solution
 * n(x,y) = (sin(theta(x,y)), cos(theta(x,y)), 0)
 * where theta(x,y) = k * log_10(Sqrt((x-a)^2 + (y-b)^2))
 * 
 * Boundary Conditions Conform to the above solution
 * 
 * Initial Guess
 * Is a random perturbation of the analytical theta between - magnitude and + magnitude
 * at each value (x, y)
 * 
 * NOTE: For k=-4.5, a=0.5, and b=-0.1, the Free energy should be 8.7174
*/

#include <math.h>

// Frank Constants
void get_frank_constants(REAL *val) {
  val[0] = 1.0;
  val[1] = 1.0;
  val[2] = 1.0;
}

void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL k = -4.5;
  REAL a = 0.5;
  REAL b = -0.1;

  REAL theta = k * log10(sqrt((x[0] - a)*(x[0] - a) + (x[1] - b)*(x[1] - b)));

  REAL magnitude = 0.1;
  REAL perturbation = (((REAL) rand())/((REAL) RAND_MAX) - 0.5) * 2.0 * magnitude;

  REAL approx_zero = 0.0000000000001;
  REAL approx_one = 0.9999999999999;

  if(x[0] <= approx_zero || x[1] <= approx_zero || x[0] >= approx_one || x[1] >= approx_one) {
    val[0] = sin(theta);
    val[1] = cos(theta);
    val[2] = 0.0;
    val[3] = 0.0;
  }
  else {
    val[0] = sin(theta + perturbation);
    val[1] = cos(theta + perturbation);
    val[2] = 0.0;
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
  return 0;
}
