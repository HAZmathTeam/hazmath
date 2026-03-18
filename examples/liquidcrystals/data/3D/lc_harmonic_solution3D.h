/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * 3D harmonic mapping solutions as used in
 * Error Estimators and Marking Strategies for Electrically Coupled Liquid Crystal Systems
 * 
 * Dirichlet boundary conditions throughout the domain
 * 
 * See the paper draft for a clean description of the solution.
 * 
 * Initial Guess
 * Is a random perturbation of the analytical solution with a magnitude set in the code below
 * 
 * NOTE: For the choice of harmonic solution in the paper, the Free energy should be 8.847
*/

#include <math.h>

// Frank Constants
void get_frank_constants(REAL *val) {
  val[0] = 1.0;
  val[1] = 1.0;
  val[2] = 1.0;
}

void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL magnitude = 0.1;
  REAL perturbation_1 = (((double) rand())/((double) RAND_MAX) - 0.5) * 2.0 * magnitude;
  REAL perturbation_2 = (((double) rand())/((double) RAND_MAX) - 0.5) * 2.0 * magnitude;
  REAL perturbation_3 = (((double) rand())/((double) RAND_MAX) - 0.5) * 2.0 * magnitude;

  // Cohen, Kinderlehrer Generalized Harmonic
  
  REAL x_ = x[0] + 0.2;
  REAL y_ = x[1] + 0.1;
  REAL z_ = x[2];
  
  // unit length domain projection
  REAL p1 = x_/(sqrt(x_*x_ + y_*y_ + z_*z_));
  REAL p2 = y_/(sqrt(x_*x_ + y_*y_ + z_*z_));
  REAL p3 = z_/(sqrt(x_*x_ + y_*y_ + z_*z_));

  // forward stereographic projection
  REAL pi1 = p1/(1.0 - p3);
  REAL pi2 = p2/(1.0 - p3);

  // Holomorphic function
  // REAL r1 = pi2*pi2 - pi1*pi1;
  // REAL r2 = 2*pi1*pi2;
  REAL r1 = pi1/(pi1*pi1 + pi2*pi2) + pi1*pi1 - pi2*pi2;
  REAL r2 = 2*pi1*pi2 - pi2/(pi1*pi1 + pi2*pi2);

  REAL approx_zero = 0.0000000000001;
  REAL approx_one = 0.9999999999999;

  if(x[0] <= approx_zero || x[1] <= approx_zero || x[2] <= approx_zero || 
                    x[0] >= approx_one || x[1] >= approx_one || x[2] >= approx_one) {
    val[0] = (2.0*r1)/(1.0 + r1*r1 + r2*r2);
    val[1] = (2.0*r2)/(1.0 + r1*r1 + r2*r2);
    val[2] = (r1*r1 + r2*r2 - 1.0)/(1.0 + r1*r1 + r2*r2);
    val[3] = 0.0;
  }
  else {
    REAL n1_p = (2.0*r1)/(1.0 + r1*r1 + r2*r2) + perturbation_1;
    REAL n2_p = (2.0*r2)/(1.0 + r1*r1 + r2*r2) + perturbation_2;
    REAL n3_p = (r1*r1 + r2*r2 - 1.0)/(1.0 + r1*r1 + r2*r2) + perturbation_3;
    REAL length = sqrt(n1_p*n1_p + n2_p*n2_p + n3_p*n3_p);
    val[0] = n1_p/length;
    val[1] = n2_p/length;
    val[2] = n3_p/length;
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
