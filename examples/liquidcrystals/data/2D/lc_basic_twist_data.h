/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the LC Elastic example.
 *
 * Started updating on 02/01/2020
 * \briefUpdate: This also has the relabeling of the boundaries
 *               for the periodic boundary cases.
 */

// Frank Constants
void get_frank_constants(REAL *val) {
  val[0] = 1.0; //1 //1 //1
  val[1] = 1.2;//420094; //1 //1.2 //0.62903
  val[2] = 1.0; //1 //1 //1.32258
}

// Initial Guesses
// n = (1,0,0) when y = 0 and n = (0,0,1) when y = 1.
// and periodic at the other points.
void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL phi=M_PI*0.5;
  REAL th = M_PI*0.25;

  REAL approx_zero = 0.0000000000001;
  REAL approx_one = 0.9999999999999;

  if(x[1] <= approx_zero) {
    val[0] = 1.0;
    val[1] = 0.0;
    val[2] = 0.0;
    val[3] = 0.0;
  }
  else if (x[1] >= approx_one) {
    val[0] = 0.0;
    val[1] = 0.0;
    val[2] = 1.0;
    val[3] = 0.0;
  }
  else {
    val[0] = cos(th)*sin(phi); // n1 = sqrt(2)/2
    val[1] = sin(th)*sin(phi); // n2 = sqrt(2)/2
    val[2] = cos(phi); // n3 = 0
    val[3] = 0.0; // lamda = 0
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
