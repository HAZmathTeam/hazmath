/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * 2D harmonic mapping solutions as used in
 * Error Estimators and Marking Strategies for Electrically Coupled Liquid Crystal Systems
 * 
 * The goal is to provide functions here that, for a given point, return the value of the true solution and the
 * gradients of the true solution at a provided point
*/

#include <math.h>


void solution_value(REAL *val, REAL *x, REAL k, REAL a, REAL b) {
  REAL theta = k * log10(sqrt((x[0] - a)*(x[0] - a) + (x[1] - b)*(x[1] - b)));
  val[0] = sin(theta);
  val[1] = cos(theta);
  val[2] = 0.0;
  return;
}

void n1_gradient(REAL *val, REAL *x, REAL k, REAL a, REAL b) {
  REAL px = x[0];
  REAL py = x[1];

  REAL theta = k * log10(sqrt((px-a)*(px-a) + (py-b)*(py-b)));
  REAL thetax = k*(px - a)/(log(10)*((px-a)*(px-a) + (py-b)*(py-b))));
  REAL thetay = k*(py - b)/(log(10)*((px-a)*(px-a) + (py-b)*(py-b))));

  val[0] = cos(theta)*thetax;
  val[1] = cos(theta)*thetay;
  val[2] = 0.0;
  return;
}

void n2_gradient(REAL *val, REAL *x, REAL k, REAL a, REAL b) {
  REAL px = x[0];
  REAL py = x[1];

  REAL theta = k * log10(sqrt((px-a)*(px-a) + (py-b)*(py-b)));
  REAL thetax = k*(px - a)/(log(10)*((px-a)*(px-a) + (py-b)*(py-b))));
  REAL thetay = k*(py - b)/(log(10)*((px-a)*(px-a) + (py-b)*(py-b))));

  val[0] = -sin(theta)*thetax;
  val[1] = -sin(theta)*thetay;
  val[2] = 0.0;
  return;
}
