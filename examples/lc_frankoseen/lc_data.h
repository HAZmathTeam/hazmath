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
  val[1] = 1.0; //1 //1.2 //0.62903
  val[2] = 1.0; //1 //1 //1.32258
}

// Initial Guesses - There are three test cases, you need to Command+/ to
// use the test you would like.

// TEST CASE 1:
// Dirichlet boundary conditions in the y bounds and
// Periodic boundary conditions in the x bounds
void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
  REAL phi=M_PI*0.5;
  REAL th = M_PI*0.25;
  if(x[1]==0 || x[1]==1) {
    val[0] = 1.0;
    val[1] = 0.0;
    val[2] = 0.0;
    val[3] = 0.0;
  } else {
    val[0] = cos(th)*sin(phi); // n1 = sqrt(2)/2
    val[1] = sin(th)*sin(phi); // n2 = sqrt(2)/2
    val[2] = cos(phi); // n3 = 0
    val[3] = 0.0; // lamda = 0
  }
  return;
}

// // TEST CASE 2:
// // Another test case where n = (1,0,0) when y = 0 and n = (0,0,1) when y = 1.
//     // and periodic at the other points.
// void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
//   REAL phi=M_PI*0.5;
//   REAL th = M_PI*0.25;
//   if(x[1]==0) {
//     val[0] = 1.0;
//     val[1] = 0.0;
//     val[2] = 0.0;
//     val[3] = 0.0;
//   }
//   else if (x[1]==1) {
//     val[0] = 0.0;
//     val[1] = 0.0;
//     val[2] = 1.0;
//     val[3] = 0.0;
//   }
//   else {
//     val[0] = cos(th)*sin(phi); // n1 = sqrt(2)/2
//     val[1] = sin(th)*sin(phi); // n2 = sqrt(2)/2
//     val[2] = cos(phi); // n3 = 0
//     val[3] = 0.0; // lamda = 0
//   }
//   return;
// }

// // TEST CASE 3:
// // Another test case where n = (1,0,0) when y = 0 and n = (0,0,1) when y = 1.
//     // and periodic at the other points.
// void initial_guess(REAL *val,REAL *x,REAL time,void *param) {
//   REAL r   = 0.25;
//   REAL s   = 0.95;
//   REAL phi = M_PI*0.5;
//   REAL th  = M_PI*0.25;
//   if(x[1] == 0 || x[1] == 1) {
//     REAL X_m = (-s*sin(2*M_PI*(x[0]+r)))/(-s*cos(2*M_PI*(x[0]+r)) - 1);
//     REAL X_p = (-s*sin(2*M_PI*(x[0]+r)))/(-s*cos(2*M_PI*(x[0]+r)) + 1);
//     val[0] = 0.0;
//     val[1] = cos(r*(M_PI + (2*atan(X_m))-(2*atan(X_p)) ) );
//     val[2] = sin(r*(M_PI + (2*atan(X_m))-(2*atan(X_p)) ) );
//     val[3] = 0.0;
//   }
//   else {
//     val[0] = cos(th)*sin(phi); // n1 = sqrt(2)/2
//     val[1] = sin(th)*sin(phi); // n2 = sqrt(2)/2
//     val[2] = cos(phi); // n3 = 0
//     val[3] = 0.0; // lamda = 0
//   }
//   return;
// }


// Boundary Conditions for the update during each newton step
void bc(REAL *val, REAL *x, REAL time,void *param) {
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
  val[3] = 0.0;
 return;
}



// ----------------------------------------------------------------------

// DO NOT ERASE ANYTHING...PLEASE :)
// INITIAL GUESSES FOR OTHER TEST CASES ---------------------------------------

// // assigning the value of the director, n.
// val[0] = cos(th)*sin(phi); // n1 = sqrt(2)/2
// val[1] = sin(th)*sin(phi); // n2 = sqrt(2)/2
// val[2] = cos(phi); // n3 = 0
// val[3] = 0.0; // lamda = 0
// // change the values of n on the boundaries
// if(x[0]==0 || x[0]==1 || x[1]==0) {
//   val[0] = 1.0;
//   val[1] = 0.0;
//   val[2] = 0.0;
// }
// if(x[1]==1) {
//   val[0] = 1.0;
//   val[1] = 0.0;
//   val[2] = 0.0;

  // REAL phi=M_PI*0.5;
  // REAL th = M_PI*0.25;
  // val[0] = cos(th)*sin(phi);
  // val[1] = sin(th)*sin(phi);
  // val[2] = cos(phi);
  // val[3] = 0.0;
  // if(x[1]==0 || x[0]==0 || x[0]==1) {
  //   val[0] = 1.0;
  //   val[1] = 0.0;
  //   val[2] = 0.0;
  // }
  // if(x[1]==1) {
  //   val[0] = 0.0;
  //   val[1] = 0.0;
  //   val[2] = 1.0;
  // }

// different test case.
//  REAL th = 3.0*log10(sqrt((x[0]+0.1)*(x[0]+0.1) + (x[1]+0.1)*(x[1]+0.1)));
//  if(x[1]==0 || x[1]==1 || x[0]==0 || x[0]==1) {
//    val[0] = sin(th);
//    val[1] = cos(th);
//    val[2] = 0.0;
//    val[3] = 0.0;
//  } else {
//    val[0] = 0.0;
//    val[1] = 0.0;
//    val[2] = 0.0;
//    val[3] = 0.0;
//  }

// ----------------------------------------------------------------------------

//// RHS (not actually used, but the function needs to exist
//void source(REAL *val, REAL *x, REAL time) {
//  val[0] = 0.0;
//  val[1] = 0.0;
//  val[2] = 0.0;
//  val[3] = 0.0;
//  return;
//}
