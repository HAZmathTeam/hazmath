/*! \file examples/heat_equation/heat_equation_data.h
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the Data parameters and coefficients
*        for the heat_equation example.  This includes exact solutions,
*        RHS functions, initial conditions, and boundary conditions.
*
*/

// PDE Coefficients
void diffusion_coeff(REAL *val,REAL* x,REAL time,void *param) {
  // a(x)
  *val = 1.0;
}

// Exact Solution (if you have one)
// Change as needed for different dimensions
void exactsol2D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*exp(-2.0*M_PI*M_PI*time);
}
void exactsol3D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2])*exp(-3*M_PI*M_PI*time);
}

// Right-hand Side
void myrhs(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.0;
}

// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time,void *param) {
  *val= 0.0;
}

// Initial Conditions
// Change as needed for different dimensions
void initial_conditions2D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1]);
}
void initial_conditions3D(REAL *val,REAL* x,REAL time,void *param) {
  *val = sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}
