/*! \file StokesData.h
 *
 *  Created by Peter Ohm on 1/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Stokes example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 */

// Exact Solutions
void exact_sol(REAL *val,REAL *x,REAL time) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];

}
void exact_sol2D(REAL *val, REAL *x, REAL time){

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];

}

// Gradients of Exact Solution
void Dexact_sol(REAL *val, REAL *x, REAL time) {

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
}
void Dexact_sol2D(REAL *val, REAL *x, REAL time){

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[4] = -1.0;
  val[5] = 0.0;

}

// RHS
void source3D(REAL *val, REAL *x, REAL time) {
  double pi = M_PI;
  val[0] = -3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  val[1] = 3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = -3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  val[3] = 0.0;
}
void source2D(REAL *val, REAL *x, REAL time) {
  double pi = M_PI;
  val[0] = 2*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]) -1.0;
  val[1] = -2*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  val[2] = 0.0;
}

// Boundary Conditions
void bc_ux(REAL *val, REAL *x, REAL time) {
  REAL myexact[4];
  exact_sol(myexact,x,time);
  *val = myexact[0];
}

void bc_uy(REAL *val, REAL *x, REAL time) {
  REAL myexact[4];
  exact_sol(myexact,x,time);
  *val = myexact[1];
}

void bc_uz(REAL *val, REAL *x, REAL time) {
  REAL myexact[4];
  exact_sol(myexact,x,time);
  *val = myexact[2];
}

void bc_p(REAL *val, REAL *x, REAL time) {
  REAL myexact[4];
  exact_sol(myexact,x,time);
  *val = myexact[3];
}

void bc2D(REAL *val, REAL *x, REAL time) {

  exact_sol2D(val,x,time);

}

void bc3D(REAL *val, REAL *x, REAL time) {

  exact_sol(val,x,time);

}

// Matlab Dump
void print_matlab_vector_field(dvector* ux, dvector* uy, dvector* uz, fespace* FE ){
 FILE *fid = fopen("output/usol_vfield.mat","w");
 INT i;
 for(i=0; i<ux->row; i++) {
  fprintf(fid,"%f\t%f\t%f\t%f\t%f\t%f\n",FE->cdof->x[i],FE->cdof->y[i],FE->cdof->z[i],ux->val[i],uy->val[i],uz->val[i]);
 }
 fclose(fid);
 return;
}
