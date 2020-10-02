/*! \file examples/stokes/stokes_data.h
 *
 *  Created by Peter Ohm on 1/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Stokes example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 * \note We include different examples for a 2D and 3D problem.
 *
 */

// Exact Solutions
void exact_sol3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];
  return;
}
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){
  /*
  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];
  */

  double xx = x[0];
  double y = x[1];
  double lame_lambda = 1.;//000000.;//000000.;//000000.;//000000.;
  val[0] = (xx*xx)*y*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*2.0+(xx*xx)*(y*y)*(y*2.0-2.0)*pow(xx-1.0,2.0); //sin(x[0]) * sin(x[1]) + (1./lame_lambda) * x[0];
  val[1] = xx*(y*y)*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*-2.0-(xx*xx)*(y*y)*(xx*2.0-2.0)*pow(y-1.0,2.0);//cos(x[0]) * cos(x[1]) + (1./lame_lambda) * x[1];
  val[2] = 0.;

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
  /*
  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[4] = -1.0;
  val[5] = 0.0;
  */

  double xx = x[0];
  double y = x[1];

  double lame_lambda = 1.;//000000.;//000000.;//000000.;//000000.;
  
  val[0] = xx*(y*y)*(y*2.0-2.0)*pow(xx-1.0,2.0)*2.0+(xx*xx)*y*(xx*2.0-2.0)*pow(y-1.0,2.0)*2.0+(xx*xx)*(y*y)*(xx*2.0-2.0)*(y*2.0-2.0)+xx*y*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*4.0;//2*x*y^2*(2*y - 2)*(x - 1)^2 + 2*x^2*y*(2*x - 2)*(y - 1)^2 + x^2*y^2*(2*x - 2)*(2*y - 2) + 4*x*y*(x - 1)^2*(y - 1)^2; //
  //cos(x[0])*sin(x[1]) + 1./lame_lambda;
  val[1] =(xx*xx)*(y*y)*pow(xx-1.0,2.0)*2.0+(xx*xx)*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*2.0+(xx*xx)*y*(y*2.0-2.0)*pow(xx-1.0,2.0)*4.0;// 2*x^2*y^2*(x - 1)^2 + 2*x^2*(x - 1)^2*(y - 1)^2 + 4*x^2*y*(2*y - 2)*(x - 1)^2;//cos(x[1])*sin(x[0]);
  val[2] = (xx*xx)*(y*y)*pow(y-1.0,2.0)*-2.0-(y*y)*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*2.0-xx*(y*y)*(xx*2.0-2.0)*pow(y-1.0,2.0)*4.0;;//- 2*x^2*y^2*(y - 1)^2 - 2*y^2*(x - 1)^2*(y - 1)^2 - 4*x*y^2*(2*x - 2)*(y - 1)^2; //-cos(x[1])*sin(x[0]);
  val[3] = xx*(y*y)*(y*2.0-2.0)*pow(xx-1.0,2.0)*-2.0-(xx*xx)*y*(xx*2.0-2.0)*pow(y-1.0,2.0)*2.0-(xx*xx)*(y*y)*(xx*2.0-2.0)*(y*2.0-2.0)-xx*y*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*4.0;//- cos(x[0])*sin(x[1]) + 1./lame_lambda;
  val[4] = 0.;
  val[5] = 0.;

  
  return;
}

// RHS
void source3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  val[0] = -3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  val[1] = 3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = -3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  val[3] = 0.0;
  return;
}
void source2D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;

  //val[0] = 2*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]);
  //val[1] = -2*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  val[2] = 0.0;
  double xx = x[0];
  double y = x[1];
  val[0] = (xx*xx)*(y*y)*(y*2.0-2.0)*-2.0-y*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*4.0-(xx*xx)*(y*2.0-2.0)*pow(xx-1.0,2.0)*6.0-(y*y)*(y*2.0-2.0)*pow(xx-1.0,2.0)*2.0-(xx*xx)*y*pow(xx-1.0,2.0)*1.2E1-(xx*xx)*y*pow(y-1.0,2.0)*4.0-xx*(y*y)*(xx*2.0-2.0)*(y*2.0-2.0)*4.0-xx*y*(xx*2.0-2.0)*pow(y-1.0,2.0)*8.0;//2. * sin(x[0]) * sin(x[1]);
    val[1] = (xx*xx)*(y*y)*(xx*2.0-2.0)*2.0+xx*pow(xx-1.0,2.0)*pow(y-1.0,2.0)*4.0+(xx*xx)*(xx*2.0-2.0)*pow(y-1.0,2.0)*2.0+(y*y)*(xx*2.0-2.0)*pow(y-1.0,2.0)*6.0+xx*(y*y)*pow(xx-1.0,2.0)*4.0+xx*(y*y)*pow(y-1.0,2.0)*1.2E1+(xx*xx)*y*(xx*2.0-2.0)*(y*2.0-2.0)*4.0+xx*y*(y*2.0-2.0)*pow(xx-1.0,2.0)*8.0;//2. * cos(x[0]) * cos(x[1]);

  return;
}

// Boundary Conditions
void bc2D(REAL *val, REAL *x, REAL time,void *param) {

  exact_sol2D(val,x,time,param);
  return;
}
void bc3D(REAL *val, REAL *x, REAL time,void *param) {

  exact_sol3D(val,x,time,param);
  return;
}
