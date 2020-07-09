/*! \file examples/elasticity/elasticity_data.h
 *
 *  Created by James Adler on 07/04/2020
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Elasticity example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 * \note We include different examples for a 2D and 3D problems.
 *
 */

// Parameters
 // Youngs Modulus
 void get_young(REAL *val,REAL* x,REAL time,void *param) {
   *val = 1.0;
 }

 // Poisson Ratio
 void get_nu(REAL *val,REAL* x,REAL time,void *param) {
   *val = 0.499;
 }

 // Lame Coefficients
 void get_lam(REAL *val,REAL* x,REAL time,void *param) {
   REAL E;
   REAL nu;
   get_nu(&nu,x,time,param);
   get_young(&E,x,time,param);
   *val = E*nu/((1-2*nu)*(1+nu));
 }
 void  get_mu(REAL *val,REAL* x,REAL time,void *param) {
   REAL E;
   REAL nu;
   get_nu(&nu,x,time,param);
   get_young(&E,x,time,param);
   *val = E/(2+2*nu);
 }

// Exact Solutions
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];
  return;
}
void exact_sol3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];
  return;
}
// Gradients of Exact Solution
void Dexact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[4] = -1.0;
  val[5] = 0.0;
  return;
}
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

// Boundary Conditions
void bc2D(REAL *val, REAL *x, REAL time,void *param) {
  exact_sol2D(val,x,time,param);
  return;
}
void bc3D(REAL *val, REAL *x, REAL time,void *param) {
  exact_sol3D(val,x,time,param);
  return;
}

// RHS
void source2D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL mu=0.0;
  get_mu(&mu,x,time,param);
  REAL lam=0.0;
  get_lam(&lam,x,time,param);
  val[0] = 2*mu*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]) -1.0;
  val[1] = -2*mu*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  val[2] = -(1.0/lam)*(0.5 - x[0]);
  return;
}
void source3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL mu=0.0;
  get_mu(&mu,x,time,param);
  REAL lam=0.0;
  get_lam(&lam,x,time,param);
  val[0] = -3*mu*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  val[1] = 3*mu*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = -3*mu*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  val[3] =  -(1.0/lam)*(0.5 - x[0]);
  return;
}

// Just grab stuff for u for primal formulation
void uexact_sol2D(REAL *val, REAL *x, REAL time,void *param){
  REAL tempval[3];
  exact_sol2D(tempval,x,time,param);
  val[0] = tempval[0];
  val[1] = tempval[1];
  return;
}
void uexact_sol3D(REAL *val,REAL *x,REAL time,void *param) {
  REAL tempval[4];
  exact_sol3D(tempval,x,time,param);
  val[0] = tempval[0];
  val[1] = tempval[1];
  val[2] = tempval[2];
  return;
}
// Gradients of Exact Solution
void uDexact_sol2D(REAL *val, REAL *x, REAL time,void *param){
  REAL tempval[6];
  Dexact_sol2D(tempval,x,time,param);
  val[0] = tempval[0];
  val[1] = tempval[1];
  val[2] = tempval[2];
  val[3] = tempval[3];
  return;
}
void uDexact_sol3D(REAL *val, REAL *x, REAL time,void *param) {
  REAL tempval[12];
  Dexact_sol3D(tempval,x,time,param);
  val[0] = tempval[0];
  val[1] = tempval[1];
  val[2] = tempval[2];
  val[3] = tempval[3];
  val[4] = tempval[4];
  val[5] = tempval[5];
  val[6] = tempval[6];
  val[7] = tempval[7];
  val[8] = tempval[8];
  return;
}

// Boundary Conditions
void ubc2D(REAL *val, REAL *x, REAL time,void *param) {
  uexact_sol2D(val,x,time,param);
  return;
}
void ubc3D(REAL *val, REAL *x, REAL time,void *param) {
  uexact_sol3D(val,x,time,param);
  return;
}

// RHS
void usource2D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL mu=0.0;
  get_mu(&mu,x,time,param);
  REAL lam=0.0;
  get_lam(&lam,x,time,param);
  val[0] = 2*mu*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]);
  val[1] = -2*mu*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  return;
}
void usource3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL mu=0.0;
  get_mu(&mu,x,time,param);
  REAL lam=0.0;
  get_lam(&lam,x,time,param);
  val[0] = -3*mu*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  val[1] = 3*mu*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = -3*mu*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  return;
}


// Relabel boundary flags:
// 1: x=0  2: x=1   3: y=0    4: y=1
void relabel_boundary2D(fespace* FE) {

  INT i;
  for(i=0;i<FE->ndof;i++) {
    if(FE->dof_flag[i]==1) {
      if(FE->cdof->x[i]==0.0) {
        FE->dof_flag[i] = 1;
      }
      if(FE->cdof->x[i]==1.0) {
        FE->dof_flag[i] = 2;
      }
      if(FE->cdof->y[i]==0.0) {
        FE->dof_flag[i] = 3;
      }
      if(FE->cdof->y[i]==1.0) {
          FE->dof_flag[i] = 4;
      }
    }
  }
  return;
}

// 1: x=0  2: x=1   3: y=0    4: y=1
void relabel_boundary3D(fespace* FE) {

  INT i;
  for(i=0;i<FE->ndof;i++) {
    if(FE->dof_flag[i]==1) {
      if(FE->cdof->x[i]==0.0) {
        FE->dof_flag[i] = 1;
      }
      if(FE->cdof->x[i]==1.0) {
        FE->dof_flag[i] = 2;
      }
      if(FE->cdof->y[i]==0.0) {
        FE->dof_flag[i] = 3;
      }
      if(FE->cdof->y[i]==1.0) {
          FE->dof_flag[i] = 4;
      }
      if(FE->cdof->z[i]==0.0) {
        FE->dof_flag[i] = 5;
      }
      if(FE->cdof->z[i]==1.0) {
          FE->dof_flag[i] = 6;
      }
    }
  }
  return;
}
