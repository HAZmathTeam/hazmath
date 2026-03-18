/*! \file ns_data.h
 *
 *  Created by James Adler on 01/30/2024
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the Data parameters and coefficients
 *        for the Navier-Stokes example.  This includes exact solutions,
 *        RHS functions, and boundary conditions.
 *
 * \note We include different examples for a 2D and 3D problem.
 *
 */

// Reynolds number (Re)
void get_reynolds_number(REAL *val){
  REAL reynolds = 1.0;
  *val = reynolds;
  return;
}

// Exact Solutions
void exact_sol3D(REAL *val,REAL *x,REAL time,void *param) {

  val[0] = -sin(M_PI*x[0])*sin(M_PI*(x[1]-x[2]));
  val[1] = sin(M_PI*x[1])*sin(M_PI*(x[0]-x[2]));
  val[2] = -sin(M_PI*x[2])*sin(M_PI*(x[0]-x[1]));
  val[3] = 0.5 - x[0];
  return;
}
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.5 - x[0];
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

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);

  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);

  val[4] = -1.0;
  val[5] = 0.0;
  return;
}

// RHS
void sourcerhs(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  REAL Re=0.0;
  get_reynolds_number(&Re);

  // 2D
  val[0] = pi * sin(pi*x[0]) * (Re*cos(pi*x[0]) + 2.0*pi*cos(pi*x[1])) -1.0;
  val[1] = pi * sin(pi*x[1]) * (Re*cos(pi*x[1]) - 2.0*pi*cos(pi*x[0]));
  val[2] = 0.0;
  
  // 3D
  // val[0] = -3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  // val[1] = 3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  // val[2] = -3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  // val[3] = 0.0;
  return;
}

// Initial Conditions
void initial_guess(REAL *val, REAL *x, REAL time,void *param) {

  INT* flag = (INT*) param; // The FE projection functions assume that param is the DoF flag
  
  // Set to zero unless it's the Dirichlet boundary
  REAL sol_exact[3]; // Exact solution
  // 2D
  exact_sol2D(sol_exact,x,0.0,NULL);
  // 3D
  // exact_sol3D(sol_exact,x,0.0,NULL);

  if(flag[0]==0) { // Interior is zero
    val[0] = 0.0;
    val[1] = 0.0;
    val[2] = 0.0;
    // val[3] = 0.0; // 3D
  } else { // boundary 
    val[0] = sol_exact[0];
    val[1] = sol_exact[1];
    val[2] = sol_exact[2];
    // val[3] = sol_exact[3]; // 3D
  }

  return;
}

// Boundary Conditions (for nonlinear so this is BC for update and should be zero)
void bc(REAL *val, REAL *x, REAL time,void *param) {

  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;
  // val[3] = 0.0; // 3D
  return;
}

// Mean value of p for mean value constraint
void pmean(REAL* val) {
  *val = 0.0;
}

// These routines will relable the mesh and the boundary flags for your FEM spaces if you like to have different conditions on different boundaries
// Relabel boundary flags we'll have the following labels:
// 1: x=0  2: x=1   3: y=0    4: y=1  5: z = 0  6: z = 1 
void relabel_boundary(fespace* FE, INT dim) {

  INT i,j;
  INT ndof = FE->ndof;
  for(i=0;i<ndof;i++) {
    if(FE->dof_flag[i]!=0) {
      for(j=0;j<dim;j++) {
        if(FE->cdof->x[j*ndof+i]==0.0) {
          FE->dof_flag[i] = 2*j+1;
        }
        if(FE->cdof->x[j*ndof+i]==1.0) {
          FE->dof_flag[i] = 2*j+2;
        }
      }
    }
  }
  return;
}

// To be safe we should relabel all the mesh flags as well.
// 1: x=0  2: x=1   3: y=0    4: y=1    5: z = 0  6: z = 1  
void relabel_mesh(mesh_struct* mesh) {

  INT i,j;

  INT dim = mesh->dim;

    // Vertices
    for(i=0;i<mesh->nv;i++) {
      if(mesh->v_flag[i]!=0) {
        for(j=0;j<dim;j++) {
          if(mesh->cv->x[j*mesh->nv+i]==0.0) {
            mesh->v_flag[i] = 2*j+1;
          }
          if(mesh->cv->x[j*mesh->nv+i]==1.0) {
            mesh->v_flag[i] = 2*j+2;
          }
        }
      }
    }
    // Edges
    for(i=0;i<mesh->nedge;i++) {
      if(mesh->ed_flag[i]!=0) {
        for(j=0;j<dim;j++) {
          if(mesh->ed_mid[i*dim+j]==0.0) {
            mesh->ed_flag[i] = 2*j+1;
          }
          if(mesh->ed_mid[i*dim+j]==1.0) {
            mesh->ed_flag[i] = 2*j+2;
          }
        }
      }
    }
    // Faces
    for(i=0;i<mesh->nface;i++) {
      if(mesh->f_flag[i]!=0) {
        for(j=0;j<dim;j++) {
          if(mesh->f_mid[i*dim+j]==0.0) {
            mesh->f_flag[i] = 2*j+1;
          }
          if(mesh->f_mid[i*dim+j]==1.0) {
            mesh->f_flag[i] = 2*j+2;
          }
        }
      }
    }
  return;
}
