//#include "elasticity_data.h"
/*! \file examples/stokes/elasticity_error.h
 *
 *  Created by SLee (20200812)
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
  val[3] = 0.0; // fake slot for eg
  val[4] = 0.5 - x[0];
  return;
}
void exact_sol2D(REAL *val, REAL *x, REAL time,void *param){

  double lame_lambda = LAME_LAMBDA_GLOBAL;

  val[0] = sin(x[0]) * sin(x[1]);
  val[1] = cos(x[0]) * cos(x[1]);
  val[2] = 0.;
  val[3] = sin(M_PI*x[0])* cos(M_PI*x[1]);

  val[0] = sin(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -cos(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = 0.;
  val[3] = 0.5 - x[0];

  return;
}


void exact_sol2D_dt(REAL *val, REAL *x, REAL time,void *param){
  double pi = M_PI;

  double lame_mu = LAME_MU_GLOBAL;
  double lame_lambda = LAME_LAMBDA_GLOBAL;


  val[0] = 0.;// -exp(-time)*(sin(2*pi*x[1])*(cos(2*pi*x[0]) - 2) + (sin(pi*x[0])*sin(pi*x[1]))/(lame_mu + lame_lambda));
  val[1] = 0.;//exp(-time)*(sin(2*pi*x[0])*(cos(2*pi*x[1]) - 2) - (sin(pi*x[0])*sin(pi*x[1]))/(lame_mu + lame_lambda));

  /*
  val[0] = x[0]/lame_lambda + x[0]*cos(time*x[0])*sin(time*x[1]) + x[1]*cos(time*x[1])*sin(time*x[0]);
  val[1] = x[1]/lame_lambda - x[0]*cos(time*x[1])*sin(time*x[0]) - x[1]*cos(time*x[0])*sin(time*x[1]);
  */
  //printf("val[0[ = %f, \n", val[0]);

  return;
}

void exact_sol3D_dt(REAL *val, REAL *x, REAL time,void *param){
  double pi = M_PI;

  double lame_mu = LAME_MU_GLOBAL;
  double lame_lambda = LAME_LAMBDA_GLOBAL;


  val[0] = 0.;// -exp(-time)*(sin(2*pi*x[1])*(cos(2*pi*x[0]) - 2) + (sin(pi*x[0])*sin(pi*x[1]))/(lame_mu + lame_lambda));
  val[1] = 0.;//exp(-time)*(sin(2*pi*x[0])*(cos(2*pi*x[1]) - 2) - (sin(pi*x[0])*sin(pi*x[1]))/(lame_mu + lame_lambda));
  val[2] = 0.;
  /*
  val[0] = x[0]/lame_lambda + x[0]*cos(time*x[0])*sin(time*x[1]) + x[1]*cos(time*x[1])*sin(time*x[0]);
  val[1] = x[1]/lame_lambda - x[0]*cos(time*x[1])*sin(time*x[0]) - x[1]*cos(time*x[0])*sin(time*x[1]);
  */
  //printf("val[0[ = %f, \n", val[0]);

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

  // fake slots for eg
  val[9] = 0.0;
  val[10] = 0.0;
  val[11] = 0.0;

  val[12] = -1.0;
  val[13] = 0.0;
  val[14] = 0.0;
  return;
}
void Dexact_sol2D(REAL *val, REAL *x, REAL time,void *param){


  double lame_mu = 1.;
  double lame_lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.;//000000.;//000000.;//000000.;
  double pi = M_PI;

  val[0] = cos(x[0])*sin(x[1]);
  val[1] = cos(x[1])*sin(x[0]);
  val[2] = -cos(x[1])*sin(x[0]);
  val[3] = -cos(x[0])*sin(x[1]);
  val[4] = 0.;
  val[5] = 0.;
  val[6] = pi*cos(pi*x[1])*cos(pi*x[0]);//M_PI * cos(M_PI*x[0]*time)*cos(M_PI*x[1]*time);
  val[7] = -pi*sin(pi*x[0])*sin(pi*x[1]);////-M_PI * sin(M_PI*x[0]*time)*sin(M_PI*x[1]*time);

  val[0] = M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[1] = -M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*x[1]);
  val[3] = -M_PI*cos(M_PI*x[0])*cos(M_PI*x[1]);
  val[4] = 0.0;
  val[5] = 0.0;
  val[6] = -1.0;
  val[7] = 0.0;

  return;
}

// RHS
void source3D(REAL *val, REAL *x, REAL time,void *param) {
  double pi = M_PI;
  val[0] = -3*pow(pi,2)*sin(pi*x[0])*sin(pi*(x[1]-x[2])) - 1.0;
  val[1] = 3*pow(pi,2)*sin(pi*x[1])*sin(pi*(x[0]-x[2]));
  val[2] = -3*pow(pi,2)*sin(pi*x[2])*sin(pi*(x[0]-x[1]));
  val[3] = 0.0; // fake slot for eg
  val[4] = 0.0;
  return;
}
void source2D(REAL *val, REAL *x, REAL time,void *param) {

  REAL pi = M_PI;
  REAL xx = x[0];
  REAL yy = x[1];
  REAL alpha = 1.;
  //printf(" ********************************** \n");
  //printf(" TIME IN Source 2d = %f \n", time);
  //printf(" ********************************** \n");
  double lame_mu = 1.;
  double lame_lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.;//000000.;//000000.;//000000.;

  val[0] = 2*lame_mu*sin(xx)*sin(yy) + pi*cos(pi*xx)*cos(pi*yy);
  val[1] = 2*lame_mu*cos(xx)*cos(yy) - pi*sin(pi*xx)*sin(pi*yy);
  val[2] = 0.;
  val[3] = 0.;//sin(M_PI*x[0])* cos(M_PI*x[1]);

  val[0] = 2*pow(pi,2) * sin(pi*x[0]) * cos(pi*x[1]) -1.0;
  val[1] = -2*pow(pi,2) * cos(pi*x[0]) * sin(pi*x[1]);
  val[2] = 0.0;
  val[3] = 0.0;
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


void zquad_face(qcoordinates *cqbdry,INT nq1d, INT dim, REAL *xf, REAL farea)
{
  // Flag for errors
  SHORT status;

  INT q,j; /* Loop indices */
  INT nq = 0;

  if(dim==2) { // face is edge in 2D
    nq = nq1d;
  } else if(dim==3) {
    nq = nq1d*nq1d;
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  /* Coordinates of vertices of face are in xf ordered by rows (every
     vertex corresponds to a row)*/

  // Gaussian points for reference element
  // (edges: [-1,1] faces: tri[(0,0),(1,0),(0,1)])
  REAL* gp = (REAL *) calloc((dim-1)*nq,sizeof(REAL));
  /* Gaussian weights for reference element */
  REAL* gw = (REAL *) calloc(nq,sizeof(REAL));

  REAL r,s;      	/* Points on Reference Face */

  REAL w = 0.0;
  if(dim==2) { // Faces are Edges in 2D
    w = 0.5*farea; /* Jacobian = 1/2 |e| */
  } else if(dim==3) {
    w = 2.*farea; /* Jacobian = 2*|f| */
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Get coordinates of vertices for given edge/face
  //  get_incidence_row(dof,dof_v,thisdof_v);
  // Get Quad Nodes and Weights
  if(dim==2) { // face is an edge
    quad1d(gp,gw,nq1d);
  } else if(dim==3) { // face is a face
    triquad_(gp,gw,nq1d);
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Map to Real Face
  if(dim==2) {
    // Edges: x = 0.5(x1*(1-r) + x2*(1+r))
    //        y = 0.5(y1*(1-r) + y2*(1+r))
    //        z = 0.5(z1*(1-r) + z2*(1+r))
    //        w = 0.5*Edge Length*wref
    for (q=0; q<nq1d; q++) {
      r = gp[q];
      /* cqbdry->x[q] = 0.5*cvdof->x[0]*(1-r) + 0.5*cvdof->x[1]*(1+r); */
      /* cqbdry->y[q] = 0.5*cvdof->y[0]*(1-r) + 0.5*cvdof->y[1]*(1+r); */
      cqbdry->x[q] =				\
	0.5*(xf[0*dim+0]*(1-r) +		\
	     xf[1*dim+0]*(1+r));
      cqbdry->y[q] =				\
	0.5*(xf[0*dim+1]*(1-r) +			\
	     xf[1*dim+1]*(1+r));
      cqbdry->w[q] = w*gw[q];
    }
  } else if(dim==3) {
    // Faces: x = x1*(1-r-s) + x2*r + x3*s
    //        y = y1*(1-r-s) + y2*r + y3*s
    //        z = z1*(1-r-s) + z2*r + z3*s
    //        w = 2*Element Area*wref
    for (q=0; q<nq1d*nq1d; q++) {
      r = gp[q];
      s = gp[nq1d*nq1d+q];
      /* cqbdry->x[q] = cvdof->x[0]*(1-r-s) + cvdof->x[1]*r + cvdof->x[2]*s; */
      /* cqbdry->y[q] = cvdof->y[0]*(1-r-s) + cvdof->y[1]*r + cvdof->y[2]*s; */
      /* cqbdry->z[q] = cvdof->z[0]*(1-r-s) + cvdof->z[1]*r + cvdof->z[2]*s; */
      cqbdry->x[q] =	      \
	xf[0*dim+0]*(1-r-s) + \
	xf[1*dim+0]*r +       \
	xf[2*dim+0]*s;
      cqbdry->y[q] =	      \
	xf[0*dim+1]*(1-r-s) + \
	xf[1*dim+1]*r +       \
	xf[2*dim+1]*s;
      cqbdry->z[q] =	      \
	xf[0*dim+2]*(1-r-s) + \
	xf[1*dim+2]*r +       \
	xf[2*dim+2]*s;
      cqbdry->w[q] = w*gw[q];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(gp) free(gp);
  if(gw) free(gw);
  return;
}
/*********************************************************************************/

void L2error_block_EG
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
{
  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  coordinates *barycenter = allocatecoords(1,dim);

  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i]=0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20)
    //Scalar Element
      ncomp[i]=1;
    else // Vector Element
      ncomp[i] = dim;
    nun += ncomp[i];
  }


  printf("nun = %d, dof_per_elm = %d \n", nun, dof_per_elm);
  //exit(0);

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  //REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  INT *local_dof_on_elm;
  INT dof;
  REAL u0_value_at_q = 0.;
  REAL u1_value_at_q = 0.;
  REAL u2_value_at_q = 0.;
  REAL eg_vals[dim];
  REAL p_value_at_q = 0.;
  REAL* u_comp;
  REAL eg_xterm,eg_yterm,eg_zterm;

  // Loop over all Elements
  for (elm=0; elm<mesh->nelm; elm++) {

    barycenter->x[0] = mesh->el_mid[elm*dim];
    barycenter->y[0] = mesh->el_mid[elm*dim + 1];
    if(dim==3) barycenter->z[0] = mesh->el_mid[elm*dim + 2];


    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
	//printf("i = %d, jcntr = %d, dof_on_elm = %d \n", i,  jcntr,   dof_on_elm[jcntr]);
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      eg_xterm = qx[0] - barycenter->x[0];
      eg_yterm = qx[1] - barycenter->y[0];
      if(dim==3) eg_zterm = qx[2] - barycenter->z[0];

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));


      //for(INT ii = 0; ii < nun; ++ii)
      //printf("val_true[%d] = %f \n", ii, val_true[ii]);
      //exit(0);

      local_dof_on_elm = dof_on_elm;
      u0_value_at_q = 0.;
      u1_value_at_q = 0.;
      u2_value_at_q = 0.;
      p_value_at_q = 0.;

      //REAL u2_value_at_q = 0.;

      u_comp = u;

      // TODO: Rewriting this with some hazmath FUNCTIONS
      // u0
      FE_Interpolation(&u0_value_at_q,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
      get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);
      local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
      u_comp += FE->var_spaces[0]->ndof;

      // u1
      FE_Interpolation(&u1_value_at_q,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      u_comp += FE->var_spaces[1]->ndof;

      // u2
      if(dim==3) {
        FE_Interpolation(&u2_value_at_q,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
        get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
        local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
        u_comp += FE->var_spaces[2]->ndof;
      }

      // eg
      get_FEM_basis(FE->var_spaces[dim]->phi,FE->var_spaces[dim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim]);
      // Interpolate
      eg_vals[0] = 0.0;
      eg_vals[1] = 0.0;
      if(dim==3) eg_vals[2] = 0.0;
      dof = local_dof_on_elm[0]; // get the dof for the last component
      //printf("--u2 dof = %d, u _comp= %f, q[0] = %f, q[1] = %f,   \n", dof, u_comp[dof],qx[0] - barycenter->x[0],qx[1] - barycenter->x[1]);
      eg_vals[0] = u_comp[dof]*eg_xterm;//   (qx[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      eg_vals[1] = u_comp[dof]*eg_yterm;//   (qx[1] - barycenter->y[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      if(dim==3) eg_vals[2] = u_comp[dof]*eg_zterm;
      //printf("--3) u2_value_at_q= %f  \n", val[0]);
      //printf("--3) u2_value_at_q= %f  \n", val[1]);

      //printf("***u2  dof = %d, val[0]= %f,  val[1] = %f, x = %f, y = %f  \n", dof,  val[0],  val[1], (qx[0] - barycenter->x[0]) ,  (qx[1] - barycenter->y[0]));
      // ************************************************* //
      // u4 /////////////////////////////////////////////////
      local_dof_on_elm += FE->var_spaces[dim]->dof_per_elm;
      u_comp += FE->var_spaces[dim]->ndof;

      // p
      FE_Interpolation(&p_value_at_q,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[dim+1],mesh);
      get_FEM_basis(FE->var_spaces[dim+1]->phi,FE->var_spaces[dim+1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim+1]);
      local_dof_on_elm += FE->var_spaces[dim+1]->dof_per_elm;
      u_comp += FE->var_spaces[dim+1]->ndof;

      // Compute Square of Error on Element for each component of FE space
      err[0] += w*(ABS(u0_value_at_q  +eg_vals[0]- val_true[0]))*(ABS(u0_value_at_q +eg_vals[0] - val_true[0]));
      err[1] += w*(ABS(u1_value_at_q  +eg_vals[1]- val_true[1]))*(ABS(u1_value_at_q +eg_vals[1] - val_true[1]));
      if(dim==3) err[2] += w*(ABS(u2_value_at_q +eg_vals[2]- val_true[2]))*(ABS(u2_value_at_q +eg_vals[2] - val_true[2]));


      //err[0] += w*(ABS(val[0]- val_true[0]))*(ABS(val[0] - val_true[0]));
      //err[1] += w*(ABS(val[1]- val_true[1]))*(ABS(val[1] - val_true[1]));


      //err[0] += w*(ABS(u0_value_at_q  - val_true[0]))*(ABS(u0_value_at_q  - val_true[0]));
      //err[1] += w*(ABS(u1_value_at_q  - val_true[1]))*(ABS(u1_value_at_q  - val_true[1]));

      err[dim] = 0;

      // For Pressure Error
      // bEG_Pressure
      err[dim+1] += w*(ABS(p_value_at_q    -  val_true[dim+1]))*(ABS(p_value_at_q - val_true[dim+1]));
      //err[3] += w*(ABS(p_value_at_q     -  val_true[3]))*(ABS(p_value_at_q  - val_true[3]));

      //exit(0);

    }//quad
  }

  //Include Pressure
  for(i=0;i<dim;i++) {
    err[i] = sqrt(err[i]);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  //if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);
  if(barycenter) {
    free_coords(barycenter);
    free(barycenter);
    barycenter=NULL;
  }

  return;
}

void HDerror_block_EG
//(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)

{

  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;


  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i] = 0.0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) /* Scalar Gradient */
      ncomp[i]=dim;
    else if(FE->var_spaces[i]->FEtype==20 && dim==2) /* Curl in 2D is scalar */
      ncomp[i] = 1;
    else if(FE->var_spaces[i]->FEtype==20 && dim==3) /* Curl in 3D is vector */
      ncomp[i] = dim;
    else /* Div is scalar */
      ncomp[i] = 1;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  //sanity checks //slee
  //printf("nspaces = %d |  ncomp[0],ncomp[1],ncomp[2]= %d, %d, %d, %d,  %d  \n", nspaces, ncomp[0],ncomp[1],ncomp[2],ncomp[3],ncomp[4]);
  //exit(0);

  // THIS NEEDS TO BE CODED BETTER // slee todo
  //double d_Lame_coeff_mu = 1.;
  //double d_Lame_coeff_lambda = LAME_LAMBDA_GLOBAL ; //000000.;

  INT *local_dof_on_elm;
  INT dof,k;
  REAL gradu0[dim];
  REAL gradu1[dim];
  REAL gradu2[dim];
  REAL gradeg0[dim];
  REAL gradeg1[dim];
  REAL gradeg2[dim];
  REAL gradp0[dim];
  REAL gradp1[dim];
  REAL gradp2[dim];
  REAL* u_comp;

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      //blockFE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      local_dof_on_elm = dof_on_elm;
      u_comp = u;

      // TODO: Again rewriting with hazmath functions

      // u0
      FE_DerivativeInterpolation(gradu0,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
      get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);
      local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
      u_comp += FE->var_spaces[0]->ndof;

      // u1
      FE_DerivativeInterpolation(gradu1,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      u_comp += FE->var_spaces[1]->ndof;

      // u2
      if(dim==3) {
        FE_DerivativeInterpolation(gradu2,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
        get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
        local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
        u_comp += FE->var_spaces[2]->ndof;
      }

      // eg
      get_FEM_basis(FE->var_spaces[dim]->phi,FE->var_spaces[dim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim]);
      // Interpolate
      gradeg0[0] = 0.0;
      gradeg0[1] = 0.0;
      if(dim==3) gradeg0[2] = 0.0;
      dof = local_dof_on_elm[0]; // get the dof for the last component
      gradeg0[0] += u_comp[dof]*1.;
      gradeg1[0] = 0.0;
      gradeg1[1] = 0.0;
      if(dim==3) gradeg1[2] = 0.0;
      gradeg1[1] += u_comp[dof]*1.;
      if(dim==3) {
        gradeg2[0] = 0.0;
        gradeg2[1] = 0.0;
        gradeg2[2] = 0.0;
        gradeg2[2] += u_comp[dof]*1.;
      }


      //val[0] = u_comp[dof]*   (qx[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      //val[1] = u_comp[dof]*   (qx[1] - barycenter->y[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point


      local_dof_on_elm += FE->var_spaces[dim]->dof_per_elm;
      u_comp += FE->var_spaces[dim]->ndof;

      // p
      // TODO: Not sure why there are so many components for p?
      gradp0[0]=0.0;
      gradp0[1]=0.0;
      gradp1[0]=0.0;
      gradp1[1]=0.0;
      if(dim==3) {
        gradp0[2] = 0.0;
        gradp1[2] = 0.0;
        gradp2[0] = 0.0;
        gradp2[1] = 0.0;
        gradp2[2] = 0.0;
      }

      // p
      get_FEM_basis(FE->var_spaces[dim+1]->phi,FE->var_spaces[dim+1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim+1]);
      // Interpolate
      dof = local_dof_on_elm[0];
      //grad_val4[0] += u_comp[dof]*FE->var_spaces[3]->dphi[0];
      //grad_val4[1] += u_comp[dof]*FE->var_spaces[3]->dphi[1];


      local_dof_on_elm += FE->var_spaces[dim+1]->dof_per_elm;
      u_comp += FE->var_spaces[dim+1]->ndof;

      /*
      err[0]+=w*(ABS(grad_val0[0] - val_true[0]))*(ABS(grad_val0[0] - val_true[0]));
      err[0]+=w*(ABS(grad_val1[0] - val_true[2]))*(ABS(grad_val1[0] - val_true[2]));

      err[1]+=w*(ABS(grad_val0[1] - val_true[1]))*(ABS(grad_val0[1] - val_true[1]));
      err[1]+=w*(ABS(grad_val1[1] - val_true[3]))*(ABS(grad_val1[1] - val_true[3]));
      */
      // bEG
      err[0]+=w*(ABS(gradu0[0]+gradeg0[0] - val_true[0]))*(ABS(gradu0[0]+gradeg0[0] - val_true[0]));
      err[0]+=w*(ABS(gradu1[0]+gradeg1[0] - val_true[dim]))*(ABS(gradu1[0]+gradeg1[0] - val_true[dim]));
      if(dim==3) err[0]+=w*(ABS(gradu2[0]+gradeg2[0] - val_true[2*dim]))*(ABS(gradu2[0]+gradeg2[0] - val_true[2*dim]));

      // bEG
      err[1]+=w*(ABS(gradu0[1] + gradeg0[1] - val_true[1]))*(ABS(gradu0[1] + gradeg0[1] - val_true[1]));
      err[1]+=w*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]))*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]));
      if(dim==3) err[1]+=w*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]))*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]));

      if(dim==3) {
        err[2]+=w*(ABS(gradu0[2] + gradeg0[2] - val_true[2]))*(ABS(gradu0[2] + gradeg0[2] - val_true[2]));
        err[2]+=w*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]))*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]));
        err[2]+=w*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]))*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]));
      }

      // bEG PRESSURE
      //err[3]+=w*(ABS(grad_val4[0]  - val_true[6]))*(ABS(grad_val4[0]  - val_true[6]));
      //err[3]+=w*(ABS(grad_val4[1]  - val_true[7]))*(ABS(grad_val4[1]  - val_true[7]));

      err[dim+1]+=w*(ABS(gradp0[0]  - val_true[(dim+1)*dim]))*(ABS(gradp0[0]  - val_true[(dim+1)*dim]));
      err[dim+1]+=w*(ABS(gradp0[1]  - val_true[(dim+1)*dim+1]))*(ABS(gradp0[1]  - val_true[(dim+1)*dim+1]));
      if(dim==3) err[dim+1]+=w*(ABS(gradp0[2] - val_true[(dim+1)*dim+2]))*(ABS(gradp0[2] - val_true[(dim+1)*dim+2]));
    }
    /*
      jcntr=0;
      for(i=0;i<nspaces;i++) {
        for(j=0;j<ncomp[i];j++) {
          err[i]+=w*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]));
        }
        jcntr+=ncomp[i];
      }
    }
      */
  }


  for(i=0;i<dim;i++) {
    err[i] = sqrt(err[i]);
  }

  //exit(0);


  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);

  return;
}

void HDsemierror_block_Stress
//(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)

{

  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;


  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i] = 0.0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) /* Scalar Gradient */
      ncomp[i]=dim;
    else if(FE->var_spaces[i]->FEtype==20 && dim==2) /* Curl in 2D is scalar */
      ncomp[i] = 1;
    else if(FE->var_spaces[i]->FEtype==20 && dim==3) /* Curl in 3D is vector */
      ncomp[i] = dim;
    else /* Div is scalar */
      ncomp[i] = 1;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  //sanity checks //slee
  //printf("nspaces = %d |  ncomp[0],ncomp[1],ncomp[2]= %d, %d, %d  \n", nspaces, ncomp[0],ncomp[1],ncomp[2]);
  //exit(0);

  // THIS NEEDS TO BE CODED BETTER // slee todo
  double d_Lame_coeff_mu = 1.;
  double d_Lame_coeff_lambda = LAME_LAMBDA_GLOBAL ;//000000.;//000000.;//000000.; //000000.;

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      blockFE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      // Compute Square of Error on Element for each component of FE space
      jcntr=0;
      for(i=0;i<dim;i++) {
        for(j=0;j<ncomp[i];j++) {
          // \mu * grad u
          err[i]+=d_Lame_coeff_mu*w*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]));
        }
        jcntr+=ncomp[i];
      }
      // if(dim == 3) // TODO: adding this now
      // {
        //printf("ERROR IN STRESS ERROR \n"); exit(0);
        err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
        err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[dim] - val_true[dim]))*(ABS(val_sol[dim] - val_true[dim]));
        if(dim==3) err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[2*dim] - val_true[2*dim]))*(ABS(val_sol[2*dim] - val_true[2*dim]));

        err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[1] - val_true[1]))*(ABS(val_sol[1] - val_true[1]));
        err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[dim+1] - val_true[dim+1]))*(ABS(val_sol[dim+1] - val_true[dim+1]));
        if(dim==3) err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[2*dim+1] - val_true[2*dim+1]))*(ABS(val_sol[2*dim+1] - val_true[2*dim+1]));

        if(dim==3) {
          err[2] += d_Lame_coeff_mu * w*(ABS(val_sol[2] - val_true[2]))*(ABS(val_sol[2]-val_true[2]));
          err[2] += d_Lame_coeff_mu * w*(ABS(val_sol[dim+2]-val_true[dim+2]))*(ABS(val_sol[dim+2]-val_true[dim+2]));
          err[2] += d_Lame_coeff_mu * w*(ABS(val_sol[2*dim+2]-val_true[2*dim+2]))*(ABS(val_sol[2*dim+2]-val_true[2*dim+2]));
        }
        // div u
        err[0] += d_Lame_coeff_lambda * w * (ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
        err[1] += d_Lame_coeff_lambda * w * (ABS(val_sol[dim+1] - val_true[dim+1]))*(ABS(val_sol[dim+1] - val_true[dim+1]));
        if(dim==3) err[2] += d_Lame_coeff_lambda * w * (ABS(val_sol[2*dim+2] - val_true[2*dim+2]))*(ABS(val_sol[2*dim+2]-val_true[2*dim+2]));
  //     }
  //     else if(dim == 2){
  //
	// // grad u^T
  //
	// err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
	// err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[2] - val_true[2]))*(ABS(val_sol[2] - val_true[2]));
  //
	// err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[1] - val_true[1]))*(ABS(val_sol[1] - val_true[1]));
	// err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[3] - val_true[3]))*(ABS(val_sol[3] - val_true[3]));
  //
	// // div u
	// err[0] += d_Lame_coeff_lambda * w * (ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
	// err[1] += d_Lame_coeff_lambda * w * (ABS(val_sol[3] - val_true[3]))*(ABS(val_sol[3] - val_true[3]));
  //     }
  //


    }
  }

  for(i=0;i<dim;i++) {
    err[i] = sqrt(err[i]);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);

  return;
}



void HDsemierror_block_Stress_EG
//(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)

{

  // Loop Indices
  INT i,elm,quad,j,rowa,rowb,jcntr;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;


  INT dof_per_elm = 0;
  INT* ncomp = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT nun=0;
  for(i=0;i<FE->nspaces;i++) {
    err[i] = 0.0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) /* Scalar Gradient */
      ncomp[i]=dim;
    else if(FE->var_spaces[i]->FEtype==20 && dim==2) /* Curl in 2D is scalar */
      ncomp[i] = 1;
    else if(FE->var_spaces[i]->FEtype==20 && dim==3) /* Curl in 3D is vector */
      ncomp[i] = dim;
    else /* Div is scalar */
      ncomp[i] = 1;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  //sanity checks //slee
  //printf("nspaces = %d |  ncomp[0],ncomp[1],ncomp[2]= %d, %d, %d  \n", nspaces, ncomp[0],ncomp[1],ncomp[2]);
  //exit(0);

  // THIS NEEDS TO BE CODED BETTER // slee todo
  double d_Lame_coeff_mu = 1.;
  double d_Lame_coeff_lambda = LAME_LAMBDA_GLOBAL ;///000000.;//000000.;//000000.; //000000.;

  INT *local_dof_on_elm;
  INT dof,k;
  REAL gradu0[dim];
  REAL gradu1[dim];
  REAL gradu2[dim];
  REAL gradeg0[dim];
  REAL gradeg1[dim];
  REAL gradeg2[dim];
  REAL* u_comp;

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }
    // Find vertices for given Element
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Get True Solution at Quadrature Nodes
      (*D_truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      //blockFE_DerivativeInterpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

      local_dof_on_elm = dof_on_elm;
      u_comp = u;

      // TODO: Again rewriting with hazmath functions

      // u0
      FE_DerivativeInterpolation(gradu0,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
      get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);
      local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
      u_comp += FE->var_spaces[0]->ndof;

      // u1
      FE_DerivativeInterpolation(gradu1,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      u_comp += FE->var_spaces[1]->ndof;

      // u2
      if(dim==3) {
        FE_DerivativeInterpolation(gradu2,u_comp,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
        get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
        local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
        u_comp += FE->var_spaces[2]->ndof;
      }

      // eg
      get_FEM_basis(FE->var_spaces[dim]->phi,FE->var_spaces[dim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim]);
      // Interpolate
      gradeg0[0] = 0.0;
      gradeg0[1] = 0.0;
      if(dim==3) gradeg0[2] = 0.0;
      dof = local_dof_on_elm[0]; // get the dof for the last component
      gradeg0[0] += u_comp[dof]*1.;
      gradeg1[0] = 0.0;
      gradeg1[1] = 0.0;
      if(dim==3) gradeg1[2] = 0.0;
      gradeg1[1] += u_comp[dof]*1.;
      if(dim==3) {
        gradeg2[0] = 0.0;
        gradeg2[1] = 0.0;
        gradeg2[2] = 0.0;
        gradeg2[2] += u_comp[dof]*1.;
      }

      //val_sol 0 -> grad_val0 [0]
      //val_sol 1 -> grad_val0 [1]
      //val_sol 2 -> grad_val1 [0]
      //val_sol 3 -> grad_val1 [1]


      // grad u
      err[0] += d_Lame_coeff_mu * w*(ABS(gradu0[0] + gradeg0[0]- val_true[0]))*(ABS(gradu0[0] + gradeg0[0] - val_true[0]));
      err[0] += d_Lame_coeff_mu * w*(ABS(gradu0[1] + gradeg0[1] - val_true[1]))*(ABS(gradu0[1] + gradeg0[1] - val_true[1]));
      if(dim==3) err[0] += d_Lame_coeff_mu * w*(ABS(gradu0[2] + gradeg0[2] - val_true[2]))*(ABS(gradu0[2] + gradeg0[2] - val_true[2]));

      err[1] += d_Lame_coeff_mu * w*(ABS(gradu1[0] + gradeg1[0] - val_true[dim]))*(ABS(gradu1[0] + gradeg1[0] - val_true[dim]));
      err[1] += d_Lame_coeff_mu * w*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]))*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]));
      if(dim==3) err[1] += d_Lame_coeff_mu * w*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]))*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]));

      if(dim==3) {
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu2[0] + gradeg2[0] - val_true[2*dim]))*(ABS(gradu2[0] + gradeg2[0] - val_true[2*dim]));
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]))*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]));
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]))*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]));
      }


      // grad u^T
      err[0] += d_Lame_coeff_mu * w*(ABS(gradu0[0] + gradeg0[0] - val_true[0]))*(ABS(gradu0[0] + gradeg0[0] - val_true[0]));
      err[0] += d_Lame_coeff_mu * w*(ABS(gradu1[0] + gradeg1[0] - val_true[dim]))*(ABS(gradu1[0] + gradeg1[0] - val_true[dim]));
      if(dim==3) err[0] += d_Lame_coeff_mu * w*(ABS(gradu2[0] + gradeg2[0] - val_true[2*dim]))*(ABS(gradu2[0] + gradeg2[0] - val_true[2*dim]));

      err[1] += d_Lame_coeff_mu * w*(ABS(gradu0[1] + gradeg0[1] - val_true[1]))*(ABS(gradu0[1] + gradeg0[1] - val_true[1]));
      err[1] += d_Lame_coeff_mu * w*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]))*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]));
      if(dim==3) err[1] += d_Lame_coeff_mu * w*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]))*(ABS(gradu2[1] + gradeg2[1] - val_true[2*dim+1]));

      if(dim==3) {
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu0[2] + gradeg0[2] - val_true[2]))*(ABS(gradu0[2] + gradeg0[2] - val_true[2]));
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]))*(ABS(gradu1[2] + gradeg1[2] - val_true[dim+2]));
        err[2] += d_Lame_coeff_mu * w*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]))*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]));
      }

      // div u
      err[0] += d_Lame_coeff_lambda * w * (ABS(gradu0[0] + gradeg0[0] - val_true[0]))*(ABS(gradu0[0] + gradeg0[0] - val_true[0]));
      err[1] += d_Lame_coeff_lambda * w * (ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]))*(ABS(gradu1[1] + gradeg1[1] - val_true[dim+1]));
      if(dim==3) err[2] += d_Lame_coeff_lambda * w * (ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]))*(ABS(gradu2[2] + gradeg2[2] - val_true[2*dim+2]));

    }//quad loop
  }//element loop

  // if(nspaces ==3)
  //exit(0);

  for(i=0;i<dim;i++) {
    err[i] = sqrt(err[i]);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(val_true) free(val_true);
  if(val_sol) free(val_sol);
  if(ncomp) free(ncomp);

  return;
}





//(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)

void HDsemierror_block_EnergyNorm_EG_FaceLoop(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time){
  // Loop Indices
  INT i,j,quad,dof,rowa,rowb;
  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));


  INT dof_per_elm = 0;

  for(i=0;i<FE->nspaces;i++) {
    err[i] = 0.0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));


  // FEM Stuff
  INT nspaces = FE->nspaces;

  //sanity checks //slee
  //printf("nspaces = %d |  ncomp[0],ncomp[1],ncomp[2]= %d, %d, %d  \n", nspaces, ncomp[0],ncomp[1],ncomp[2]);
  //exit(0);

  // THIS NEEDS TO BE CODED BETTER // slee todo
  double d_Lame_coeff_mu = 1.;
  double d_Lame_coeff_lambda = LAME_LAMBDA_GLOBAL ;//000000.;//000000.;//000000.; //000000.;
  /*
    INT dof_per_face = 0;
    INT FEtype;
    for(i=0;i<FE->nspaces;i++) {
    FEtype = FE->var_spaces[i]->FEtype;
    if(FEtype>=1 && FEtype<10) { // PX Elements
    dof_per_face += dim + (FEtype-1)*(2*dim-3);
    }
    else if (FEtype==0) { // P0
    // Questionable handling of P0 for the face integral
    dof_per_face += 1;
    } else {
    printf("Block face integration isn't set up for the FEM space you chose\n");
    exit(0);
    }
    }
  */

  INT nq1d_face=3;
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face.
  REAL *xfi=data_face; //coords of vertices on face i.
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);

  INT jk,k,face,quad_face, jcntr,ed, jkl;
  INT maxdim=dim;
  REAL qx_face[maxdim];
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;
  iCSRmat *f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);

  INT* local_dof_on_elm_face = NULL;
  INT *local_dof_on_elm_face_interface = NULL;
  INT *local_dof_on_elm_neighbor = NULL;
  REAL w_face;

  REAL* u_comp_neighbor=NULL;
  INT dof_neighbor;
  REAL u0_value_at_q_neighbor;
  REAL u1_value_at_q_neighbor;
  REAL u2_value_at_q_neighbor;
  REAL p_value_at_q_neighbor;

  // penalty
  double penalty_term;
  //
  REAL* u_comp_face;
  INT dof_face;
  REAL u0_value_at_q;
  REAL u1_value_at_q;
  REAL u2_value_at_q;
  REAL p_value_at_q;

  REAL* u_comp;
  //
  REAL val[dim];
  REAL val_face[dim];
  REAL val_neighbor[dim];
  /////////////////////////////////////////////////

  INT iu0 = 0;
  INT iu1 = 1;
  INT iu2 = 2;
  INT ieg = dim;
  INT ip = dim+1;
  INT u0dofpelm = FE->var_spaces[iu0]->dof_per_elm;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dofpelm;
  if(dim==3) u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  INT egdofpelm = FE->var_spaces[ieg]->dof_per_elm;
  INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;


  // Now, get the basis for the neighbor
  REAL* neighbor_basis_u0_phi  = (REAL *) calloc(u0dofpelm,sizeof(REAL));
  REAL* neighbor_basis_u0_dphi = (REAL *) calloc(u0dofpelm*mesh->dim,sizeof(REAL));

  REAL* neighbor_basis_u1_phi  = (REAL *) calloc(u1dofpelm,sizeof(REAL));
  REAL* neighbor_basis_u1_dphi = (REAL *) calloc(u1dofpelm*mesh->dim,sizeof(REAL));

  REAL* neighbor_basis_u2_phi;
  REAL* neighbor_basis_u2_dphi;
  if(dim==3) {
    neighbor_basis_u2_phi = (REAL *) calloc(u2dofpelm,sizeof(REAL));
    neighbor_basis_u2_dphi = (REAL *) calloc(u2dofpelm*dim,sizeof(REAL));
  }

  REAL* neighbor_basis_p_phi  = (REAL *) calloc(pdofpelm,sizeof(REAL));
  REAL* neighbor_basis_p_dphi = (REAL *) calloc(pdofpelm*mesh->dim,sizeof(REAL));
  // Get the BD values
  fprintf(stdout,"\nNQUAD=%d\n",nquad_face);fflush(stdout);
  REAL* val_true_face = (REAL *) calloc(mesh->dim+2,sizeof(REAL));
  REAL* val_true_face_neighbor = (REAL *) calloc(mesh->dim+2,sizeof(REAL));

  INT rowa_neighbor,rowb_neighbor,jcntr_neighbor;
  INT k_neighbor, j_neighbor;
  double fiarea;
  coordinates *barycenter = allocatecoords(1,dim);
  coordinates *barycenter_neighbor = allocatecoords(1,dim);
  //printf("===== FACE = %d  ===== \n",face);
  INT neighbor_index[2];
  INT counter = 0,pq,nbr0;
  /* Loop over all Faces */
  for (face=0; face<mesh->nface; face++) {

    neighbor_index[0] = -1;
    neighbor_index[1] = -1;
    counter = 0;

    for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){

      //printf("-- pq = %d, f_el->IA[face] = %d, f_el->IA[face+1] = %d \n", pq, f_el->IA[face], f_el->IA[face+1]);
      nbr0=f_el->JA[pq];
      //printf("-- nbr0 = %d  \n", nbr0);

      neighbor_index[counter] = nbr0;
      counter++;

    }

  //   for(pq=0;pq<2;++pq)
  //     {
	// if(counter == 2)
	//   //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
  //
	//   if(counter == 1){
	//     if(pq == 0)
	//       {
	// 	//printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
	//       }
	//     else{
	//       //printf("-\n");
	//     }
	//   }
  //     }

    // Saniyy Check
    // counter == 1 means the BD
    if(counter == 1 && mesh->f_flag[face] == 0)
      {printf("check-error FALSE\n"); exit(0);}


    if(neighbor_index[0] == -1)
      {
	printf("UHIH \n");
	exit(0);
      }

    // Find DOF on element
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[neighbor_index[0]];
      rowb = FE->var_spaces[i]->el_dof->IA[neighbor_index[0]+1];
      for (j=rowa; j<rowb; j++) {
	dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
	jcntr++;
      }
    }

    // Find vertices for given Element
    get_incidence_row(neighbor_index[0],mesh->el_v,v_on_elm);

    // Find barycenter for given element
    barycenter->x[0] = mesh->el_mid[neighbor_index[0]*dim];
    barycenter->y[0] = mesh->el_mid[neighbor_index[0]*dim + 1];
    if(dim==3) barycenter->z[0] = mesh->el_mid[neighbor_index[0]*dim+2];


    // Find DOF on the neighbor Element
    if(counter == 2)// in the interface
      {
	if(neighbor_index[1] == -1)
	  {
	    printf("!#!@$#@!$ \n");
	    exit(0);
	  }

	jcntr_neighbor = 0;
	for(k_neighbor=0;k_neighbor<FE->nspaces;k_neighbor++) {
	  rowa_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[neighbor_index[1]];
	  rowb_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[neighbor_index[1]+1];
	  for (j_neighbor=rowa_neighbor; j_neighbor<rowb_neighbor; j_neighbor++) {
	    dof_on_elm_neighbor[jcntr_neighbor] = FE->var_spaces[k_neighbor]->el_dof->JA[j_neighbor];
	    jcntr_neighbor++;
	  }
	}
	// Find vertices for given Element
	get_incidence_row(neighbor_index[1],mesh->el_v,v_on_elm_neighbor);

	// Find barycenter for given element
	barycenter_neighbor->x[0] = mesh->el_mid[neighbor_index[1]*dim];
	barycenter_neighbor->y[0] = mesh->el_mid[neighbor_index[1]*dim + 1];
  if(dim==3) barycenter_neighbor->z[0] = mesh->el_mid[neighbor_index[1]*dim + 2];
      }
    /*
      for(INT zz=0; zz<v_per_elm; ++zz)
      {
      printf("V_on_elm[%d] = %d  ||  dof_on_elm[%d] =%d \n", zz,v_on_elm[zz],zz,dof_on_elm[zz]);
      }

      if(counter ==2) // if it has neighbor
      for(INT zz=0; zz<v_per_elm; ++zz)
      {
      printf("V_on_elm_neighbor[%d] = %d  ||  dof_on_elm_neighbor[%d] =%d \n", zz,v_on_elm_neighbor[zz],zz,dof_on_elm_neighbor[zz]);
      }

      printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[0], barycenter->x[0], barycenter->y[0]);
      if(counter ==2) // if it has neighbor
      printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[1], barycenter_neighbor->x[0], barycenter_neighbor->y[0]);
    */
    // Get normal vector values.
    finrm[0]=mesh->f_norm[face*dim+0];
    finrm[1]=mesh->f_norm[face*dim+1];
    if(dim==3) finrm[2]=mesh->f_norm[face*dim+2];
    //printf("FACE = %d, ELM = %d, (Nx,Ny) = (%f,%f) \n", face, neighbor_index[0], finrm[0], finrm[1]);

    ///////////////////////////////////////HEREHERE
    // Get xfi points..
    for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
      //printf("M** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);


      j=jkl-mesh->f_v->IA[face];
      k=mesh->f_v->JA[jkl];

      //printf("M** j = %d, k = %d \n", j, k );


      xfi[j*dim+0]=mesh->cv->x[k];
      xfi[j*dim+1]=mesh->cv->y[k];
      if(dim==3) xfi[j*dim+2]=mesh->cv->z[k];

      //printf("M** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);

    }

    //define the length of the face-
    fiarea=mesh->f_area[face];
    //printf("FACE = %d, ELM = %d, fiarea = %f \n", face, neighbor_index[0], fiarea);

    // get the quadrautre
    zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

    // get the penalty (EnergyNorm)
    penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
    //penalty_term*=LAME_LAMBDA_GLOBAL;
    // SLEE
    //
    for (quad_face=0;quad_face<nquad_face;quad_face++) {

      qx_face[0] = cq_face->x[quad_face];
      qx_face[1] = cq_face->y[quad_face];
      if(dim==3) qx_face[2] = cq_face->z[quad_face];

      w_face = cq_face->w[quad_face];

      // TODO: Don't think this needs to be hardcoded?
      (*truesol)(val_true_face,qx_face,time,
		     &(mesh->el_flag[neighbor_index[0]])); //DUMMY VALUE at the end // ???


      ///////////////////////////////////////////////////
      if(counter == 1){ //only at boundary


	//printf("== BD = %d \n", face);

	//Sanity Check
	if(mesh->f_flag[face] == 0)
	  {printf("check-error0\n"); exit(0);}

	u_comp = u;
	u0_value_at_q = 0.;
	u1_value_at_q = 0.;
  u2_value_at_q = 0.;
	p_value_at_q = 0.;

	//local_row_index=0;
	//unknown_index=0;
	local_dof_on_elm_face = dof_on_elm;

  // TODO: Rewriting this with some hazmath FUNCTIONS
  // u0
  FE_Interpolation(&u0_value_at_q,u_comp,qx_face,local_dof_on_elm_face,v_on_elm,FE->var_spaces[0],mesh);
  // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
  local_dof_on_elm_face += FE->var_spaces[0]->dof_per_elm;
  u_comp += FE->var_spaces[0]->ndof;

  // u1
  FE_Interpolation(&u1_value_at_q,u_comp,qx_face,local_dof_on_elm_face,v_on_elm,FE->var_spaces[1],mesh);
  // get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[1]);
  local_dof_on_elm_face += FE->var_spaces[1]->dof_per_elm;
  u_comp += FE->var_spaces[1]->ndof;

  // u2
  if(dim==3) {
    FE_Interpolation(&u2_value_at_q,u_comp,qx_face,local_dof_on_elm_face,v_on_elm,FE->var_spaces[2],mesh);
    // get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[2]);
    local_dof_on_elm_face += FE->var_spaces[2]->dof_per_elm;
    u_comp += FE->var_spaces[2]->ndof;
  }

  // eg
  // Interpolate
	dof = local_dof_on_elm_face[0]; // get the dof for the last component
	val[0] = u_comp[dof]*   (qx_face[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
	val[1] = u_comp[dof]*   (qx_face[1] - barycenter->y[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
  if(dim==3) val[2] = u_comp[dof]* (qx_face[2] - barycenter->z[0]);
  local_dof_on_elm_face += FE->var_spaces[dim]->dof_per_elm;
  u_comp += FE->var_spaces[dim]->ndof;

  // p
  FE_Interpolation(&p_value_at_q,u_comp,qx_face,local_dof_on_elm_face,v_on_elm,FE->var_spaces[dim+1],mesh);
  // get_FEM_basis(FE->var_spaces[dim+1]->phi,FE->var_spaces[dim+1]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[dim+1]);
  local_dof_on_elm_face += FE->var_spaces[dim+1]->dof_per_elm;
  u_comp += FE->var_spaces[dim+1]->ndof;

	err[0] += penalty_term * w_face * (ABS(u0_value_at_q  +val[0]- val_true_face[0]))*(ABS(u0_value_at_q +val[0] - val_true_face[0]));
	err[1] += penalty_term * w_face * (ABS(u1_value_at_q  +val[1]- val_true_face[1]))*(ABS(u1_value_at_q +val[1] - val_true_face[1]));
  if(dim==3) err[2] += penalty_term * w_face * (ABS(u2_value_at_q + val[2] - val_true_face[2]))*(ABS(u2_value_at_q + val[2] - val_true_face[2]));

	//printf("** EG-u0  ERROR = %f \n", err[0] );
	//printf("** EG-u1  ERROR = %f \n", err[1] );

      } // if on Bonudary conuter == 1
      else if(counter ==2){ //on interface

	//printf("== INTERFACE = %d \n", face);

	if(mesh->f_flag[face] != 0)
	  {printf("check-error --- \n"); exit(0);}



	///////////////////////////////////////////////////////////////////////////////////////
	u_comp_face = u;
	u0_value_at_q = 0.;
	u1_value_at_q = 0.;
  u2_value_at_q = 0.;
	p_value_at_q = 0.;

	local_dof_on_elm_face_interface = dof_on_elm;

  FE_Interpolation(&u0_value_at_q,u_comp_face,qx_face,local_dof_on_elm_face_interface,v_on_elm,FE->var_spaces[0],mesh);
  // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
  local_dof_on_elm_face_interface += FE->var_spaces[0]->dof_per_elm;
  u_comp_face += FE->var_spaces[0]->ndof;

  // u1
  FE_Interpolation(&u1_value_at_q,u_comp_face,qx_face,local_dof_on_elm_face_interface,v_on_elm,FE->var_spaces[1],mesh);
  // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
  local_dof_on_elm_face_interface += FE->var_spaces[1]->dof_per_elm;
  u_comp_face += FE->var_spaces[1]->ndof;

  if(dim==3) {
    FE_Interpolation(&u2_value_at_q,u_comp_face,qx_face,local_dof_on_elm_face_interface,v_on_elm,FE->var_spaces[2],mesh);
    // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
    local_dof_on_elm_face_interface += FE->var_spaces[2]->dof_per_elm;
    u_comp_face += FE->var_spaces[2]->ndof;
  }

	val_face[0] = 0.0;
	val_face[1] = 0.0;
  if(dim==3) val_face[2] = 0.0;
	dof_face = local_dof_on_elm_face_interface[0]; // get the dof for the last component
	val_face[0] = u_comp_face[dof_face]*   (qx_face[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
	val_face[1] = u_comp_face[dof_face]*   (qx_face[1] - barycenter->y[0]);
  if(dim==3) val_face[2] = u_comp_face[dof_face]* (qx_face[2] - barycenter->z[0]);

	//DEBUG
	local_dof_on_elm_face_interface += FE->var_spaces[dim]->dof_per_elm;
	u_comp_face += FE->var_spaces[dim]->ndof;

  // p
  FE_Interpolation(&p_value_at_q,u_comp_face,qx_face,local_dof_on_elm_face_interface,v_on_elm,FE->var_spaces[dim+1],mesh);
  local_dof_on_elm_face_interface += FE->var_spaces[dim+1]->dof_per_elm;
  u_comp_face += FE->var_spaces[dim+1]->ndof;

	///////////////////////////////////////////////////////////////////////////////////////
	u_comp_neighbor = u;
	u0_value_at_q_neighbor = 0.;
	u1_value_at_q_neighbor = 0.;
  u2_value_at_q_neighbor = 0.;
	p_value_at_q_neighbor = 0.;
	local_dof_on_elm_neighbor = dof_on_elm_neighbor;

  FE_Interpolation(&u0_value_at_q_neighbor,u_comp_neighbor,qx_face,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[0],mesh);
  // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
  local_dof_on_elm_neighbor += FE->var_spaces[0]->dof_per_elm;
  u_comp_neighbor += FE->var_spaces[0]->ndof;

  // u1
  FE_Interpolation(&u1_value_at_q_neighbor,u_comp_neighbor,qx_face,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[1],mesh);
  // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
  local_dof_on_elm_neighbor += FE->var_spaces[1]->dof_per_elm;
  u_comp_neighbor += FE->var_spaces[1]->ndof;

  if(dim==3) {
    FE_Interpolation(&u2_value_at_q_neighbor,u_comp_neighbor,qx_face,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[2],mesh);
    // get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[0]);
    local_dof_on_elm_neighbor += FE->var_spaces[2]->dof_per_elm;
    u_comp_neighbor += FE->var_spaces[2]->ndof;
  }

	val_neighbor[0] = 0.0;
	val_neighbor[1] = 0.0;
  if(dim==3) val_neighbor[2] = 0.0;
	dof_neighbor = local_dof_on_elm_neighbor[0]; // get the dof for the last component
	val_neighbor[0] = u_comp_neighbor[dof_neighbor]*   (qx_face[0] - barycenter_neighbor->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
	val_neighbor[1] = u_comp_neighbor[dof_neighbor]*   (qx_face[1] - barycenter_neighbor->y[0]);
  if(dim==3) val_neighbor[2] = u_comp_neighbor[dof_neighbor]* (qx_face[2] - barycenter_neighbor->z[0]);

	//DEBUG
	local_dof_on_elm_face_interface += FE->var_spaces[dim]->dof_per_elm;
	u_comp_face += FE->var_spaces[dim]->ndof;

  // p
  FE_Interpolation(&p_value_at_q_neighbor,u_comp_neighbor,qx_face,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[dim+1],mesh);
  local_dof_on_elm_neighbor += FE->var_spaces[dim+1]->dof_per_elm;
  u_comp_neighbor += FE->var_spaces[dim+1]->ndof;

	if(ABS(u0_value_at_q - u0_value_at_q_neighbor) > 0.0000001)
	  {
	    printf("u0_value_at_q_face - u0_value_at_q_face_neighbor = %f \n", u0_value_at_q - u0_value_at_q_neighbor);
	    //exit(0);
	  }

	if(ABS(u1_value_at_q - u1_value_at_q_neighbor) > 0.0000001)
	  {
	    printf("u1_value_at_q_face - u1_value_at_q_face_neighbor = %f \n", u1_value_at_q - u1_value_at_q_neighbor);

	    //exit(0);
	  }

	//printf("penalty = %f \n", penalty_term);
	err[0] += penalty_term * w_face * ABS( ( val_face[0]- val_neighbor[0] ) )*ABS( val_face[0]- val_neighbor[0] );
	err[1] += penalty_term * w_face * ABS( ( val_face[1]- val_neighbor[1] ) )*ABS( val_face[1]- val_neighbor[1] );
  if(dim==3) err[2] += penalty_term * w_face * ABS( (val_face[2]- val_neighbor[2] ) )*ABS( val_face[2]- val_neighbor[2] );

      }// on interface

    }//quad
  }//face

  if(neighbor_basis_u0_phi) free(neighbor_basis_u0_phi);
  if(neighbor_basis_u0_dphi) free(neighbor_basis_u0_dphi);
  if(neighbor_basis_u1_phi) free(neighbor_basis_u1_phi);
  if(neighbor_basis_u1_dphi) free(neighbor_basis_u1_dphi);
  if(dim==3) {
    if(neighbor_basis_u2_phi) free(neighbor_basis_u2_phi);
    if(neighbor_basis_u2_dphi) free(neighbor_basis_u2_dphi);
  }
  if(neighbor_basis_p_phi) free(neighbor_basis_p_phi);
  if(neighbor_basis_p_dphi) free(neighbor_basis_p_dphi);
  if(val_true_face) free(val_true_face);
  if(val_true_face_neighbor) free(val_true_face_neighbor);
  icsr_free(f_el);
  if(f_el) free(f_el);
  if(barycenter) {
    free_coords(barycenter);
    free(barycenter);
    barycenter = NULL;
  }
  if(barycenter_neighbor) {
    free_coords(barycenter_neighbor);
    free(barycenter_neighbor);
    barycenter_neighbor = NULL;
  }
  if(cq) {
    free_qcoords(cq_face);
    free(cq_face);
    cq_face=NULL;
  }
  if(data_face) free(data_face);
  if(dof_on_elm) free(dof_on_elm);
  if(dof_on_elm_neighbor) free(dof_on_elm_neighbor);
  if(v_on_elm) free(v_on_elm);
  if(v_on_elm_neighbor) free(v_on_elm_neighbor);
  // Get the BD values
  // for(i=0;i<2;i++) {
  //err[i] = sqrt(err[i]);
  //}

  //exit(0);
}
