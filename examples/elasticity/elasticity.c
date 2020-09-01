/*! \file examples/Elasticity/elasticity.c
 *
 *  This code is extended version of stokes.c  
 *  created by Peter Ohm on 2/5/17. (S. Lee 8/12/20)
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Elasticity PDE using finite elements
 *
 *      -2* mu * div(eps(u))  - lambda div (div (u) )  = f
 *        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient,
 *        mu and lambda are lame coefficients
 *        in 2D or 3D.
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed for u and Neumann for p.  P2-P1 or P2-P0 can be used,
 *        though others can be implemented.
 *
 * \note This example shows how to build your own bilinear form for a system.
 *       The forms are found in elasticity_system.h and all Problem Data is found in
 *       elasticity_data.h to simplify the code.  This example also illustrates how to
 *       construct block versions of the finite-element spaces, linear systems, and
 *       how to use block solvers.
 *
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
//#include "elasticity_data_2.h"
#include "elasticity_data.h"
//#include "elasticity_system.h"
/*********************************************************************************/

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

/******************************************************************************************************/
/*!
 * \fn void FEM_Block_RHS_Local(REAL* bLoc,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local Right hand side vector for a block FEM system
 *        b_i  = <f,phi_i>
 *
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific Vertices on element
 * \param elm           Current element
 * \param rhs           Function that gives RHS (in FEM block ordering
 * \param time          Physical Time if time dependent
 *
 * \return bLoc         Local RHS Vector
 *
 *
 */

// ISSUE : this doesn't allow to get the exact solution out...  
void FEM_Block_RHS_Local_Elasticity(REAL* bLoc,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Loop Indices
  INT i,quad,test;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_elm = 0;
  INT nun=0;

  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) /* Scalar Element */
      nun++;
    else /* Vector Element */
      nun += dim;
  }
  INT* local_dof_on_elm = NULL;
  INT local_row_index=0;
  INT unknown_index=0;

  // Quadrature Weights and Nodes
  REAL w;
  INT maxdim=4;
  REAL qx[maxdim];

  // Right-hand side function at Quadrature Nodes
  REAL rhs_val[nun];
  // SLEE
  // Get the BD values 
  REAL* val_true_face = (REAL *) calloc(nun,sizeof(REAL));

  //printf("nun = %d \n", nun );
  //exit(0);
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    (*rhs)(rhs_val,qx,time,&(mesh->el_flag[elm]));

    local_row_index=0;
    unknown_index=0;
    local_dof_on_elm=dof_on_elm;
    
    for(i=0;i<FE->nspaces;i++) {

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);
      
      // Loop over test functions and integrate rhs
      if(FE->var_spaces[i]->FEtype<20) { // Scalar Element

	//printf("1 = %d\n",i);
	for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {

	  //sanity check
	

	  //SLEE
	  if(i == 0 || i == 1)
	    {
	      //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      
	      bLoc[(local_row_index+test)] += w*rhs_val[unknown_index]*FE->var_spaces[i]->phi[test];
	    }
	    else if(i == 2)
	    {
	      //SLEE
	      // Note that new DOF is \phi^3 = [x ; y] 
	      //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      bLoc[(local_row_index+test)] +=w*(rhs_val[0]*  qx[0]  +rhs_val[1]* qx[1]);
	    }
	  
	  //sanity check
	  if(unknown_index != i)
	    {
	      printf("unknown index != i \n");
	      exit(0);
	    }
	  //SLEE

	  
        }
	//SLEE : postpone this for later ?? No
        unknown_index++;
	
      } else if (FE->var_spaces[i]->FEtype==61) { // bubble
        for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
          bLoc[(local_row_index+test)] += w*(rhs_val[unknown_index]*FE->var_spaces[i]->phi[test*dim] +
					     rhs_val[unknown_index+1]*FE->var_spaces[i]->phi[test*dim+1]);
          if(dim==3) bLoc[(local_row_index+test)] += w*rhs_val[unknown_index+2]*FE->var_spaces[i]->phi[test*dim+2];
        }
        //
        // no update of unknown_index.
        // Assuming that bubble is the first FE space and that the space it is enriching
        // follow immediately
        //
      }
      else { // Vector Element
        for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
	  //SLEE
	  printf("NOT HERE : %d\n",i);
	  exit(0);
          bLoc[(local_row_index+test)] +=
	    w*(rhs_val[unknown_index]*FE->var_spaces[i]->phi[test*dim] +
	       rhs_val[unknown_index+1]*FE->var_spaces[i]->phi[test*dim+1]);
	  
	  if(dim==3) bLoc[(local_row_index+test)] += w*rhs_val[unknown_index+2]*FE->var_spaces[i]->phi[test*dim+2];
        }
        unknown_index += dim;
      }

      //SLEE : postpone this for later ?? No
      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
      local_row_index += FE->var_spaces[i]->dof_per_elm;
    }
  }

  //DEBUG
  //exit(0);

  REAL penalty_term = 0.;

  
  for (INT face=0; face<mesh->nface; face++) {

    //NEED TO DEBUG!!!!
    //SLEE // SEEK BD flag...?????

    //printf("mesh->nface = %d, mesh->f_flag[face] = %d \n", mesh->nface,  mesh->f_flag[0]);
    //printf("mesh->nface = %d, mesh->f_flag[face] = %d \n", mesh->nface, mesh->f_flag[1]);
    //printf("mesh->nface = %d, mesh->f_flag[face] = %d \n", mesh->nface, mesh->f_flag[2]);
    //exit(0);

    // TODO : Ludmil
    // IS THIS THE RIGHT WAY TO CHECK THE BD??
    if(mesh->f_flag[face]>=0 && mesh->f_flag[face]<=1) {

    INT i,j,rowa, rowb, jcntr,ed;
    // Quadrature Weights and Nodes
    INT nq_face = 2*dim-3; // = ed_per_face
    
    REAL* qx_face = (REAL *) calloc(nq_face,sizeof(REAL));

    // 3D: Using triangle midpoint rule, so qx is midpoint of edges and w is |F|/3
    REAL w_face = mesh->f_area[face]/3.0;

    // Find the edges for the given face
    // NEED TO CHECK ?? 
    INT* ed_on_f = (INT *) calloc(nq_face,sizeof(INT));
    rowa = mesh->f_ed->IA[face];
    rowb = mesh->f_ed->IA[face+1];
    jcntr = 0;
    for(i=rowa;i<rowb;i++) {
      ed_on_f[jcntr] = mesh->f_ed->JA[i];
      jcntr++;
    }

    local_row_index=0;
    unknown_index=0;
    
    local_dof_on_elm = NULL;
    local_dof_on_elm=dof_on_elm;
    
    // Get normal vector components on face
    REAL nx = mesh->f_norm[face*dim];
    REAL ny = mesh->f_norm[face*dim+1];
    REAL nz = 0.0;
    if(dim==3) nz = mesh->f_norm[face*dim+2];
    
    // Loop over each block and get dof_per_face
    INT dof_per_face = 0;
    INT* dof_per_face_blk = (INT*)calloc(FE->nspaces,sizeof(INT));
    INT FEtype;
    for(i=0;i<FE->nspaces;i++) {
      FEtype = FE->var_spaces[i]->FEtype;
      if(FEtype>=1 && FEtype<10) { // PX Elements
	dof_per_face_blk[i] = dim + (FEtype-1)*(2*dim-3);
	dof_per_face += dim + (FEtype-1)*(2*dim-3);
      } else if (FEtype==20) { // Nedelec Elements
	dof_per_face_blk[i] = 2*dim - 3;
	dof_per_face += 2*dim - 3;
      } else if (FEtype==30) { // Raviart-Thomas Elements
	dof_per_face_blk[i] = 1;
	dof_per_face += 1;
      } else if (FEtype==61) { // Bubbles
	dof_per_face_blk[i] = 1;
	dof_per_face += 1;
      } else if (FEtype==0) { // P0
	// Questionable handling of P0 for the face integral
	dof_per_face_blk[i] = 1;
	dof_per_face += 1;
      } else {
	printf("Block face integration isn't set up for the FEM space you chose\n");
	exit(0);
      }
    }

    //printf("DOF per FACE = %d \n", dof_per_face);
    //exit(0);

      
    //  Sum over midpoints of edges
    for (quad=0;quad<nq_face;quad++) {

      ed = ed_on_f[quad];

      qx_face[0] = mesh->ed_mid[ed*dim];
      qx_face[1] = mesh->ed_mid[ed*dim+1];

      
      // TODO : Ludmil
      // SANITY CHECK ??
      // CHECK :: CENTER OF MESH
      // 
      //! midpoint of face
      //REAL* f_mid;
      //qx_face[0] = mesh->f_mid[face*dim];
      //qx_face[1] = mesh->f_mid[face*dim+1];
      
      if(fabs(qx_face[0] -mesh->f_mid[face*dim])>1e-10){
	  printf("NOT SAME ! X \n");
      }      
      if(fabs(qx_face[1] -mesh->f_mid[face*dim+1])>1e-10){
	printf("NOT SAME ! Y \n");
      }      
      // SLEE
      // Get the BD values 
      // Get True Solution at Quadrature Nodes
      (*exact_sol2D)(val_true_face,qx_face,time,&(mesh->el_flag[elm]));      
      //printf("qx_face[0] = %f\n",qx_face[0]);
      //printf("qx_face[1] = %f\n",qx_face[1]);
      //printf("local_dof_on_face = %d \n", local_dof_on_elm[0]);
      //printf("local_dof_on_face = %d \n", local_dof_on_elm[1]);
      //printf("local_dof_on_face = %d \n", local_dof_on_elm[2]);
      
      for(i=0;i<FE->nspaces;i++) {
	
	//FE->var_spaces[i]->phi = NULL;
	//FE->var_spaces[i]->dphi = NULL;
	PX_H1_basis(FE->var_spaces[i]->phi,
		    FE->var_spaces[i]->dphi,
		    qx_face,
		    local_dof_on_elm,
		    FE->var_spaces[i]->FEtype,
		    mesh);

	for (INT test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
	  //SLEE
	  if(i==0)
	    {
	      // printf("(test=0 / 2) i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      bLoc[(local_row_index+test)] += penalty_term * w_face*
		(val_true_face[unknown_index]* nx * FE->var_spaces[i]->phi[test] * nx);
	    }
	  
	  else if(i == 1)
	    {
	      //printf("(test=1 / 3) i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      bLoc[(local_row_index+test)] += penalty_term *  w_face* val_true_face[unknown_index]* ny
		* FE->var_spaces[i]->phi[test] * ny;
	    }
	  else if(i == 2)
	    {
	      //SLEE
	      // Note that new DOF is \phi^3 = [x ; y] 
	      //printf("(i=2?) i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      
	      bLoc[(local_row_index+test)] +=
		penalty_term *  w_face * (val_true_face[0]*  nx  + val_true_face[1]*  ny)
		* (qx_face[0]*  nx  + qx_face[1]*  ny);
	    }
	  //bLoc[(local_row_index+test)] += w * val_true[unknown_index] *
	  //* FE->var_spaces[i]->phi[test] * ;

	    
	  //sanity check
	  if(unknown_index != i)
	    {
	      printf("unknown index != i \n");
	      exit(0);
	    }
	  //SLEE

	}
	unknown_index++;
	
	//printf(" i = %d\n", i);
	//printf("FE->var_spaces[i]->dof_per_elm = %d \n", FE->var_spaces[i]->dof_per_elm);
	//printf("local_dof_on_face = %d \n", local_dof_on_elm[0]);
	//printf("local_dof_on_face = %d \n", local_dof_on_elm[1]);
	//printf("local_dof_on_face = %d \n", local_dof_on_elm[2]);
	local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;	
	local_row_index += FE->var_spaces[i]->dof_per_elm; 
	
      }
      //exit(0);
      
	
	 }
    
    }
  }

  return;
}


void L2error_block_EG
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
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
    err[i]=0;
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    if(FE->var_spaces[i]->FEtype<20) 
    //Scalar Element 
      ncomp[i]=1;
    else // Vector Element 
      ncomp[i] = dim;
    nun += ncomp[i];
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  // Loop over all Elements 
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
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      // Interpolate FE solution to quadrature point
      blockFE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);


      // in blockFE_Interopation for val_sol[dim]
      
      //for k==FE->nspaces
      //copied and modified from FE_Interpolation
      INT i,j,dof;
      REAL val[dim];

      //printf(" qx[0] = %f \n " ,  qx[0]);
      
      val[0] = 0.0;
      val[1] = 0.0;

      //printf("dof_per_elm = %d \n", dof_per_elm);
      dof = dof_on_elm[dof_per_elm-1]; // get the dof for the last component
      val[0] += u[dof]*   qx[0]; //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      val[1] += u[dof]*   qx[1]; //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
     
      // Compute Square of Error on Element for each component of FE space
      jcntr=0;
      for(i=0;i<nspaces;i++) {
        for(j=0;j<ncomp[i];j++) {

	  //  printf("nspaces = %d \n" , nspaces);
	  //  printf("ncomp[0] = %d ,  ncomp[1] = %d, \n", ncomp[0],ncomp[1]);
	  //  printf("(i,j,jcntr) = %d, %d, %d \n", i,j,jcntr);
	 
          err[i]+=w*(ABS(val_sol[jcntr+j] + val[jcntr+j]  - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] +val[jcntr+j] - val_true[jcntr+j]));
        }
        jcntr+=ncomp[i];
      }

      //exit(0);
      
    }
  }


  printf("---------------------1\n");
  for(i=0;i<nspaces;i++) {
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
  double d_Lame_coeff_lambda = 1.; //000000.;
  
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
      for(i=0;i<nspaces;i++) {
        for(j=0;j<ncomp[i];j++) {
	  // \mu * grad u 
	  err[i]+=d_Lame_coeff_mu*w*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]))*(ABS(val_sol[jcntr+j] - val_true[jcntr+j]));
	}
        jcntr+=ncomp[i];
      }
      if(dim == 3)
	{
	  printf("ERROR IN STRESS ERROR \n"); exit(0);
	}
      else if(dim == 2){

	// grad u^T
	
	err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
	err[0] += d_Lame_coeff_mu * w*(ABS(val_sol[2] - val_true[2]))*(ABS(val_sol[2] - val_true[2]));	
	
	err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[1] - val_true[1]))*(ABS(val_sol[1] - val_true[1]));
	err[1] += d_Lame_coeff_mu * w*(ABS(val_sol[3] - val_true[3]))*(ABS(val_sol[3] - val_true[3]));

	// div u
	err[0] += d_Lame_coeff_lambda * w * (ABS(val_sol[0] - val_true[0]))*(ABS(val_sol[0] - val_true[0]));
	err[1] += d_Lame_coeff_lambda * w * (ABS(val_sol[3] - val_true[3]))*(ABS(val_sol[3] - val_true[3]));
      }


      
    }
  }

  for(i=0;i<nspaces;i++) {
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

void local_assembly_Elasticity(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  // Loop indices
  INT i,j,idim,quad,test,trial;

  // Mesh and FE data
  INT dof_per_elm = 0;
  
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;

  // SANITY CHECK SLEE
  //printf("------\n");
  
  //printf("CHECK POINT: FE->nspaces = %d, dof_per_elem = %d, FE->var_spaces[0]->dof_per_elm = %d,  FE->var_spaces[1]->dof_per_elm = %d ,  FE->var_spaces[1]->dof_per_elm = %d \n", FE->nspaces, dof_per_elm, FE->var_spaces[0]->dof_per_elm,FE->var_spaces[1]->dof_per_elm, FE->var_spaces[2]->dof_per_elm);

  //printf("-------\n");
  //exit(0);	 
 
  INT* local_dof_on_elm;
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  double lambda = 1.; //000000.0;

  
  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];

    if(mesh->dim==3)
      qx[2] = cq->z[elm*cq->nq_per_elm+quad];

    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u = (u1,u2,u3) and v = (v1,v2,v3)
    get_FEM_basis(FE->var_spaces[0]->phi,
		  FE->var_spaces[0]->dphi,
		  qx,
		  v_on_elm,
		  dof_on_elm,
		  mesh,
		  FE->var_spaces[0]);
    
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;

    //printf(" local_dof_on_elm = %d,  dof_on_elm = %d, FE->var_spaces[0]->dof_per_elm = %d, \n",
    //	   *local_dof_on_elm,
    //	   *dof_on_elm,
    //	   FE->var_spaces[0]->dof_per_elm);
    
    get_FEM_basis(FE->var_spaces[1]->phi,
		  FE->var_spaces[1]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[1]);
    if(dim==3){
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;

      get_FEM_basis(FE->var_spaces[2]->phi,
		    FE->var_spaces[2]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[2]);
    }

    // p
    local_dof_on_elm += FE->var_spaces[dim-1]->dof_per_elm;
    
    get_FEM_basis(FE->var_spaces[dim]->phi,
		  FE->var_spaces[dim]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[dim]);

    //printf("* local_dof_on_elm = %d,  dof_on_elm = %d, FE->var_spaces[0]->dof_per_elm = %d, \n",
    //	   local_dof_on_elm[1],
    //	   *dof_on_elm,
    //	   FE->var_spaces[0]->dof_per_elm);

    //exit(0);
    
    // u1-v1 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	kij = 2.0*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[0]->dphi[trial*dim];
	kij += 1.0*FE->var_spaces[0]->dphi[test*dim+1]*FE->var_spaces[0]->dphi[trial*dim+1];

	// NEW SLEE : elasticity div part
	// u1-v1 block : <dx(u1),dx(v1)>
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[0]->dphi[trial*dim];
	
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    // u1-v2 block <dy(u1),dx(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[1]->dphi[test*dim+0]*FE->var_spaces[0]->dphi[trial*dim+1];
	
	// NEW SLEE : elasticity div part
	// u1-v2 block : <dx(u1),dx(v2)> 
	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[0]->dphi[trial*dim];

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    // u2-v1 block : <dx(u2),dy(v1)>>
    local_row_index = 0;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[0]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+0];

	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)> 
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[1]->dphi[trial*dim+1];
	
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }



    // u2-v2 block <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	kij = 2.0*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+1];
	kij += 1.0*FE->var_spaces[1]->dphi[test*dim]*FE->var_spaces[1]->dphi[trial*dim];

	// NEW SLEE : elasticity div part
	// u2-v2 block : <dy(u2),dx(v2)> 
   	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+1];
     	
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }


    //EG PART
    
    //NEW SLEE :
    // u1- q block < 2mu e(u1) * e(q) > + lambda <div u1 : e(q) >
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	kij = 2.0*(FE->var_spaces[0]->dphi[trial*dim]) * 1.; //FE->var_spaces[dim]->dphi[test*dim];
	
	// Divergence
	// u1-q block : <dx(u1),dx(v1)>
	kij += lambda*FE->var_spaces[0]->dphi[trial*dim]*1.;//FE->var_spaces[dim]->dphi[test*dim];	
	kij += lambda*FE->var_spaces[0]->dphi[trial*dim]*1.;


	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    // u2-q block
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	kij = 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;

	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)> 
	kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;	
	kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    
    // q-v1 block: 
    local_row_index = 0;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	kij = 2.0*(FE->var_spaces[0]->dphi[test*dim]) * 1.; //FE->var_spaces[dim]->dphi[test*dim];
	
	// Divergence
	// u1-q block : <dx(u1),dx(v1)>
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*1.;//FE->var_spaces[dim]->dphi[test*dim];
	kij += lambda*FE->var_spaces[0]->dphi[test*dim]*1.;

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    
    // q-v2 block: 
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	kij = 2.0*FE->var_spaces[1]->dphi[test*dim+1] * 1.;

	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)> 
	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;	
	kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    
    //q-q block: 
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	kij = 2.0 * 1. * 1.;
	kij += 2.0 * 1. * 1.;
	
	// NEW SLEE : elasticity div part
	// u2-v1 block : <dy(u2),dx(v1)> 
	kij += lambda * 1. * 1.;	
	kij += lambda * 1. * 1.;

	kij += lambda * 1. * 1.;
	kij += lambda * 1. * 1.;
	
	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }
    
    
    


    if(dim==3){
      // u3-v3 block <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij=FE->var_spaces[2]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+2];
          for(idim=0;idim<dim;idim++){
            kij += FE->var_spaces[2]->dphi[test*dim+idim]*FE->var_spaces[2]->dphi[trial*dim+idim];
          }
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // p-v3 block: -<p, dz(u3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
          kij = -FE->var_spaces[dim]->phi[trial]*(FE->var_spaces[2]->dphi[test*dim+2]);
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u3-q block: -<dz(u3), q>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = -(FE->var_spaces[2]->dphi[trial*dim+2])*FE->var_spaces[dim]->phi[test];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u1-v3 block <dz(u1),dx(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
          kij = FE->var_spaces[2]->dphi[test*dim+0]*FE->var_spaces[0]->dphi[trial*dim+2];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u3-v1 block <dx(u3),dz(v1)>
      local_row_index = 0;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = FE->var_spaces[0]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+0];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u2-v3 block <dz(u2),dy(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
          kij = FE->var_spaces[2]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+2];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u3-v2 block <dy(u3),dz(v2)>
      local_row_index = FE->var_spaces[0]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = FE->var_spaces[1]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+1];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }
    }
  }

  
  //SLEE TODO
  // 1. Need Face Quadrature
  // 2. Need FE_Basis on the Faces / Boundaries
  //    -- this is not developed in HAZmath.
  //    -- currently only works in P1 ? 
  //printf("----------$\n");
  REAL penalty_term = 0.;
  
  // Loop over each block and get dof_per_face
  INT dof_per_face = 0;
  INT* dof_per_face_blk = (INT *) calloc(FE->nspaces,sizeof(INT));
  INT FEtype;
  for(i=0;i<FE->nspaces;i++) {
    FEtype = FE->var_spaces[i]->FEtype;
    if(FEtype>=1 && FEtype<10) { // PX Elements
      dof_per_face_blk[i] = dim + (FEtype-1)*(2*dim-3);
      dof_per_face += dim + (FEtype-1)*(2*dim-3);
    } else if (FEtype==20) { // Nedelec Elements
      dof_per_face_blk[i] = 2*dim - 3;
      dof_per_face += 2*dim - 3;
    } else if (FEtype==30) { // Raviart-Thomas Elements
      dof_per_face_blk[i] = 1;
      dof_per_face += 1;
    } else if (FEtype==61) { // Bubbles
      dof_per_face_blk[i] = 1;
      dof_per_face += 1;
    } else if (FEtype==0) { // P0
      // Questionable handling of P0 for the face integral
      dof_per_face_blk[i] = 1;
      dof_per_face += 1;
    } else {
      printf("Block face integration isn't set up for the FEM space you chose\n");
      exit(0);
    }
  }
  //LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLl  
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  REAL *data_face_end=finrm+dim; //end
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
  REAL fiarea=0e0;
  INT jk,k,face,quad_face,rowa, rowb, jcntr,ed;
  INT maxdim=4;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face=quad<cq_face->nq_per_elm;
  //Ludmil: this loops over all faces in the mesh!
  for (face=0; face<mesh->nface; face++) {
    for(jk=mesh->f_v->IA[face];jk<mesh->f_v->IA[face+1];jk++){
      j=jk-mesh->f_v->IA[face];
      k=mesh->f_v->JA[jk];
      xfi[j*dim+0]=mesh->cv->x[k];
      if(dim>1)
	xfi[j*dim+1]=mesh->cv->y[k];
      if(dim>2)
	xfi[j*dim+2]=mesh->cv->z[k];	  
    }
    finrm[0]=mesh->f_norm[face*dim+0];
    if(dim>1)
      finrm[1]=mesh->f_norm[face*dim+1];
    if(dim>2)
      finrm[2]=mesh->f_norm[face*dim+2];
    // Only grab the faces on the flagged boundary
    //NEED TO DEBUG!!!!
    //SLEE // SEEK BD flag...?????
    fiarea=mesh->f_area[face];
    zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);
    // elements
    if(mesh->f_flag[face]>=0 && mesh->f_flag[face]<=1) { 
      //CHECK.ARE WE USING THE SAME QUAD??
      // Quadrature Weights and Nodes  
      //  Sum over midpoints of edges
      for (quad_face=0;quad_face<nquad_face;quad_face++) {	
	// TODO : Ludmil  // This doesn't work!
	qx_face[0] = cq_face->x[quad_face];
	qx_face[1] = cq_face->y[quad_face];
	if(dim==3) qx_face[2] = cq_face->z[quad_face]; 
	w_face = cq_face->w[quad_face];
	get_FEM_basis(FE->var_spaces[0]->phi,
		      FE->var_spaces[0]->dphi,
		      qx_face,
		      v_on_elm,
		      dof_on_elm,
		      mesh,
		      FE->var_spaces[0]);
	
	local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
	
	
	get_FEM_basis(FE->var_spaces[1]->phi,
		      FE->var_spaces[1]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm,
		      mesh,
		      FE->var_spaces[1]);
	
	// p
	//local_dof_on_elm += FE->var_spaces[dim-1]->dof_per_elm;
	
	/*get_FEM_basis(FE->var_spaces[dim]->phi,
		      FE->var_spaces[dim]->dphi,
		      qx,
		      v_on_elm,
		      local_dof_on_elm,
		      mesh,
		      FE->var_spaces[dim]);
	*/
	
	//SLEE
	//WE DON"T USE THIS BUT MANUALLY 
	// FE->var_spaces[dim]->phi  =  X
	//FE->var_spaces[dim+1]->phi   = Y

	
	// u1-v1 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
	local_row_index = 0;
	local_col_index = 0;


	double penalty_term = 0.;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	    
	    kij = penalty_term * FE->var_spaces[0]->phi[test]* finrm[0] * FE->var_spaces[0]->dphi[trial] * finrm[1];
	    
	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	  }
	}
	
      }// quad face
      
      
      
    }
  }
  

  return;
}
/*****************************************************************************************************/


/****** MAIN DRIVER **************************************************************/
int main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve Elasticity Equation.\n");
  printf("===========================================================================\n");

  // Aug.3.2020 SLEE
  // Define variables forthe error convergence test
  int total_num_cycle = 3; 
  // SLEE initialize the vectors to save the errors for each cycle
  double L2_error_per_cycle[total_num_cycle];
  double L2_error_p_per_cycle[total_num_cycle];
  double H1_error_per_cycle[total_num_cycle];
  double H1_stress_error_per_cycle[total_num_cycle];

  double L2_EG_error_per_cycle[total_num_cycle];

  // SLEE vector to save the DOF
  int dof_per_cycle[total_num_cycle];
  // SLEE vector to save the convergence rate
  double L2_conv_rate_per_cycle[total_num_cycle];
  double L2_p_conv_rate_per_cycle[total_num_cycle];
  double H1_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_conv_rate_per_cycle[total_num_cycle];

  double L2_EG_conv_rate_per_cycle[total_num_cycle];
  
  int global_dim_space = 0;
  
  for(int cycle=0; cycle<total_num_cycle; ++cycle)
    {
      //Aug.3.2020 SLEE initilize
      L2_error_per_cycle[cycle] = 0.;
      L2_error_p_per_cycle[cycle] = 0.;
      H1_error_per_cycle[cycle] = 0.;
      H1_stress_error_per_cycle[cycle] = 0.;

      L2_EG_error_per_cycle[cycle] = 0.;
      
      dof_per_cycle[cycle]=0;
      L2_conv_rate_per_cycle[cycle]=0.;
      L2_p_conv_rate_per_cycle[cycle]=0.;
      H1_conv_rate_per_cycle[cycle]=0.;
      H1_stress_conv_rate_per_cycle[cycle]=0.;

      L2_EG_conv_rate_per_cycle[cycle]=0.;
      
      printf("************ CYCLE   %d  /   %d  ************** \n", cycle, total_num_cycle);
    
  /****** INITIALIZE PARAMETERS **************************************************/
  INT i;

  // Overall CPU Timing
  clock_t clk_overall_start = clock();

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  // Open gridfile for reading
  printf("\nCreating mesh and FEM spaces:\n");
  //FILE* gfid = HAZ_fopen(inparam.gridfile,"r");
  //SLEE
  FILE* gfid;

  //Jul.10.2020 SLEE: setup the code to read the different mesh files for each cycle 
  //gfid = HAZ_fopen(inparam.gridfile,"r");
  char filename_per_cycle[100]={'\0'};
  sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle);


  
  gfid = HAZ_fopen(filename_per_cycle,"r");
  if(gfid == NULL){
    perror("Could not find and open the file !!!! ");
  }

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  INT mesh_type = 0;
  clock_t clk_mesh_start = clock(); // Time mesh generation FE setup
  mesh_struct mesh;

  //printf(" --> loading grid from file: %s\n",inparam.gridfile);
  //Jul.10. 2020 SLEE
  printf(" --> loading grid from file: %s\n",filename_per_cycle);
  
  initialize_mesh(&mesh);   // Q1. Why only here? 
  creategrid_fread(gfid,mesh_type,&mesh);
  fclose(gfid);
  
  INT dim = mesh.dim;
  // Jul.10.2020 SLEE : for further use in the convergence test
  global_dim_space = dim;

  // Get Quadrature Nodes for the Mesh
  INT nq1d = inparam.nquad; /* Quadrature points per dimension */
  qcoordinates *cq = get_quadrature(&mesh,nq1d);

  // Get info for and create FEM spaces
  // Order of elements: 0 - P0; 1 - P1; 2 - P2; 20 - Nedlec; 30 - Raviart-Thomas
  INT order_u = 1;
  INT order_p = 0;

  // Need Spaces for each component of the velocity plus pressure
  fespace FE_ux; // Velocity in x direction
  create_fespace(&FE_ux,&mesh,order_u);
  fespace FE_uy; // Velocity in y direction
  create_fespace(&FE_uy,&mesh,order_u);
  fespace FE_uz; // Velocity in z direction
  if(dim==3) create_fespace(&FE_uz,&mesh,order_u);
  fespace FE_p; // Pressure
  create_fespace(&FE_p,&mesh,order_p);

  // Set Dirichlet Boundaries
  set_dirichlet_bdry(&FE_ux,&mesh,1,1);
  set_dirichlet_bdry(&FE_uy,&mesh,1,1);
  if(dim==3) set_dirichlet_bdry(&FE_uz,&mesh,1,1);
  set_dirichlet_bdry(&FE_p,&mesh,1,1);
  for(i=0;i<FE_p.ndof;i++) {
    FE_p.dirichlet[i] = 0;
  }

  // Create Block System with ordering (u,p)
  INT ndof = FE_ux.ndof + FE_uy.ndof + FE_p.ndof;
  if(dim==3) ndof += FE_uz.ndof;
  // Get Global FE Space
  block_fespace FE;
  FE.nun = dim+1;  // SLEE? 
  FE.ndof = ndof;
  FE.nbdof = FE_ux.nbdof + FE_uy.nbdof + FE_p.nbdof;
  if(dim==3) FE.nbdof += FE_uz.nbdof;
  FE.nspaces = dim+1; // SLEE?
  FE.var_spaces = (fespace **) calloc(dim+1,sizeof(fespace *));

  FE.var_spaces[0] = &FE_ux;
  FE.var_spaces[1] = &FE_uy;
  if(dim==3) FE.var_spaces[2] = &FE_uz;
  FE.var_spaces[dim] = &FE_p;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(&FE,&mesh);


  clock_t clk_mesh_end = clock(); // End of timing for mesh and FE setup
  printf(" --> elapsed CPU time for mesh and FEM space construction = %f seconds.\n\n",
         (REAL) (clk_mesh_end - clk_mesh_start)/CLOCKS_PER_SEC);
  /*******************************************************************************/

  printf("***********************************************************************************\n");
  printf("Number of Elements = %d\tOrder of Quadrature = %d\n",mesh.nelm,2*nq1d-1);
  printf("\n\t--- Element Type ---\n");
  printf("Velocity Element Type = %d\tPressure Element Type = %d\n",order_u,order_p);
  printf("\n\t--- Degrees of Freedom ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nv,mesh.nedge,mesh.nface);
  printf("\t--> DOF: %d\n",FE.ndof);
  printf("\n\t--- Boundaries ---\n");
  printf("Vertices: %-7d\tEdges: %-7d\tFaces: %-7d",mesh.nbv,mesh.nbedge,mesh.nbface);
  printf("\t--> Boundary DOF: %d\n",FE.nbdof);
  printf("***********************************************************************************\n\n");

  //Jul.10.2020 SLEE: insert the total number of DOF for convergnece computation 
  dof_per_cycle[cycle] = FE.ndof + FE.nbdof;
  
  /*** Assemble the matrix and right hand side *******************************/
  /* Here we assemble the discrete system:
   *  The weak form is:
   *
   *  <2*eps(u), eps(v)> - <p, div v> = <f, v>
   *                   - <div u, q> = 0
   */
  printf("Assembling the matrix and right-hand side:\n");
  clock_t clk_assembly_start = clock();

  // Allocate the right-hand and declare the csr matrix
  dvector b;
  //SLEE Aug 27 2020
  // WHAT DO I NEED ? 
  //dvector b_bdry;// = dvec_create(FE_q.ndof);
  //dvec_set(FE_q.ndof,&b_bdry,0.0);
  
  // Put into Block Form
  block_dCSRmat A;
  bdcsr_alloc(dim+1,dim+1,&A);

  // Assemble the matricies without BC first
  if(dim==2) assemble_global_block(&A,&b,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity,&FE,&mesh,cq,source2D,0.0);
  //if(dim==3) assemble_global_block(&A,&b,local_assembly_Elasticity,FEM_Block_RHS_Local,&FE,&mesh,cq,source3D,0.0);

  // assemble_global_block(&A,&b,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity_Face,&FE,&mesh,cq,source2D,0.0);
  
  // Boundary Integral <g,r*n>_boundary
  // Flag is which boundary you want to compute this
  // SLEE
  //INT flag0 = 1;
  //INT flag1 = 15;
  //assemble_global_RHS_face_block(&b_bdry,NULL,local_assemble_Elasticity_bdryRHS,&FE,&mesh,cq,source2D,0.0,flag0,flag1);

  //SLEE I don't like this, let me add. 
  // Add RHS vectors together
  //for(i=0;i<FE.ndof;i++) {
  //b.val[i] += b_bdry.val[i];
  //}

  printf("------ 6\n");
  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) {
    eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  }
  if(dim==3) {
    eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);
  }

  /**************************************************/
  //  Apply Pressure "BCs" (removes singularity)
  REAL pressureval =0.5;
  INT pressureloc = 0;
  /**************************************************/

  clock_t clk_assembly_end = clock();
  printf(" --> elapsed CPU time for assembly = %f seconds.\n\n",(REAL)
         (clk_assembly_end-clk_assembly_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /************ Prepare Preconditioners **************************************************/

  // Prepare diagonal blocks
  dCSRmat *A_diag;
  A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));

  for(i=0;i<dim;i++){ // copy block diagonal to A_diag
    dcsr_alloc(A.blocks[i*(dim+2)]->row, A.blocks[i*(dim+2)]->col, A.blocks[i*(dim+2)]->nnz, &A_diag[i]);
    dcsr_cp(A.blocks[i*(dim+2)], &A_diag[i]);
  }

  // Get Mass Matrix for p
  dCSRmat Mp;
  assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  dcsr_cp(&Mp, &A_diag[dim]);
  /*******************************************************************************************/

  /***************** Solve *******************************************************************/
  printf("Solving the System:\n");
  clock_t clk_solve_start = clock();

  INT solver_flag = -20;

  // Allocate solution
  dvector sol = dvec_create(ndof);
  dvector v_ux = dvec_create(FE_ux.ndof);
  dvector v_uy = dvec_create(FE_uy.ndof);
  dvector v_uz;
  if(dim==3) v_uz = dvec_create(FE_uz.ndof);
  dvector v_p = dvec_create(FE_p.ndof);

  // Set initial guess to be all zero
  dvec_set(sol.row, &sol, 0.0);
  // Set initial guess for pressure to match the known "boundary condition" for pressure
  if(dim==2) sol.val[FE_ux.ndof + FE_uy.ndof + pressureloc]  = pressureval;
  if(dim==3) sol.val[FE_ux.ndof + FE_uy.ndof + FE_uz.ndof + pressureloc]  = pressureval;

  // Set parameters for linear iterative methods
  linear_itsolver_param linear_itparam;
  param_linear_solver_init(&linear_itparam);
  param_linear_solver_set(&linear_itparam,&inparam);

  // Solve
  void *numeric=NULL;  // prepare for direct solve:
  numeric=(void *)block_factorize_UMF(&A,0);//inparam.print_level);
  solver_flag=(INT )block_solve_UMF(&A,
				    &sol, // this is the rhs here. 
				    &b,  // this is the solution here. 
				    numeric,
				    0);//     inparam.print_level);
  /* if(dim==2){ */
  /*   if (linear_itparam.linear_precond_type == PREC_NULL) { */
  /*     solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam); */
  /*   } else { */
  /*     solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag); */
  /*   } */
  /* } else if (dim==3) { */
  /*   if (linear_itparam.linear_precond_type == PREC_NULL) { */
  /*     solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam); */
  /*   } else { */
  /*     solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag); */
  /*   } */
  /* } */
  
  // Error Check
  if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

  clock_t clk_solve_end = clock();
  printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
         (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /********************* Compute Errors if you have exact solution ****************************/
  clock_t clk_error_start = clock();
  REAL* solerrL2 = (REAL *) calloc(dim+1, sizeof(REAL));
  REAL* solerrH1 = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  if(dim==2){
    L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
    //if(order_p > 0)
    HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  }
  else if(dim==3){
    L2error_block(solerrL2, sol.val, exact_sol3D, &FE, &mesh, cq, 0.0);
    if(order_p > 0) HDerror_block(solerrH1, sol.val, exact_sol3D, Dexact_sol3D, &FE, &mesh, cq, 0.0);
  }

  REAL uerrL2 = 0;
  REAL uerrH1 = 0;
  for(i=0;i<dim;i++) uerrL2 += solerrL2[i]*solerrL2[i];
  for(i=0;i<dim;i++) uerrH1 += solerrH1[i]*solerrH1[i];
  uerrL2 = sqrt(uerrL2);
  uerrH1 = sqrt(uerrH1);

  REAL perrL2 = solerrL2[dim];
  REAL perrH1 = solerrH1[dim];

  printf("*******************************************************\n");
  printf("L2 Norm of u error    = %26.13e\n",uerrL2);
  //printf("L2 Norm of p error    = %26.13e\n",perrL2);
  printf("H1 Norm of u error    = %26.13e\n",uerrH1);
  //printf("H1 Norm of p error    = %26.13e\n",perrH1);
  printf("*******************************************************\n\n");

  //Jul. 10. 2020 SLEE save the errors for convergence computation
  L2_error_per_cycle[cycle] = uerrL2; 
  H1_error_per_cycle[cycle] = uerrH1; 
  L2_error_p_per_cycle[cycle] = perrL2; 

  //NEW SLEE Aug 23 2020
  //NEW ERROR FOR EG

  REAL* solerrL2_EG = (REAL *) calloc(dim+1, sizeof(REAL));
  //REAL* solerrH1_EG = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  L2error_block_EG(solerrL2_EG, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
  REAL uerrL2_EG = 0;
  //REAL uerrH1_EG = 0;
  for(i=0;i<dim;i++)
    uerrL2_EG += solerrL2_EG[i]*solerrL2_EG[i];
  //for(i=0;i<dim;i++)
  //uerrH1_EG += solerrH1_EG[i]*solerrH1_EG[i];
  uerrL2_EG = sqrt(uerrL2_EG);
  //uerrH1_EG = sqrt(uerrH1_EG);

  L2_EG_error_per_cycle[cycle] = uerrL2_EG;
  printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG);
  
  //NEW SLEE Aug 17 2020
  //NEW ERROR FOR STRESS, \mu \epsilon(u) + \lambda \div u
  REAL* solerr_stress = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  HDsemierror_block_Stress(solerr_stress, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  REAL uerr_stressH1 = 0;
  for(i=0;i<dim;i++)
    uerr_stressH1 += solerr_stress[i]*solerr_stress[i];
  uerr_stressH1 = sqrt(uerr_stressH1);

  H1_stress_error_per_cycle[cycle] = uerr_stressH1; 
  
  FILE *fptmp = NULL;
  /////////////////////////////////////
  fptmp=fopen("output/a.dat","w");
  /* fprintf(stdout,"\n%d,%d\n\n",A.brow,A.bcol); */
  /* for(i=0;i<A.brow;i++){ */
  /*   for(j=0;j<A.brow;j++){ */
  /*     fprintf(stdout,"\n(%d,%d):::%d,%d,%d\n\n",i,j,	\ */
  /* 	      A.blocks[j*A.bcol+i]->row,		\ */
  /* 	      A.blocks[j*A.bcol+i]->col,		\ */
  /* 	      A.blocks[j*A.bcol+i]->nnz); */
  /*   } */
  /* } */
  /* fflush(stdout); */
  bdcsr_print_matlab(fptmp,&A);
  fclose(fptmp);
  /////////////////////////////////////
  
  
  clock_t clk_error_end = clock();
  printf("Elapsed CPU time for getting errors = %lf seconds.\n\n",(REAL)
         (clk_error_end-clk_error_start)/CLOCKS_PER_SEC);


 




  
  /*******************************************************************************************/

  // Plotting
  get_unknown_component(&v_ux,&sol,&FE,0);
  get_unknown_component(&v_uy,&sol,&FE,1);
  if(dim==3) get_unknown_component(&v_uz,&sol,&FE,2);
  get_unknown_component(&v_p,&sol,&FE,dim);

  char** varname;
  if(inparam.print_level > 3){
    char* soldump = "output/solution.vtu";
    varname = malloc(10*FE.nspaces*sizeof(char *));
    varname[0] = "ux";
    varname[1] = "uy";
    if(dim==3) varname[2] = "uz";
    varname[dim] = "p ";
    dump_blocksol_vtk(soldump,varname,&mesh,&FE,sol.val);

    // Print in Matlab format to show vector field in nice way.
    if(dim==3) print_matlab_vector_field(&v_ux,&v_uy,&v_uz,&FE_ux);
  }

  
  /************ Free All the Arrays ***********************************************************/
  // CSR
  bdcsr_free( &A );
  dcsr_free( &Mp);
  for(i=0;i<dim+1;i++)
    dcsr_free( &A_diag[i] );
  if(A_diag) free(A_diag);
  
  // Vectors
  if(solerrL2) free(solerrL2);
  if(solerrL2_EG) free(solerrL2_EG);
  if( solerr_stress ) free( solerr_stress);
  if(solerrH1) free(solerrH1);
  dvec_free( &b );
  dvec_free( &sol );
  dvec_free( &v_ux );
  dvec_free( &v_uy );
  if(dim==3) dvec_free( &v_uz );
  dvec_free( &v_p );
  
  // FE Spaces
  free_fespace(&FE_ux);
  free_fespace(&FE_uy);
  if(dim==3) free_fespace(&FE_uz);
  free_fespace(&FE_p);
  free_blockfespace(&FE);
  
  // Quadrature
  if(cq){
    free_qcoords(cq);
    free(cq);
    cq=NULL;
  }
  
  // Mesh
  free_mesh(&mesh);

  // Strings

  if(inparam.print_level > 3){
    if(varname) free(varname);
  }

 
  /*******************************************************************************************/
  clock_t clk_overall_end = clock();
  printf("\nEnd of Program: Total CPU Time = %f seconds.\n\n",
         (REAL) (clk_overall_end-clk_overall_start)/CLOCKS_PER_SEC);


  }//SLEE :: cycle loop end

  
  // Jul.10.2020 SLEE: Error Computation
  //SLEE compute the convergence rate and print
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
      if(tmp == 0){
	L2_conv_rate_per_cycle[tmp] = 0;
	L2_p_conv_rate_per_cycle[tmp] = 0;
	H1_conv_rate_per_cycle[tmp] = 0;
	H1_stress_conv_rate_per_cycle[tmp] = 0;

      }
      else{
	// multiplied dim since we use DOFs not h here.
	L2_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_per_cycle[tmp]) -log(L2_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) );  
	L2_p_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_p_per_cycle[tmp]) -log(L2_error_p_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) ); 
	H1_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_per_cycle[tmp]) -log(H1_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) ); 
	H1_stress_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_error_per_cycle[tmp]) -log(H1_stress_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) );

	L2_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_EG_error_per_cycle[tmp]) -log(L2_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) );  
	
      }
      
      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_per_cycle[tmp], dof_per_cycle[tmp],L2_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_per_cycle[tmp], dof_per_cycle[tmp],H1_conv_rate_per_cycle[tmp]);
      printf("H1 Stress Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_error_per_cycle[tmp], dof_per_cycle[tmp],H1_stress_conv_rate_per_cycle[tmp]);
      
      printf("L2 EG_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_EG_error_per_cycle[tmp], dof_per_cycle[tmp],
	     L2_EG_conv_rate_per_cycle[tmp]);
 
      printf("----------------------------------------------------\n");      
    }

  
  //for Latex Print Table
  printf("** LATEX TABLE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
      printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\ hline \n",  dof_per_cycle[tmp],  L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_stress_error_per_cycle[tmp], H1_stress_conv_rate_per_cycle[tmp] ); 
    }

  
    return 0;
  
  
}  /* End of Program */
/*********************************************************************************************/
