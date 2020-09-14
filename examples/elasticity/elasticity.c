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

/******/
void block_LocaltoGlobal_neighbor(INT *dof_on_elm,INT *dof_on_elm_neighbor,block_fespace *FE,block_dCSRmat *A,REAL *ALoc)
{
  INT i,j,k,col_a,col_b,acol,block_row,block_col;
  INT local_row,local_col;

  // Loop over all the blocks
  INT nblocks = FE->nspaces;
  INT dof_per_elm_test = 0;
  INT dof_per_elm_trial = 0;
  INT local_row_index = 0;
  INT local_col_index = 0;
  INT global_row_index = 0;
  INT block_dof_per_elm = 0;

  // Get total dof_per_elm for indexing
  for(block_row=0;block_row<nblocks;block_row++) {
    block_dof_per_elm += FE->var_spaces[block_row]->dof_per_elm;
  }

  // Loop through all the blocks
  for(block_row=0;block_row<nblocks;block_row++) {
    dof_per_elm_test = FE->var_spaces[block_row]->dof_per_elm;

    for(block_col=0;block_col<nblocks;block_col++) {
      dof_per_elm_trial = FE->var_spaces[block_col]->dof_per_elm;


      /* Rows of Local Stiffness (test space)*/
      for (i=0; i<dof_per_elm_test; i++) {
        local_row = dof_on_elm[local_row_index+i];
      
        /* Columns of Local Stiffness (trial space)*/
        for (j=0; j<dof_per_elm_trial; j++) {
          local_col = dof_on_elm_neighbor[local_col_index + j];
          /* Columns of A */
          if(A->blocks[block_row*nblocks+block_col]) {
            col_a = A->blocks[block_row*nblocks+block_col]->IA[local_row];
            col_b = A->blocks[block_row*nblocks+block_col]->IA[local_row+1];
            for (k=col_a; k<col_b; k++) {
              acol = A->blocks[block_row*nblocks+block_col]->JA[k];
              if (acol==local_col) {	/* If they match, put it in the global matrix */
                A->blocks[block_row*nblocks+block_col]->val[k]
		  += ALoc[(local_row_index+i)*block_dof_per_elm+(local_col_index+j)];
              }
            }
          }
        }
      }
      local_col_index += dof_per_elm_trial;
    }
    local_col_index = 0;
    global_row_index += FE->var_spaces[block_row]->ndof;
    local_row_index += dof_per_elm_test;
  }
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

/******************************************************************************************************/
/*!
* \fn assemble_global_block(block_dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,REAL),void (*local_rhs_assembly)(REAL *,block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the global stiffness BLOCK matrix and rhs for any a(u,v) = <f,v> bilinear form using various element
*        types (eg. P0, P1, P2, Nedelec, and Raviart-Thomas).
*        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
*        This allows for several matrices to be assembled then added or concatenated together.
*
*        For this problem we compute:
*
*        Lu = f  ---->   a(u,v) = <f,v>
*
*        which gives Ax = b,
*
*        A_ij = a( phi_j, psi_i)
*        b_i  = <f,psi_i>
*
* \note All matrices are assumed to be blocks and indexed at 0 in the CSR formatting.
*
* \param local_assembly     Routine to get local matrices
* \param local_rhs_assembly Routine to get local rhs vectors
* \param FE                 block FE Space
* \param mesh               Mesh Data
* \param cq                 Quadrature Nodes
* \param rhs                Routine to get RHS function (NULL if only assembling matrix)
* \param time               Physical Time if time dependent
*
* \return A                 Global stiffness BLOCK CSR matrix
* \return b                 Global RHS vector
*
*/
void assemble_global_block_neighbor(block_dCSRmat* A,dvector *b,void (*local_assembly)(block_dCSRmat* A,dvector *b,REAL *,block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,REAL,INT*),void (*local_rhs_assembly)(REAL *,block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dof_per_elm = 0;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,k,testdof,trialdof;

  // Get block data first
  INT nblocks = A->brow;
  // Check for errors
  if(nblocks!=A->bcol) {
    printf("Your block matrix is not square.  It is an %d x %d matrix.\n\n",A->brow,A->bcol);
    exit(0);
  }
  if(nblocks!=FE->nspaces) {
    printf("You have %d FEM spaces, but only %dx%d blocks.  They must be consistent.\n\n",FE->nspaces,A->brow,A->bcol);
    exit(0);
  }
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
    dvec_set(b->row,b,0.0);
  }

  // Loop over each block and build sparsity structure of matrices
  for(i=0;i<nblocks;i++) {
    for(j=0;j<nblocks;j++) {
      testdof = FE->var_spaces[i]->ndof;
      trialdof = FE->var_spaces[j]->ndof;
      if(A->blocks[i*nblocks+j]) {
        A->blocks[i*nblocks+j]->row = testdof; // test functions
        A->blocks[i*nblocks+j]->col = trialdof; // trial functions
        A->blocks[i*nblocks+j]->IA = (INT *) calloc(testdof+1,sizeof(INT));

        // Get Sparsity Structure First
        // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
        create_CSR_rows_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);

        // Columns of A -> JA
        A->blocks[i*nblocks+j]->JA = (INT *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(INT));
        create_CSR_cols_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);

        // Set values
        A->blocks[i*nblocks+j]->val = (REAL *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(REAL));
        for (k=0; k<A->blocks[i*nblocks+j]->nnz; k++) {
          A->blocks[i*nblocks+j]->val[k] = 0;
        }
      }
    }
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  // Now Build Global Matrix entries

  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL)
  bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  // Loop over elements

  INT* switch_on_face = (INT *) calloc(mesh->nface,sizeof(INT));
  for (INT mm=0; mm<mesh->nface; mm++) {
    switch_on_face[mm] =0.;
  }

  
  for (i=0; i<mesh->nelm; i++) {
    // Zero out local matrices
    for (j=0; j<local_size; j++) {
      ALoc[j]=0;
    }
    if(rhs!=NULL) {
      for (j=0; j<dof_per_elm; j++) {
        bLoc[j]=0;
      }
    }

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(k=0;k<nblocks;k++) {
      rowa = FE->var_spaces[k]->el_dof->IA[i];
      rowb = FE->var_spaces[k]->el_dof->IA[i+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find vertices for given Element
    get_incidence_row(i,mesh->el_v,v_on_elm);


    
    
    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(A,b,ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,time,switch_on_face);
    if(rhs!=NULL)
    (*local_rhs_assembly)(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    block_LocaltoGlobal(dof_on_elm,FE,b,A,ALoc,bLoc);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}


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
 
  coordinates *barycenter = allocatecoords(1,dim);
  barycenter->x[0] = mesh->el_mid[elm*dim];
  barycenter->y[0] = mesh->el_mid[elm*dim + 1];

  
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

    /*
    printf("----  RHS  ------------------ \n");
    for(INT mm=0; mm <  sizeof(dof_on_elm); ++mm)
      {
	printf("0. :: %d :: dof_on_elm[%d] = %u\n", elm,  mm, dof_on_elm[mm]);
	//printf("3. :: %d :: local_dof_on_elm[%d] = %u\n", elm,  mm, *local_dof_on_elm);
	//local_dof_on_elm++;
      }    	
    printf("----------------------\n");
    */
    
    for(i=0;i<FE->nspaces;i++) {

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);
      
      // Loop over test functions and integrate rhs
      if(FE->var_spaces[i]->FEtype<20) { // Scalar Element

	for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {

	  //SLEE
	  if(i == 0 || i == 1)
	    {
	      //This is for  u0 and u1 (CG part)
	      //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      bLoc[(local_row_index+test)] += w*rhs_val[unknown_index]*FE->var_spaces[i]->phi[test];
	    }
	    else if(i == 2)
	    {
	      //SLEE
	      // This is for  u2:: (EG part)
	      //Note that new DOF is \phi^3 = [x ; y] 
	      //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	      bLoc[(local_row_index+test)] +=w*(rhs_val[0]*  (qx[0]-barycenter->x[0])  +rhs_val[1]* (qx[1]-barycenter->y[0]));
	    }

	  //sanity check
	  if(unknown_index != i)
	    {
	      printf("unknown index != i \n");
	      exit(0);
	    }
        }
        unknown_index++;
	
      }
      else { 
       	  //SLEE
	  printf("NOT HERE : %d\n",i);
	  exit(0);
      }

      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
      local_row_index  += FE->var_spaces[i]->dof_per_elm;


      
    }
  }

  // SLEE
  // NOW FOR FACES (at BOUNDARY)
  
  INT* local_dof_on_elm_face = NULL;
 
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

  //Sanity Check // SLEE
  //for(i=0;i<FE->nspaces;i++) {
  //printf("dof_per_face = %d,   dof_per_face_blk[%d] = %d \n", dof_per_face, i, dof_per_face_blk[i]);
  //}
  //exit(0);
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  REAL *data_face_end=finrm + dim; //end
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
  REAL fiarea=0e0;
  INT jk,k,face,quad_face,rowa, rowb, jcntr,ed, jkl, j;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;
  
  // SLEE
  // Get the BD values 
  REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));

  bool bWeakBC_RHS = true; //false;
  if(bWeakBC_RHS){
    
    for(jk=mesh->el_f->IA[elm];jk<mesh->el_f->IA[elm+1];jk++){
      
      //printf("jk = %d, mesh->el_f->IA[element] = %d, mesh->el_f->IA[element+1] = %d  \n",
      //     jk, mesh->el_f->IA[elm], mesh->el_f->IA[elm+1]);
      //  j = face0, face1, face2
      j=jk - mesh->el_f->IA[elm];
      // face is the global number
      face=mesh->el_f->JA[jk];
      
      // elements
      // 0 -> is the interface
      if(mesh->f_flag[face]>0) { 	
	//printf("j = %d, jk = %d,  mesh->el_f->IA[element] = %d \n",j, jk, mesh->el_f->IA[elm]);
	//printf("face = %d, mesh->el_f->JA[jk] = %d \n", face, mesh->el_f->JA[jk]);
	
	// Get normal vector values. 
	finrm[0]=mesh->f_norm[face*dim+0];
	if(dim>1)
	  finrm[1]=mesh->f_norm[face*dim+1];
	if(dim>2)
	  finrm[2]=mesh->f_norm[face*dim+2];
	
	for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
	  
	  //printf("** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);
	  
	  j=jkl-mesh->f_v->IA[face];
	  k=mesh->f_v->JA[jkl];
	  
	  //printf("** j = %d, k = %d \n", j, k );
	  
	  xfi[j*dim+0]=mesh->cv->x[k];
	  if(dim>1)
	    xfi[j*dim+1]=mesh->cv->y[k];
	  if(dim>2)
	    xfi[j*dim+2]=mesh->cv->z[k];	  
	  
	  //printf("** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);
	  
	}

	// Only grab the faces on the flagged boundary
	//NEED TO DEBUG!!!!
	//SLEE // SEEK BD flag...?????
	fiarea=mesh->f_area[face];
	//nq1d_face == 3, 3 quad points. 
	zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

	REAL edge_length = mesh->ed_len[face];
	//if(edge_length != fiarea)
	  {
	    printf("face=%d, elm=%d: length = %f, area = %f \n", face,elm,edge_length, fiarea );
	    //   exit(0);
	  }
	double penalty_term =  100./fiarea;
	
      	
	for (quad_face=0;quad_face<nquad_face;quad_face++) {
	  
	  qx_face[0] = cq_face->x[quad_face];
	  qx_face[1] = cq_face->y[quad_face];	  

	  if(dim==3) qx_face[2] = cq_face->z[quad_face];
	  
	  w_face = cq_face->w[quad_face];
	  
	  (*exact_sol2D)(val_true_face,qx_face,time,
			 &(mesh->el_flag[elm]));  // ???      

	  local_row_index=0;
	  unknown_index=0;  
	  local_dof_on_elm_face = dof_on_elm;

	  /*
	  printf("----  RHS  FACE ------------------ \n");
	  for(INT mm=0; mm <  sizeof(local_dof_on_elm_face); ++mm)
	    {
	      printf("0. :: %d :: dof_on_elm[%d] = %u\n", elm,  mm, local_dof_on_elm_face[mm]);
	      //printf("3. :: %d :: local_dof_on_elm[%d] = %u\n", elm,  mm, *local_dof_on_elm);
	      //local_dof_on_elm++;
	    }    	
	  printf("----------------------\n");
	  */
	  
	  for(i=0;i<FE->nspaces;i++) {

	    /*
	    PX_H1_basis(FE->var_spaces[i]->phi,
			FE->var_spaces[i]->dphi,
			qx_face,
			local_dof_on_elm_face,
			FE->var_spaces[i]->FEtype,
			mesh);
	    */
	    get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,
			  qx_face,
			  v_on_elm,
			  local_dof_on_elm_face,
			  mesh,FE->var_spaces[i]);
	    
	    for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
	      //SLEE
	      if(i==0 || i ==1)
		{
		  //printf("quad_face = %d :: i = %d,  local_row_index = %d,  qx_face[0] = %f, qx_face[1] = %f, val_true_face[%d] = %f, test = %d \n",
		  //	 quad_face, i,  local_row_index,  qx_face[0], qx_face[1] , unknown_index,  val_true_face[unknown_index], test);
		  bLoc[(local_row_index+test)] += penalty_term * w_face*
		    (val_true_face[unknown_index]*  FE->var_spaces[i]->phi[test]);
		}
	      
	      else if(i == 2)
		{
		  //SLEE
		  // Note that new DOF is \phi^3 = [x ; y] 
		  //printf("i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
		  bLoc[(local_row_index+test)] +=
		    penalty_term *  w_face * (val_true_face[0] * (qx_face[0] - barycenter->x[0])
					      + val_true_face[1]*  (qx_face[1] - barycenter->y[0]));
		}
	      //sanity check
	      if(unknown_index != i)
		{
		  printf("unknown index != i \n");
		  exit(0);
		}
	      //SLEE
	      
	    }
	    unknown_index++;
	    
	    local_dof_on_elm_face += FE->var_spaces[i]->dof_per_elm;	
	    local_row_index += FE->var_spaces[i]->dof_per_elm; 
	    
	  } // i = 0
	   
	} //face
	
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
  
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));

  // Loop over all Elements 
  for (elm=0; elm<mesh->nelm; elm++) {

    barycenter->x[0] = mesh->el_mid[elm*dim];
    barycenter->y[0] = mesh->el_mid[elm*dim + 1];


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

      // Get True Solution at Quadrature Nodes
      (*truesol)(val_true,qx,time,&(mesh->el_flag[elm]));

      
      INT *local_dof_on_elm = dof_on_elm;
      INT i,j,dof;
      REAL u0_value_at_q = 0.;
      REAL u1_value_at_q = 0.;

      //REAL u2_value_at_q = 0.;

      REAL* u_comp = u;
      
      // Interpolate FE solution to quadrature point
      //blockFE_Interpolation(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);
      get_FEM_basis(FE->var_spaces[0]->phi,
		    FE->var_spaces[0]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[0]);

      for(j=0; j<3; j++) {
	dof = local_dof_on_elm[j];
	//printf("--u1  dof = %d, u _comp= %f, u0 = %f  \n", dof, u_comp[dof],FE->var_spaces[0]->phi[j]);
	u0_value_at_q += u_comp[dof]*FE->var_spaces[0]->phi[j];
      }


      //printf("--1) local_dof_on_elm[0]= %d  \n", local_dof_on_elm[0]);

      
      local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
      u_comp += FE->var_spaces[0]->ndof;

      
      //printf("--2) local_dof_on_elm[0]= %d  \n", local_dof_on_elm[0]);

      
      get_FEM_basis(FE->var_spaces[1]->phi,
		    FE->var_spaces[1]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[1]);
      
      
      for(j=0; j<3; j++) {
	dof =  local_dof_on_elm[j];
       	//printf("--u1  dof = %d, u _comp= %f, u1 = %f  \n", dof, u_comp[dof],FE->var_spaces[1]->phi[j]);
	u1_value_at_q += u_comp[dof]*FE->var_spaces[1]->phi[j];
      }

      // EG
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      u_comp += FE->var_spaces[1]->ndof;

      REAL val[dim];      
      val[0] = 0.0;
      val[1] = 0.0;

      //printf("--3) local_dof_on_elm[0]= %d  \n", local_dof_on_elm[0]);
      //printf("----\n");
   
      dof = local_dof_on_elm[0]; // get the dof for the last component
      //printf("--u2 dof = %d, u _comp= %f, q[0] = %f, q[1] = %f,   \n", dof, u_comp[dof],qx[0] - barycenter->x[0],qx[1] - barycenter->x[1]);
      //printf("----\n");
   
      val[0] = u_comp[dof]*   (qx[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      val[1] = u_comp[dof]*   (qx[1] - barycenter->y[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point

      
      local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
      u_comp += FE->var_spaces[2]->ndof;
      
      // Compute Square of Error on Element for each component of FE space


      err[0] += w*(ABS(u0_value_at_q  +val[0]- val_true[0]))*(ABS(u0_value_at_q +val[0] - val_true[0]));
      err[1] += w*(ABS(u1_value_at_q  +val[1]- val_true[1]))*(ABS(u1_value_at_q +val[1] - val_true[1]));

      //err[0] += w*(ABS(val[0]- val_true[0]))*(ABS(val[0] - val_true[0]));
      //err[1] += w*(ABS(val[1]- val_true[1]))*(ABS(val[1] - val_true[1]));


      //err[0] += w*(ABS(u0_value_at_q  - val_true[0]))*(ABS(u0_value_at_q  - val_true[0]));
      //err[1] += w*(ABS(u1_value_at_q  - val_true[1]))*(ABS(u1_value_at_q  - val_true[1]));
    
      //blockFE_Interpolation_new(val_sol,u,qx,dof_on_elm,v_on_elm,FE,mesh);

    }//quad
  }


  printf("---------------------1\n");
  for(i=0;i<nspaces-1;i++) {
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


void HDerror_block_EG
//(REAL *err,REAL *u,void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)
(REAL *err,REAL *u,void (*truesol)(REAL *,REAL *,REAL,void *),void (*D_truesol)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,REAL time)

{
	printf("============== HI\n");

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

  
   
      INT* local_dof_on_elm = dof_on_elm;
      REAL* u_comp = u;
      INT dof,j,k,i;
      
      get_FEM_basis(FE->var_spaces[0]->phi,
		    FE->var_spaces[0]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[0]);

      REAL grad_val0[dim];
      grad_val0[0]=0;
      grad_val0[1]=0;
      REAL grad_val1[dim];
      grad_val1[0]=0;
      grad_val1[1]=0;

      ////
      /// GET the values 
      dof = local_dof_on_elm[0];
      grad_val0[0] += u_comp[dof]*FE->var_spaces[0]->dphi[0];
      
      dof = local_dof_on_elm[1];
      grad_val0[0] += u_comp[dof]*FE->var_spaces[0]->dphi[2];
      
      dof = local_dof_on_elm[2];
      grad_val0[0] += u_comp[dof]*+FE->var_spaces[0]->dphi[4];
      
      
      dof = local_dof_on_elm[0];
      grad_val0[1] += u_comp[dof]*FE->var_spaces[0]->dphi[1];
      
      dof = local_dof_on_elm[1];
      grad_val0[1] += u_comp[dof]* FE->var_spaces[0]->dphi[3];
      
      dof = local_dof_on_elm[2];
      grad_val0[1] += u_comp[dof]*FE->var_spaces[0]->dphi[5];
      //////////
      

      /// Increments	
      local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
      u_comp += FE->var_spaces[0]->ndof;

      get_FEM_basis(FE->var_spaces[1]->phi,
		    FE->var_spaces[1]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[1]);
      
      ////
      /// GET the values 
      dof = local_dof_on_elm[0];
      grad_val1[0] += u_comp[dof]*FE->var_spaces[1]->dphi[0];

      dof = local_dof_on_elm[1];
      grad_val1[0] += u_comp[dof]*FE->var_spaces[1]->dphi[2];

      dof = local_dof_on_elm[2];
      grad_val1[0] += u_comp[dof]*+FE->var_spaces[1]->dphi[4];
      
      dof = local_dof_on_elm[0];
      grad_val1[1] += u_comp[dof]*FE->var_spaces[1]->dphi[1];

      dof = local_dof_on_elm[1];
      grad_val1[1] += u_comp[dof]* FE->var_spaces[1]->dphi[3];

      dof = local_dof_on_elm[2];
      grad_val1[1] += u_comp[dof]*FE->var_spaces[1]->dphi[5];

      /// Increments
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      u_comp += FE->var_spaces[1]->ndof;


      ////
      /// GET the values 
      // bEG
      REAL grad_val2[dim];      
      grad_val2[0] = 0.0;
      grad_val2[1] = 0.0;

      REAL grad_val3[dim];      
      grad_val3[0] = 0.0;
      grad_val3[1] = 0.0;
    
      dof = local_dof_on_elm[0]; // get the dof for the last component

      //printf("-- dof = %d, u _comp= %f  \n", dof, u_comp[dof]);

      grad_val2[0] = u_comp[dof]*1.;
      grad_val2[1] = 0.;

      
      grad_val3[0] = 0.;
      grad_val3[1] = u_comp[dof]*1.;

      
      //val[0] = u_comp[dof]*   (qx[0] - barycenter->x[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      //val[1] = u_comp[dof]*   (qx[1] - barycenter->y[0]); //FE->phi[j]; -> this basis function is a vector <x,y> - the quadrature point
      
      
      local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
      u_comp += FE->var_spaces[2]->ndof;


      
      err[0]+=w*(ABS(grad_val0[0] - val_true[0]))*(ABS(grad_val0[0] - val_true[0]));
      err[0]+=w*(ABS(grad_val1[0] - val_true[2]))*(ABS(grad_val1[0] - val_true[2]));

      // bEG
      err[0]+=w*(ABS(grad_val2[0] - val_true[0]))*(ABS(grad_val2[0] - val_true[0]));
      err[0]+=w*(ABS(grad_val3[0] - val_true[2]))*(ABS(grad_val3[0] - val_true[2]));
    
      
      err[1]+=w*(ABS(grad_val0[1] - val_true[1]))*(ABS(grad_val0[1] - val_true[1]));
      err[1]+=w*(ABS(grad_val1[1] - val_true[3]))*(ABS(grad_val1[1] - val_true[3]));

      // bEG
      err[1]+=w*(ABS(grad_val2[1] - val_true[1]))*(ABS(grad_val2[1] - val_true[1]));
      err[1]+=w*(ABS(grad_val3[1] - val_true[3]))*(ABS(grad_val3[1] - val_true[3]));
    
      
    }
  }

  //exit(0);

  for(i=0;i<nspaces;i++) {
    if(err[i] < 0)
      {
	printf("minus ? \n");
      }
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

void local_assembly_Elasticity(block_dCSRmat* A,dvector *b,REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time, INT* switch_on_face)
//void local_assembly_Elasticity(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)

{

  bool bEG = true;//true;//true;//true;//false;//true;//false;//false;//true;//false;//true;//false;//true;//false;
 
  // Loop indices
  INT i,j,idim,quad,test,trial;
  INT dim = mesh->dim;
  // Mesh and FE data
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  
  //printf("CHECK POINT: FE->nspaces = %d, dof_per_elem = %d, FE->var_spaces[0]->dof_per_elm = %d,  FE->var_spaces[1]->dof_per_elm = %d ,  FE->var_spaces[1]->dof_per_elm = %d \n",
  //	 FE->nspaces, dof_per_elm, FE->var_spaces[0]->dof_per_elm,FE->var_spaces[1]->dof_per_elm, FE->var_spaces[2]->dof_per_elm);
  //exit(0);	 
 
  INT *local_dof_on_elm = NULL; 
  INT *local_dof_on_elm_face = NULL;
  INT *local_dof_on_elm_face_interface = NULL;
  INT *local_dof_on_elm_neighbor = NULL;
  //printf("dof per elm = %d" , dof_per_elm);
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));
       
  //SLEE
  coordinates *barycenter = allocatecoords(1,dim);
  coordinates *barycenter_neighbor = allocatecoords(1,dim);
  barycenter->x[0] = mesh->el_mid[elm*dim];
  barycenter->y[0] = mesh->el_mid[elm*dim + 1];

  printf("============================================================\n");
  printf("ELEMENT = %d, x = %f , y = %f \n", elm,  barycenter->x[0],   barycenter->y[0]);
  INT local_size = dof_per_elm*dof_per_elm;  
  REAL* ALoc_neighbor = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  
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

    get_FEM_basis(FE->var_spaces[1]->phi,
		  FE->var_spaces[1]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[1]);
  
    // p
    local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;

    //DEBUG
    local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;

     
    // u1-v1 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	//printf("quad = %d ::   local_row_index = %d,  qx[0] = %f, qx[1] = %f,  test = %d \n",
	//     quad,  local_row_index,  qx[0], qx[1] ,   test);

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
    //Local Matrix 
    if(bEG){
      // exit(0);
      //EG PART
      
      //NEW SLEE :
      // u0- q block < 2mu e(u1) * e(q) > + lambda <div u1 : e(q) >
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	  //u0 - q0
	  kij  = 2.0*FE->var_spaces[0]->dphi[trial*dim] * 1.; //FE->var_spaces[dim]->dphi[test*dim];
	  //u0 - q1
	  kij += 0.;
	  
	  // Divergence
	  // u0-q0 block : <dx(u1),dx(v1)>
	  kij += lambda*FE->var_spaces[0]->dphi[trial*dim]* 1.;//FE->var_spaces[dim]->dphi[test*dim];	

	  // u0-q1
	  kij += lambda*FE->var_spaces[0]->dphi[trial*dim]* 1.;
	  
	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
      
      // u1-q block
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	  //u1 - q0
	  // = 0 
	  //u1 - q1
	  kij = 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;
	  
	  // NEW SLEE : elasticity div part
	  // u2-v1 block : <dy(u2),dx(v1)> 
	  kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;	
	  kij += lambda*FE->var_spaces[1]->dphi[trial*dim+1] * 1.;
	  
	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
      
      
      // q-v0 block: 
      local_row_index = 0;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	  //q0 - v0
	  kij = 2.0* 1. * FE->var_spaces[0]->dphi[test*dim]; //FE->var_spaces[dim]->dphi[test*dim];
	  //q1 - v0
	  // = 0.;
	  
	  // Divergence
	  // u1-q block : <dx(u1),dx(v1)>
	  kij += lambda * 1. * FE->var_spaces[0]->dphi[test*dim];//FE->var_spaces[dim]->dphi[test*dim];
	  kij += lambda * 1. * FE->var_spaces[0]->dphi[test*dim];
	  
	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
      
      
      // q-v1 block: 
      local_row_index = FE->var_spaces[0]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

	  //q0 - v1
	  //q1 - v1
	  kij = 2.0* 1. * FE->var_spaces[1]->dphi[test*dim+1];
	  
	  // NEW SLEE : elasticity div part
	  // u2-v1 block : <dy(u2),dx(v1)> 
	  kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;	
	  kij += lambda*FE->var_spaces[1]->dphi[test*dim+1] * 1.;
	  
	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
      
      
    
      
      
    } // bEG



    
  }//QUAD
  
  INT  dof_per_face = 0;
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
  
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  REAL *data_face_end=finrm + dim; //end
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
  qcoordinates *cq_face_neighbor = allocateqcoords_bdry(nq1d_face,1,dim,2);

  REAL fiarea=0e0;
  INT jk,k,face,quad_face,rowa, rowb, jcntr,ed, jkl;
  INT maxdim=4;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;
  
  iCSRmat *f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);
  INT pq, nbr, nbr0;
  INT facein=0; //num faces interior
  INT facebd=0;// num faces on boundary
  

  //LEE : WE only need to loop over faces only on th mesh
  
  for(jk=mesh->el_f->IA[elm];jk<mesh->el_f->IA[elm+1];jk++){
      
    //printf("M* jk = %d, mesh->el_f->IA[element] = %d, mesh->el_f->IA[element+1] = %d  \n",jk, mesh->el_f->IA[elm], mesh->el_f->IA[elm+1]);

 local_dof_on_elm_face = NULL;
 local_dof_on_elm_face_interface = NULL;
 local_dof_on_elm_neighbor = NULL;

    
    //  j = face0, face1, face2
    j=jk - mesh->el_f->IA[elm];
    // face is the global number
    face=mesh->el_f->JA[jk];

    printf("** \n");
    printf("switch_on_face[%d] = %d \n", face,  switch_on_face[face]); 
    switch_on_face[face] += 1.;
    printf("switch_on_face[%d] = %d \n", face,  switch_on_face[face]);      
    // elements
    // 0 -> is the interface
    bool bWeakBC_Matrix = true;
    
    if(bWeakBC_Matrix){

      //neighbor indicator
      nbr=-1;
      
      for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){	
	
	nbr0=f_el->JA[pq];
	
	if(nbr0==elm)
	  continue;
	
	nbr=nbr0;
      }
      
      // Get normal vector values. 
      finrm[0]=mesh->f_norm[face*dim+0];
      if(dim>1)
	finrm[1]=mesh->f_norm[face*dim+1];
      if(dim>2)
	finrm[2]=mesh->f_norm[face*dim+2];


      //DEBUG
      if(switch_on_face[face] == 2)
	{
	  printf("switch the sign of the normal vector \n"); 
	  finrm[0] *= -1;
	  finrm[1] *= -1;
	}
      
      //Get the area of face
      fiarea=mesh->f_area[face];
      double penalty_term = 100./fiarea;


      // To get xfi 'coordinates' for quadrature points
      for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
	
	//printf("M** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);
	
	j=jkl-mesh->f_v->IA[face];
	k=mesh->f_v->JA[jkl];
	
	//printf("M** j = %d, k = %d \n", j, k );
	
	xfi[j*dim+0]=mesh->cv->x[k];
	if(dim>1)
	  xfi[j*dim+1]=mesh->cv->y[k];
	if(dim>2)
	  xfi[j*dim+2]=mesh->cv->z[k];	  
	
	//printf("M** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);
      }
      
      
      //nq1d_face == 3, 3 quad points. 
      zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

      /////////////////////////////////////////////////
      if(nbr<0) {

	printf("-----------------\n");
	printf("(BD) ELEMENT = %d - Face = %d - Neighbor = %d,  \n", elm,  face, nbr);
	printf("BOUNDARY = %d,  \n", nbr);
	printf("N[0] = %f, N[1] = %f \n", finrm[0], finrm[1]);
	printf("-----------------\n");

	
	// boundary
	facebd++;
	//printf("%d: face %d in element %d is on the boundary\n",
	
	//Sanity Check
	if(mesh->f_flag[face] == 0)
	  {printf("check-error0\n"); exit(0);}
	
	
	for (quad_face=0;quad_face<nquad_face;quad_face++) {
	  
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
	  
	  local_dof_on_elm_face = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
	  
	  
	  get_FEM_basis(FE->var_spaces[1]->phi,
			FE->var_spaces[1]->dphi,
			qx_face,
			v_on_elm,
			local_dof_on_elm_face,
			mesh,
			FE->var_spaces[1]);
	  
	  local_dof_on_elm_face += FE->var_spaces[1]->dof_per_elm;


	  //DEBUG
	  local_dof_on_elm_face += FE->var_spaces[2]->dof_per_elm;

	  
	  // u0-v0 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
	  local_row_index = 0;
	  local_col_index = 0;
	  
	  // u0-v0 block: 
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      kij = -2.0 * FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * FE->var_spaces[0]->phi[test];
	      kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] *    FE->var_spaces[0]->phi[test];
	      
	      kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[0] * FE->var_spaces[0]->phi[test];
	      
	      //penalty term
	      kij += penalty_term * FE->var_spaces[0]->phi[trial] * FE->var_spaces[0]->phi[test];
	      
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face *kij;
	    }
	  }
	  
	  // u0-v1 block <dy(u1),dx(v2)>
	  local_row_index = FE->var_spaces[0]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      kij = -FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[1]->phi[test];
	      
	      kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[1] * FE->var_spaces[1]->phi[test];
	      
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-v0 block : 
	  local_row_index = 0;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      
	      kij = - FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * FE->var_spaces[0]->phi[test];
	      
	      kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[0]->phi[test];
	      
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-v1 block 
	  local_row_index = FE->var_spaces[0]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      
	      kij = -2.0 * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test] ;
	      kij -= 1.0 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0]   * FE->var_spaces[1]->phi[test];
	      
	      // NEW SLEE : elasticity div part
	      // u2-v2 block : <dy(u2),dx(v2)> 
	      kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test];
	      
	      //penalty term
	      kij += penalty_term *  FE->var_spaces[1]->phi[trial]  * FE->var_spaces[1]->phi[test];
	      
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
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
	      
	      kij = -2.0 * 1. * finrm[0] * (qx_face[0] - barycenter->x[0]);
	      kij -= 2.0 * 1. * finrm[1] * (qx_face[1] - barycenter->y[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);	
	      kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	      
	      kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	      
	      
	      //penalty term
	      kij += penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	      kij += penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);
	      
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // Matrix Face Boundary
	  if(bEG){
	    //EG PART 
	    //NEW SLEE :
	    // u0- q block 
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	    local_col_index = 0;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

		//u0-q0
		//u0-q1
		kij = -2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
		kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
		
		// Divergence
		// u1-q block : <dx(u1),dx(v1)>
		kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
		kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
		
		//penalty term
		kij += penalty_term * FE->var_spaces[0]->phi[trial] * (qx_face[0]- barycenter->x[0]);
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    // u1-q block
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	    local_col_index = FE->var_spaces[0]->dof_per_elm;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

		//u1-q0
		//u1-q1
		kij = -2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
		kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
	     
		kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
		kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);
		
		
		//penalty term
		kij += penalty_term * FE->var_spaces[1]->phi[trial] * (qx_face[1]- barycenter->y[0]);
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    
	    // q-v0 block: 
	    local_row_index = 0;
	    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

		//q0 - v0
		//q1 - v0
		kij = -2. * 1. * finrm[0] * FE->var_spaces[0]->phi[test];
		
		// Divergence
		// u1-q block : <dx(u1),dx(v1)>
		kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];//FE->var_spaces[dim]->dphi[test*dim];
		kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];
		
		
		//penalty term
		kij += penalty_term * FE->var_spaces[0]->phi[test] * (qx_face[0]- barycenter->x[0]);
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    
	    // q-v1 block: 
	    local_row_index = FE->var_spaces[0]->dof_per_elm;
	    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){

		//q0-v1
		//q1-v1
		kij = -2.0* 1. * finrm[1] * FE->var_spaces[1]->phi[test];
		
		// NEW SLEE : elasticity div part
		// u2-v1 block : <dy(u2),dx(v1)> 
		kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];	
		kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];
		
		//penalty term
		kij += penalty_term * FE->var_spaces[1]->phi[test] * (qx_face[1]- barycenter->y[0]);
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	  }//bEG

	  
	  
	  
	}//quadface
	
      } // if boundary	if(nbr<0) {
      
      //else if(mesh->f_flag[face] == 0)
      // INTERFACE
      else {

	//ALoc_neighbor = NULL;
	for (j=0; j<local_size; j++) {
	  ALoc_neighbor[j]=0;
	}

	local_dof_on_elm_face_interface = NULL;
	local_dof_on_elm_neighbor = NULL;
	
	if(nbr>elm)
	  facein++; // only count each face once, just for counting
	
	//Sanity Cehck
	if(mesh->f_flag[face] != 0)
	  {printf("well? \n"); exit(0);}

	barycenter_neighbor->x[0] = mesh->el_mid[nbr*dim];
	barycenter_neighbor->y[0] = mesh->el_mid[nbr*dim + 1];

	printf("-----------------\n");
	printf("ELEMENT = %d - Face = %d - Neighbor = %d,  \n", elm,  face, nbr);
	printf("NEIGHBOR = %d, x = %f , y = %f \n", nbr,  barycenter_neighbor->x[0],   barycenter_neighbor->y[0]);
	printf("N[0] = %f, N[1] = %f \n", finrm[0], finrm[1]);
	printf("-----------------\n");
	
	for (quad_face=0;quad_face<nquad_face;quad_face++) {

	  //printf("--------------\n");
	  //printf("quad_face = %d, nquad-face = %d, \n", quad_face, nquad_face);
	  
	  //neighbor, we can stil share the quadrature points, where the funtions are evaluated..
	  qx_face[0] = cq_face->x[quad_face];
	  qx_face[1] = cq_face->y[quad_face];
	  if(dim==3) qx_face[2] = cq_face->z[quad_face];
	  
	  w_face = cq_face->w[quad_face];


	  // Again, get the basis functions for the element
	  // Not for the neighbor yet
	  get_FEM_basis(FE->var_spaces[0]->phi,
			FE->var_spaces[0]->dphi,
			qx_face,
			v_on_elm,
			dof_on_elm,
			mesh,
			FE->var_spaces[0]);
	  
	  local_dof_on_elm_face_interface = dof_on_elm + FE->var_spaces[0]->dof_per_elm;

	  
	  get_FEM_basis(FE->var_spaces[1]->phi,
			FE->var_spaces[1]->dphi,
			qx_face,
			v_on_elm,
			local_dof_on_elm_face_interface,
			mesh,
			FE->var_spaces[1]);
	  
	  local_dof_on_elm_face_interface += FE->var_spaces[1]->dof_per_elm;

	  //DEBUG
	  local_dof_on_elm_face_interface += FE->var_spaces[2]->dof_per_elm;

	  /////////////
	  /// SLEE
	  // FOR NEIGHBOR
	  // Find DOF for given Element
	  //get_incidence_row(nbr,FE->var_spaces[i]->el_dof,dof_on_elm_neighbor);
	  // Find DOF for given Element
	  // Note this is "local" ordering for the given FE space of the block
	  // Not global ordering of all DOF
	  INT rowa_neighbor,rowb_neighbor,jcntr_neighbor;
	  //NT nblocks = 3;
	  jcntr_neighbor = 0;
	  
	  for(INT k_neighbor=0;k_neighbor<FE->nspaces;k_neighbor++) {
	    
	    rowa_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[nbr];
	    rowb_neighbor = FE->var_spaces[k_neighbor]->el_dof->IA[nbr+1];
	    
	    //printf(" rowa_neighbor = %d,  rowb_neighbor = %d \n",  rowa_neighbor,  rowb_neighbor);
	    
	    for (INT j_neighbor=rowa_neighbor; j_neighbor<rowb_neighbor; j_neighbor++) {
	      
	      dof_on_elm_neighbor[jcntr_neighbor] = FE->var_spaces[k_neighbor]->el_dof->JA[j_neighbor];
	      
	      jcntr_neighbor++;
	     
	      
	    }
	  }
	  
	  local_dof_on_elm_neighbor = dof_on_elm_neighbor;


	  for(INT zz=0; zz<dof_per_elm; ++zz)
	    {
	      printf("dof_on_elm[%d] = %d \n", zz,dof_on_elm[zz]);
	      printf("dof_on_elm_neighbor[%d] = %d \n", zz,dof_on_elm_neighbor[zz]);

	      // if(dof_on_elm[zz]!= dof_on_elm_neighbor[zz])
	      //{
	      //  printf("Uh oh ? \n");
	      //}
	    }
	  
	  
	  INT v_per_elm = mesh->v_per_elm;
	  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));
	  // Find vertices for given Element
	  
	  get_incidence_row(nbr,mesh->el_v,v_on_elm_neighbor);
	  
	  printf("v_per_elm = %d \n", v_per_elm);
	  
	  for(INT zz=0; zz<dof_per_elm; ++zz)
	    {
	      printf("V_on_elm[%d] = %d \n", zz,v_on_elm[zz]);
	      printf("V_on_elm_neighbor[%d] = %d \n", zz,v_on_elm_neighbor[zz]);
	      
	      if(v_on_elm[zz]!= v_on_elm_neighbor[zz])
		{
		  printf("Uh oh ? \n");
		}
	    }
		  
	  // Now, get the basis for the neighbor
	  REAL* neighbor_basis_0_phi  = (REAL *) calloc(3,sizeof(REAL));
	  REAL* neighbor_basis_0_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));
	  
	  REAL* neighbor_basis_1_phi  = (REAL *) calloc(3,sizeof(REAL));
	  REAL* neighbor_basis_1_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));
	  
	  
	  get_FEM_basis( neighbor_basis_0_phi ,
			 neighbor_basis_0_dphi ,
			 qx_face,
			 v_on_elm_neighbor,
			 local_dof_on_elm_neighbor,
			 mesh,
			 FE->var_spaces[1]);
	    
	  local_dof_on_elm_neighbor = dof_on_elm_neighbor + FE->var_spaces[0]->dof_per_elm;
	  
	  
	  get_FEM_basis( neighbor_basis_1_phi ,
			 neighbor_basis_1_dphi ,
			 qx_face,
			 v_on_elm_neighbor,
			 local_dof_on_elm_neighbor,
			 mesh,
			 FE->var_spaces[1]);
	  
	  // p
	  local_dof_on_elm_neighbor += FE->var_spaces[1]->dof_per_elm;

	  //DEBUG
	  local_dof_on_elm_neighbor += FE->var_spaces[2]->dof_per_elm;


	  //Sanity Check
	    
	    for(INT zz=0; zz<FE->var_spaces[0]->dof_per_elm; ++zz)
	      {
		printf("FE->var_spaces[0]->phi[%d] = %f \n", zz,   FE->var_spaces[0]->phi[zz]);
		
	      }
	    printf("=====\n");
	 
	    for(INT zz=0; zz<FE->var_spaces[1]->dof_per_elm; ++zz)
	      {
		printf("FE->var_spaces[1]->phi[%d] = %f \n", zz,   FE->var_spaces[1]->phi[zz]);
		
	      }


	    printf("=====\n");
	    
	    for(INT zz=0; zz<3; ++zz)
	      {
		printf("neighbor_basis_0_phi[%d]  = %f \n", zz,   neighbor_basis_0_phi[zz]);
		
	      }

	    printf("=====\n");
	 
	    
	    for(INT zz=0; zz<3; ++zz)
	      {
		printf("neighbor_basis_1_phi[%d]  = %f \n", zz,   neighbor_basis_1_phi[zz]);
		
	      }

	      for(INT zz=0; zz<3; ++zz)
	      {
		if(FE->var_spaces[0]->phi[zz] != neighbor_basis_0_phi[zz])
		  {
		    printf("Uh-oh 00\n");
		    printf("FE->var_spaces[0]->phi[zz] = %f, neighbor_basis_0_phi[zz] = %f \n",
			   FE->var_spaces[0]->phi[zz],neighbor_basis_0_phi[zz]);
		    //	    exit(0);
		  }

		
		if(FE->var_spaces[1]->phi[zz] != neighbor_basis_1_phi[zz])
		  {
		    printf("Uh-oh  11\n");
		    //exit(0);
		  }
		//printf("neighbor_basis_1_phi[%d]  = %f \n", zz,   neighbor_basis_1_phi[zz]);
		
	      }
	    
	    printf("new DOF : phi[0] = %f (q[0] = %f), phi[1] = %f (q[1]=%f) \n", qx_face[0]- barycenter->x[0],qx_face[0],qx_face[1]- barycenter->y[0],qx_face[1]);
	    printf("elm barycenter->x[0] = %f, barycenter->y[0]=%f \n", barycenter->x[0], barycenter->y[0]);
	    
	    printf("new DOF : phi_neighb[0] = %f (q[0] = %f), q_neighb[1] = %f (q[1] = %f)\n", qx_face[0]- barycenter_neighbor->x[0], qx_face[0], qx_face[1]- barycenter_neighbor->y[0],qx_face[1]);
	    printf("elm barycenter_neighb->x[0] = %f, barycenter_neighb->y[0]=%f \n", barycenter_neighbor->x[0], barycenter_neighbor->y[0]);
	    
	    //Matrix Face Interface
	  if(bEG){
	    // u0 - q block
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	    local_col_index = 0;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

		//CG-DG
		// u0 - q0
		// u0 - q1
		kij = -0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
		kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
		
		// Divergence
		// u1-q block : 
		kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
		kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    // u1-q block
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	    local_col_index = FE->var_spaces[0]->dof_per_elm;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
		
		kij = -0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
		kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
		
		kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		// NEW SLEE : elasticity div part
		// u2-v1 block : <dy(u2),dx(v1)> 
		kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
		kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);
		
		
		
		ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }

	    // SLEE TODO LUDMIL
	    // NEIGHBOR CELL ACCESS
	    // NEIGHBOR QUAD -maybe not.
	    // NEIGHBOR BASIS FUNCTIONS
	    // NEIGHBOR DOF INDEX
	    // NEIGHBOR  ALOC_NEIGHBOR ?

	    // u0 - q (where u0 is from neighbor)
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
	    local_col_index = 0;
	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){


		//DEBUG
		
		kij = -0.5 * 2.0*neighbor_basis_0_dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
		kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
		
		// Divergence
		// u1-q block : 
		kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
		kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
		
		/*
		kij = -0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
		kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
		
		// Divergence
		// u1-q block : 
		kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
		kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
		*/
		ALoc_neighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    // u1-q (where u1 is from neighbor)
	    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;

	    // DEBUG? I think this is fine...? 
	    local_col_index = FE->var_spaces[0]->dof_per_elm;

	    // Loop over Test Functions (Rows)
	    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	      // Loop over Trial Functions (Columns)
	      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

		//DEBUG
		
		kij = -0.5 * 2.0*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
		kij -= 0.5 * neighbor_basis_1_dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
		
		kij -= 0.5 * neighbor_basis_1_dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		// NEW SLEE : elasticity div part
		// u2-v1 block : <dy(u2),dx(v1)> 
		kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
		kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->x[1]);
		/*
		kij = -0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
		kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
		
		kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
		
		// NEW SLEE : elasticity div part
		// u2-v1 block : <dy(u2),dx(v1)> 
		kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
		kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);
		*/
		
		
		ALoc_neighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	      }
	    }
	    
	    
	  }// bEG

	  
	  
	  //q-q block: 
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	      
	      kij = -0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);	
	      kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	      
	      kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);    
	      
	      //penalty term
	      kij += penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	      kij += penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);
	         
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  //sanity
	  if(finrm[0] == 0 && finrm[1] ==0)
	    {printf("!!\n");
	      exit(0);
	    }
	  // SLEE TODO LUDMIL
	  // NEIGHBOR CELL ACCESS
	  // NEIGHBOR QUAD -maybe not.
	  // NEIGHBOR BASIS FUNCTIONS
	  // NEIGHBOR DOF INDEX
	  // NEIGHBOR  ALOC_NEIGHBOR ?
	  
	  //q-q block: where q is from neighbor 
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	      
	      kij = 0.5* -2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);	
	      kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	      
	      kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);    
	       
	      
	      //penalty term
	      kij += -penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter->x[0]);
	      kij += -penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter->y[0]);
	    
	      ALoc_neighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }


	} //for face quad

	// add A loc neighbor for each 'faces' to global 
	block_LocaltoGlobal_neighbor(dof_on_elm, dof_on_elm_neighbor,FE,A,ALoc_neighbor);


      } // if INTERFACE
	    
	  
      
    }// if bWEAK
      

  } //face
  
  printf("*********** {DONE} @ END  ************************ \n");
  
  
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
  int total_num_cycle = 4; 
  // SLEE initialize the vectors to save the errors for each cycle
  double L2_error_per_cycle[total_num_cycle];
  double L2_error_p_per_cycle[total_num_cycle];
  double H1_error_per_cycle[total_num_cycle];
  double H1_stress_error_per_cycle[total_num_cycle];

  double L2_EG_error_per_cycle[total_num_cycle];
  double H1_stress_EG_error_per_cycle[total_num_cycle];
  
  // SLEE vector to save the DOF
  int dof_per_cycle[total_num_cycle];
  // SLEE vector to save the convergence rate
  double L2_conv_rate_per_cycle[total_num_cycle];
  double L2_p_conv_rate_per_cycle[total_num_cycle];
  double H1_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_conv_rate_per_cycle[total_num_cycle];

  double L2_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_EG_conv_rate_per_cycle[total_num_cycle];
  
  int global_dim_space = 0;
  
  for(int cycle=0; cycle<total_num_cycle; ++cycle)
    {
      //Aug.3.2020 SLEE initilize
      L2_error_per_cycle[cycle] = 0.;
      L2_error_p_per_cycle[cycle] = 0.;
      H1_error_per_cycle[cycle] = 0.;
      H1_stress_error_per_cycle[cycle] = 0.;

      L2_EG_error_per_cycle[cycle] = 0.;
      H1_stress_EG_error_per_cycle[cycle] = 0.;
      
      dof_per_cycle[cycle]=0;
      L2_conv_rate_per_cycle[cycle]=0.;
      L2_p_conv_rate_per_cycle[cycle]=0.;
      H1_conv_rate_per_cycle[cycle]=0.;
      H1_stress_conv_rate_per_cycle[cycle]=0.;

      L2_EG_conv_rate_per_cycle[cycle]=0.;
      H1_stress_EG_conv_rate_per_cycle[cycle]=0.;
      
      printf("************ CYCLE   %d  /   %d  ************** \n", cycle, total_num_cycle);
      fflush(stdout);
    
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
  //sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle);
  
  //DEBUG SIMPLE MESH
  sprintf(filename_per_cycle, "%s%d0.haz", inparam.gridfile,cycle);

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
  
  for(i=0;i<FE_ux.ndof;i++) {
    FE_ux.dirichlet[i] = 0;
  }
  for(i=0;i<FE_uy.ndof;i++) {
    FE_uy.dirichlet[i] = 0;
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
  //set_dirichlet_bdry_block(&FE,&mesh);


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
  //assemble_global_block(&A,&b,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity,&FE,&mesh,cq,source2D,0.0);
 
  assemble_global_block_neighbor(&A,&b,local_assembly_Elasticity,FEM_Block_RHS_Local_Elasticity,&FE,&mesh,cq,source2D,0.0);
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

  /*
  printf("------ 6\n");
  // Eliminate boundary conditions in matrix and rhs
  if(dim==2) {
    eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
  }
  if(dim==3) {
    eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);
  }
  */
  
  /**************************************************/
  //  Apply Pressure "BCs" (removes singularity)
  REAL pressureval =0.;
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
  //dCSRmat Mp;
  //assemble_global(&Mp,NULL,assemble_mass_local,&FE_p,&mesh,cq,NULL,one_coeff_scal,0.0);
  //dcsr_alloc(Mp.row, Mp.col, Mp.nnz, &A_diag[dim]);
  //dcsr_cp(&Mp, &A_diag[dim]);
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
  if(dim==2){
    if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
      solver_flag = linear_solver_bdcsr_krylov_block_3(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
  } else if (dim==3) {
    if (linear_itparam.linear_precond_type == PREC_NULL) {
      solver_flag = linear_solver_bdcsr_krylov(&A, &b, &sol, &linear_itparam);
    } else {
      solver_flag = linear_solver_bdcsr_krylov_block_4(&A, &b, &sol, &linear_itparam, NULL, A_diag);
    }
  }

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

  L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
  //HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  HDsemierror_block(solerrH1, sol.val, Dexact_sol2D, &FE, &mesh, cq, 0.0);
 
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
  printf("H1 Norm of u error    = %26.13e\n",uerrH1);
  printf("*******************************************************\n\n");

  //Jul. 10. 2020 SLEE save the errors for convergence computation
  L2_error_per_cycle[cycle] = uerrL2; 
  H1_error_per_cycle[cycle] = uerrH1; 
  L2_error_p_per_cycle[cycle] = perrL2; 

  //NEW SLEE Aug 23 2020
  //NEW ERROR FOR EG

  REAL* solerrL2_EG = (REAL *) calloc(dim+1, sizeof(REAL));
  REAL* solerrH1_EG = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
  L2error_block_EG(solerrL2_EG, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
  HDerror_block_EG(solerrH1_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
  //HDsemierror_block(solerrH1_EG, sol.val, Dexact_sol2D, &FE, &mesh, cq, 0.0);

  REAL uerrL2_EG = 0;
  REAL uerrH1_EG = 0;
  for(i=0;i<dim;i++)
    uerrL2_EG += solerrL2_EG[i]*solerrL2_EG[i];
  for(i=0;i<dim;i++)
    uerrH1_EG += solerrH1_EG[i]*solerrH1_EG[i];

  uerrL2_EG = sqrt(uerrL2_EG);
  uerrH1_EG = sqrt(uerrH1_EG);

  L2_EG_error_per_cycle[cycle] = uerrL2_EG;
  H1_stress_EG_error_per_cycle[cycle] = uerrH1_EG;
  
  printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG);
  printf("H1-stres Norm of u (EG) error    = %26.13e\n",uerrH1_EG);
  
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
  //dcsr_free( &Mp);
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
	H1_stress_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_EG_error_per_cycle[tmp]) -log(H1_stress_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle[tmp-1]) -log(dof_per_cycle[tmp]) );
	
      }
      
      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_per_cycle[tmp], dof_per_cycle[tmp],L2_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_per_cycle[tmp], dof_per_cycle[tmp],H1_conv_rate_per_cycle[tmp]);
      printf("H1 Stress Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_error_per_cycle[tmp], dof_per_cycle[tmp],H1_stress_conv_rate_per_cycle[tmp]);
      
      printf("EG - L2 EG_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_EG_error_per_cycle[tmp], dof_per_cycle[tmp],
	     L2_EG_conv_rate_per_cycle[tmp]);
      printf("EG - H1 EG_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_EG_error_per_cycle[tmp], dof_per_cycle[tmp],
	   H1_stress_EG_conv_rate_per_cycle[tmp]);
 
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
