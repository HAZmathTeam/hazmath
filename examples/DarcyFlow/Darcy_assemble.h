/*! \file DarcySystem.h
 *
 *  Created by Adler, Hu, Zikatanov on 8/30/16.
 *  Copyright 2016_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routines (LHS, RHS, and boundary terms)
 *        for the DarcyFlow example.
 *
 */

/*!
 * \fn void steady_state_Darcy(REAL* ALoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time)
 *
 * \brief Computes the local stiffness matrix for the Darcy Flow system
 *      (steady-state part)
 *
 *      For this problem we compute LHS of:
 *
 *       <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
 *       <div q, v> = -<W,v>
 *
 *
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific vertices on element
 * \param elm           Current element
 * \param time          Physical Time if time dependent
 *
 * \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (q,h)
 *
 * \note Assumes 3D only
 *
 */
void steady_state_Darcy(REAL* ALoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time) 
{

  // Loop indices
  INT i,j,quad,test,trial;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;
  INT dim = mesh->dim;
  
  // Quadrature Weights and Nodes
  REAL w;
  // Stiffness Matrix Entry
  REAL kij = 0.0,alocpq=-10.;

  REAL* qx = (REAL *) calloc(dim,sizeof(REAL)); 
  // Porosity Coefficient Data
  REAL *K=(REAL *)calloc(dim*dim,sizeof(REAL)); //ltz1
  REAL *Kinv=(REAL *)calloc(dim*dim,sizeof(REAL)); //ltz1
  void *wrk=(void *)calloc(dim*(dim+1)*sizeof(REAL)+dim*sizeof(INT),sizeof(char));//ltz1
  // Keep track of local indexing
  INT local_row_index, local_col_index,di,pdim,qdim;
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions at each quadrature node
    // q
    get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,dof_on_elm,mesh,FE->var_spaces[0]);

    // h (need to shift local index)
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
    
    // Data
    porosity(K,qx,time,NULL); //ltz1
    invfull(Kinv,dim, K, wrk); //ltz1
    //       prtmat(dim,dim,Kinv,"Kinv");
    //    REAL Kinv1=Kinv[0];
    //    REAL Kinv2=Kinv[4];
    //    REAL Kinv3=Kinv[8];
    //    prtmat(dim,dim,Kinv,"Kinv");
    //    Kinv1 = 1.0/K[0];
    //    Kinv2 = 1.0/K[4];
    //    Kinv3 = 1.0/K[8];

    // q-q block: <K^(-1) q, r>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      pdim=test*dim;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	qdim=trial*dim;
	alocpq=0.;
	for(i=0;i<dim;i++){
	  di=dim*i;
	  for(j=0;j<dim;j++){	  
	    alocpq += Kinv[i*dim+j]*FE->var_spaces[0]->phi[qdim+j]*FE->var_spaces[0]->phi[pdim+i];
	  }
	}
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*alocpq;
      }
    }
/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx     */

/*     // q-q block: <K^(-1) q, r> */
/*     local_row_index = 0; */
/*     local_col_index = 0; */
/*     // Loop over Test Functions (Rows) */
/*     for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) { */
/*       // Loop over Trial Functions (Columns) */
/*       for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) { */
/*         kij = Kinv1*FE->var_spaces[0]->phi[trial*dim]*FE->var_spaces[0]->phi[test*dim] + */
/*             Kinv2*FE->var_spaces[0]->phi[trial*dim+1]*FE->var_spaces[0]->phi[test*dim+1]; */
/*         if(dim==3) kij+=Kinv3*FE->var_spaces[0]->phi[trial*dim+2]*FE->var_spaces[0]->phi[test*dim+2]; */
/*         ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij; */
/*       } */
/*     } */
    
    // q-h block:  <h, div r>
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
        kij = FE->var_spaces[1]->phi[trial]*FE->var_spaces[0]->dphi[test];
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }

    // h-q block: <div q, v>
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
        kij = FE->var_spaces[0]->dphi[trial]*FE->var_spaces[1]->phi[test];
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
  }

  // Free stuff
  if(K) free(K);
  if(qx) free(qx);
  return;
}

/*!
 * \fn void steady_state_Darcy_RHS(REAL* bLoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local volume rhs for the Darcy Flow system
 *        (steady-state part)
 *
 *        For this problem we compute the non-boundary RHS of:
 *
 *       <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
 *       <div q, v> = -<W,v>
 *
 *
 * \param FE            Block FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific vertices on element
 * \param elm           Current element
 * \param rhs           Function for RHS (steady-state part)
 * \param time          Physical Time if time dependent
 *
 * \return bLoc         Local RHS
 *
 * \note Assumes 3D only
 *
 */
void steady_state_Darcy_RHS(REAL* bLoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop Indices
  INT quad,test,i,j;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;
  INT dim = mesh->dim;
  
  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // RHS function value at quadrature
  REAL rhs_val;

  // Local index counter
  INT local_row_index;

  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions at each quadrature node
    // q
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    // h
    get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);

    (*rhs)(&rhs_val,qx,time,NULL);

    // q block: 0
    local_row_index = 0;

    // h block: -<W,v>
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      kij = -rhs_val*FE->var_spaces[1]->phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    }
  }

  // Free Stuff
  if(qx) free(qx);

  return;
}

/*!
 * \fn steady_state_Darcy_bdryRHS(REAL* bLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_f,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local boundary rhs for the Darcy Flow system
 *        (steady-state part)
 *
 *        For this problem we compute the boundary RHS of:
 *
 *       <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
 *       <div q, v> = -<W,v>
 *
 *
 * \param old_sol       Solution at previous Newton step (NULL if linear)
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_f      Specific DOF on face
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific vertices on element
 * \param dof_per_f     # DOF per face
 * \param face          Current face
 * \param elm           Current element
 * \param rhs           Function for Boundary RHS (steady-state part)
 * \param time          Physical Time if time dependent
 *
 * \return bLoc         Local RHS
 *
 * \note Assumes 3D only
 *
 */
void steady_state_Darcy_bdryRHS(REAL* bLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_f,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop Indices
  INT i,j,quad,test,doft,rowa,rowb,jcntr,ed;
  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Quadrature Weights and Nodes
  INT nq = 2*dim-3; // = ed_per_face
  REAL* qx = (REAL *) calloc(nq,sizeof(REAL));
  // 3D: Using triangle midpoint rule, so qx is midpoint of edges and w is |F|/3
  REAL w = mesh->f_area[face]/3.0;
  // Find the edges for the given face
  INT* ed_on_f = (INT *) calloc(nq,sizeof(INT));
  rowa = mesh->f_ed->IA[face]-1;
  rowb = mesh->f_ed->IA[face+1]-1;
  jcntr = 0;
  for(i=rowa;i<rowb;i++) {
    ed_on_f[jcntr] = mesh->f_ed->JA[i];
    jcntr++;
  }

  // Get normal vector components on face
  REAL nx = mesh->f_norm[face*dim];
  REAL ny = mesh->f_norm[face*dim+1];
  REAL nz = 0.0;
  if(dim==3) nz = mesh->f_norm[face*dim+2];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Value of g
  REAL rhs_val;

  //  Sum over midpoints of edges
  for (quad=0;quad<nq;quad++) {
    ed = ed_on_f[quad]-1;
    qx[0] = mesh->ed_mid[ed*dim];
    qx[1] = mesh->ed_mid[ed*dim+1];
    if(mesh->dim==3) qx[2] = mesh->ed_mid[ed*dim+2];
    
    //  Get the Basis Functions at each quadrature node
    get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

    // Get RHS at quadrature points
    (*rhs)(&rhs_val,qx,time,NULL);

    // Loop over Test Functions (Rows)
    for (test=0; test<dof_per_f;test++) {
      // Make sure ordering for global matrix is right
      for(j=0;j<FE->dof_per_elm;j++) {
        if(dof_on_f[test]==dof_on_elm[j]) {
          doft = j;
        }
      }
      kij = rhs_val*(nx*FE->phi[doft*dim] + ny*FE->phi[doft*dim+1]);
      if(dim==3) kij +=rhs_val*nz*FE->phi[doft*dim+2];
      bLoc[test] += w*kij;
    }
  }

  // Free Stuff
  if(qx) free(qx);
  if(ed_on_f) free(ed_on_f);
  return;
}
/******************************************************************************************************/
