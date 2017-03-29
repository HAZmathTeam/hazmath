/********************** Local Assembly **************************************************/
void steady_state_Darcy(REAL* ALoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time) 
{
	
  /* Computes the local stiffness matrix for the Darcy Flow system 
   * (steady-state part)
   *
   * For this problem we compute LHS of:
   *
   *  <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
   *  <div q, v> = -<W,v>
   *
   * 
   *    INPUT:
   *            uprev               Solution at previous Newton step (NULL if Linear)
   *            FE		    Finite-Element Space Struct in block form
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            v_on_elm            Vertices on element
   *            elm                 Current element
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	ALoc                Local Stiffness Matrix (Full Matrix)
   */
	
  INT i,j;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) 
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;

  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test,trial;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Data
  REAL* K = (REAL *) calloc(dim*dim,sizeof(REAL));
  REAL Kinv1=0.0;
  REAL Kinv2=0.0;
  REAL Kinv3=0.0;
  
  // Basis Functions and its derivatives if necessary
  // q
  REAL* q_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm*dim,sizeof(REAL));
  REAL* divq_phi = (REAL *) calloc(FE->var_spaces[0]->dof_per_elm,sizeof(REAL));
  // h
  REAL* h_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));

  INT local_row_index, local_col_index;
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions and old solutions at each quadrature node
    // q
    rt_basis(q_phi,divq_phi,qx,v_on_elm,dof_on_elm,mesh);
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    
    // h
    PX_H1_basis(h_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[1]->FEtype,mesh);
    
    // Data
    porosity(K,qx,0.0);
    Kinv1 = 1.0/K[0];
    Kinv2 = 1.0/K[4];
    Kinv3 = 1.0/K[8];

    // q-q block: <K^(-1) q, r> 
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm; trial++) {
	kij = Kinv1*q_phi[trial*dim]*q_phi[test*dim] + Kinv2*q_phi[trial*dim+1]*q_phi[test*dim+1];
	if(dim==3) kij+=Kinv3*q_phi[trial*dim+2]*q_phi[test*dim+2];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
    
    // q-h block:  <h, div r>
    local_col_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm; trial++) {
	kij = h_phi[trial]*divq_phi[test];
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
	kij = divq_phi[trial]*h_phi[test];
	ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    }
  }
    
  if (h_phi) free(h_phi);
  if (q_phi) free(q_phi);
  if(divq_phi) free(divq_phi);
  if(K) free(K);

  if(qx) free(qx);
  return;
}

void steady_state_Darcy_RHS(REAL* bLoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
	
  /* Computes the local volume rhs for the Darcy Flow system 
   * (steady-state part)
   *
   * For this problem we compute the non-boundary RHS of:
   *
   *  <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
   *  <div q, v> = -<W,v> 
   * 
   * 
   *    INPUT:
   *            uprev               Solution at previous Newton step (NULL if linear)
   *            FE		    Finite-Element Space Struct in block form
   *	      	mesh                Mesh Struct
   *            cq                  Quadrature Coordinates and Weights
   *            dof_on_elm          Specific DOF on element
   *            v_on_elm            Vertices on element
   *            elm                 Current element
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	bLoc                Local RHS
   */
	
  INT i,j;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) 
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm = NULL;

  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test;
  // Quadrature Weights and Nodes
  REAL w; 
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij = 0.0;  
	
  // Basis Functions and its derivatives if necessary
  // h
  REAL* h_phi = (REAL *) calloc(FE->var_spaces[1]->dof_per_elm,sizeof(REAL));

  REAL rhs_val;

  INT local_row_index;
      
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {        
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    //  Get the Basis Functions and old solutions at each quadrature node
    // q
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    // h
    PX_H1_basis(h_phi,NULL,qx,local_dof_on_elm,FE->var_spaces[1]->FEtype,mesh);

    (*rhs)(&rhs_val,qx,0.0);

    // q block: 0
    local_row_index = 0;

    // h block: -<W,v>
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      kij = -rhs_val*h_phi[test];
      bLoc[(local_row_index+test)] += w*kij;
    } 
  }
    
  if (h_phi) free(h_phi);

  if(qx) free(qx);
  return;
}

void steady_state_Darcy_bdryRHS(REAL* bLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_f,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
	
  /* Computes the local boundary rhs for the Darcy Flow system 
   * (steady-state part)
   *
   * For this problem we compute the boundary RHS of:
   *
   *  <K^(-1) q, r> + <h, div r> = <g,r*n>_boundary
   *  <div q, v> = -<W,v> 
   * 
   * 
   *    INPUT:
   *            old_sol               Solution at previous Newton step (NULL if linear)
   *            FE		    Finite-Element Space Struct in block form
   *	      	mesh                Mesh Struct
   *            dof_on_f            Specific DOF on face
   *            dof_per_f           # DOF per face
   *            dof_on_elm          Specific DOF on element
   *            v_on_elm            Vertices on element
   *            face                Current face
   *            elm                 Current element
   *            rhs                 RHS function
   *            time                Time if time dependent
   *
   *    OUTPUT:
   *	       	bLoc                Local RHS
   */
	
  // Mesh and FE data
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  
  // Loop Indices
  INT i,j,quad,test,doft,rowa,rowb,jcntr,ed;
  
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
	
  // Basis Functions and its derivatives if necessary
  // q
  REAL* q_phi = (REAL *) calloc(dim*dof_per_elm,sizeof(REAL));
  REAL* divq_phi = (REAL *) calloc(dof_per_elm,sizeof(REAL));
  
  // Value of g
  REAL rhs_val;
      
  //  Sum over midpoints of edges
  for (quad=0;quad<nq;quad++) {
    ed = ed_on_f[quad]-1;
    qx[0] = mesh->ed_mid[ed*dim];
    qx[1] = mesh->ed_mid[ed*dim+1];
    if(mesh->dim==3) qx[2] = mesh->ed_mid[ed*dim+2];
    
    //  Get the Basis Functions at each quadrature node
    rt_basis(q_phi,divq_phi,qx,v_on_elm,dof_on_elm,mesh);

    (*rhs)(&rhs_val,qx,0.0);

    // Loop over Test Functions (Rows)
    for (test=0; test<dof_per_f;test++) {
      // Make sure ordering for global matrix is right
      for(j=0;j<FE->dof_per_elm;j++) {
	if(dof_on_f[test]==dof_on_elm[j]) {
	  doft = j;
	}
      }
      kij = rhs_val*(nx*q_phi[doft*dim] + ny*q_phi[doft*dim+1]);
      if(dim==3) kij +=rhs_val*nz*q_phi[doft*dim+2];
      bLoc[test] += w*kij;
    } 
  }
    
  if (q_phi) free(q_phi);
  if(divq_phi) free(divq_phi);

  if(qx) free(qx);
  if(ed_on_f) free(ed_on_f);
  return;
}
/******************************************************************************************************/
