/*! \file examples/Elasticity/elasticity.c
 *
 *  This code is EG for elasticity
 *  created by SLee (20200812)
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

#include "elasticity_params.h"
#include "elasticity_error.h"
#include "elasticity_system.h"
/*********************************************************************************/



void local_assembly_Elasticity_FACE(block_dCSRmat* A, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, \
				    INT face, iCSRmat *f_el)
{
  bool bEG = true;//false;//false;//true;//false;//true;//false;//true;

  // INDEX
  INT i,j,k;
  INT dim = mesh->dim;
  //INT face = 0;
  INT rowa,rowb;
  
  // Mesh Stuff 
  INT nspaces = FE->nspaces;


  INT v_per_elm = mesh->v_per_elm;
  
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));

  //
  INT test, trial;
  INT local_row_index, local_col_index;
  REAL kij = 0.0;
  
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));

  REAL* ALoc_u_vneighbor = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_u_v = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_uneighbor_v = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* ALoc_uneighbor_vneighbor = (REAL *) calloc(local_size,sizeof(REAL));


  for (j=0; j<local_size; j++) {
    
    ALoc[j]=0;
    ALoc_u_vneighbor[j]=0;
    ALoc_u_v[j]=0;
    ALoc_uneighbor_v[j]=0;
    ALoc_uneighbor_vneighbor[j]=0;
  }

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));


  // THIS NEEDS TO BE CODED BETTER // slee todo
  double mu = 1.;
  double lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.; //000000.;

  
  INT quad,quad_face;
  INT jk, jkl,jcntr,ed;

  //iCSRmat *f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  //icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);
  
  /* Loop over all Faces */
  //for (face=0; face<mesh->nface; face++) {
    
  INT rowa_neighbor,rowb_neighbor,jcntr_neighbor;
  INT k_neighbor, j_neighbor;
  double fiarea;
  
  coordinates *barycenter = allocatecoords(1,dim);
  coordinates *barycenter_neighbor = allocatecoords(1,dim);
  
  
  // printf("===== FACE = %d  ===== \n",face);
  INT neighbor_index[2];
  neighbor_index[0] = -1;
  neighbor_index[1] = -1;
  INT counter = 0;


  //printf("=============================\n");
  //printf("** FACE = %d \n", face );
  //printf("=============================\n");
  int pq,nbr0;
  for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){	
    
    //printf("-- pq = %d, f_el->IA[face] = %d, f_el->IA[face+1] = %d \n", pq, f_el->IA[face], f_el->IA[face+1]);
    nbr0=f_el->JA[pq];
    //printf("-- nbr0 = %d  \n", nbr0);
    
    neighbor_index[counter] = nbr0;
    counter++;	
    
  }
  /// print out
  for(pq=0;pq<2;++pq)
    {
      if(counter == 2)
	{
	  //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
	}
      else if(counter == 1){
	if(pq == 0)
	  {
	    //printf("neighbor_index[%d]= %d || counter  = %d\n", pq, neighbor_index[pq],counter);
	  }
	else{
	  //printf("-\n");
	}
      }
    }
  
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

  INT* local_dof_on_elm_face = NULL;	
  INT *local_dof_on_elm_face_interface = NULL;
  INT *local_dof_on_elm_neighbor = NULL;
  
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  //REAL *data_face_end=finrm + dim; //end
    
  // Get normal vector values. 
  finrm[0]=mesh->f_norm[face*dim+0];
  finrm[1]=mesh->f_norm[face*dim+1];
  //printf("FACE = %d, ELM = %d, (Nx,Ny) = (%f,%f) \n", face, neighbor_index[0], finrm[0], finrm[1]);
    
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
    
  INT maxdim=2;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;
    
  // Get xfi points..
  for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
    //printf("M** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);
      
      
    j=jkl-mesh->f_v->IA[face];
    k=mesh->f_v->JA[jkl];
      
    //printf("M** j = %d, k = %d \n", j, k );
      
      
    xfi[j*dim+0]=mesh->cv->x[k];
    xfi[j*dim+1]=mesh->cv->y[k];
      
    //printf("M** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);
      
  }

  //define the length of the face-
  fiarea=mesh->f_area[face];
  //printf("FACE = %d, ELM = %d, fiarea = %f \n", face, neighbor_index[0], fiarea);
    
  // get the quadrautre 
  zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);
    
  // get the penalty (FACE ASSEMBLE)
  double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,(REAL )(dim-1)));
  double penalty_term_new =  0.;//1000000./fiarea;
    
  // SLEE
  // Get the BD values 
  //REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));
    
  for (quad_face=0;quad_face<nquad_face;quad_face++) {
      
    qx_face[0] = cq_face->x[quad_face];
    qx_face[1] = cq_face->y[quad_face];	  
      
    if(dim==3) qx_face[2] = cq_face->z[quad_face];
      
    w_face = cq_face->w[quad_face];
      
    //(*exact_sol2D)(val_true_face,qx_face,time,
    //	     &(mesh->el_flag[neighbor_index[0]])); //DUMMY VALUE at the end // ???      
      
      
    ///////////////////////////////////////////////////
    if(counter == 1){ //only at boundary 
      //printf("== BD = %d \n", face);
      //Sanity Check
      if(mesh->f_flag[face] == 0)
	{printf("check-error0\n"); exit(0);}

      bool bWeakBC_FACE = true;//false;//true;//FACE
	
      if(bWeakBC_FACE){
	  
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
	  
	//for the EG part
	local_dof_on_elm_face += FE->var_spaces[2]->dof_per_elm;
	  
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

	    kij += penalty_term_new * FE->var_spaces[0]->phi[trial] * finrm[0]* FE->var_spaces[0]->phi[test]* finrm[0];
	     
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

	      
	    kij += penalty_term_new * FE->var_spaces[0]->phi[trial] * finrm[0]* FE->var_spaces[1]->phi[test]* finrm[1];

	      
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

	    kij += penalty_term_new * FE->var_spaces[1]->phi[trial] * finrm[1]* FE->var_spaces[0]->phi[test]* finrm[0];

	      
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

	    kij += penalty_term_new *  FE->var_spaces[1]->phi[trial] * finrm[1] *  FE->var_spaces[1]->phi[test] * finrm[1];
	     
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

	      
	    kij += penalty_term_new * (qx_face[0]- barycenter->x[0]) * finrm[0] * (qx_face[0]- barycenter->x[0]) * finrm[0];
	    kij += penalty_term_new * (qx_face[1]- barycenter->y[0]) * finrm[1] * (qx_face[1]- barycenter->y[0]) * finrm[1];

	    kij += penalty_term_new * (qx_face[0]- barycenter->x[0]) * finrm[0] * (qx_face[1]- barycenter->y[0]) * finrm[1];
	    kij += penalty_term_new * (qx_face[1]- barycenter->y[0]) * finrm[1] * (qx_face[0]- barycenter->x[0]) * finrm[0];
	    
	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Face Boundary
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

	      kij += penalty_term_new * FE->var_spaces[0]->phi[trial] * finrm[0] * ( (qx_face[0]- barycenter->x[0]) *finrm[0] + (qx_face[1]- barycenter->y[0]) *finrm[1]);
	
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

	      kij += penalty_term_new * FE->var_spaces[1]->phi[trial] * finrm[1] * ((qx_face[0]- barycenter->x[0]) * finrm[0] + (qx_face[1]- barycenter->y[0]) * finrm[1]);

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

	      kij += penalty_term_new * FE->var_spaces[0]->phi[test] * finrm[0] * ((qx_face[0]- barycenter->x[0]) * finrm[0] + (qx_face[1]- barycenter->y[0]) * finrm[1]);
	
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

	      kij += penalty_term_new * FE->var_spaces[1]->phi[test] * finrm[1] *  ((qx_face[0]- barycenter->x[0])* finrm[0] + (qx_face[1]- barycenter->y[0])* finrm[1]);
	
	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	    
	} // bEG
      }
    }
    else if(counter ==2)
      { //on interface

	local_dof_on_elm_face_interface = NULL;
	local_dof_on_elm_neighbor = NULL;
	// HERE IN THE INTERFACE, WE HAVE FOUR DIFFERENT LOCAL MATRICES
	// CASE I   u v
	// CASE II  u v_neighbor
	// CASE III u_neighbor v
	// CASE IV  u_neighbor v_neighbor
	//Sanity Cehck
	if(mesh->f_flag[face] != 0)
	  {printf("well? \n"); exit(0);}

	
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
	
	//FOR EG
	local_dof_on_elm_face_interface += FE->var_spaces[2]->dof_per_elm;
	
	// Now, get the basis for the neighbor
	REAL* neighbor_basis_0_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_0_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));
	
	REAL* neighbor_basis_1_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_1_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));
	
	get_FEM_basis( neighbor_basis_0_phi ,
		       neighbor_basis_0_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[0]);
	
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




	if(bEG){
	  // u0 - q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = -0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1]    * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]    * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
	      
	      // Divergence
	      // u1-q block : 
	      kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
	      kij -= 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
	      
	      ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-q block
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      
	      kij = -0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0]       * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1]       * (qx_face[0]- barycenter->x[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
	      kij -= 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);
	      
	      ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  
	  // u0 - q (where u0 is from neighbor)
	  // column (trial is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      //DEBUG
	      kij = -0.5 * 2.0*neighbor_basis_0_dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[1]    * (qx_face[0]- barycenter->x[0]);
	      kij -= 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[0]    * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
	      
	      // Divergence
	      // u1-q block : 
	      kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
	      kij -= 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
	      
	      ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-q (where u1 is from neighbor)
	  // column (trial is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      //DEBUG
	      kij = -0.5 * 2.0*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * neighbor_basis_1_dphi[trial*dim]       * finrm[0] * (qx_face[1]- barycenter->y[0]);
	      kij -= 0.5 * neighbor_basis_1_dphi[trial*dim]       * finrm[1] * (qx_face[0]- barycenter->x[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);	
	      kij -= 0.5 *lambda*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);

	      ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u0 - q_neighbor block
	  // (q has minus sign)
	  // row (test is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = 0.5 * 2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter_neighbor->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij += 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);
	      kij += 0.5 * FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter_neighbor->y[0]); //FE->var_spaces[dim]->phi[test];
	      
	      // Divergence
	      // u1-q block : 
	      kij += 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
	      kij += 0.5 * lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      
	      ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-q_neighbor block
	  // (q has minus sign)
	  // row (test is from neighbor)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      
	      kij = 0.5 * 2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij += 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]);	
	      kij += 0.5 *lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]);
	      
	      
	      ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	  // u0-neighbor - q_neighbor block
	  // (q has minus sign)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = 0;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
	      
	      //CG-DG
	      // u0 - q0
	      // u0 - q1
	      kij = 0.5 * 2.0*neighbor_basis_0_dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter_neighbor->x[0]); //FE->var_spaces[dim]->phi[test];
	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);
	      
	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter_neighbor->y[0]); //FE->var_spaces[dim]->phi[test];
	      
	      // Divergence
	      // u1-q block : 
	      kij += 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);//FE->var_spaces[dim]->dphi[test*dim];	
	      kij += 0.5 * lambda*neighbor_basis_0_dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      
	      ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	  
	  // u1-neighbor - q_neighbor block
	  // (q has minus sign)
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  local_col_index = FE->var_spaces[0]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
	      
	      kij = 0.5 * 2.0*neighbor_basis_1_dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * neighbor_basis_1_dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter_neighbor->y[0]);
	      kij += 0.5 * neighbor_basis_0_dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter_neighbor->x[0]);
	      
	      // NEW SLEE : elasticity div part
	      // u2-v1 block : <dy(u2),dx(v1)> 
	      kij += 0.5 *lambda*neighbor_basis_0_dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]);	
	      kij += 0.5 *lambda*neighbor_basis_0_dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]);
	      
	      
	      ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }
	}// bEG TERMS.

	//q-q block: 
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    
	    kij = -0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5* 2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);
	    
	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)> 
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);	
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	    
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);    
	    
	    //penalty term
	    kij += penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	    kij += penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);
	    
	    
	    kij += penalty_term_new * (qx_face[0]- barycenter->x[0]) *  finrm[0] * (qx_face[0]- barycenter->x[0]) *  finrm[0];
	    kij += penalty_term_new * (qx_face[1]- barycenter->y[0]) *  finrm[1] * (qx_face[1]- barycenter->y[0]) *  finrm[1];

	    kij += penalty_term_new * (qx_face[0]- barycenter->x[0]) *  finrm[0] * (qx_face[1]- barycenter->y[0]) *  finrm[1];
	    kij += penalty_term_new * (qx_face[1]- barycenter->y[0]) *  finrm[1] * (qx_face[0]- barycenter->x[0]) *  finrm[0];

	    
	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	
	//q-neighbor q block: where q is from neighbor 
	// column (trial is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }
	    
	    kij = 0.5* -2.0 * 1. * finrm[0] * (qx_face[0]- barycenter->x[0]);
	    kij -= 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter->y[0]);
	    
	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)> 
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);	
	    kij -= 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
	    
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
	    kij -= 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);    
	    
	    //penalty_term
	    kij += -penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter->x[0]) ;
	    kij += -penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter->y[0]) ;
	    
	    
	    //penalty term
	    
	    kij += -penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[0]- barycenter->x[0]) *  finrm[0];
	    kij += -penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[1]- barycenter->y[0]) *  finrm[1];

	    kij += -penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[1]- barycenter->y[0]) *  finrm[1];
	    kij += -penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[0]- barycenter->x[0]) *  finrm[0];
	 
	    
	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

		
	//q q-neighbor block: where q (Test) is from neighbor => we need minus sign 
	// row (test is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }
	    
	    kij = 0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    
	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)> 
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);	
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);    
	    
	    //penalty_term
	    kij += -penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter->x[0]) ;
	    kij += -penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter->y[0]) ;
	    
	    
	    //penalty term    
	    kij += -penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[0]- barycenter->x[0]) *  finrm[0];
	    kij += -penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[1]- barycenter->y[0]) *  finrm[1];

	    kij += -penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[1]- barycenter->y[0]) *  finrm[1];
	    kij += -penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[0]- barycenter->x[0]) *  finrm[0];
	  
	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	//q-neighobr q-neighbor block: where q (Test) is from neighbor => we need minus sign 
	// column (trial  is from neighbor)
	// row (test is from neighbor)
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
	    if(test == 1)
	      {
		printf("oops !? \n");
		exit(0);
	      }
	    
	    kij = 0.5* 2.0 * 1. * finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    kij += 0.5*2.0 * 1. * finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    
	    // NEW SLEE : elasticity div part
	    // u2-v1 block : <dy(u2),dx(v1)> 
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);	
	    kij += 0.5*lambda * 1. *  finrm[0] * (qx_face[0]- barycenter_neighbor->x[0]);
	    
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);
	    kij += 0.5*lambda * 1. *  finrm[1] * (qx_face[1]- barycenter_neighbor->y[0]);    
	    
	    //penalty_term
	    kij += penalty_term * (qx_face[0]- barycenter_neighbor->x[0]) * (qx_face[0]- barycenter_neighbor->x[0]) ;
	    kij += penalty_term * (qx_face[1]- barycenter_neighbor->y[0]) * (qx_face[1]- barycenter_neighbor->y[0]) ;
	    
	    
	    //penalty term
	    kij += penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[0]- barycenter_neighbor->x[0]) *  finrm[0];
	    kij += penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[1]- barycenter_neighbor->y[0]) *  finrm[1];

	    kij += penalty_term_new * (qx_face[0]- barycenter_neighbor->x[0])  *  finrm[0]* (qx_face[1]- barycenter_neighbor->y[0]) *  finrm[1];
	    kij += penalty_term_new * (qx_face[1]- barycenter_neighbor->y[0])  *  finrm[1]* (qx_face[0]- barycenter_neighbor->x[0]) *  finrm[0];
	  
	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	
      }
  }// face q loop
    
  if(counter == 1)
    block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);
  else if(counter == 2)
    {
      //printf("=============== 111111111111 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc_u_v);
      // row is neighbor
	
      //printf("=============== 222222222222 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm,FE,A,ALoc_u_vneighbor);
      //block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm_neighbor,FE,A,ALoc_u_vneighbor);

      //printf("=============== 333333333333 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm_neighbor,FE,A,ALoc_uneighbor_v);
      //block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm,FE,A,ALoc_uneighbor_v); 

      //printf("=============== 444444444444 =================== \n");
      block_LocaltoGlobal_neighbor(dof_on_elm_neighbor,dof_on_elm_neighbor,FE,A,ALoc_uneighbor_vneighbor);
    }
  else{
    printf("!@#$\n");exit(0);
  }

     
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
void FEM_Block_RHS_Local_Elasticity(dvector *b,REAL* bLoc,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
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
 
  //Sanity Check // SLEE
  //for(i=0;i<FE->nspaces;i++) {
  //printf("dof_per_face = %d,   dof_per_face_blk[%d] = %d \n", dof_per_face, i, dof_per_face_blk[i]);
  //}
  //exit(0);
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal vector on face i
  //REAL *data_face_end=finrm + dim; //end
  INT nq1d_face=3;
  qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);
  //REAL fiarea=0e0;
  INT jk,k,face,quad_face,rowa, rowb, jcntr,ed, jkl, j;
  REAL qx_face[maxdim];
  REAL w_face;
  INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;
  
  // SLEE
  // Get the BD values 
  REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));

  bool bWeakBC_RHS = true;//false;//true; //false;
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
	double fiarea=mesh->f_area[face];
	double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,(REAL )(dim-1)));
	double lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.; //000000.0;
	double penalty_term_new = 0.;//1000000./fiarea;
	//nq1d_face == 3, 3 quad points. 
	zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

	REAL edge_length = mesh->ed_len[face];
	//if(edge_length != fiarea)
	//{
	//  printf("length = %f, area = %f \n", edge_length, fiarea );
	//   exit(0);
	//}

	/*
	  double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,(REAL )(dim-1)));
	  
	  printf("Penlaty RHS= %f \n", penalty_term);
	*/
	  
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


		  bLoc[(local_row_index+test)] += penalty_term_new * w_face*
		    (val_true_face[0]* finrm[0] + val_true_face[1]* finrm[1]) * FE->var_spaces[i]->phi[test]* finrm[i];
		}
	      
	      else if(i == 2)
		{
		  //SLEE
		  // Note that new DOF is \phi^3 = [x ; y] 
		  //printf("i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
		  bLoc[(local_row_index+test)] +=
		    penalty_term *  w_face * (val_true_face[0] * (qx_face[0] - barycenter->x[0])
					      + val_true_face[1]*  (qx_face[1] - barycenter->y[0]));
		  bLoc[(local_row_index+test)] +=
		    penalty_term_new *  w_face * (val_true_face[0] * finrm[0] * (qx_face[0] - barycenter->x[0]) * finrm[0]
						  + val_true_face[1] * finrm[1] * (qx_face[1] - barycenter->y[0]) * finrm[1]
						  + val_true_face[0] * finrm[0] * (qx_face[1] - barycenter->y[0]) * finrm[1]
						  + val_true_face[1] * finrm[1] * (qx_face[0] - barycenter->x[0]) * finrm[0]);
	
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

  block_LocaltoGlobal_RHS(dof_on_elm,FE,b,bLoc);

  
  return;
}





void local_assembly_Elasticity(block_dCSRmat* A,dvector *b,REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time, INT* switch_on_face)
//void local_assembly_Elasticity(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)

{
  /*
    for(INT kk=0;kk<7;++kk)
    printf("dof_on_elm[%d]=%d \n", kk, dof_on_elm[kk]);
    exit(0);
  */
  bool bEG =true;//false;//true;// false;//true;//false;//true;//false;//true;//false;//true;//false;//true;// true;//false;//true;//false;//true;//false;//true;//true;//true;//false;//true;//false;//false;//true;//false;//true;//false;//true;//false;
  ///////////
  
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

  //printf("============================================================\n");
  //printf("ELEMENT = %D, x = %f , y = %f \n", elm,  barycenter->x[0],   barycenter->y[0]);
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

  double lambda = LAME_LAMBDA_GLOBAL ;//000000.;//000000.;//000000.; //000000.0;

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


  block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);

  //printf("*********** {DONE} @ END  ************************ \n");
  
  
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
  int total_num_cycle = TOTAL_NUM_CYCLES_GLOBAL ; 
  // SLEE initialize the vectors to save the errors for each cycle
  double L2_error_per_cycle[total_num_cycle];
  double L2_error_p_per_cycle[total_num_cycle];
  double H1_error_per_cycle[total_num_cycle];
  double H1_stress_error_per_cycle[total_num_cycle];

  double L2_EG_error_per_cycle[total_num_cycle];
  double H1_EG_error_per_cycle[total_num_cycle];
  double H1_stress_EG_error_per_cycle[total_num_cycle];
  double H1_energy_EG_error_per_cycle[total_num_cycle];
 
  
  // SLEE vector to save the DOF
  int dof_per_cycle_CG[total_num_cycle];
  int dof_per_cycle_EG[total_num_cycle];

  double mesh_size_per_cycle[total_num_cycle];
 
  // SLEE vector to save the convergence rate
  double L2_conv_rate_per_cycle[total_num_cycle];
  double L2_p_conv_rate_per_cycle[total_num_cycle];
  double H1_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_conv_rate_per_cycle[total_num_cycle];

  double L2_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_stress_EG_conv_rate_per_cycle[total_num_cycle];
  double H1_energy_EG_conv_rate_per_cycle[total_num_cycle];
 
  
  int global_dim_space = 0;
  
  for(int cycle=0; cycle<total_num_cycle; ++cycle)
    {
      //Aug.3.2020 SLEE initilize
      L2_error_per_cycle[cycle] = 0.;
      //L2_error_p_per_cycle[cycle] = 0.;
      H1_error_per_cycle[cycle] = 0.;
      H1_stress_error_per_cycle[cycle] = 0.;

      L2_EG_error_per_cycle[cycle] = 0.;
      H1_EG_error_per_cycle[cycle] = 0.;
      H1_stress_EG_error_per_cycle[cycle] = 0.;
      H1_energy_EG_error_per_cycle[cycle] = 0.;
     
      
      dof_per_cycle_CG[cycle]=0;
      dof_per_cycle_EG[cycle]=0;
      mesh_size_per_cycle[cycle]=0.;
      
      L2_conv_rate_per_cycle[cycle]=0.;
      L2_p_conv_rate_per_cycle[cycle]=0.;
      H1_conv_rate_per_cycle[cycle]=0.;
      H1_stress_conv_rate_per_cycle[cycle]=0.;

      L2_EG_conv_rate_per_cycle[cycle]=0.;      
      H1_EG_conv_rate_per_cycle[cycle]=0.;
      H1_stress_EG_conv_rate_per_cycle[cycle]=0.;
      H1_energy_EG_conv_rate_per_cycle[cycle]=0.;
     
      
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
      char filename_per_cycle[512]={'\0'};
      //sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle);
    
      //DEBUG SIMPLE MESH
      sprintf(filename_per_cycle, "%s%d.haz", inparam.gridfile,cycle+1); 
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
      /*
	for(i=0;i<FE_ux.ndof;i++) {
	FE_ux.dirichlet[i] = 0;
	}
	for(i=0;i<FE_uy.ndof;i++) {
	FE_uy.dirichlet[i] = 0;
	}
      */    
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
      dof_per_cycle_CG[cycle]  = FE_ux.ndof + FE_uy.ndof;////FE.ndof + FE.nbdof;
      dof_per_cycle_EG[cycle] =  FE_ux.ndof + FE_uy.ndof + FE_p.ndof;////FE.ndof + FE.nbdof;
  
  
      printf("FE.ndof = %d | ndof = %d | FE_ux.ndof = %d | FE_uy.ndof = %d |  FE_p.ndof = %d \n",
	     FE.ndof , ndof, FE_ux.ndof  , FE_uy.ndof , FE_p.ndof );
      printf("FE.nbdof = %d \n", FE.nbdof);

      printf("###########\n");
      printf("CG DOF = %d\n",   dof_per_cycle_CG[cycle] );
      printf("EG DOF = %d\n",   dof_per_cycle_EG[cycle] );
      printf("###########\n");
      /*** Assemble the matrix and right hand side *******************************/
      /* Here we assemble the discrete system:
       *  The weak form is:
       *
       *  <2*eps(u), eps(v)> - <p, div v> = <f, v>
       *                   - <div u, q> = 0
       */
      printf("Assembling the matrix and right-hand side:\n");fflush(stdout);
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

      // Get the minimum mesh size
      int zz=0;
      double min_mesh_size = 10000000.;
      double tmp_size =0;
      mesh_struct *mesh_2 = &mesh;
  
      for (zz=0; zz<mesh_2->nface; zz++) {

	tmp_size=mesh_2->f_area[zz];

	if(tmp_size < min_mesh_size)
	  min_mesh_size = tmp_size;    
      }
  
      mesh_size_per_cycle[cycle] = min_mesh_size;
      //
      //      fprintf(stdout,"\n------ CHECK: 5\n");fflush(stdout);
      //
      assemble_global_block_neighbor(&A,&b,				\
				     local_assembly_Elasticity_FACE,	\
				     local_assembly_Elasticity,		\
				     FEM_Block_RHS_Local_Elasticity,\
				     &FE,			    \
				     &mesh,			    \
				     cq,source2D,0.0);
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

  
      //      fprintf(stdout,"\n------ CHECK: 6\n");fflush(stdout);
      // Eliminate boundary conditions in matrix and rhs
      if(dim==2) {
	eliminate_DirichletBC_blockFE_blockA(bc2D,&FE,&mesh,&b,&A,0.0);
      }
      if(dim==3) {
	eliminate_DirichletBC_blockFE_blockA(bc3D,&FE,&mesh,&b,&A,0.0);
      }
  
  
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
      /*******************************************************************************************/

      /***************** Solve *******************************************************************/
      printf("Solving the System:\n");fflush(stdout);
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
      /***********************************************************************/
      AMG_param amgparam;
      param_amg_init(&amgparam);
      //  param_amg_set(&amgparam, &inparam);
      //  param_amg_print(&amgparam);
      // get preconditioner diagonal blocks
      // Prepare diagonal blocks
      dCSRmat *A_diag;
      A_diag = (dCSRmat *)calloc(dim+1, sizeof(dCSRmat));
      for(i=0;i<dim+1;i++){ // copy block diagonal to A_diag
	dcsr_alloc(A.blocks[i*(dim+1)+i]->row, \
		   A.blocks[i*(dim+1)+i]->col, \
		   A.blocks[i*(dim+1)+i]->nnz, \
		   &A_diag[i]);
	dcsr_cp(A.blocks[i*(dim+1)+i], &A_diag[i]);
      }
      // for the last block, we only take its diagonal:
      INT nzd=0,j,k;
      for(i=0;i<A_diag[dim].row;i++){
	for(k=A_diag[dim].IA[i];k<A_diag[dim].IA[i+1];k++){
	  j=A_diag[dim].JA[k];
	  if(i!=j) continue;
	  A_diag[dim].JA[nzd]=A_diag[dim].JA[k];//=i
	  A_diag[dim].val[nzd]=A_diag[dim].val[k];//aii
	  nzd++;
	}
      }
      for(i=0;i<A_diag[dim].row;i++) A_diag[dim].IA[i+1]=A_diag[dim].IA[i]+1;
      // realloc them:
      A_diag[dim].JA=realloc(A_diag[dim].JA,nzd*sizeof(INT));
      A_diag[dim].val=realloc(A_diag[dim].val,nzd*sizeof(REAL));
      // finally: set nnz:
      A_diag[dim].nnz=nzd;
  
      ////////////////////////////
      //PRINTING
      ////////////////////////////
      /* FILE *fptmp; */
      /* fptmp=fopen("output/a.dat","w"); */
      /* bdcsr_print_matlab(fptmp,&A); */
      /* fclose(fptmp); */
      /* //////////////////////////// */
      /* fptmp=fopen("output/d1.dat","w"); */
      /* csr_print_matlab(fptmp,&A_diag[0]); */
      /* fclose(fptmp); */
      /* fptmp=fopen("output/d2.dat","w"); */
      /* csr_print_matlab(fptmp,&A_diag[1]); */
      /* fclose(fptmp); */
      /* fptmp=fopen("output/d3.dat","w"); */
      /* csr_print_matlab(fptmp,&A_diag[2]); */
      /* fclose(fptmp); */
      /* fptmp=fopen("output/ndofs.m","w"); */
      /* fprintf(fptmp,							\ */
      /* 	      "n_dofs=[%d %d %d];\n",FE_ux.ndof,FE_uy.ndof,FE_p.ndof); */
      /* fprintf(fptmp,							\ */
      /* 	      "n_row=[%d %d %d];\n",A_diag[0].row,A_diag[1].row,A_diag[2].row); */
      /* fprintf(fptmp,							\ */
      /* 	      "n_col=[%d %d %d];\n",A_diag[0].col,A_diag[1].col,A_diag[2].col); */
      /* fprintf(fptmp,							\ */
      /* 	      "n_nnz=[%d %d %d];\n",A_diag[0].nnz,A_diag[1].nnz,A_diag[2].nnz); */
      /* fclose(fptmp); */
      /////////////////////////////////////
      //END OF PRINTING  
      /////////////////////////////////////
  
      // SOLVING
      solver_flag = linear_solver_bdcsr_krylov_block_3(&A,&b,&sol,	\
						       &linear_itparam, \
						       &amgparam, A_diag);  
      /***********************************************************************/
      //////////////////////////////////////////////////////////////
      /* void *numeric=NULL;  // prepare for direct solve: */
      /* numeric=(void *)block_factorize_UMF(&A,0);//inparam.print_level); */
      /* solver_flag=(INT )block_solve_UMF(&A, */
      /* 				    &b, // this is the rhs here.  */
      /* 				    &sol,  // this is the solution here.  */
      /* 				    numeric, */
      /* 				    0);//     inparam.print_level); */
      /////////////////////////////////////////////////////////////////
      /*
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
      */
      // Error Check
      if (solver_flag < 0) printf("### ERROR: Solver does not converge with error code = %d!\n",solver_flag);

      clock_t clk_solve_end = clock();
      printf("Elapsed CPU Time for Solve = %f seconds.\n\n",
	     (REAL) (clk_solve_end-clk_solve_start)/CLOCKS_PER_SEC);
      /*******************************************************************************************/

      //////////////////////////////////////////////////////////////////////////////////////////////
      /********************* Compute Errors if you have exact solution ****************************/
      clock_t clk_error_start = clock();
      REAL* solerrL2 = (REAL *) calloc(dim+1, sizeof(REAL));
      REAL* solerrH1 = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_stress = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
 
      L2error_block(solerrL2, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
      //HDerror_block(solerrH1, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      HDsemierror_block(solerrH1, sol.val, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      //NEW SLEE Aug 17 2020
      HDsemierror_block_Stress(solerr_stress, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);

      REAL uerrL2 = 0;
      REAL uerrH1 = 0;
      REAL uerr_stressH1 = 0;
  
      for(i=0;i<dim;i++)
	uerrL2 += solerrL2[i]*solerrL2[i];
      for(i=0;i<dim;i++)
	uerrH1 += solerrH1[i]*solerrH1[i];
      for(i=0;i<dim;i++)
	uerr_stressH1 += solerr_stress[i]*solerr_stress[i];

      uerrL2 = sqrt(uerrL2);
      uerrH1 = sqrt(uerrH1);
      uerr_stressH1 = sqrt(uerr_stressH1);

      //REAL perrL2 = solerrL2[dim];
      //REAL perrH1 = solerrH1[dim];

      printf("*******************************************************\n");
      printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2);
      printf("[CG] H1 Norm of u error    = %26.13e\n",uerrH1);
      printf("[CG] Stress Norm of u error    = %26.13e\n",uerrH1);
      printf("*******************************************************\n\n");

      //Jul. 10. 2020 SLEE save the errors for convergence computation
      L2_error_per_cycle[cycle] = uerrL2; 
      H1_error_per_cycle[cycle] = uerrH1;
      H1_stress_error_per_cycle[cycle] = uerr_stressH1; 
  
      //L2_error_p_per_cycle[cycle] = perrL2; 
      //NEW SLEE Aug 23 2020
      //NEW ERROR FOR EG

      REAL* solerrL2_EG = (REAL *) calloc(dim+1, sizeof(REAL));
      REAL* solerrH1_EG = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_stress_EG = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements
      REAL* solerr_energy_EG = (REAL *) calloc(dim+1, sizeof(REAL)); // Note: No H1 error for P0 elements

      printf("---1\n");
      L2error_block_EG(solerrL2_EG, sol.val, exact_sol2D, &FE, &mesh, cq, 0.0);
      printf("---2\n");
      HDerror_block_EG(solerrH1_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      printf("---3\n");
      HDsemierror_block_Stress_EG(solerr_stress_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      printf("---4\n");
      //HDsemierror_block_EnergyNorm_EG(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      HDsemierror_block_EnergyNorm_EG_FaceLoop(solerr_energy_EG, sol.val, exact_sol2D, Dexact_sol2D, &FE, &mesh, cq, 0.0);
      printf("---5\n");
  
      //NEW SLEE Aug 17 2020
      //NEW ERROR FOR STRESS, \mu \epsilon(u) + \lambda \div u
      REAL uerrL2_EG = 0;
      REAL uerrH1_EG = 0;
      REAL uerr_stress_EG = 0;
      REAL uerr_energy_EG = 0;

  
      for(i=0;i<dim;i++)
	uerrL2_EG += solerrL2_EG[i]*solerrL2_EG[i];
      for(i=0;i<dim;i++)
	uerrH1_EG += solerrH1_EG[i]*solerrH1_EG[i]; 
      for(i=0;i<dim;i++)
	uerr_stress_EG += solerr_stress_EG[i]*solerr_stress_EG[i]; 
      //for(i=0;i<dim;i++)
      //uerr_energy_EG += solerr_energy_EG[i]*solerr_energy_EG[i]; 

      uerr_energy_EG = solerr_energy_EG[0]+solerr_energy_EG[1];
  

      //DEBUG
      uerr_energy_EG += uerrH1_EG;
      //
 
      uerrL2_EG = sqrt(uerrL2_EG);
      uerrH1_EG = sqrt(uerrH1_EG);
      uerr_stress_EG = sqrt(uerr_stress_EG);
      uerr_energy_EG = sqrt(uerr_energy_EG);
  
      L2_EG_error_per_cycle[cycle] = uerrL2_EG;
      H1_EG_error_per_cycle[cycle] = uerrH1_EG;
      H1_stress_EG_error_per_cycle[cycle] = uerr_stress_EG;
      H1_energy_EG_error_per_cycle[cycle] = uerr_energy_EG;

      printf("# of elements = %d \n", mesh.nelm); 
      printf("*******************************************************\n");
      printf("L2 Norm of u (EG) error    = %26.13e\n",uerrL2_EG);
      printf("H1 Norm of u (EG) error    = %26.13e\n",uerrH1_EG);
      printf("Stress Norm of u (EG) error    = %26.13e\n",uerr_stress_EG);
      printf("Energy Norm of u (EG) error    = %26.13e\n",uerr_energy_EG);
      printf("*******************************************************\n");
  
  
      /////////////////////////////////////
      /*fprintf(stdout,"\n%d,%d\n\n",A.brow,A.bcol);
	INT j;
	for(i=0;i<A.brow;i++){ 
	for(j=0;j<A.brow;j++){ 
	fprintf(stdout,"\n(%d,%d):::%d,%d,%d\n\n",i,j,	 
	A.blocks[j*A.bcol+i]->row,		 
	A.blocks[j*A.bcol+i]->col,		 
	A.blocks[j*A.bcol+i]->nnz); 
	} 
	}
      */ 
      /* fflush(stdout); */
  
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

      printf("OUTPUT?\n");

  
      if(inparam.print_level > 3){

	char output_filename_per_cycle[512]={'\0'};
	sprintf( output_filename_per_cycle, "output/solution_%d.vtu", cycle);
	char* soldump = output_filename_per_cycle;//"output/solution.vtu";
	//char* soldump = "output/solution.vtu";
   
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
	H1_conv_rate_per_cycle[tmp] = 0;
	H1_stress_conv_rate_per_cycle[tmp] = 0;

	L2_EG_conv_rate_per_cycle[tmp] = 0;
	H1_EG_conv_rate_per_cycle[tmp] = 0;
	H1_stress_EG_conv_rate_per_cycle[tmp] = 0;
	H1_energy_EG_conv_rate_per_cycle[tmp] = 0;

      }
      else{
	// multiplied dim since we use DOFs not h here.
	L2_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_error_per_cycle[tmp]) -log(L2_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );  
	H1_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_error_per_cycle[tmp]) -log(H1_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) ); 
	H1_stress_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_error_per_cycle[tmp]) -log(H1_stress_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_CG[tmp-1]) -log(dof_per_cycle_CG[tmp]) );

	L2_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(L2_EG_error_per_cycle[tmp]) -log(L2_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );
	H1_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_EG_error_per_cycle[tmp]) -log(H1_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	H1_stress_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_stress_EG_error_per_cycle[tmp]) -log(H1_stress_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

	H1_energy_EG_conv_rate_per_cycle[tmp] = global_dim_space * (log(H1_energy_EG_error_per_cycle[tmp]) -log(H1_energy_EG_error_per_cycle[tmp-1]) ) /  (log(dof_per_cycle_EG[tmp-1]) -log(dof_per_cycle_EG[tmp]) );

      }
      
      printf("L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_error_per_cycle[tmp], dof_per_cycle_CG[tmp],L2_conv_rate_per_cycle[tmp]);
      printf("H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_conv_rate_per_cycle[tmp]);
      printf("H1 Stress Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_error_per_cycle[tmp], dof_per_cycle_CG[tmp],H1_stress_conv_rate_per_cycle[tmp]);
      
      printf("EG - L2 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, L2_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     L2_EG_conv_rate_per_cycle[tmp]);
      printf("EG - H1 Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_EG_conv_rate_per_cycle[tmp]);
   
      printf("EG - Stress_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_stress_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_stress_EG_conv_rate_per_cycle[tmp]);

      printf("EG - EnergyNorm_Error for cycle %d = %f  with %d DOFs  -- convergence rate  =  %f\n", tmp, H1_energy_EG_error_per_cycle[tmp], dof_per_cycle_EG[tmp],
	     H1_energy_EG_conv_rate_per_cycle[tmp]);
 
      printf("----------------------------------------------------\n");      
    }

  
  //for Latex Print Table
  printf("** LATEX TABLE ** \n");
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
      printf("%d & %f & %f &  %f &  %f   &  %f &  %f \\\\ \\hline \n", dof_per_cycle_CG[tmp], L2_error_per_cycle[tmp], L2_conv_rate_per_cycle[tmp],H1_error_per_cycle[tmp], H1_conv_rate_per_cycle[tmp],H1_stress_error_per_cycle[tmp], H1_stress_conv_rate_per_cycle[tmp] ); 
    }
  
  for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
      printf("%f \n",  mesh_size_per_cycle[tmp]);
    }

  for(int tmp=0; tmp<total_num_cycle; ++tmp)
    {
      printf("%d & %f & %f &  %f &  %f   &  %f &  %f  & %f & %f \\\\ \\hline \n",  dof_per_cycle_EG[tmp],L2_EG_error_per_cycle[tmp], L2_EG_conv_rate_per_cycle[tmp],H1_EG_error_per_cycle[tmp], H1_EG_conv_rate_per_cycle[tmp],
	     H1_stress_EG_error_per_cycle[tmp], H1_stress_EG_conv_rate_per_cycle[tmp], H1_energy_EG_error_per_cycle[tmp], H1_energy_EG_conv_rate_per_cycle[tmp] ); 
    }
  

  //printf("%f\n",H1_EG_error_per_cycle[total_num_cycle - 1]);
  //printf("%f\n",H1_EG_conv_rate_per_cycle[total_num_cycle - 1]);

 
  //printf("%f\n",H1_energy_EG_error_per_cycle[total_num_cycle - 1]);
  //printf("%f\n",H1_energy_EG_conv_rate_per_cycle[total_num_cycle - 1]);
  
  
  return 0;
  
  
}  /* End of Program */
/*********************************************************************************************/
