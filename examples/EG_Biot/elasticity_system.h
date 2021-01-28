/*********************************************************************************/
iCSRmat *extend_el_dof(fespace *FE,mesh_struct *mesh, INT *ix)
{
  // extends the DOF so that every DOF in element i interacts with all
  // DOFs of the elements sharing face with i. It is DG like stiffness
  // matrix. ix is working vector of at least FE->el_dofs->col elements
  iCSRmat *el_dof=NULL; // output
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;
  INT pq,nbr,nel;
  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
  iCSRmat *f_el=malloc(1*sizeof(iCSRmat)); // face_to_element;
  iCSRmat *el_f=mesh->el_f;
  for(i=0;i<el_f->IA[mesh->nelm];i++)
    el_f->val[i]=1;
  icsr_trans(el_f,f_el); // f_el=transpose(el_f);
  iCSRmat *el2el=malloc(1*sizeof(iCSRmat));
  icsr_mxm(el_f,f_el,el2el);
  //  icsr_tri(el2el,'l');
  icsr_nodiag(el2el);// remove diagonal.
  icsr_free(f_el);
  ///////////////////////////////////////////////
  INT nrows = FE->ndof;
  INT ncols = FE->ndof;
  /////////////
  for(i=0;i<ncols;i++) ix[i]=-1;
  ////////
  icp=0;
  for (nel=0; nel<el2el->row; nel++) {
    j_a=FE->el_dof->IA[nel];
    j_b=FE->el_dof->IA[nel+1];
    for(j=j_a;j<j_b;j++){
      k=FE->el_dof->JA[j];
      icp++;
      ix[k]=nel;
    }
    for(pq=el2el->IA[nel];pq<el2el->IA[nel+1];pq++){	
      nbr=el2el->JA[pq];
      j_a=FE->el_dof->IA[nbr];
      j_b=FE->el_dof->IA[nbr+1];
      for(j=j_a;j<j_b;j++){
	k=FE->el_dof->JA[j];
	if(ix[k] != nel){
	  icp++;
	  ix[k]=nel;
	}
      }
    }
  }
  j_a=FE->el_dof->row;
  j_b=FE->el_dof->col;
  ////////////////
  el_dof=(iCSRmat *)icsr_create_p(j_a,j_b,icp);
  ////////////////
  for(i=0;i<ncols;i++) ix[i]=-1;
  icp=0;
  // REPEAT THE LOOP but now store nonzeroes
  el_dof->IA[0]=icp;
  for (nel=0; nel<el2el->row; nel++) {
    j_a=FE->el_dof->IA[nel];
    j_b=FE->el_dof->IA[nel+1];
    for(j=j_a;j<j_b;j++){
      k=FE->el_dof->JA[j];
      el_dof->JA[icp]=k;
      icp++;
      ix[k]=nel;
    }
    //    fprintf(stdout,"\nneighbors[%d]=(", nel);      
    for(pq=el2el->IA[nel];pq<el2el->IA[nel+1];pq++){	
      nbr=el2el->JA[pq];
      //      fprintf(stdout,"%d ", nbr);
      j_a=FE->el_dof->IA[nbr];
      j_b=FE->el_dof->IA[nbr+1];
      for(j=j_a;j<j_b;j++){
	k=FE->el_dof->JA[j];
	if(ix[k] != nel){
	  el_dof->JA[icp]=k;
	  icp++;	  
	  ix[k]=nel;
	}
      }
    }
    el_dof->IA[nel+1]=icp;
    //    fprintf(stdout,")");      
  }  
  //  fprintf(stdout,"\n\nnnz[el2el]=%d;nnz[el_dof]=%d;check=%d\n",el2el->nnz,el_dof->nnz,el_dof->IA[el_dof->row]);
  // printing only
  /* for(nel=0;nel<mesh->nelm;nel++){ */
  /*   fprintf(stdout,"\nold-dofs[%d]=(", nel);       */
  /*   j_a=FE->el_dof->IA[nel]; */
  /*   j_b=FE->el_dof->IA[nel+1]; */
  /*   for(j=j_a;j<j_b;j++){ */
  /*     k=FE->el_dof->JA[j]; */
  /*     fprintf(stdout,"%d ", k); */
  /*   } */
  /*   fprintf(stdout,")");       */
  /*   fprintf(stdout,"\nNEW-dofs[%d]=(", nel);       */
  /*   j_a=el_dof->IA[nel]; */
  /*   j_b=el_dof->IA[nel+1]; */
  /*   for(j=j_a;j<j_b;j++){ */
  /*     k=el_dof->JA[j]; */
  /*     fprintf(stdout,"%d ", k); */
  /*   } */
  /*   fprintf(stdout,")");       */
  /* } */
  /* fprintf(stdout,"\n\n");       */
  /* exit(255); */
  return el_dof;
}
/*********************************************************************************/
void sl_create_CSR_cols_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2, \
			       mesh_struct *mesh)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;
  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2,v_f;
  //ltz  icsr_trans(FE2->el_dof,&dof_el_2);
  /*YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY*/
  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;
  INT *ix = NULL;
  if(nrows<ncols)
    ix=(INT *) calloc(ncols,sizeof(INT));
  else
    ix=(INT *) calloc(nrows,sizeof(INT));
    
  iCSRmat *el_dof=extend_el_dof(FE2,mesh,ix);
  ////////////////////////////////////////
  icsr_trans(el_dof,&dof_el_2);
  //////////////////////////////////
  for (i=0; i<ncols; i++) {
    ix[i] = -1;
  }
  // Loop over all DOF of test space and build A->JA
  icp=0;
  for (i=0; i<nrows; i++) {
    // Loop over all Elements connected to particular edge
    j_a = dof_el_2.IA[i];
    j_b = dof_el_2.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el_2.JA[j];
      k_a = FE1->el_dof->IA[if1];
      k_b = FE1->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE1->el_dof->JA[k];
        if (ix[mydof]!=i) { /* We haven't been here */
	  A->JA[icp] = mydof;
        icp++;
        ix[mydof] = i;
	}
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el_2);
  free(el_dof); // no need of icsr_free here;
return;
}
/********************************************************************************************/
/********************************************************************************************/
/****************************************************************************/
void sl_create_CSR_rows_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2,	\
			       mesh_struct *mesh)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2;
  //ltz  icsr_trans(FE2->el_dof,&dof_el_2);
  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;
  INT *ix = NULL;
  if(nrows<ncols)
    ix=(INT *) calloc(ncols,sizeof(INT));
  else
    ix=(INT *) calloc(nrows,sizeof(INT));
  iCSRmat *el_dof=extend_el_dof(FE2,mesh,ix);
  /////////////////////////////////////
  icsr_trans(el_dof,&dof_el_2);
  /////////////////////////////////////
  for (i=0; i<ncols; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF of test space and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=0;
  for (i=0; i<nrows; i++) {
    A->IA[i] = icp;
    // Loop over all Elements connected to particular DOF of test space
    j_a = dof_el_2.IA[i];
    j_b = dof_el_2.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el_2.JA[j];
      // For this given element grab the DOF in the trial space
      k_a = FE1->el_dof->IA[if1];
      k_b = FE1->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE1->el_dof->JA[k];
        if (ix[mydof]!=i) { /* We haven't been here  */
        icp++;
	//	fprintf(stdout,"\n%d:%d;mydof=%d;;;r=%d;c=%d",i,if1,mydof,nrows,ncols); fflush(stdout);
        ix[mydof] = i;
      }
    }
  }
}
  A->IA[nrows] = icp;
  A->nnz = icp;

 if(ix) free(ix);
 icsr_free(&dof_el_2);
 free(el_dof); // no need of icsr_free here;
 
return;
}
/******/
/******/
void block_LocaltoGlobal_RHS(INT *dof_on_elm,block_fespace *FE,dvector *b,REAL *bLoc)
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
        // Adjust Right-hand side globally
        if(bLoc!=NULL && block_col==0)
        b->val[local_row+global_row_index] += bLoc[local_row_index+i];
      }
      local_col_index += dof_per_elm_trial;
    }
    local_col_index = 0;
    global_row_index += FE->var_spaces[block_row]->ndof;
    local_row_index += dof_per_elm_test;
  }
  return;
}
/******************************************************************************************************/

void block_LocaltoGlobal_neighbor(INT *dof_on_elm,INT *dof_on_elm_neighbor,block_fespace *FE,block_dCSRmat *A,REAL *ALoc)
{
  INT i,j,k,col_a,col_b,acol,arow,block_row,block_col,row_a, row_b;
  INT local_row,local_col;

  // Loop over all the blocks
  INT nblocks = FE->nspaces;
  INT dof_per_elm_test = 0;
  INT dof_per_elm_trial = 0;
  INT local_row_index = 0;
  INT local_col_index = 0;
  INT global_row_index = 0;
  INT block_dof_per_elm = 0;

  /*
  for(i = 0; i<7; i++)
    printf("dof_on_elm[i] =  %d \n", dof_on_elm[i]);

  for(i = 0; i<7; i++)
    printf("dof_on_elm_neighbor[i] =  %d \n", dof_on_elm_neighbor[i]);
  */
  
  // Get total dof_per_elm for indexing
  for(block_row=0;block_row<nblocks;block_row++) {
    block_dof_per_elm += FE->var_spaces[block_row]->dof_per_elm;
  }


  //printf("*** ALoc[48] = %f \n", ALoc[48]);
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
	  /*
	  printf("\n ===================================\n");
	  printf("-[blockrow, blockcol] =[%d, %d] - (i,j) = (%d,%d) local_row_index = %d, local_row = %d || local_col_index = %d, local_col = %d \n",
	  	 block_row, block_col, i,j, local_row_index, local_row, local_col_index, local_col); 
	  printf("===================================\n");
	  */
          /* Columns of A */
          if(A->blocks[block_row*nblocks+block_col]) {
	    
	    //printf("NOW in the GLOBAL MATRIX = A->blocks[%d], \n", block_row*nblocks+block_col);
	    
	    col_a = A->blocks[block_row*nblocks+block_col]->IA[local_row];
	    col_b = A->blocks[block_row*nblocks+block_col]->IA[local_row+1];

	    for (k=col_a; k<col_b; k++) {
	      	     
	      acol = A->blocks[block_row*nblocks+block_col]->JA[k];

	      //arow = A->blocks[block_row*nblocks+block_col]->IA[k];
	      //if(arow == local_row)
	      //{
	      //  printf("######################## ROW :: [k=%d] --  arow = %d || local_row = %d\n ",k,arow,local_row);       
	      //}
	      //printf("$$ COL :: [k=%d] --  acol = %d || local_col = %d\n ",k,acol,local_col); 
		
	      if (acol==local_col)
		{
		  /* If they match, put it in the global matrix */
		  A->blocks[block_row*nblocks+block_col]->val[k]
		    += ALoc[(local_row_index+i)*block_dof_per_elm+(local_col_index+j)];
		  
		  //printf("A[%d](%d,%d) = %f \n", block_row*nblocks+block_col, local_row, acol, A->blocks[block_row*nblocks+block_col]->val[k]);
		  //printf("ALoc[%d] = %f \n",(local_row_index+i)*block_dof_per_elm+(local_col_index+j), ALoc[(local_row_index+i)*block_dof_per_elm+(local_col_index+j)]);
		  //printf("$$ COL :: [k=%d] --  acol = %d || local_col = %d\n ",k,acol,local_col); 
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
  /*
  FILE *fptmp = NULL;
 
  fptmp=fopen("output/a.dat","w");
  bdcsr_print_matlab(fptmp,A);
  fclose(fptmp);
  
  exit(0);
  */
return;
}

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
void assemble_global_block_neighbor(block_dCSRmat* A,dvector *b, REAL *solution, void (*local_assembly_face)(block_dCSRmat* A,block_fespace *,mesh_struct *,qcoordinates *,INT, iCSRmat*, REAL,REAL), \
				    void (*local_assembly)(block_dCSRmat* A,dvector *b,REAL *,block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT,REAL,REAL,INT*), \
				    void (*local_rhs_assembly)(dvector *b,REAL *, REAL *, block_fespace *,mesh_struct *,qcoordinates *,INT *,INT *,INT, \
							       void (*)(REAL *,REAL *,REAL,void *),REAL,REAL, \
							       void (*)(REAL *,REAL *,REAL,void *),void (*)(REAL *,REAL *,REAL,void *),void (*)(REAL *,REAL *,REAL,void *), REAL *,INT *), \
				    block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, \
				    void (*rhs)(REAL *,REAL *,REAL,void *), \
				    void (*truesol)(REAL *,REAL *,REAL,void *),\
				    void (*D_truesol)(REAL *,REAL *,REAL,void *), \
				    void (*D_truesol_dt)(REAL *,REAL *,REAL,void *), \
				    REAL time,		\
				    REAL timestep
				   )
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
        //create_CSR_rows_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);
	sl_create_CSR_rows_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i],mesh);

        // Columns of A -> JA
        A->blocks[i*nblocks+j]->JA = (INT *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(INT));
        //create_CSR_cols_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);
	sl_create_CSR_cols_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i],mesh);

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

  
  REAL* tmperrL2 = (REAL *) calloc(2+1+1+1, sizeof(REAL));

  for(i=0;i<FE->nspaces;i++) {
    tmperrL2[i]=0;
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
    //printf("====== Ele = %d \n", i);    
    //for(INT kk=0;kk<7;++kk)
    //printf("dof_on_elm[%d]=%d \n", kk, dof_on_elm[kk]);
    // Find vertices for given Element
    get_incidence_row(i,mesh->el_v,v_on_elm);


    
    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(A,b,ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,time,timestep,switch_on_face);
    if(rhs!=NULL)
      (*local_rhs_assembly)(b,bLoc,solution,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time,timestep,truesol,D_truesol, D_truesol_dt, tmperrL2,switch_on_face);


    //printf("local_size = %d \n", local_size);
    //for (j=0; j<local_size; j++) {
    //  printf(" ALoc[%d]=%f \n", j, ALoc[j]);
    //}
    //exit(0);

    // Loop over DOF and place in appropriate slot globally
    //block_LocaltoGlobal(dof_on_elm,FE,b,A,ALoc,bLoc);
    // --> this is now embedded for each assembly. 
  }// loop over elements

  
  iCSRmat *f_el=NULL;
  bool bEG = true;//false;
  // NOW ASSEMBLE FACES
  if(bEG){
    // First generate global face-to-element map:
    f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
    icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);
    //LOOP OVER FACES    
    for (i=0; i<mesh->nface; i++) {
      (*local_assembly_face)(A,FE,mesh,cq,i,f_el,time,timestep);
      //      fprintf(stdout,"\n------ FACE: %7d\n",mesh->nface-i);fflush(stdout);   
      //if(i == 3)
      //exit(0);
      
    }
    icsr_free(f_el);
  }


  /*
  REAL uerrL2 = 0;
  REAL uerrL2_p = 0;


  for(i=0;i<2;i++) {
    tmperrL2[i] = sqrt(tmperrL2[i]);
  }
  
  for(i=0;i<2;i++)
    uerrL2 += tmperrL2[i]*tmperrL2[i];
  uerrL2_p = tmperrL2[3];
	      
  
  printf("$$$$$$$$$$$$$$$ MECHANCIS   *****************************\n");
  printf("[CG] L2 Norm of u error    = %26.13e\n",uerrL2);	      
  printf("$$$$$$$$$$$$$$$ PRESSURE    *****************************\n");
  printf("[CG] L2 Norm of p error    = %26.13e\n",uerrL2_p);
  printf("$$$$$$$$$$$$$$$****************************************\n\n");
  */	      

	      
  
  //printf("--end  \n");
  
  /*
  FILE *fptmp = NULL;
 
  fptmp=fopen("output/a.dat","w");
  bdcsr_print_matlab(fptmp,A);
  fclose(fptmp);
  */
  //exit(0);
  //  
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}
