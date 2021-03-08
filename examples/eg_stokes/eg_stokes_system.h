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
/******************************************************************************/

/*************** LOCAL STUFF **************************************************/

/**************************************************************/
void local_assembly_Elasticity_FACE(block_dCSRmat* A, block_fespace *FE, \
				    mesh_struct *mesh, qcoordinates *cq, \
				    INT face, iCSRmat *f_el,REAL time,	\
				    REAL timestep)
{

  //printf("Local Assemble Face - start\n");
  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;
  REAL alpha = 1.;

  // INDEX
  INT i,j,k;
  INT dim = mesh->dim;
  INT rowa,rowb;

  // Mesh Stuff
  INT nspaces = FE->nspaces;

  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));

  INT test, trial;
  INT local_row_index, local_col_index;
  REAL kij = 0.0;

  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  INT local_size = dof_per_elm*dof_per_elm;

  //Local Matrix for Each Cases on Faces
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

  double lambda = LAME_LAMBDA_GLOBAL;

  INT quad,quad_face;
  INT jk, jkl,jcntr,ed;

  INT rowa_neighbor,rowb_neighbor,jcntr_neighbor;
  INT k_neighbor, j_neighbor;
  double fiarea;

  // For Locking-Free EG
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
  // Sanity Check
  /*
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
  */

  // Saniyy Check
  // counter == 1 means the BD
  if(counter == 1 && mesh->f_flag[face] == 0)
    {printf("check-error FALSE\n"); exit(0);}

  if(neighbor_index[0] == -1)
    {
      printf("Something is Wrong! \n");
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

  // Sanity Check
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

  //Sanity Check; Print the normal vectors..
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
  double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
  double BC_penalty_term = BC_PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

  double penalty_term_pressure = PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
  double BC_penalty_term_pressure = BC_PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));


  // DEBUG mode
  // Get the BD values
  //REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));

  for (quad_face=0;quad_face<nquad_face;quad_face++) {

    qx_face[0] = cq_face->x[quad_face];
    qx_face[1] = cq_face->y[quad_face];

    if(dim==3) qx_face[2] = cq_face->z[quad_face];

    w_face = cq_face->w[quad_face];

    //DEBUG mode
    //(*exact_sol2D)(val_true_face,qx_face,time,
    //	     &(mesh->el_flag[neighbor_index[0]])); //DUMMY VALUE at the end // ???


    ///////////////////////////////////////////////////
    if(counter == 1){ //only at boundary
      //printf("== BD = %d \n", face);
      //Sanity Check
      if(mesh->f_flag[face] == 0)
	{printf("check-error0\n"); exit(0);}

      bool bWeakBC_FACE = BOOL_WEAKLY_IMPOSED_BC;

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


	get_FEM_basis(FE->var_spaces[3]->phi,
		      FE->var_spaces[3]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face,
		      mesh,
		      FE->var_spaces[3]);


	//for the P part
	local_dof_on_elm_face += FE->var_spaces[3]->dof_per_elm;

	//---- #2 BOUNDARY ---------------------------- START //
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
	    kij += BC_penalty_term * FE->var_spaces[0]->phi[trial] * FE->var_spaces[0]->phi[test];

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
	    kij += BC_penalty_term *  FE->var_spaces[1]->phi[trial]  * FE->var_spaces[1]->phi[test];

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
	    kij += BC_penalty_term * (qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
	    kij += BC_penalty_term * (qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);


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

	} // bEG
	//----- #2 BOUNDARY ------------------------------------------- END //


	//---- Poro P2 Boundary -------------------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// porop2
	// On BOUNDARY
	// PoroElasticity   + ({p},  [u] dot n)
	// Mechanics coupeld with pressure
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG u0
	//       trial (column) 3  test (row) 0
	// test - u0 cg
	local_row_index = 0;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - u0 cg
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	if(bEG){

	  //      ({p},  [u] dot n)
	  //       CG  EG
	  //       trial (column) 3  test (row) 2
	  // test - u0 eg
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  // trial - p cg
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	      kij = alpha * FE->var_spaces[3]->phi[trial] *
		(finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	}

	//---- Poro P2 Boundary -------------------------------------- END //


	//---- Stokes P2 Boundary -------------------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// On BOUNDARY
	// STokes  + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG u0
	//       trial (column) 3  test (row) 0
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// trial - u0 cg
	local_col_index = 0;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}
	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// trial - u0 cg
	local_col_index = FE->var_spaces[0]->dof_per_elm;

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * FE->var_spaces[3]->phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	if(bEG){

	  //      ({p},  [u] dot n)
	  //       CG  EG
	  //       trial (column) 3  test (row) 2
	  // test - p cg
	  local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	  // trial - u0 cg
	  local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
	  // Loop over Test Functions (Rows)
	  for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	    // Loop over Trial Functions (Columns)
	    for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	      kij = alpha * FE->var_spaces[3]->phi[test] *
		(finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	    }
	  }

	}

	//---- Stokes P2 Boundary -------------------------------------- END //

      }// Weakly Face

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
	  {printf("well? something is wrong \n"); exit(0);}


	//printf("counter 2 \n");

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

	get_FEM_basis(FE->var_spaces[3]->phi,
		      FE->var_spaces[3]->dphi,
		      qx_face,
		      v_on_elm,
		      local_dof_on_elm_face_interface,
		      mesh,
		      FE->var_spaces[3]);

	//FOR P
	local_dof_on_elm_face_interface += FE->var_spaces[3]->dof_per_elm;

	// Now, get the basis for the neighbor
	REAL* neighbor_basis_0_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_0_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));

	REAL* neighbor_basis_1_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_1_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));

	REAL* neighbor_basis_3_phi  = (REAL *) calloc(3,sizeof(REAL));
	REAL* neighbor_basis_3_dphi = (REAL *) calloc(3*mesh->dim,sizeof(REAL));


	get_FEM_basis( neighbor_basis_0_phi ,
		       neighbor_basis_0_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[0]);

	local_dof_on_elm_neighbor = dof_on_elm_neighbor + FE->var_spaces[0]->dof_per_elm;

	//printf("counter 111 \n");

	get_FEM_basis( neighbor_basis_1_phi ,
		       neighbor_basis_1_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       local_dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[1]);

	// u2
	local_dof_on_elm_neighbor += FE->var_spaces[1]->dof_per_elm;

	//printf("counter 222 \n");
	// u eg
	local_dof_on_elm_neighbor += FE->var_spaces[2]->dof_per_elm;


	get_FEM_basis( neighbor_basis_3_phi ,
		       neighbor_basis_3_dphi ,
		       qx_face,
		       v_on_elm_neighbor,
		       local_dof_on_elm_neighbor,
		       mesh,
		       FE->var_spaces[3]);

	// p
	local_dof_on_elm_neighbor += FE->var_spaces[3]->dof_per_elm;


	// ------ #2 INTERFACE -------------------------- START //
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


	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


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
	// ------ #2 INTERFACE -------------------------- END //


	// ------ PORO P1 INTERFACE -------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	//porop1
	//PoroElasticity   + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG
	//       trial (column) 3  test (row) 0
	// test - u0 cg

	local_row_index = 0;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  neighbor_basis_0_phi[test]*finrm[0];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *  neighbor_basis_0_phi[test]*finrm[0];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// test - u1 cg
	local_row_index = FE->var_spaces[0]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *  neighbor_basis_1_phi[test]*finrm[1];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *  neighbor_basis_1_phi[test]*finrm[1];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}




	//porop1
	//      ({p},  [u] dot n)
	//       CG  DG
	//       trial (column) 3  test (row) 2
	// test - u eg
	local_row_index = FE->var_spaces[0]->dof_per_elm +FE->var_spaces[1]->dof_per_elm;
	// trial - p cg
	local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));
	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[trial] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}



	// ------ PORO P1 INTERFACE -------------------------- END//


	// ------ STOKES P1 INTERFACE -------------------------- START //
	///////////////////////////////////////////////////////////////////////////
	// Stokes   + ({p},  [u] dot n)
	//////////////////////////////////////////////////////////////////////////
	//      ({p},  [u] dot n)
	//       CG  CG
	//       trial (column) 3  test (row) 0
	// trial - u0 cg
	local_col_index = 0;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *  neighbor_basis_0_phi[trial]*finrm[0];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *  FE->var_spaces[0]->phi[trial]*finrm[0];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *  neighbor_basis_0_phi[trial]*finrm[0];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//      ({p},  [u] dot n)
	//       CG  CG u1
	//       trial (column) 3  test (row) 1
	// trial - u1 cg
	local_col_index = FE->var_spaces[0]->dof_per_elm;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *  neighbor_basis_1_phi[trial]*finrm[1];

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *  FE->var_spaces[1]->phi[trial]*finrm[1];

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *  neighbor_basis_1_phi[trial]*finrm[1];

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//Stokes1
	//      ({p},  [u] dot n)
	//       CG  DG
	//       test (column) 3  trial (row) 2
	// trial - u eg
	local_col_index = FE->var_spaces[0]->dof_per_elm +FE->var_spaces[1]->dof_per_elm;
	// test - p cg
	local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * FE->var_spaces[3]->phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));
	    ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * FE->var_spaces[3]->phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = alpha * 0.5 * neighbor_basis_3_phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));

	    ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}

	// Loop over Test Functions (Rows)
	for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	  // Loop over Trial Functions (Columns)
	  for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	    kij = -alpha * 0.5 * neighbor_basis_3_phi[test] *
	      (finrm[0] *  (qx_face[0]- barycenter_neighbor->x[0]) + finrm[1] *  (qx_face[1]- barycenter_neighbor->y[0]));

	    ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
	  }
	}


	//Stokes1


	// Stokes Interafec -- end //



	free(neighbor_basis_0_phi);
	free(neighbor_basis_0_dphi);

	free(neighbor_basis_1_phi);
	free(neighbor_basis_1_dphi);

	free(neighbor_basis_3_phi);
	free(neighbor_basis_3_dphi);
      }
  }// face q loop


  // printf("ok\n");
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





  //printf("Local Assemble Face - end\n");
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
 \return bLoc         Local RHS Vector
 *
 *
 */

// ISSUE : this doesn't allow to get the exact solution out...
void FEM_Block_RHS_Local_Elasticity(dvector *b,REAL* bLoc, REAL *solution, \
				    block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,
				    void (*rhs)(REAL *,REAL *,REAL,void *),REAL time, REAL timestep, \
				    void (*truesol)(REAL *,REAL *,REAL,void *),
				    void (*D_truesol)(REAL *,REAL *,REAL,void *),
				    void (*truesol_dt)(REAL *,REAL *,REAL,void *),
				    REAL *err, INT* switch_on_face)
{

  //printf("RHS ASSEMBLE  \n");
  //fflush(stdout);

  // Loop Indices
  INT i,quad,test;
  REAL C_0 = 0.;
  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_elm = 0;
  INT nun=0;

  REAL alpha = 1.;

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
  INT maxdim=nun;
  REAL qx[maxdim];

  // Right-hand side function at Quadrature Nodes
  REAL rhs_val[nun];

  // For initial condition
  REAL* val_true = (REAL *) calloc(nun,sizeof(REAL));
  REAL* val_true_D = (REAL *) calloc(nun*2,sizeof(REAL));

  // for the preivous time step value
  //REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));


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

    //if(time == timestep){

    //printf("********\n");
    //printf("********\n");
    //printf("***** time = %f ***\n", time);
    (*truesol)(val_true,qx,time - timestep,&(mesh->el_flag[elm]));
    (*D_truesol)(val_true_D,qx,time - timestep,&(mesh->el_flag[elm]));


    //if(time != timestep){
    //printf("qx[0] = %f, qx[1] = %f, time = %f, val_true_D[0] = %f, val_true_D[3] = %f \n",
    //	   qx[0], qx[1], time-timestep, val_true_D[0], val_true_D[3]);
    //exit(0);
    //}

    local_row_index=0;
    unknown_index=0;
    local_dof_on_elm=dof_on_elm;

    REAL* u_comp = solution;
    REAL u0_value_at_q = 0.;
    REAL u1_value_at_q = 0.;
    REAL val[dim];
    val[0] = 0.0;
    val[1] = 0.0;

    REAL p_value_at_q = 0.;
    REAL p_eg_value_at_q = 0.;
    int j,dof;

    REAL grad_val0[dim];
    grad_val0[0]=0;
    grad_val0[1]=0;
    REAL grad_val1[dim];
    grad_val1[0]=0;
    grad_val1[1]=0;

    // U  bEG
    REAL grad_val2[dim];
    grad_val2[0] = 0.0;
    grad_val2[1] = 0.0;

    REAL grad_val3[dim];
    grad_val3[0] = 0.0;
    grad_val3[1] = 0.0;

    // Pressure
    REAL grad_val4[dim];
    grad_val4[0]=0;
    grad_val4[1]=0;


    //printf("time = %f, timestep = %f, val_true[3] = %f, val_true[4] = %f \n", time, timestep, val_true[3], val_true[4]);

    for(i=0;i<FE->nspaces;i++) {

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->var_spaces[i]->phi,
		    FE->var_spaces[i]->dphi,
		    qx,
		    v_on_elm,
		    local_dof_on_elm,
		    mesh,
		    FE->var_spaces[i]);




      for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {

	//SLEE
	if(i == 0 || i == 1)
	  {
	    //This is for  u0 and u1 (CG part)
	    //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	    bLoc[(local_row_index+test)] += w*rhs_val[i]*FE->var_spaces[i]->phi[test];
	  }
	else if(i == 2)
	  {
	    //SLEE
	    // This is for  u2:: (EG part)
	    //Note that new DOF is \phi^3 = [x ; y]
	    //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
	    bLoc[(local_row_index+test)] += w*(rhs_val[0]*  (qx[0]-barycenter->x[0])  +rhs_val[1]* (qx[1]-barycenter->y[0]));
	  }
	else if(i == 3)
	  {

	    bLoc[(local_row_index+test)] += w*rhs_val[3]*FE->var_spaces[3]->phi[test];

	    //printf(" FE->var_spaces[3]->phi[test] = %f\n", FE->var_spaces[3]->phi[test]);


	  }



      }

      local_row_index  += FE->var_spaces[i]->dof_per_elm;



    }
  }// end quad


  INT v_per_elm = mesh->v_per_elm;
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

  bool bWeakBC_RHS = BOOL_WEAKLY_IMPOSED_BC;

  if(bWeakBC_RHS){

    //printf("=============================\n");
    //printf("** ELEMENT = %d \n", elm );
    //printf("=============================\n");
    for(jk=mesh->el_f->IA[elm];jk<mesh->el_f->IA[elm+1];jk++){

      //printf("jk = %d, mesh->el_f->IA[element] = %d, mesh->el_f->IA[element+1] = %d  \n",
      //     jk, mesh->el_f->IA[elm], mesh->el_f->IA[elm+1]);
      //  j = face0, face1, face2
      j=jk - mesh->el_f->IA[elm];
      // face is the global number
      face=mesh->el_f->JA[jk];

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

      // Get the BD values
      REAL* val_true_face = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_face_n = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_face_n_neighbor = (REAL *) calloc(nun,sizeof(REAL));
      REAL* val_true_dt_face = (REAL *) calloc(nun,sizeof(REAL));


      // FOR NEIGHBOR..
      INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));
      INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));
      INT *local_dof_on_elm_face_interface = NULL;
      INT *local_dof_on_elm_neighbor = NULL;


      // NOW FOR FACES (at BOUNDARY)
      INT* local_dof_on_elm_face = NULL;
      //Neighbor
      INT neighbor_index[2];
      neighbor_index[0] = -1;
      neighbor_index[1] = -1;
      INT counter = 0;

      iCSRmat *f_el=NULL;
      f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
      icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f);

      // printf("=============================\n");
      //printf("** ELM  = %d   FACE = %d \n",  elm, face );
      //printf("=============================\n");
      int pq,nbr0;
      for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){

	//printf("-- pq = %d, f_el->IA[face] = %d, f_el->IA[face+1] = %d \n", pq, f_el->IA[face], f_el->IA[face+1]);
	nbr0=f_el->JA[pq];
	//printf("-- nbr0 = %d  \n", nbr0);

	neighbor_index[counter] = nbr0;
	counter++;

      }

      //Sanity Check
      /* print out
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
      */

      double fiarea=mesh->f_area[face];
      double lambda = LAME_LAMBDA_GLOBAL;//000000.;//000000.;//000000.; //000000.0;
      double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
      double BC_penalty_term = BC_PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

      //penalty_term*=lambda;

      double penalty_term_pressure =  PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
      double BC_penalty_term_pressure =  BC_PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

      //nq1d_face == 3, 3 quad points.
      zquad_face(cq_face,nq1d_face,dim,xfi,fiarea);

      REAL edge_length = mesh->ed_len[face];

      if(mesh->f_flag[face]>0) {


	for (quad_face=0;quad_face<nquad_face;quad_face++) {

	  qx_face[0] = cq_face->x[quad_face];
	  qx_face[1] = cq_face->y[quad_face];

	  if(dim==3) qx_face[2] = cq_face->z[quad_face];

	  w_face = cq_face->w[quad_face];

	  (*exact_sol2D)(val_true_face,qx_face,time,
			 &(mesh->el_flag[elm]));  // ???


	  (*exact_sol2D)(val_true_face_n,qx_face,time-timestep,
			 &(mesh->el_flag[elm]));  // ???


	  //true solution is in n+1 time
	  (*truesol_dt)(val_true_dt_face,
			qx_face,
			time,
			&(mesh->el_flag[elm]));

	  local_row_index=0;
	  unknown_index=0;
	  local_dof_on_elm_face = dof_on_elm;

	  for(i=0;i<FE->nspaces;i++) {

	    get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,
			  qx_face,
			  v_on_elm,
			  local_dof_on_elm_face,
			  mesh,FE->var_spaces[i]);

	    for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
	      //SLEE
	      if(i==0 || i ==1)
		{
		  bLoc[(local_row_index+test)] += BC_penalty_term * w_face*
		    (val_true_face[unknown_index]*  FE->var_spaces[i]->phi[test]);

		  //DEBUG100
		  //bLoc[(local_row_index+test)] += penalty_term * w_face*
		  //(val_true_face[3]*  FE->var_spaces[i]->phi[test] *finrm[i]);

		}

	      else if(i == 2)
		{
		  //SLEE
		  // Note that new DOF is \phi^3 = [x ; y]
		  //printf("i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
		  bLoc[(local_row_index+test)] +=
		    BC_penalty_term *  w_face * (val_true_face[0] * (qx_face[0] - barycenter->x[0])
						 + val_true_face[1]*  (qx_face[1] - barycenter->y[0]));

		}
	      else if(i == 3)
		{

		  //Stokes 3
		  bLoc[(local_row_index+test)] +=
		    w_face * (val_true_face[0] * finrm[0] + val_true_face[1] * finrm[1])
		    * FE->var_spaces[3]->phi[test];

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

      }// else if


      //printf("FEREE \n");
      icsr_free(f_el);
      //printf("FEREE DONE\n");
    }// for each face


  }


  block_LocaltoGlobal_RHS(dof_on_elm,FE,b,bLoc);

  //if(val_true) free(val_true);
  //if(val_sol) free(val_sol);

  //printf("RHS ASSEMBLE END \n");
  return;
}





void local_assembly_Elasticity(block_dCSRmat* A,dvector *b,REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time, REAL timestep,  INT* switch_on_face)
//void local_assembly_Elasticity(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  //printf("Local Assembly -- start\n");

  bool bEG = BOOL_EG_MECHANICS;
  bool bEG_Pressure = BOOL_EG_PRESSURE;

  REAL alpha = 1.;

  REAL C_0 = 0.001;
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

  double lambda = LAME_LAMBDA_GLOBAL ;

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
    /////////

    // u2
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;

    get_FEM_basis(FE->var_spaces[1]->phi,
		  FE->var_spaces[1]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[1]);

    // u_eg
    local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
    get_FEM_basis(FE->var_spaces[2]->phi,
		  FE->var_spaces[2]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[2]);


    // p
    local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;

    get_FEM_basis(FE->var_spaces[3]->phi,
		  FE->var_spaces[3]->dphi,
		  qx,
		  v_on_elm,
		  local_dof_on_elm,
		  mesh,
		  FE->var_spaces[3]);

    // p_eg
    local_dof_on_elm += FE->var_spaces[3]->dof_per_elm;



    //------#1---------------------------------------------------- START//
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
      //EG PART
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
    //------#1---------------------------------------------------- END//


    //------PORO #0---------------------------------------------------- START//
    //////////////////////////////////////////////////////////////////////////
    // porop0
    // PoroElasticity   -alpha(p, div v)
    // Add Pressure to Mechanics
    //////////////////////////////////////////////////////////////////////////
    //   -alpha(p, div v)
    //          CG  CG
    //          trial (column) -3  test (row) -0
    // test - u0 cg
    local_row_index = 0.;
    // trial - p
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[trial] *  (FE->var_spaces[0]->dphi[test*dim]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // test - u1 cg
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    // trial - p
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[trial] *  (FE->var_spaces[1]->dphi[test*dim+1]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    if(bEG){
      //   -alpha(p, div v)
      //          CG  EG
      //        trial (column) 3  test (row) 2
      // test - u eg
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // trial - p
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

	  kij = -alpha * FE->var_spaces[3]->phi[trial] * 2.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }
    }


    //------PORO #0---------------------------------------------------- END//


    //------STOKES  #0---------------------------------------------------- START//
    //////////////////////////////////////////////////////////////////////////
    // Transpose of Stokes #0
    // PoroElasticity   -(div v, p)
    // Add Pressure to Mechanics
    //////////////////////////////////////////////////////////////////////////
    //   -alpha(div v, p)
    //          CG  CG
    //          trial (column) -3  test (row) -0
    // test - p cg
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // trial - u0 cg
    local_col_index = 0.;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[test] *  (FE->var_spaces[0]->dphi[trial*dim]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // test - p cg
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
    // trial - u1 cg
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){

	kij = -alpha * FE->var_spaces[3]->phi[test] *  (FE->var_spaces[1]->dphi[trial*dim+1]);

	ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    if(bEG){
      //   -alpha(div v, p)
      //          EG  CG
      //        test (column) 3  trial (row) 2
      // test - p cg
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // trial - u eg
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
	// Loop over Trial Functions (Columns)
	for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){

	  kij = -alpha * FE->var_spaces[3]->phi[test] * 2.;

	  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
	}
      }

    }

    //------STOKES #0---------------------------------------------------- END//



    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // PRESSURE BLOCK
    // p-p block:
    /*
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){

      kij =  FE->var_spaces[3]->phi[test] *  FE->var_spaces[3]->phi[trial];

      ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      }
    */



  }//QUAD


  block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);

  //printf("*********** {DONE} @ END  ************************ \n");


  //printf("Local Assembly --end\n");

  return;
}
/*********************************************************************/
