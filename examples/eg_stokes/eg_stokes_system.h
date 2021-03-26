/*********************************************************************************/
iCSRmat *extend_el_dof(fespace *FE,mesh_struct *mesh, INT *ix) {
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
  if(f_el) free(f_el);
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

  // Frees
  icsr_free(el2el);
  if(el2el) free(el2el);
  return el_dof;
}
/*********************************************************************************/
void sl_create_CSR_cols_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2, mesh_struct *mesh) {
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
        // We haven't been here
        if (ix[mydof]!=i) {
          A->JA[icp] = mydof;
          icp++;
          ix[mydof] = i;
        }
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el_2);
  if(el_dof) {
    free(el_dof); // no need of icsr_free here;
    el_dof=NULL;
  }
  return;
}
/********************************************************************************************/
/********************************************************************************************/
/****************************************************************************/
void sl_create_CSR_rows_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2,mesh_struct *mesh) {
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
        // We haven't been here
        if (ix[mydof]!=i) {
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
  if(el_dof) {
    free(el_dof); // no need of icsr_free here;
    el_dof=NULL;
  }

  return;
}
/******/
/******/
void block_LocaltoGlobal_RHS(INT *dof_on_elm,block_fespace *FE,dvector *b,REAL *bLoc) {
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

void block_LocaltoGlobal_neighbor(INT *dof_on_elm,INT *dof_on_elm_neighbor,block_fespace *FE,block_dCSRmat *A,REAL *ALoc) {
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


  REAL* tmperrL2 = (REAL *) calloc(FE->nspaces,sizeof(REAL));

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
    if(f_el) free(f_el);
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
if(switch_on_face) free(switch_on_face);
if(tmperrL2) free(tmperrL2);

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
    INT i,j,k,idim;
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
    REAL kij_unv = 0.0;
    REAL kij_uvn = 0.0;
    REAL kij_unvn = 0.0;

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
      if(dim==3) barycenter_neighbor->z[0] = mesh->el_mid[neighbor_index[1]*dim+2];
    }

    // Sanity Check
    //
    //   for(INT zz=0; zz<v_per_elm; ++zz)
    //   {
    //   printf("V_on_elm[%d] = %d  ||  dof_on_elm[%d] =%d \n", zz,v_on_elm[zz],zz,dof_on_elm[zz]);
    // }
    //
    // if(counter ==2) // if it has neighbor
    // for(INT zz=0; zz<v_per_elm; ++zz)
    // {
    // printf("V_on_elm_neighbor[%d] = %d  ||  dof_on_elm_neighbor[%d] =%d \n", zz,v_on_elm_neighbor[zz],zz,dof_on_elm_neighbor[zz]);
    // }
    //
    // printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[0], barycenter->x[0], barycenter->y[0]);
    // if(counter ==2) // if it has neighbor
    // printf("FACE = %d, ELM = %d, (X,Y) = (%f,%f) \n", face, neighbor_index[1], barycenter_neighbor->x[0], barycenter_neighbor->y[0]);


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
    if(dim==3) finrm[2]=mesh->f_norm[face*dim+2];

    //Sanity Check; Print the normal vectors..
    //printf("FACE = %d, ELM = %d, (Nx,Ny) = (%f,%f) \n", face, neighbor_index[0], finrm[0], finrm[1]);

    INT nq1d_face=3;
    qcoordinates *cq_face = allocateqcoords_bdry(nq1d_face,1,dim,2);

    INT maxdim=3;
    REAL qx_face[maxdim];
    REAL w_face;
    INT nquad_face= cq_face->nq_per_elm; //quad<cq_face->nq_per_elm;

    // Make life easier for eg TERMS
    REAL eg_xterm, eg_yterm, eg_zterm, eg_xterm_neighbor, eg_yterm_neighbor, eg_zterm_neighbor;


    // Get xfi points..
    for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){
      //printf("M** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);
      j=jkl-mesh->f_v->IA[face];
      k=mesh->f_v->JA[jkl];

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

    // get the penalty (FACE ASSEMBLE)
    double penalty_term = PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
    double BC_penalty_term = BC_PENALTY_PARAMETER_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));

    double penalty_term_pressure = PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));
    double BC_penalty_term_pressure = BC_PENALTY_PARAMETER_PRESSURE_GLOBAL / (pow(fiarea,1e0/(REAL )(dim-1)));


    // DEBUG mode
    // Get the BD values
    //REAL* val_true_face = (REAL *) calloc(nquad_face,sizeof(REAL));
    // Easier to keep track of test and trial FUNCTIONS
    INT iu0 = 0;
    INT iu1 = 1;
    INT iu2 = 2;
    INT ieg = dim;
    INT ip = dim+1;
    INT u0dofpelm = FE->var_spaces[iu0]->dof_per_elm;
    INT u0dofpface = FE->var_spaces[iu0]->dof_per_face;
    INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
    INT u1dofpface = FE->var_spaces[iu1]->dof_per_face;
    INT u2dofpelm,u2dofpface;
    if(dim==3) {
      u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
      u2dofpface = FE->var_spaces[iu2]->dof_per_face;
    }
    INT egdofpelm = FE->var_spaces[ieg]->dof_per_elm;
    INT egdofpface = FE->var_spaces[ieg]->dof_per_face;
    INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;
    INT pdofpface = FE->var_spaces[ip]->dof_per_face;
    REAL u0,u1,u2,p,v0,v1,v2,q,u0x,u0y,u0z,u1x,u1y,u1z,u2x,u2y,u2z,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z;
    REAL un0,un1,un2,pn,vn0,vn1,vn2,qn,un0x,un0y,un0z,un1x,un1y,un1z,un2x,un2y,un2z,vn0x,vn0y,vn0z,vn1x,vn1y,vn1z,vn2x,vn2y,vn2z;

    REAL* neighbor_basis_u0_phi  = (REAL *) calloc(u0dofpelm,sizeof(REAL));
    REAL* neighbor_basis_u0_dphi = (REAL *) calloc(u0dofpelm*dim,sizeof(REAL));

    REAL* neighbor_basis_u1_phi  = (REAL *) calloc(u1dofpelm,sizeof(REAL));
    REAL* neighbor_basis_u1_dphi = (REAL *) calloc(u1dofpelm*dim,sizeof(REAL));

    REAL* neighbor_basis_u2_phi;
    REAL* neighbor_basis_u2_dphi;
    if(dim==3) {
      neighbor_basis_u2_phi = (REAL *) calloc(u2dofpelm,sizeof(REAL));
      neighbor_basis_u2_dphi = (REAL *) calloc(u2dofpelm*dim,sizeof(REAL));
    }

    REAL* neighbor_basis_p_phi  = (REAL *) calloc(pdofpelm,sizeof(REAL));
    REAL* neighbor_basis_p_dphi = (REAL *) calloc(pdofpelm*dim,sizeof(REAL));

    for (quad_face=0;quad_face<nquad_face;quad_face++) {

      qx_face[0] = cq_face->x[quad_face];
      qx_face[1] = cq_face->y[quad_face];
      if(dim==3) qx_face[2] = cq_face->z[quad_face];

      w_face = cq_face->w[quad_face];

      // EG Stuff
      eg_xterm = qx_face[0] - barycenter->x[0];
      eg_yterm = qx_face[1] - barycenter->y[0];
      eg_xterm_neighbor = qx_face[0] - barycenter_neighbor->x[0];
      eg_yterm_neighbor = qx_face[1] - barycenter_neighbor->y[0];
      if(dim==3) {
        eg_zterm = qx_face[2] - barycenter->z[0];
        eg_zterm_neighbor = qx_face[2] - barycenter_neighbor->z[0];
      }

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

          // get basis functions
          local_dof_on_elm_face = dof_on_elm;
          // u,v
          for(idim=0;idim<dim;idim++) {
            get_FEM_basis(FE->var_spaces[idim]->phi,FE->var_spaces[idim]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[idim]);
            local_dof_on_elm_face += FE->var_spaces[idim]->dof_per_elm;
          }
          // eg skip
          local_dof_on_elm_face += egdofpelm;
          // p
          get_FEM_basis(FE->var_spaces[ip]->phi,FE->var_spaces[ip]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[ip]);
          local_dof_on_elm_face += pdofpelm;

          //---- #2 BOUNDARY ---------------------------- START //
          // TODO: Is this theta=0 only? IIPG?
          local_row_index = 0;
          // v0 test rows
          for (test=0; test<u0dofpelm;test++){
            v0 = FE->var_spaces[iu0]->phi[test];
            v0x = FE->var_spaces[iu0]->dphi[test*dim];
            v0y = FE->var_spaces[iu0]->dphi[test*dim+1];
            if(dim==3) v0z = FE->var_spaces[iu0]->dphi[test*dim+2];

            // u0-v0 block: (-2*u0x*n0 - u0y*n1 - u0z*n2)*v0 + theta*(2*v0x*n0 + v0y*n1 + v0z*n2)*u0 + pen*u0*v0
            local_col_index = 0;
            for (trial=0; trial<u0dofpelm;trial++){
              u0 = FE->var_spaces[iu0]->phi[trial];
              u0x = FE->var_spaces[iu0]->dphi[trial*dim];
              u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
              if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

              kij = -2.0*u0x*finrm[0]*v0 - u0y*finrm[1]*v0;
              if(dim==3) kij -= u0z*finrm[2]*v0;
              //  kij = -2.0 * FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * FE->var_spaces[0]->phi[test];
              //  kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] *    FE->var_spaces[0]->phi[test];
              // TODO: Below is 0 for Stokes right?  So I'm skipping this term
              //kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[0] * FE->var_spaces[0]->phi[test];
              //penalty term
              kij += BC_penalty_term*u0*v0;
              //  kij += BC_penalty_term * FE->var_spaces[0]->phi[trial] * FE->var_spaces[0]->phi[test];
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face *kij;
            } // End u0 col loop
            local_col_index += u0dofpelm;

            // u1-v0 block : -u1x*n1*v0 + th*v0y*n0*u1
            for (trial=0; trial<u1dofpelm;trial++){
              u1 = FE->var_spaces[iu1]->phi[trial];
              u1x = FE->var_spaces[iu1]->dphi[trial*dim];
              u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
              if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

              kij = -u1x*finrm[1]*v0;
              // kij = - FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * FE->var_spaces[0]->phi[test];
              // kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[0]->phi[test];

              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // End u1 loop
            local_col_index+= u1dofpelm;

            if(dim==3) {
              // u2-v0 block: -u2x*n2*v0 + th*(v0z*n0*u2)
              for (trial=0; trial<u2dofpelm;trial++){
                u2 = FE->var_spaces[iu2]->phi[trial];
                u2x = FE->var_spaces[iu2]->dphi[trial*dim];
                u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
                u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];
                kij = -u2x*finrm[2]*v0;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u2 loop
              local_col_index+= u2dofpelm;
            }

            // eg-v0 block: -2*(n0*v0) + th*(messy) + pen*(x-xT)*v0
            if(bEG) {
              for (trial=0; trial<egdofpelm;trial++){
                kij = -2.0* 1. * finrm[0] * v0;

                // TODO: Gone?
                // // Divergence
                // // u1-q block : <dx(u1),dx(v1)>
                // kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];//FE->var_spaces[dim]->dphi[test*dim];
                // kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];

                //penalty term
                kij += penalty_term * v0 * eg_xterm;//(qx_face[0]-barycenter->x[0]);
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // eg col loop
            } // end if bEG
            local_col_index += egdofpelm;

            // p-v0 bloc: p*v0*n0
            for (trial=0; trial<pdofpelm;trial++){
              p = FE->var_spaces[ip]->phi[trial];
              kij = alpha * p * v0 * finrm[0];
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // end p col
          } // End v0 loop
          local_row_index = u0dofpelm;

          // v1 test rows
          for (test=0; test<u1dofpelm;test++){
            v1 = FE->var_spaces[iu1]->phi[test];
            v1x = FE->var_spaces[iu1]->dphi[test*dim];
            v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
            if(dim==3) v1z = FE->var_spaces[iu1]->dphi[test*dim+2];

            local_col_index = 0;
            // u0-v1 block: -u0y*n0*v1 + th*v1x*n1*u0
            for (trial=0; trial<u0dofpelm;trial++){
              u0 = FE->var_spaces[iu0]->phi[trial];
              u0x = FE->var_spaces[iu0]->dphi[trial*dim];
              u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
              if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

              kij = -u0y*finrm[0]*v1;
              // kij = -FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0] * FE->var_spaces[1]->phi[test];
              // kij -= lambda * FE->var_spaces[0]->dphi[trial*dim] * finrm[1] * FE->var_spaces[1]->phi[test];
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // End u0 Loop
            local_col_index += u0dofpelm;

            // u1-v1: -u1x*n0*v1 - 2*u1y*n1*v1 - u1z*n2*v1 + th*(v1x*n0*u1 + 2*v1y*n1*u1 + v1z*n2*u1) + pen*u1*v1
            for (trial=0; trial<u1dofpelm;trial++){
              u1 = FE->var_spaces[iu1]->phi[trial];
              u1x = FE->var_spaces[iu1]->dphi[trial*dim];
              u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
              if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

              kij = -2.0*u1y*finrm[1]*v1 - u1x*finrm[0]*v1;
              if(dim==3) kij -= u1z*finrm[2]*v1;
              // kij = -2.0 * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test] ;
              // kij -= 1.0 * FE->var_spaces[1]->dphi[trial*dim] * finrm[0]   * FE->var_spaces[1]->phi[test];
              // // NEW SLEE : elasticity div part
              // // u2-v2 block : <dy(u2),dx(v2)>
              // kij -= lambda * FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * FE->var_spaces[1]->phi[test];
              //penalty term
              kij += BC_penalty_term*u1*v1;
              //  kij += BC_penalty_term *  FE->var_spaces[1]->phi[trial]  * FE->var_spaces[1]->phi[test];
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // End u1 Loop
            local_col_index+= u1dofpelm;

            if(dim==3) {
              // u2-v1 block: -u2y*n2*v1 + th*(v1z*n1*u2)
              for (trial=0; trial<u2dofpelm;trial++){
                u2 = FE->var_spaces[iu2]->phi[trial];
                u2x = FE->var_spaces[iu2]->dphi[trial*dim];
                u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
                u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];
                kij = -u2y*finrm[2]*v1;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u2 loop
              local_col_index+= u2dofpelm;
            }

            // eg-v1 block: -2*(n1*v1) + th*(messy) + pen*(y-yT)*v1
            if(bEG) {
              for (trial=0; trial<egdofpelm;trial++){
                kij = -2.0* 1. * finrm[1] * v1;

                // TODO: Gone?
                // // NEW SLEE : elasticity div part
                // // u2-v1 block : <dy(u2),dx(v1)>
                // kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];
                // kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];

                //penalty term
                kij += penalty_term * v1 * eg_yterm;//(qx_face[1]-barycenter->y[0]);
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // eg col loop
            } // end if bEG
            local_col_index += egdofpelm;

            // p-v1 bloc: p*v1*n1
            for (trial=0; trial<pdofpelm;trial++){
              p = FE->var_spaces[ip]->phi[trial];
              kij = alpha * p * v1 * finrm[1];
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // end p col

          } // End v1 Loop
          local_row_index += u1dofpelm;

          // v2 test rows
          if(dim==3) {
            for (test=0; test<u2dofpelm;test++){
              v2 = FE->var_spaces[iu2]->phi[test];
              v2x = FE->var_spaces[iu2]->dphi[test*dim];
              v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
              v2z = FE->var_spaces[iu2]->dphi[test*dim+2];

              local_col_index = 0;
              // u0-v2 block: -u0z*n0*v2 + th*v2x*n2*u0
              for (trial=0; trial<u0dofpelm;trial++){
                u0 = FE->var_spaces[iu0]->phi[trial];
                u0x = FE->var_spaces[iu0]->dphi[trial*dim];
                u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
                u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];
                kij = -u0z*finrm[0]*v2;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u0 Loop
              local_col_index += u0dofpelm;

              // u1-v2: -u1z*n1*v2 + th*(v2y*n2*u1)
              for (trial=0; trial<u1dofpelm;trial++){
                u1 = FE->var_spaces[iu1]->phi[trial];
                u1x = FE->var_spaces[iu1]->dphi[trial*dim];
                u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
                u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
                kij = -u1z*finrm[1]*v2;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u1 Loop
              local_col_index+= u1dofpelm;

              // u2-v2 block: -u2x*n0*v2 - u2y*n1*v2 - 2*u2z*n2*v2 + th*(v2x*n0*u2 + v2y*n1*u2 + 2*v2z*n2*u2) + pen*u2*v2
              for (trial=0; trial<u2dofpelm;trial++){
                u2 = FE->var_spaces[iu2]->phi[trial];
                u2x = FE->var_spaces[iu2]->dphi[trial*dim];
                u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
                u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];
                kij = -u2x*finrm[0]*v2 - u2y*finrm[1]*v2 - 2*u2z*finrm[2]*v2;
                // penalty term
                kij += BC_penalty_term*u2*v2;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u2 loop
              local_col_index+= u2dofpelm;

              // eg-v2 block: -2*(n2*v2) + th*(messy) + pen*(z-zT)*v2
              if(bEG) {
                for (trial=0; trial<egdofpelm;trial++){
                  kij = -2.0* 1. * finrm[2] * v2;
                  //penalty term
                  kij += penalty_term * v2* eg_zterm;//(qx_face[2]-barycenter->z[0]);
                  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
                } // eg col loop
              } // end if bEG
              local_col_index += egdofpelm;

              // p-v2 bloc: p*v2*n2
              for (trial=0; trial<pdofpelm;trial++){
                p = FE->var_spaces[ip]->phi[trial];
                kij = alpha * p * v2 * finrm[2];
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // end p col
            } // End v2 Loop
            local_row_index += u2dofpelm;
          } // end if dim==3

          //eg block rows
          for (test=0; test<egdofpelm;test++){

            //EG PART
            //NEW SLEE :
            local_col_index = 0;
            if(bEG) {
              // u0-eg block: -( (2*u0x*n0 + u0y*n1 + u0z*n2)*(x-xT) + (u0y*n0)*(y-yT) + (u0z*n0)*(z-zT)) + th*u0*n0 + pen*u0*(x-xT)
              for (trial=0; trial<u0dofpelm;trial++){
                u0 = FE->var_spaces[iu0]->phi[trial];
                u0x = FE->var_spaces[iu0]->dphi[trial*dim];
                u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
                if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

                kij = (-2.0*u0x*finrm[0] - u0y*finrm[1])*eg_xterm;
                kij -= u0y*finrm[0]*eg_yterm;
                if(dim==3) kij += -u0z*finrm[2]*eg_xterm - u0z*finrm[0]*eg_zterm;

                // TODO: Gone?
                // // Divergence
                // // u1-q block : <dx(u1),dx(v1)>
                // kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
                // kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);

                //penalty term
                kij += penalty_term * u0 * eg_xterm;//(qx_face[0]- barycenter->x[0]);
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // End u0 col
              local_col_index += u0dofpelm;

              // u1-eg block: -( (u1x*n1)*(x-xT) + (u1x*n0+2*u1y*n1+u1z*n2)*(y-yT) + (u1z*n1)*(z-zT)) + th*u1*n1 + pen*u1*(y-yT)
              for (trial=0; trial<u1dofpelm;trial++){
                u1 = FE->var_spaces[iu1]->phi[trial];
                u1x = FE->var_spaces[iu1]->dphi[trial*dim];
                u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
                if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

                kij = -u1x*finrm[1]*eg_xterm;
                kij += (-u1x*finrm[0] - 2.0*u1y*finrm[1])*eg_yterm;
                if(dim==3) kij += -u1z*finrm[2]*eg_yterm - u1z*finrm[1]*eg_zterm;
                // kij = -2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
                // kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
                //
                // kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
                //
                // kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);
                // kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);

                //penalty term
                kij += penalty_term * u1*eg_yterm;//FE->var_spaces[1]->phi[trial] * (qx_face[1]- barycenter->y[0]);
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // end u1 col
              local_col_index += u1dofpelm;

              if(dim==3) {
                // u2-eg block: -( (u2x*n2)*(x-xT) + (u2y*n2)*(y-yT) + (u2x*n0+u2y*n1+2*u2z*n2)*(z-zT)) + th*u2*n2 + pen*u2*(z-zT)
                for (trial=0; trial<u2dofpelm;trial++){
                  u2 = FE->var_spaces[iu2]->phi[trial];
                  u2x = FE->var_spaces[iu2]->dphi[trial*dim];
                  u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
                  u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

                  kij = -u2x*finrm[2]*eg_xterm;
                  kij += -u2y*finrm[2]*eg_yterm;
                  kij += (-u2x*finrm[0]-u2y*finrm[1]-2.0*u2z*finrm[2])*eg_zterm;

                  //penalty term
                  kij += penalty_term * u2*eg_zterm;
                  ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
                } // end u2 col
                local_col_index += u2dofpelm;
              } //end if dim==3
            } // end if bEG

            // TODO: ISN"T THIS EG?? -> why isn't it under if(bEG)??
            // eg-eg: (2*n0(x-xT) + 2*n1*(y-yT) + 2*n2*(z-zT))*(-1+th) + pen*( (x-xT)^2 + (y-yT)^2 + (z-zT)^2)
            local_col_index = u0dofpelm+u1dofpelm;
            if(dim==3) local_col_index += u2dofpelm;
            for (trial=0; trial<egdofpelm;trial++){

              kij = -2.0 * 1. * finrm[0] * eg_xterm;//(qx_face[0] - barycenter->x[0]);
              kij -= 2.0 * 1. * finrm[1] * eg_yterm; //(qx_face[1] - barycenter->y[0]);
              if(dim==3) kij -= 2.0 * 1. * finrm[2] * eg_zterm; //(qx_face[2] - barycenter->z[0]);

              // TODO: This is gone right?
              // // NEW SLEE : elasticity div part
              // // u2-v1 block : <dy(u2),dx(v1)>
              // kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
              // kij -= lambda * 1. *  finrm[0] * (qx_face[0]- barycenter->x[0]);
              //
              // kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
              // kij -= lambda * 1. *  finrm[1] * (qx_face[1]- barycenter->y[0]);
              //

              //penalty term
              kij += BC_penalty_term * eg_xterm*eg_xterm;//(qx_face[0]- barycenter->x[0]) * (qx_face[0]- barycenter->x[0]);
              kij += BC_penalty_term * eg_yterm*eg_yterm;//(qx_face[1]- barycenter->y[0]) * (qx_face[1]- barycenter->y[0]);
              if(dim==3) kij += BC_penalty_term * eg_zterm*eg_zterm;//(qx_face[2]- barycenter->z[0]) * (qx_face[2]- barycenter->z[0]);
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // end eg col loop
            local_col_index += egdofpelm;

            if(bEG) {
              // p-eg block: p( (x-xT)*n0 + (y-yT)*n1 + (z-zT)*n2 )
              for (trial=0; trial<pdofpelm;trial++){
                p = FE->var_spaces[ip]->phi[trial];
                kij = alpha * p * (eg_xterm*finrm[0] + eg_yterm*finrm[1]);
                if(dim==3) kij += alpha * p * eg_zterm*finrm[2];
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // end p col
            } // if bEG
          } // end eg row loop
          local_row_index += egdofpelm;

          // TODO: I am moving this up inside the other loops but will add the if statements
          // Face Boundary
          // if(bEG){
          // //EG PART
          // //NEW SLEE :
          // // u0- q block
          // local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
          // if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
          // local_col_index = 0;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
          //   // Loop over Trial Functions (Columns)
          //   for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
          //
          //     //u0-q0
          //     //u0-q1
          //     kij = -2.0*FE->var_spaces[0]->dphi[trial*dim] * finrm[0]  * (qx_face[0]- barycenter->x[0]); //FE->var_spaces[dim]->phi[test];
          //     kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[1] * (qx_face[0]- barycenter->x[0]);
          //
          //     kij -= FE->var_spaces[0]->dphi[trial*dim+1] * finrm[0]  * (qx_face[1]- barycenter->y[0]); //FE->var_spaces[dim]->phi[test];
          //
          //     // Divergence
          //     // u1-q block : <dx(u1),dx(v1)>
          //     kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[0] * (qx_face[0]- barycenter->x[0]);//FE->var_spaces[dim]->dphi[test*dim];
          //     kij -= lambda*FE->var_spaces[0]->dphi[trial*dim]* finrm[1] * (qx_face[1]- barycenter->y[0]);
          //
          //     //penalty term
          //     kij += penalty_term * FE->var_spaces[0]->phi[trial] * (qx_face[0]- barycenter->x[0]);
          //
          //     ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //   }
          // }
          //
          // // u1-q block
          // local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
          // if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
          // local_col_index = FE->var_spaces[0]->dof_per_elm;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
          //   // Loop over Trial Functions (Columns)
          //   for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
          //
          //     //u1-q0
          //     //u1-q1
          //     kij = -2.0*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] * (qx_face[1]- barycenter->y[0]);
          //     kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[0] * (qx_face[1]- barycenter->y[0]);
          //
          //     kij -= FE->var_spaces[1]->dphi[trial*dim] * finrm[1] * (qx_face[0]- barycenter->x[0]);
          //
          //     kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[0] *  (qx_face[0]- barycenter->x[0]);
          //     kij -= lambda*FE->var_spaces[1]->dphi[trial*dim+1] * finrm[1] *  (qx_face[1]- barycenter->y[0]);
          //
          //
          //     //penalty term
          //     kij += penalty_term * FE->var_spaces[1]->phi[trial] * (qx_face[1]- barycenter->y[0]);
          //
          //     ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //   }
          // }

          // q-v0 block:
          // local_row_index = 0;
          // local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
          // if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
          // Loop over Trial Functions (Columns)
          // for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
          //
          //   //q0 - v0
          //   //q1 - v0
          //   kij = -2. * 1. * finrm[0] * FE->var_spaces[0]->phi[test];
          //
          //   // Divergence
          //   // u1-q block : <dx(u1),dx(v1)>
          //   kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];//FE->var_spaces[dim]->dphi[test*dim];
          //   kij -= lambda* 1. * finrm[0] * FE->var_spaces[0]->phi[test];
          //
          //
          //   //penalty term
          //   kij += penalty_term * FE->var_spaces[0]->phi[test] * (qx_face[0]- barycenter->x[0]);
          //
          //   ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          // }
          // }


          // // q-v1 block:
          // local_row_index = FE->var_spaces[0]->dof_per_elm;
          // local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
          // if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
          //   // Loop over Trial Functions (Columns)
          //   for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
          //
          //     //q0-v1
          //     //q1-v1
          //     kij = -2.0* 1. * finrm[1] * FE->var_spaces[1]->phi[test];
          //
          //     // NEW SLEE : elasticity div part
          //     // u2-v1 block : <dy(u2),dx(v1)>
          //     kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];
          //     kij -= lambda * 1. * finrm[1] * FE->var_spaces[1]->phi[test];
          //
          //     //penalty term
          //     kij += penalty_term * FE->var_spaces[1]->phi[test] * (qx_face[1]- barycenter->y[0]);
          //
          //     ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //   }
          // }

          // } // bEG
          //----- #2 BOUNDARY ------------------------------------------- END //

          // TODO: Moving this up as well (I believe this is the pressure terms?)
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
          // local_row_index = 0;
          // // trial - p cg
          // local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
          //   // Loop over Trial Functions (Columns)
          //   for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
          //
          //     kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[0]->phi[test]*finrm[0];
          //
          //     ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //   }
          // }
          // //      ({p},  [u] dot n)
          // //       CG  CG u1
          // //       trial (column) 3  test (row) 1
          // // test - u0 cg
          // local_row_index = FE->var_spaces[0]->dof_per_elm;
          // // trial - p cg
          // local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
          // // Loop over Test Functions (Rows)
          // for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
          //   // Loop over Trial Functions (Columns)
          //   for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
          //
          //     kij = alpha * FE->var_spaces[3]->phi[trial] *  FE->var_spaces[1]->phi[test]*finrm[1];
          //
          //     ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //   }
          // }
          //
          // if(bEG){
          //
          //   //      ({p},  [u] dot n)
          //   //       CG  EG
          //   //       trial (column) 3  test (row) 2
          //   // test - u0 eg
          //   local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
          //   // trial - p cg
          //   local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
          //   // Loop over Test Functions (Rows)
          //   for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
          //     // Loop over Trial Functions (Columns)
          //     for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
          //
          //       kij = alpha * FE->var_spaces[3]->phi[trial] *
          //       (finrm[0] *  (qx_face[0]- barycenter->x[0]) + finrm[1] *  (qx_face[1]- barycenter->y[0]));
          //
          //       ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          //     }
          //   }
          //
          // }

          //---- Poro P2 Boundary -------------------------------------- END //


          //---- Stokes P2 Boundary -------------------------------------- START //
          ///////////////////////////////////////////////////////////////////////////
          // On BOUNDARY
          // STokes  + ({p},  [u] dot n)
          //////////////////////////////////////////////////////////////////////////
          //      ({p},  [u] dot n)
          //       CG  CG u0
          //       trial (column) 3  test (row) 0neighbor_basis_1_
          // test - p cg

          // q block rows: reset counters to be safe
          local_row_index = u0dofpelm + u1dofpelm + egdofpelm;
          if(dim==3) local_row_index += u2dofpelm;
          for (test=0; test<pdofpelm;test++){
            q = FE->var_spaces[ip]->phi[test];

            local_col_index = 0;
            // u0-q: u0*n0*q
            for (trial=0; trial<u0dofpelm;trial++){
              u0 = FE->var_spaces[iu0]->phi[trial];
              kij = alpha * u0*finrm[0]*q;
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // u0 cols
            local_col_index += u0dofpelm;

            // u1-q: u1*n1*q
            for (trial=0; trial<u1dofpelm;trial++){
              u1 = FE->var_spaces[iu1]->phi[trial];
              kij = alpha * u1*finrm[1]*q;
              ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // u1 cols
            local_col_index += u1dofpelm;

            // u2-q: u2*n2*q
            if(dim==3) {
              for (trial=0; trial<u2dofpelm;trial++){
                u2 = FE->var_spaces[iu2]->phi[trial];
                kij = alpha * u2*finrm[2]*q;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // u2 cols
              local_col_index += u2dofpelm;
            } // dim==3

            // eg-q: ((x-xT)*n0 + (y-yT)*n1 + (z-zT)*n2)*q
            if(bEG){
              for (trial=0; trial<egdofpelm;trial++){
                kij = alpha * (eg_xterm*finrm[0] + eg_yterm*finrm[1]) * q;
                if(dim==3) kij += alpha * eg_zterm*finrm[2] * q;
                ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // eg cols
            } // if bEG
          } // q rows
        }// Weakly Face
      } else if(counter ==2) { //on interface
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
        // get basis functions
        local_dof_on_elm_face_interface = dof_on_elm;
        // u,v
        for(idim=0;idim<dim;idim++) {
          get_FEM_basis(FE->var_spaces[idim]->phi,FE->var_spaces[idim]->dphi,qx_face,v_on_elm,local_dof_on_elm_face_interface,mesh,FE->var_spaces[idim]);
          local_dof_on_elm_face_interface += FE->var_spaces[idim]->dof_per_elm;
        }
        // eg skip
        local_dof_on_elm_face_interface += egdofpelm;
        // p
        get_FEM_basis(FE->var_spaces[ip]->phi,FE->var_spaces[ip]->dphi,qx_face,v_on_elm,local_dof_on_elm_face_interface,mesh,FE->var_spaces[ip]);
        local_dof_on_elm_face_interface += pdofpelm;

        // Now, get the basis for the neighbor
        // u0
        get_FEM_basis(neighbor_basis_u0_phi,neighbor_basis_u0_dphi,qx_face,v_on_elm_neighbor,dof_on_elm_neighbor,mesh,FE->var_spaces[iu0]);
        local_dof_on_elm_neighbor = dof_on_elm_neighbor + u0dofpelm;

        // u1
        get_FEM_basis(neighbor_basis_u1_phi,neighbor_basis_u1_dphi,qx_face,v_on_elm_neighbor,local_dof_on_elm_neighbor,mesh,FE->var_spaces[iu1]);
        local_dof_on_elm_neighbor += u1dofpelm;

        // u2
        if(dim==3) {
          get_FEM_basis(neighbor_basis_u2_phi,neighbor_basis_u2_dphi,qx_face,v_on_elm_neighbor,local_dof_on_elm_neighbor,mesh,FE->var_spaces[iu2]);
          local_dof_on_elm_neighbor += u2dofpelm;
        }

        // u eg
        local_dof_on_elm_neighbor += egdofpelm;

        // p
        get_FEM_basis(neighbor_basis_p_phi,neighbor_basis_p_dphi,qx_face,v_on_elm_neighbor,local_dof_on_elm_neighbor,mesh,FE->var_spaces[ip]);
        local_dof_on_elm_neighbor += pdofpelm;

        // TODO: reordering here and elminating any poro stuff
        local_row_index = 0;
        // v0 test rows
        for (test=0; test<u0dofpelm;test++){
          v0 = FE->var_spaces[iu0]->phi[test];
          v0x = FE->var_spaces[iu0]->dphi[test*dim];
          v0y = FE->var_spaces[iu0]->dphi[test*dim+1];
          vn0 = neighbor_basis_u0_phi[test];
          vn0x = neighbor_basis_u0_dphi[test*dim];
          vn0y = neighbor_basis_u0_dphi[test*dim+1];
          if(dim==3) {
            v0z = FE->var_spaces[iu0]->dphi[test*dim+2];
            vn0z = neighbor_basis_u0_dphi[test*dim+2];
          }

          local_col_index = 0;
          // u0-v0 block: 0
          local_col_index += u0dofpelm;

          // u1-v0 block : 0
          local_col_index+= u1dofpelm;

          if(dim==3) {
            // u2-v0 block: 0
            local_col_index+= u2dofpelm;
          }

          // eg-v0 block: 0
          local_col_index += egdofpelm;

          // p-v0 bloc: // TODO: why doesn't this vanish?
          for (trial=0; trial<pdofpelm;trial++){
            p = FE->var_spaces[ip]->phi[trial];
            pn = neighbor_basis_p_phi[trial];
            // p-v0
            kij = alpha * 0.5 * p * v0 * finrm[0];
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-v0
            kij = alpha * 0.5 * pn * v0*finrm[0];
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // p-v0neighbor
            kij = -alpha * 0.5 * p * vn0 * finrm[0];
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-v0neighbor
            kij = -alpha * 0.5 * pn *  vn0*finrm[0];
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // end p col
        } // End v0 loop
        local_row_index = u0dofpelm;

        // v1 test rows
        for (test=0; test<u1dofpelm;test++){
          v1 = FE->var_spaces[iu1]->phi[test];
          v1x = FE->var_spaces[iu1]->dphi[test*dim];
          v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
          vn1 = neighbor_basis_u1_phi[test];
          vn1x = neighbor_basis_u1_dphi[test*dim];
          vn1y = neighbor_basis_u1_dphi[test*dim+1];
          if(dim==3) {
            v1z = FE->var_spaces[iu1]->dphi[test*dim+2];
            vn1z = neighbor_basis_u1_dphi[test*dim+2];
          }

          local_col_index = 0;
          // u0-v1 block: 0
          local_col_index += u0dofpelm;

          // u1-v1 block : 0
          local_col_index+= u1dofpelm;

          if(dim==3) {
            // u2-v1 block: 0
            local_col_index+= u2dofpelm;
          }

          // eg-v1 block:0
          local_col_index += egdofpelm;

          // p-v1 bloc: // TODO: why doesn't this vanish?
          for (trial=0; trial<pdofpelm;trial++){
            p = FE->var_spaces[ip]->phi[trial];
            pn = neighbor_basis_p_phi[trial];
            // p-v1
            kij = alpha * 0.5 * p * v1 * finrm[1];
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-v0
            kij = alpha * 0.5 * pn * v1*finrm[1];
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // p-v0neighbor
            kij = -alpha * 0.5 * p * vn1 * finrm[1];
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-v0neighbor
            kij = -alpha * 0.5 * pn *  vn1*finrm[1];
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // end p col
        } // End v1 loop
        local_row_index += u1dofpelm;

        // v2 test rows
        if(dim==3) {
          for (test=0; test<u2dofpelm;test++){
            v2 = FE->var_spaces[iu2]->phi[test];
            v2x = FE->var_spaces[iu2]->dphi[test*dim];
            v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
            v2z = FE->var_spaces[iu2]->dphi[test*dim+2];
            vn2 = neighbor_basis_u2_phi[test];
            vn2x = neighbor_basis_u2_dphi[test*dim];
            vn2y = neighbor_basis_u2_dphi[test*dim+1];
            vn2z = neighbor_basis_u2_dphi[test*dim+2];

            local_col_index = 0;
            // u0-v2 block: 0
            local_col_index += u0dofpelm;

            // u1-v2 block : 0
            local_col_index+= u1dofpelm;

            // u2-v2 block: 0
            local_col_index+= u2dofpelm;

            // eg-v2 block:0
            local_col_index += egdofpelm;

            // p-v2 bloc: // TODO: why doesn't this vanish?
            for (trial=0; trial<pdofpelm;trial++){
              p = FE->var_spaces[ip]->phi[trial];
              pn = neighbor_basis_p_phi[trial];
              // p-v2
              kij = alpha * 0.5 * p * v2 * finrm[2];
              ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // pneighbor-v0
              kij = alpha * 0.5 * pn * v2*finrm[2];
              ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // p-v0neighbor
              kij = -alpha * 0.5 * p * vn2 * finrm[2];
              ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // pneighbor-v0neighbor
              kij = -alpha * 0.5 * pn *  vn2*finrm[2];
              ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // end p col
          } // End v2 loop
          local_row_index += u2dofpelm;
        } // end if dim==3

        // eg test rows
        for (test=0; test<egdofpelm;test++){

          local_col_index = 0;
          // u0-eg block:
          if(bEG) {
            for(trial=0; trial<u0dofpelm; trial++) {
              u0 = FE->var_spaces[iu0]->phi[trial];
              u0x = FE->var_spaces[iu0]->dphi[trial*dim];
              u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
              un0 = neighbor_basis_u0_phi[trial];
              un0x = neighbor_basis_u0_dphi[trial*dim];
              un0y = neighbor_basis_u0_dphi[trial*dim+1];
              if(dim==3) {
                u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];
                un0z = neighbor_basis_u0_dphi[trial*dim+2];
              }

              // u0-eg
              kij = (-0.5 * 2.0*u0x*finrm[0] - 0.5 * u0y*finrm[1])*eg_xterm;
              kij -= 0.5 * u0y*finrm[0]*eg_yterm;
              if(dim==3) kij -= 0.5 * (u0z*finrm[2]*eg_xterm + u0z*finrm[0]*eg_zterm);
              ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u0neighbor-eg
              kij = (-0.5 * 2.0*un0x*finrm[0] - 0.5 * un0y*finrm[1])*eg_xterm;
              kij -= 0.5 * un0y*finrm[0]*eg_yterm;
              if(dim==3) kij -= 0.5 * (un0z*finrm[2]*eg_xterm + un0z*finrm[0]*eg_zterm);
              ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u-egneighbor
              kij = (0.5 * 2.0*u0x*finrm[0] + 0.5 * u0y*finrm[1])*eg_xterm_neighbor;
              kij += 0.5 * u0y*finrm[0]*eg_yterm_neighbor;
              if(dim==3) kij += 0.5 * (u0z*finrm[2]*eg_xterm_neighbor + u0z*finrm[0]*eg_zterm_neighbor);
              ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // uneighbor-egneighbor
              kij = (0.5 * 2.0*un0x*finrm[0] + 0.5 * un0y*finrm[1])*eg_xterm_neighbor;
              kij += 0.5 * un0y*finrm[0]*eg_yterm_neighbor;
              if(dim==3) kij += 0.5 * (un0z*finrm[2]*eg_xterm_neighbor + un0z*finrm[0]*eg_zterm_neighbor);
              ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            } // u0 col
            local_col_index += u0dofpelm;

            // u1-eg block :
            for(trial=0; trial<u1dofpelm; trial++) {
              u1 = FE->var_spaces[iu1]->phi[trial];
              u1x = FE->var_spaces[iu1]->dphi[trial*dim];
              u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
              un1 = neighbor_basis_u1_phi[trial];
              un1x = neighbor_basis_u1_dphi[trial*dim];
              un1y = neighbor_basis_u1_dphi[trial*dim+1];
              if(dim==3) {
                u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
                un1z = neighbor_basis_u1_dphi[trial*dim+2];
              }

              // u1-eg
              kij = (-0.5 * 2.0*u1y*finrm[1] - 0.5 * u1x*finrm[0])*eg_yterm;
              kij -= 0.5 * u1x*finrm[1]*eg_xterm;
              if(dim==3) kij -= 0.5 * (u1z*finrm[2]*eg_yterm + u1z*finrm[1]*eg_zterm);
              ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u1neighbor-eg
              kij = (-0.5 * 2.0*un1y*finrm[1] - 0.5 * un1x*finrm[0])*eg_yterm;
              kij -= 0.5 * un1x*finrm[1]*eg_xterm;
              if(dim==3) kij -= 0.5 * (un1z*finrm[2]*eg_yterm + un1z*finrm[1]*eg_zterm);
              ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u-egneighbor
              kij = (0.5 * 2.0*u1y*finrm[1] + 0.5 * u1x*finrm[0])*eg_yterm_neighbor;
              kij += 0.5 * u1x*finrm[1]*eg_xterm_neighbor;
              if(dim==3) kij += 0.5 * (u1z*finrm[2]*eg_yterm_neighbor + u1z*finrm[1]*eg_zterm_neighbor);
              ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // uneighbor-egneighbor
              kij = (0.5 * 2.0*un1y*finrm[1] + 0.5 * un1x*finrm[0])*eg_yterm_neighbor;
              kij += 0.5 * un1x*finrm[1]*eg_xterm_neighbor;
              if(dim==3) kij += 0.5 * (un1z*finrm[2]*eg_yterm_neighbor + un1z*finrm[1]*eg_zterm_neighbor);
              ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // u1 col
            local_col_index+= u1dofpelm;

            if(dim==3) {
              // u2-eg block:
              for(trial=0; trial<u1dofpelm; trial++) {
                u2 = FE->var_spaces[iu2]->phi[trial];
                u2x = FE->var_spaces[iu2]->dphi[trial*dim];
                u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
                un2 = neighbor_basis_u2_phi[trial];
                un2x = neighbor_basis_u2_dphi[trial*dim];
                un2y = neighbor_basis_u2_dphi[trial*dim+1];
                if(dim==3) {
                  u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];
                  un2z = neighbor_basis_u2_dphi[trial*dim+2];
                }

                // u1-eg
                kij = (-0.5 * 2.0*u2z*finrm[2] - 0.5 * u2x*finrm[0] - 0.5 * u2y*finrm[1])*eg_zterm;
                kij -= 0.5 * u2x*finrm[2]*eg_xterm;
                kij -= 0.5 * u2y*finrm[2]*eg_yterm;
                ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

                // u1neighbor-eg
                kij = (-0.5 * 2.0*un2z*finrm[2] - 0.5 * un2x*finrm[0] - 0.5 * un2y*finrm[1])*eg_zterm;
                kij -= 0.5 * un2x*finrm[2]*eg_xterm;
                kij -= 0.5 * un2y*finrm[2]*eg_yterm;
                ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

                // u-egneighbor
                kij = (0.5 * 2.0*u2z*finrm[2] + 0.5 * u2x*finrm[0] + 0.5 *u2y*finrm[1])*eg_zterm_neighbor;
                kij += 0.5 * u2x*finrm[2]*eg_xterm_neighbor;
                kij += 0.5 * u2y*finrm[2]*eg_yterm_neighbor;
                ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

                // uneighbor-egneighbor
                kij = (0.5 * 2.0*un2x*finrm[2] + 0.5 * un2x*finrm[0] + 0.5*un2y*finrm[1])*eg_zterm_neighbor;
                kij += 0.5 * un2x*finrm[2]*eg_xterm_neighbor;
                kij += 0.5 * un2y*finrm[2]*eg_yterm_neighbor;
                ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
              } // u1 col
              local_col_index+= u2dofpelm;
            } // if dimm=3
          } // if bEG

          // TODO: Don't know why this isn't under if(bEG)?
          // eg-eg block:
          // Reset local_col_index to be safe
          local_col_index = u0dofpelm+u1dofpelm;
          if(dim==3) local_col_index += u2dofpelm;
          for (trial=0; trial<egdofpelm;trial++){

            // eg-eg
            kij = -0.5* 2.0 * 1. * finrm[0] * eg_xterm;
            kij -= 0.5* 2.0 * 1. * finrm[1] * eg_yterm;
            if(dim==3) kij -= 0.5* 2.0 * 1. * finrm[2] * eg_zterm;
            //penalty term
            kij += penalty_term * eg_xterm*eg_xterm;
            kij += penalty_term * eg_yterm*eg_yterm;
            if(dim==3) kij += penalty_term * eg_zterm*eg_zterm;
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            //egneighbor-eg
            kij_unv = 0.5* -2.0 * 1. * finrm[0] * eg_xterm;
            kij_unv -= 0.5*2.0 * 1. * finrm[1] * eg_yterm;
            if(dim==3) kij_unv -= 0.5*2.0 * 1. * finrm[2] * eg_zterm;
            //penalty_term
            kij_unv += -penalty_term * eg_xterm_neighbor*eg_xterm;
            kij_unv += -penalty_term * eg_yterm_neighbor*eg_yterm;
            if(dim==3) kij_unv+= -penalty_term * eg_zterm_neighbor*eg_zterm;
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij_unv;

            // eg-egneighbor
            kij = 0.5* 2.0 * 1. * finrm[0] * eg_xterm_neighbor;
            kij += 0.5*2.0 * 1. * finrm[1] * eg_yterm_neighbor;
            if(dim==3) kij += 0.5*2.0 * 1. * finrm[2] * eg_zterm_neighbor;
            //penalty_term
            kij += -penalty_term * eg_xterm*eg_xterm_neighbor;
            kij += -penalty_term * eg_yterm*eg_yterm_neighbor;
            if(dim==3) kij += -penalty_term * eg_zterm*eg_zterm_neighbor;
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // egneighbor-egneighbor
            kij = 0.5* 2.0 * 1. * finrm[0] * eg_xterm_neighbor;
            kij += 0.5*2.0 * 1. * finrm[1] * eg_yterm_neighbor;
            if(dim==3) kij += 0.5*2.0 * 1. + finrm[2] * eg_zterm_neighbor;
            //penalty_term
            kij += penalty_term * eg_xterm_neighbor*eg_xterm_neighbor;
            kij += penalty_term * eg_yterm_neighbor*eg_yterm_neighbor;
            if(dim==3) kij += penalty_term * eg_zterm_neighbor*eg_zterm_neighbor;
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

          } // eg col loop
          local_col_index += egdofpelm;

          // p-eg bloc:
          for (trial=0; trial<pdofpelm;trial++){
            p = FE->var_spaces[ip]->phi[trial];
            pn = neighbor_basis_p_phi[trial];

            // p-eg
            kij = alpha * 0.5 * p * (finrm[0]*eg_xterm + finrm[1]*eg_yterm);
            if(dim==3) kij += alpha * 0.5 * p * finrm[2]*eg_zterm;
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-eg
            kij = alpha * 0.5 * pn * (finrm[0]*eg_xterm + finrm[1]*eg_yterm);
            if(dim==3) kij += alpha * 0.5 * pn * finrm[2]*eg_zterm;
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // p-egneighbor
            kij = -alpha * 0.5 * p * (finrm[0]*eg_xterm_neighbor + finrm[1]*eg_yterm_neighbor);
            if(dim==3) kij -= alpha * 0.5 * p * finrm[2]*eg_zterm_neighbor;
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // pneighbor-egneighbor
            kij = -alpha * 0.5 * pn * (finrm[0]*eg_xterm_neighbor + finrm[1]*eg_yterm_neighbor);
            if(dim==3) kij -= alpha * 0.5 * pn * finrm[2]*eg_zterm_neighbor;
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // end p col
        } // End eg loop
        local_row_index += egdofpelm;

        // p test rows
        for (test=0; test<pdofpelm;test++){
          q = FE->var_spaces[ip]->phi[test];
          qn = neighbor_basis_p_phi[test];

          local_col_index = 0;
          // u0-q block: // TODO: Also not sure why this isn't 0
          for(trial=0; trial<u0dofpelm; trial++) {
            u0 = FE->var_spaces[iu0]->phi[trial];
            u0x = FE->var_spaces[iu0]->dphi[trial*dim];
            u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
            un0 = neighbor_basis_u0_phi[trial];
            un0x = neighbor_basis_u0_dphi[trial*dim];
            un0y = neighbor_basis_u0_dphi[trial*dim+1];
            if(dim==3) {
              u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];
              un0z = neighbor_basis_u0_dphi[trial*dim+2];
            }

            // u0-q
            kij = alpha * 0.5 * q * u0 * finrm[0];
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u0neighbor-q
            kij = -alpha * 0.5 * q * un0 * finrm[0];
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u0-q_neighbor
            kij = alpha * 0.5 * qn * u0 * finrm[0];
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u0neighbor-q_neighbor
            kij = -alpha * 0.5 * qn * un0 * finrm[0];
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // end u0 col
          local_col_index += u0dofpelm;

          // u1-q block : // TODO: Also not sure why this isn't 0
          for(trial=0; trial<u1dofpelm; trial++) {
            u1 = FE->var_spaces[iu1]->phi[trial];
            un1 = neighbor_basis_u1_phi[trial];

            // u1-q
            kij = alpha * 0.5 * q * u1 * finrm[1];
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u1neighbor-q
            kij = -alpha * 0.5 * q * un1 * finrm[1];
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u1-q_neighbor
            kij = alpha * 0.5 * qn * u1 * finrm[1];
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // u1neighbor-q_neighbor
            kij = -alpha * 0.5 * qn * un1 * finrm[1];
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // end u1 col
          local_col_index+= u1dofpelm;

          if(dim==3) {
            // u2-q block : // TODO: Also not sure why this isn't 0
            for(trial=0; trial<u2dofpelm; trial++) {
              u2 = FE->var_spaces[iu2]->phi[trial];
              un2 = neighbor_basis_u2_phi[trial];

              // u2-q
              kij = alpha * 0.5 * q * u2 * finrm[2];
              ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u0neighbor-q
              kij = -alpha * 0.5 * q * un2 * finrm[2];
              ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u0-q_neighbor
              kij = alpha * 0.5 * qn * u2 * finrm[2];
              ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

              // u0neighbor-q_neighbor
              kij = -alpha * 0.5 * qn * un2 * finrm[2];
              ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
            } // end u2 col
            local_col_index+= u2dofpelm;
          } // end if dim==3

          // eg-q block: // TODO: Should this be under if(bEG)??
          for (trial=0; trial<egdofpelm;trial++){
            // eg-q
            kij = alpha * 0.5 * q * (eg_xterm*finrm[0] + eg_yterm*finrm[1]);
            if(dim==3) kij += alpha * 0.5 * q * eg_zterm*finrm[2];
            ALoc_u_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // egneighbor-q
            kij = -alpha * 0.5 * q * (eg_xterm_neighbor*finrm[0] + eg_yterm_neighbor*finrm[1]);
            if(dim==3) kij -= alpha * 0.5 * q * eg_zterm_neighbor*finrm[2];
            ALoc_uneighbor_v[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // eg-q_neighbor
            kij = alpha * 0.5 * qn * (eg_xterm*finrm[0] + eg_yterm*finrm[1]);
            if(dim==3) kij += alpha * 0.5 * qn * eg_zterm*finrm[2];
            ALoc_u_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;

            // egneighbor-q_neighbor
            kij = -alpha * 0.5 * qn * (eg_xterm_neighbor*finrm[0] + eg_yterm_neighbor*finrm[1]);
            if(dim==3) kij -= alpha * 0.5 * qn * eg_zterm_neighbor*finrm[2];
            ALoc_uneighbor_vneighbor[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w_face*kij;
          } // eg col loop
          local_col_index += egdofpelm;

          // p-q bloc: 0
        } // End q loop
        local_row_index = pdofpelm;
      } // end if(counter==2)
    }// end quadrature loop

    // Map local to global
    if(counter == 1) {
      block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);
    } else if(counter == 2) {
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


    // Frees
    if(ALoc) free(ALoc);
    if(ALoc_u_v) free(ALoc_u_v);
    if(ALoc_u_vneighbor) free(ALoc_u_vneighbor);
    if(ALoc_uneighbor_v) free(ALoc_uneighbor_v);
    if(ALoc_uneighbor_vneighbor) free(ALoc_uneighbor_vneighbor);

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

    if(cq_face) {
      free_qcoords(cq_face);
      free(cq_face);
      cq_face=NULL;
    }

    if(data_face) free(data_face);

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
    if(dof_on_elm) free(dof_on_elm);
    if(dof_on_elm_neighbor) free(dof_on_elm_neighbor);
    if(v_on_elm) free(v_on_elm);
    if(v_on_elm_neighbor) free(v_on_elm_neighbor);


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
  // TODO: I am here!
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
      INT i,quad,test,idim;
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
      REAL* val_true_D = (REAL *) calloc(nun*dim,sizeof(REAL));

      // for the preivous time step value
      //REAL* val_sol = (REAL *) calloc(nun,sizeof(REAL));


      coordinates *barycenter = allocatecoords(1,dim);
      barycenter->x[0] = mesh->el_mid[elm*dim];
      barycenter->y[0] = mesh->el_mid[elm*dim + 1];
      if(dim==3) barycenter->z[0] = mesh->el_mid[elm*dim +2];
      REAL eg_xterm,eg_yterm,eg_zterm;


      //  Sum over quadrature points
      for (quad=0;quad<cq->nq_per_elm;quad++) {
        qx[0] = cq->x[elm*cq->nq_per_elm+quad];
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
        if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
        w = cq->w[elm*cq->nq_per_elm+quad];

        (*rhs)(rhs_val,qx,time,&(mesh->el_flag[elm]));

        // EG stuff
        eg_xterm = qx[0]-barycenter->x[0];
        eg_yterm = qx[1]-barycenter->y[0];
        if(dim==3) eg_zterm = qx[2]-barycenter->z[0];

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
          get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);

          for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
            //SLEE
            if(i<dim) {
              //This is for  u0 and u1 (CG part)
              //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
              bLoc[(local_row_index+test)] += w*rhs_val[i]*FE->var_spaces[i]->phi[test];
            } else if(i == dim) {
              //SLEE
              // This is for  u2:: (EG part)
              //Note that new DOF is \phi^3 = [x ; y]
              //printf(" i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
              bLoc[(local_row_index+test)] += w*(rhs_val[0]*eg_xterm + rhs_val[1]*eg_yterm);
              if(dim==3) bLoc[(local_row_index+test)] += w*(rhs_val[2]*eg_zterm);
            } else if(i == dim+1) {
              bLoc[(local_row_index+test)] += w*rhs_val[dim+1]*FE->var_spaces[dim+1]->phi[test];
            }
          } // end test loop

          local_row_index  += FE->var_spaces[i]->dof_per_elm;
          // TODO: The following line was missing...
          local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
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
          if(dim>1) finrm[1]=mesh->f_norm[face*dim+1];
          if(dim>2) finrm[2]=mesh->f_norm[face*dim+2];

          for(jkl=mesh->f_v->IA[face];jkl<mesh->f_v->IA[face+1];jkl++){

            //printf("** jkl = %d, mesh->f_v->IA[face] = %d, mesh->f_v->IA[face+1] = %d \n", jkl, mesh->f_v->IA[face], mesh->f_v->IA[face+1]);

            j=jkl-mesh->f_v->IA[face];
            k=mesh->f_v->JA[jkl];

            //printf("** j = %d, k = %d \n", j, k );

            xfi[j*dim+0]=mesh->cv->x[k];
            if(dim>1) xfi[j*dim+1]=mesh->cv->y[k];
            if(dim>2) xfi[j*dim+2]=mesh->cv->z[k];

            //printf("** xfi[j*dim+0] = %f,  xfi[j*dim+1] = %f \n",  xfi[j*dim+0] ,  xfi[j*dim+1]);
          }

          // FOR NEIGHBOR..
          local_dof_on_elm_face_interface = NULL;
          local_dof_on_elm_neighbor = NULL;


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

// TODO: Do you really mean edge here?? -> This isn't used so I'm commenting it out
//REAL edge_length = mesh->ed_len[face];

if(mesh->f_flag[face]>0) {

  for (quad_face=0;quad_face<nquad_face;quad_face++) {

    qx_face[0] = cq_face->x[quad_face];
    qx_face[1] = cq_face->y[quad_face];
    if(dim==3) qx_face[2] = cq_face->z[quad_face];

    w_face = cq_face->w[quad_face];

    // EG stuff
    eg_xterm = qx_face[0]-barycenter->x[0];
    eg_yterm = qx_face[1]-barycenter->y[0];
    if(dim==3) eg_zterm = qx_face[2]-barycenter->z[0];

    // TODO: why was this hard coded for 2D?
    (*truesol)(val_true_face,qx_face,time,&(mesh->el_flag[elm]));
    //(*exact_sol2D)(val_true_face,qx_face,time, &(mesh->el_flag[elm]));  // ???

    (*truesol)(val_true_face_n,qx_face,time-timestep,&(mesh->el_flag[elm]));
    //(*exact_sol2D)(val_true_face_n,qx_face,time-timestep,&(mesh->el_flag[elm]));  // ???

    //true solution is in n+1 time
    (*truesol_dt)(val_true_dt_face,qx_face,time,&(mesh->el_flag[elm]));

    local_row_index=0;
    unknown_index=0;
    local_dof_on_elm_face = dof_on_elm;

    for(i=0;i<FE->nspaces;i++) {

      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx_face,v_on_elm,local_dof_on_elm_face,mesh,FE->var_spaces[i]);

      for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
        //SLEE
        if(i<dim) {
          bLoc[(local_row_index+test)] += BC_penalty_term * w_face*(val_true_face[unknown_index]*  FE->var_spaces[i]->phi[test]);

          //DEBUG100
          //bLoc[(local_row_index+test)] += penalty_term * w_face*
          //(val_true_face[3]*  FE->var_spaces[i]->phi[test] *finrm[i]);

        } else if(i == dim) {
          //SLEE
          // Note that new DOF is \phi^3 = [x ; y]
          //printf("i = %d (FE->nspaces = %d), unknown_index = %d, test = %d \n", i, FE->var_spaces[i]->dof_per_elm, unknown_index, test);
          bLoc[(local_row_index+test)] += BC_penalty_term *  w_face * (val_true_face[0] * eg_xterm + val_true_face[1]*eg_yterm);
          if(dim==3) bLoc[(local_row_index+test)] += BC_penalty_term * w_face * val_true_face[2] * eg_zterm;

        } else if(i == dim+1) {
          //Stokes 3
          bLoc[(local_row_index+test)] += w_face * (val_true_face[0] * finrm[0] + val_true_face[1] * finrm[1])* FE->var_spaces[dim+1]->phi[test];
          if(dim==3) bLoc[(local_row_index+test)] += w_face * (val_true_face[2] * finrm[2])* FE->var_spaces[dim+1]->phi[test];
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
if(f_el) free(f_el);
//printf("FEREE DONE\n");
}// for each face

}

block_LocaltoGlobal_RHS(dof_on_elm,FE,b,bLoc);

// Frees
if(cq_face) {
  free_qcoords(cq_face);
  free(cq_face);
  cq_face=NULL;
}
if(val_true_face) free(val_true_face);
if(val_true_face_n) free(val_true_face_n);
if(val_true_face_n_neighbor) free(val_true_face_n_neighbor);
if(val_true_dt_face) free(val_true_dt_face);
if(v_on_elm_neighbor) free(v_on_elm_neighbor);
if(dof_on_elm_neighbor) free(dof_on_elm_neighbor);
if(val_true) free(val_true);
if(val_true_D) free(val_true_D);
if(barycenter) {
  free_coords(barycenter);
  free(barycenter);
  barycenter = NULL;
}
if(data_face) free(data_face);
//printf("RHS ASSEMBLE END \n");
return;
}

// TODO: Here

void local_assembly_Elasticity(block_dCSRmat* A,dvector *b,REAL* ALoc, block_fespace *FE,mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time,REAL timestep,  INT* switch_on_face)
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
  for (i=0; i<FE->nspaces;i++) dof_per_elm += FE->var_spaces[i]->dof_per_elm;

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
  if(dim==3) barycenter->z[0] = mesh->el_mid[elm*dim+2];

  //printf("ELEMENT = %D, x = %f , y = %f \n", elm,  barycenter->x[0],   barycenter->y[0]);
  INT local_size = dof_per_elm*dof_per_elm;
  //REAL* ALoc_neighbor = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  double lambda = LAME_LAMBDA_GLOBAL ;

  INT iu0 = 0;
  INT iu1 = 1;
  INT iu2 = 2;
  INT ieg = dim;
  INT ip = dim+1;
  INT u0dofpelm = FE->var_spaces[iu0]->dof_per_elm;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dofpelm;
  if(dim==3) {
    u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  }
  INT egdofpelm = FE->var_spaces[ieg]->dof_per_elm;
  INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;
  REAL u0,u1,u2,p,v0,v1,v2,q,u0x,u0y,u0z,u1x,u1y,u1z,u2x,u2y,u2z,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {

    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];

    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u = (u1,u2,u3) and v = (v1,v2,v3)
    local_dof_on_elm = dof_on_elm;
    for(idim=0;idim<FE->nspaces;idim++) {
      get_FEM_basis(FE->var_spaces[idim]->phi,FE->var_spaces[idim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[idim]);
      local_dof_on_elm += FE->var_spaces[idim]->dof_per_elm;
    }

    // v0 test rows
    local_row_index = 0;
    for(test=0;test<u0dofpelm;test++) {
      v0 = FE->var_spaces[iu0]->phi[test];
      v0x = FE->var_spaces[iu0]->dphi[test*dim];
      v0y = FE->var_spaces[iu0]->dphi[test*dim+1];
      if(dim==3) v0z = FE->var_spaces[iu0]->dphi[test*dim+2];
      local_col_index = 0;

      // u0-v0 XX
      for (trial=0; trial<u0dofpelm;trial++){
        u0 = FE->var_spaces[iu0]->phi[trial];
        u0x = FE->var_spaces[iu0]->dphi[trial*dim];
        u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
        if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

        kij = 2.0*u0x*v0x + u0y*v0y;
        if(dim==3) kij += u0z*v0z;

        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      } // end u0 cols
      local_col_index += u0dofpelm;

      // u1-v0 block : XX
      for (trial=0; trial<u1dofpelm;trial++){
        u1 = FE->var_spaces[iu1]->phi[trial];
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

        kij = u1x*v0y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // End u1 loop
      local_col_index+= u1dofpelm;

      if(dim==3) {
        // u2-v0 block: XX
        for (trial=0; trial<u2dofpelm;trial++){
          u2 = FE->var_spaces[iu2]->phi[trial];
          u2x = FE->var_spaces[iu2]->dphi[trial*dim];
          u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

          kij = u2x*v0z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

        } // End u2 loop
        local_col_index+= u2dofpelm;
      } // end if dim==3

      if(bEG) {
        // eg-v0: XX
        for (trial=0; trial<egdofpelm;trial++){

          kij = 2.0*v0x;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // eg cols
      }
      local_col_index += egdofpelm;

      // p-v0 bloc: XX
      for (trial=0; trial<pdofpelm;trial++){
        p = FE->var_spaces[ip]->phi[trial];

        kij = -alpha*p*v0x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      } // pcols

    } // end v0
    local_row_index += u0dofpelm;

    // v1 test rows
    for(test=0;test<u1dofpelm;test++) {
      v1 = FE->var_spaces[iu1]->phi[test];
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
      if(dim==3) v1z = FE->var_spaces[iu1]->dphi[test*dim+2];
      local_col_index = 0;

      // u0-v1 XX
      for (trial=0; trial<u0dofpelm;trial++){
        u0 = FE->var_spaces[iu0]->phi[trial];
        u0x = FE->var_spaces[iu0]->dphi[trial*dim];
        u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
        if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

        kij = u0y*v1x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // end u0 cols
      local_col_index += u0dofpelm;

      // u1-v1 block : XX
      for (trial=0; trial<u1dofpelm;trial++){
        u1 = FE->var_spaces[iu1]->phi[trial];
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

        kij = u1x*v1x + 2.0*u1y*v1y;
        if(dim==3) kij += u1z*v1z;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // End u1 loop
      local_col_index+= u1dofpelm;

      if(dim==3) {
        // u2-v1 block:XX
        for (trial=0; trial<u2dofpelm;trial++){
          u2 = FE->var_spaces[iu2]->phi[trial];
          u2x = FE->var_spaces[iu2]->dphi[trial*dim];
          u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

          kij = u2y*v1z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

        } // End u2 loop
        local_col_index+= u2dofpelm;
      } // end if dim==3

      if(bEG) {
        // eg-v1: XX
        for (trial=0; trial<egdofpelm;trial++){

          kij = 2.0*v1y;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // eg cols
      }
      local_col_index += egdofpelm;

      // p-v1 bloc: XX
      for (trial=0; trial<pdofpelm;trial++){
        p = FE->var_spaces[ip]->phi[trial];

        kij = -alpha * p*v1y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      } // pcols

    } // end v1
    local_row_index += u1dofpelm;

    if(dim==3) {
      // v2 test rows
      for(test=0;test<u2dofpelm;test++) {
        v2 = FE->var_spaces[iu2]->phi[test];
        v2x = FE->var_spaces[iu2]->dphi[test*dim];
        v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
        v2z = FE->var_spaces[iu2]->dphi[test*dim+2];
        local_col_index = 0;

        // u0-v2
        for (trial=0; trial<u0dofpelm;trial++){
          u0 = FE->var_spaces[iu0]->phi[trial];
          u0x = FE->var_spaces[iu0]->dphi[trial*dim];
          u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
          u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

          kij = u0z*v2x;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

        } // end u0 cols
        local_col_index += u0dofpelm;

        // u1-v2 block :
        for (trial=0; trial<u1dofpelm;trial++){
          u1 = FE->var_spaces[iu1]->phi[trial];
          u1x = FE->var_spaces[iu1]->dphi[trial*dim];
          u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
          u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

          kij = u1z*v2y;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // End u1 loop
        local_col_index+= u1dofpelm;

        // u2-v2 block:
        for (trial=0; trial<u2dofpelm;trial++){
          u2 = FE->var_spaces[iu2]->phi[trial];
          u2x = FE->var_spaces[iu2]->dphi[trial*dim];
          u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

          kij = u2x*v2x + u2y*v2y + 2.0*u2z*v2z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // End u2 loop
        local_col_index+= u2dofpelm;

        // eg-v2:
        if(bEG) {
          for (trial=0; trial<egdofpelm;trial++){

            kij = 2.0*v2z;
            ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
          } // eg cols
        }
        local_col_index += egdofpelm;

        // p-v2 bloc:
        for (trial=0; trial<pdofpelm;trial++){
          p = FE->var_spaces[ip]->phi[trial];

          kij = -alpha * p*v2z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

        } // pcols

      } // end v2
      local_row_index += u2dofpelm;
    } // end dim==3

    // eg test rows
    for(test=0;test<egdofpelm;test++) {

      local_col_index = 0;

      if(bEG) {
        // u0-eg XX
        for (trial=0; trial<u0dofpelm;trial++){
          u0 = FE->var_spaces[iu0]->phi[trial];
          u0x = FE->var_spaces[iu0]->dphi[trial*dim];

          kij = 2.0*u0x;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

        } // end u0 cols
        local_col_index += u0dofpelm;

        // u1-eg block : XX
        for (trial=0; trial<u1dofpelm;trial++){
          u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];

          kij = 2.0*u1y;          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // End u1 loop
        local_col_index+= u1dofpelm;

        if(dim==3) {
          // u2-eg block: XX
          for (trial=0; trial<u2dofpelm;trial++){
            u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

            kij = 2.0*u2z;
            ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
          } // End u2 loop
          local_col_index+= u2dofpelm;
        } // end if dim==3
      } // end if bEG

      // eg-eg: XX
      for (trial=0; trial<egdofpelm;trial++){

        kij = 2.0*dim; // (2.0* (1,1,1)*(1,1,1))
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // eg cols
      local_col_index += egdofpelm;

      // p-eg bloc: XX
      if(bEG) {
        for (trial=0; trial<pdofpelm;trial++){
          p = FE->var_spaces[ip]->phi[trial];

          kij = -alpha*p*dim;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // pcols
      }

    } // end eg
    local_row_index += egdofpelm;

    // q test rows
    for(test=0;test<pdofpelm;test++) {
      q = FE->var_spaces[ip]->phi[test];
      local_col_index = 0;

      // u0-q XX
      for (trial=0; trial<u0dofpelm;trial++){
        u0 = FE->var_spaces[iu0]->phi[trial];
        u0x = FE->var_spaces[iu0]->dphi[trial*dim];
        u0y = FE->var_spaces[iu0]->dphi[trial*dim+1];
        if(dim==3) u0z = FE->var_spaces[iu0]->dphi[trial*dim+2];

        kij = -alpha*u0x*q;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // end u0 cols
      local_col_index += u0dofpelm;

      // u1-q block : XX
      for (trial=0; trial<u1dofpelm;trial++){
        u1 = FE->var_spaces[iu1]->phi[trial];
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        if(dim==3) u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];

        kij = -alpha*u1y*q;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;

      } // End u1 loop
      local_col_index+= u1dofpelm;

      if(dim==3) {
        // u2-q block:
        for (trial=0; trial<u2dofpelm;trial++){
          u2 = FE->var_spaces[iu2]->phi[trial];
          u2x = FE->var_spaces[iu2]->dphi[trial*dim];
          u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z = FE->var_spaces[iu2]->dphi[trial*dim+2];

          kij = -alpha*u2z*q;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // End u2 loop
        local_col_index+= u2dofpelm;
      } // end if dim==3

      // eg-q: XX
      if(bEG) {
        for (trial=0; trial<egdofpelm;trial++){

          kij = -alpha*dim*q;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        } // eg cols
      } // end bEG
      local_col_index += egdofpelm;
    } // end q

  }//QUAD

  block_LocaltoGlobal_neighbor(dof_on_elm,dof_on_elm,FE,A,ALoc);

  // Frees
  if(dof_on_elm_neighbor) free(dof_on_elm_neighbor);
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

  return;
}

        /*********************************************************************/
