/*! \file src/assemble/assemble_utils.c   
 *
 * \brief This code will contain all the tools needed to build stiffness matrices
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 4/22/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/11/2016
 */

#include "hazmath.h"

/******************************************************************************************************/
void create_CSR_rows(dCSRmat *A, fespace *FE)
{
  /*!
   * \fn void create_CSR_rows(dCSRmat *A, fespace *FE)
   *
   * \brief Computes the "possible" number of nonzeros for the global stiffness matrix
   *        Ignores Dirichlet boundary conditions, which can be eliminated later.
   *        Entry will be nonzero if degree of freedom belongs to the element (i.e. row)
   *        Builds the IA array for A and assigns nnz.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.IA         Row structure of CSR matrix
   * \return A.nnz        Number of nonzeros of A
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=1;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Loop over all Elements connected to particular DOF
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el.JA[j-1];
      k_a = FE->el_dof->IA[if1-1];
      k_b = FE->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
        mydof = FE->el_dof->JA[k-1];
        if (ix[mydof-1]!=i+1) { /* We haven't been here  */
          icp++;
          ix[mydof-1] = i+1;
        }
      }
    }
  }
  A->IA[ndof] = icp;
  A->nnz = icp-1;

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_rows_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
{
  /*!
   * \fn void create_CSR_rows_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
   *
   * \brief Computes the "possible" number of nonzeros for the global stiffness matrix
   *        Ignores Dirichlet boundary conditions, which can be eliminated later.
   *        Entry will be nonzero if degree of freedom belongs to the element (i.e. row)
   *        Builds the IA array for A and assigns nnz.
   *
   * \note Assumes non-square system where test and trial are from different FE spaces
   *
   * \param FE1           FE Space for trial functions (u)
   * \param FE2           FE Space for test functions (v)
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.IA         Row structure of CSR matrix
   * \return A.nnz        Number of nonzeros of A
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2;
  icsr_trans_1(FE2->el_dof,&dof_el_2);

  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;

  INT* ix = (INT *) calloc(ncols,sizeof(INT));
  for (i=0; i<ncols; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF of test space and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=1;
  for (i=0; i<nrows; i++) {
    A->IA[i] = icp;
    // Loop over all Elements connected to particular DOF of test space
    j_a = dof_el_2.IA[i];
    j_b = dof_el_2.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el_2.JA[j-1];
      // For this given element grab the DOF in the trial space
      k_a = FE1->el_dof->IA[if1-1];
      k_b = FE1->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
        mydof = FE1->el_dof->JA[k-1];
        if (ix[mydof-1]!=i+1) { /* We haven't been here  */
          icp++;
          ix[mydof-1] = i+1;
        }
      }
    }
  }
  A->IA[nrows] = icp;
  A->nnz = icp-1;

  if(ix) free(ix);
  icsr_free(&dof_el_2);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_rows_withBC(dCSRmat *A, fespace *FE)
{
  /*!
   * \fn void create_CSR_rows_withBC(dCSRmat *A, fespace *FE)
   *
   * \brief Computes the "possible" number of nonzeros for the global stiffness matrix
   *        Also takes into account Dirichlet boundary conditions.  If the DOF of the corresponding row
   *        is a boundary, only include 1 nonzero for that row (it will be identity!)
   *        Entry will be nonzero if degree of freedom belongs to the element (i.e. row)
   *        Builds the IA array for A and assigns nnz.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.IA         Row structure of CSR matrix
   * \return A.nnz        Number of nonzeros of A
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=1;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Check if this is a boundary edge
    if (FE->dirichlet[i]==1) {	/* Only 1 nonzero this row */
      icp++;
    } else {			/* It's not a boundary and compute as normal */
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
        if1 = dof_el.JA[j-1];
        k_a = FE->el_dof->IA[if1-1];
        k_b = FE->el_dof->IA[if1]-1;
        for (k=k_a; k<=k_b; k++) {
          mydof = FE->el_dof->JA[k-1];
          if (ix[mydof-1]!=i+1 && FE->dirichlet[mydof-1]==0) { /* We haven't been here AND it's not a boundary */
            icp++;
            ix[mydof-1] = i+1;
          }
        }
      }
    }
  }
  A->IA[ndof] = icp;
  A->nnz = icp-1;

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_rows_flag(dCSRmat *A, fespace *FE,INT flag)
{
  /*!
   * \fn void create_CSR_rows_flag(dCSRmat *A,fespace *FE,INT flag)
   *
   * \brief Computes the "possible" number of nonzeros for the global stiffness matrix
   *        Also takes into account special boundary conditions.  Here if the boundary
   *        of the DOF is equal to "flag" (whatever that corresponds to)
   *        then we add something to the matrix.
   *        Builds the IA array for A and assigns nnz.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   * \param flag          Indicates boundaries we care about
   *
   * \return A.IA         Row structure of CSR matrix
   * \return A.nnz        Number of nonzeros of A
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=1;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Check if this is a boundary dof we want
    if (FE->dof_flag[i]==flag) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
        if1 = dof_el.JA[j-1];
        k_a = FE->el_dof->IA[if1-1];
        k_b = FE->el_dof->IA[if1]-1;
        for (k=k_a; k<=k_b; k++) {
          mydof = FE->el_dof->JA[k-1];
          if (ix[mydof-1]!=i+1 && FE->dof_flag[mydof-1]==flag) { /* We haven't been here AND it is a boundary */
            icp++;
            ix[mydof-1] = i+1;
          }
        }
      }
    }
  }
  A->IA[ndof] = icp;
  A->nnz = icp-1;

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_cols(dCSRmat *A, fespace *FE)
{
  /*!
   * \fn void create_CSR_cols(dCSRmat *A, fespace *FE)
   *
   * \brief Finds the column sparsity structure of the Global Stiffness Matrix
   *        Ignores Dirichlet boundary conditions, which can be eliminated later.
   *        Builds the JA array for A.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.JA         Columns of CSR matrix
   *
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and build A->JA
  icp=1;
  for (i=0; i<ndof; i++) {
    // Loop over all Elements connected to particular edge
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el.JA[j-1];
      k_a = FE->el_dof->IA[if1-1];
      k_b = FE->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
        mydof = FE->el_dof->JA[k-1];
        if (ix[mydof-1]!=i+1) { /* We haven't been here */
          A->JA[icp-1] = mydof;
          icp++;
          ix[mydof-1] = i+1;
        }
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_cols_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
{
  /*!
   * \fn void create_CSR_cols_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
   *
   * \brief Finds the column sparsity structure of the Global Stiffness Matrix
   *        Ignores Dirichlet boundary conditions, which can be eliminated later.
   *        Builds the JA array for A.
   *
   * \note Assumes non-square system where test and trial are from different FE spaces
   *
   * \param FE1           FE Space for trial functions (u)
   * \param FE2           FE Space for test functions (v)
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.JA         Columns of CSR matrix
   *
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2;
  icsr_trans_1(FE2->el_dof,&dof_el_2);

  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;

  INT* ix = (INT *) calloc(ncols,sizeof(INT));
  for (i=0; i<ncols; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF of test space and build A->JA
  icp=1;
  for (i=0; i<nrows; i++) {
    // Loop over all Elements connected to particular edge
    j_a = dof_el_2.IA[i];
    j_b = dof_el_2.IA[i+1]-1;
    for (j=j_a; j<=j_b; j++) {
      if1 = dof_el_2.JA[j-1];
      k_a = FE1->el_dof->IA[if1-1];
      k_b = FE1->el_dof->IA[if1]-1;
      for (k=k_a; k<=k_b; k++) {
        mydof = FE1->el_dof->JA[k-1];
        if (ix[mydof-1]!=i+1) { /* We haven't been here */
          A->JA[icp-1] = mydof;
          icp++;
          ix[mydof-1] = i+1;
        }
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el_2);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_cols_withBC(dCSRmat *A, fespace *FE)
{
  /*!
   * \fn void create_CSR_cols_withBC(dCSRmat *A, fespace *FE)
   *
   * \brief Finds the column sparsity structure of the Global Stiffness Matrix
   *        Also takes into account Dirichlet boundary conditions.  If the DOF of the corresponding row
   *        is a boundary, only include 1 nonzero for that row (it will be identity!)
   *        Builds the JA array for A.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   *
   * \return A.JA         Columns of CSR matrix
   *
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and build A->JA
  icp=1;
  for (i=0; i<ndof; i++) {
    // Check if this is a boundary edge
    if (FE->dirichlet[i]==1) {	/* Only 1 nonzero this row */
      A->JA[icp-1]=i+1;
      icp++;
    } else {			/* It's not a boundary and compute as normal */
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
        if1 = dof_el.JA[j-1];
        k_a = FE->el_dof->IA[if1-1];
        k_b = FE->el_dof->IA[if1]-1;
        for (k=k_a; k<=k_b; k++) {
          mydof = FE->el_dof->JA[k-1];
          if (ix[mydof-1]!=i+1 && FE->dirichlet[mydof-1]==0) { /* We haven't been here AND it's not a boundary */
            A->JA[icp-1] = mydof;
            icp++;
            ix[mydof-1] = i+1;
          }
        }
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void create_CSR_cols_flag(dCSRmat *A, fespace *FE,INT flag)
{
  /*!
   * \fn void create_CSR_cols_flag(dCSRmat *A, fespace *FE,INT flag)
   *
   * \brief Finds the column sparsity structure of the Global Stiffness Matrix
   *        Also takes into account special boundary conditions.  Here if the boundary
   *        of the DOF is equal to "flag" (whatever that corresponds to)
   *        then we add something to the matrix.
   *        Builds the JA array for A.
   *
   * \param FE            FE Space
   * \param A             dCSRmat Stiffness Matrix
   * \param flag          Indicates which boundary DOF to grab
   *
   * \return A.JA         Columns of CSR matrix
   *
   */

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans_1(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = 0;
  }

  // Loop over all DOF and build A->JA
  icp=1;
  for (i=0; i<ndof; i++) {
    // Check if this is a boundary edge
    if (FE->dof_flag[i]==flag) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1]-1;
      for (j=j_a; j<=j_b; j++) {
        if1 = dof_el.JA[j-1];
        k_a = FE->el_dof->IA[if1-1];
        k_b = FE->el_dof->IA[if1]-1;
        for (k=k_a; k<=k_b; k++) {
          mydof = FE->el_dof->JA[k-1];
          if (ix[mydof-1]!=i+1 && FE->dof_flag[mydof-1]==flag) { /* We haven't been here AND it is a boundary */
            A->JA[icp-1] = mydof;
            icp++;
            ix[mydof-1] = i+1;
          }
        }
      }
    }
  }

  if(ix) free(ix);
  icsr_free(&dof_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void LocaltoGlobal(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc) 
{
  /*!
   * \fn LocaltoGlobal(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
   *
   * \brief Maps the local matrix to global matrix NOT considering boundaries
   *
   * \param dof_on_elm    Array of DOF on current element
   * \param FE            FE Space
   * \param ALoc          Local stiffness matrix (full matrix)
   * \param bLoc          Local RHS vector
   *
   * \return A            Global CSR matrix
   * \return b            Global RHS vector
   *
   */

  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE->dof_per_elm; i++) { /* Rows of Local Stiffness */
    row = dof_on_elm[i]-1;
    // Adjust Right-hand side globally
    if(bLoc!=NULL)
      b->val[row] = b->val[row] + bLoc[i];

    for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
      col = dof_on_elm[j]-1;

      col_a = A->IA[row]-1;
      col_b = A->IA[row+1]-1;
      for (k=col_a; k<col_b; k++) { /* Columns of A */
        acol = A->JA[k]-1;
        if (acol==col) {	/* If they match, put it in the global matrix */
          A->val[k] = A->val[k] + ALoc[i*FE->dof_per_elm+j];
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void LocaltoGlobal_FE1FE2(INT *dof_on_elm1,fespace *FE1,INT *dof_on_elm2,fespace *FE2,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc) 
{
  /*!
   * \fn LocaltoGlobal_FE1FE2(INT *dof_on_elm1,fespace *FE1,INT *dof_on_elm2,fespace *FE2,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
   *
   * \brief Maps the local matrix to global matrix NOT considering boundaries
   *
   * \note Assumes non-square system where test and trial are from different FE spaces
   *
   * \param dof_on_elm1   Array of DOF for trial FE space on current element
   * \param FE1           FE Space for trial functions (u)
   * \param dof_on_elm2   Array of DOF for test FE space on current element
   * \param FE2           FE Space for test functions (v)
   * \param ALoc          Local stiffness matrix (full matrix)
   * \param bLoc          Local RHS vector
   *
   * \return A            Global CSR matrix
   * \return b            Global RHS vector
   *
   */

  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE2->dof_per_elm; i++) { /* Rows of Local Stiffness (test space)*/
    row = dof_on_elm2[i]-1;
    // Adjust Right-hand side globally
    if(bLoc!=NULL)
      b->val[row] = b->val[row] + bLoc[i];

    for (j=0; j<FE1->dof_per_elm; j++) { /* Columns of Local Stiffness (trial space)*/
      col = dof_on_elm1[j]-1;

      col_a = A->IA[row]-1;
      col_b = A->IA[row+1]-1;
      for (k=col_a; k<col_b; k++) { /* Columns of A */
        acol = A->JA[k]-1;
        if (acol==col) {	/* If they match, put it in the global matrix */
          A->val[k] = A->val[k] + ALoc[i*FE1->dof_per_elm+j];
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void block_LocaltoGlobal(INT *dof_on_elm,block_fespace *FE,dvector *b,block_dCSRmat *A,REAL *ALoc,REAL *bLoc)
{
  /*!
   * \fn block_LocaltoGlobal(INT *dof_on_elm,block_fespace *FE,dvector *b,block_dCSRmat *A,REAL *ALoc,REAL *bLoc)
   *
   * \brief Maps the local matrix to global block matrix NOT considering boundaries
   *
   * \note Assumes a block structure of FE spaces and matrices
   *
   * \param dof_on_elm    Array of DOF on current element (ordered by block fespace)
   * \param FE            block FE Space
   * \param ALoc          Local stiffness matrix (full matrix)
   * \param bLoc          Local RHS vector
   *
   * \return A            Global block_CSR matrix
   * \return b            Global RHS vector (ordered by block structure of FE space)
   *
   */

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
        local_row = dof_on_elm[local_row_index+i]-1;
        // Adjust Right-hand side globally
        if(bLoc!=NULL && block_col==0)
          b->val[local_row+global_row_index] += bLoc[local_row_index+i];

        /* Columns of Local Stiffness (trial space)*/
        for (j=0; j<dof_per_elm_trial; j++) {
          local_col = dof_on_elm[local_col_index + j]-1;
          /* Columns of A */
          if(A->blocks[block_row*nblocks+block_col]) {
            col_a = A->blocks[block_row*nblocks+block_col]->IA[local_row]-1;
            col_b = A->blocks[block_row*nblocks+block_col]->IA[local_row+1]-1;
            for (k=col_a; k<col_b; k++) {
              acol = A->blocks[block_row*nblocks+block_col]->JA[k]-1;
              if (acol==local_col) {	/* If they match, put it in the global matrix */
                A->blocks[block_row*nblocks+block_col]->val[k] += ALoc[(local_row_index+i)*block_dof_per_elm+(local_col_index+j)];
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
/******************************************************************************************************/

/******************************************************************************************************/
void LocaltoGlobal_withBC(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
{
  /*!
   * \fn LocaltoGlobal_withBC(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
   *
   * \brief Maps the local matrix to global matrix considering boundaries
   *
   * \param dof_on_elm    Array of DOF on current element
   * \param FE            FE Space
   * \param ALoc          Local stiffness matrix (full matrix)
   * \param bLoc          Local RHS vector
   *
   * \return A            Global CSR matrix
   * \return b            Global RHS vector
   *
   */

  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE->dof_per_elm; i++) { /* Rows of Local Stiffness */
    row = dof_on_elm[i]-1;
    if (FE->dirichlet[row]==0) { /* Only if not on a boundary */
      // Adjust Right-hand side globally
      if(bLoc!=NULL)
        b->val[row] = b->val[row] + bLoc[i];

      for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
        col = dof_on_elm[j]-1;
        if (FE->dirichlet[col]==0) { /* Only do stuff if hit a non-boundary edge */
          col_a = A->IA[row]-1;
          col_b = A->IA[row+1]-1;
          for (k=col_a; k<col_b; k++) { /* Columns of A */
            acol = A->JA[k]-1;
            if (acol==col) {	/* If they match, put it in the global matrix */
              A->val[k] = A->val[k] + ALoc[i*FE->dof_per_elm+j];
            }
          }
        } else { /* If boundary adjust Right hand side */
          if(bLoc!=NULL)
            b->val[row] = b->val[row]-(ALoc[i*FE->dof_per_elm+j])*b->val[col];
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag) 
{
  /*!
   * \fn LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag)
   *
   * \brief Maps the local matrix to global matrix considering "special" boundaries
   *        Flag indicates which types of boundaries to consider as Dirichlet
   *        We assume we are considering DOF on a face
   *
   * \param dof_on_f      Array of DOF on current face
   * \param dof_per_f     Number of DOF on current face
   * \param FE            FE Space
   * \param ALoc          Local stiffness matrix (full matrix)
   * \param bLoc          Local RHS vector
   * \param flag          Indicates which boundary DOF to grab
   *
   * \return A            Global CSR matrix
   * \return b            Global RHS vector
   *
   */

  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<dof_per_f; i++) { /* Rows of Local Stiffness */
    row = dof_on_f[i]-1;
    if (FE->dof_flag[row]==flag) { /* Only if on special boundary */
      if(bLoc!=NULL)
        b->val[row] = b->val[row] + bLoc[i];

      for (j=0; j<dof_per_f; j++) { /* Columns of Local Stiffness */
        col = dof_on_f[j]-1;
        if (FE->dof_flag[col]==flag) { /* Only do stuff if hit a special boundary */
          col_a = A->IA[row]-1;
          col_b = A->IA[row+1]-1;
          for (k=col_a; k<col_b; k++) { /* Columns of A */
            acol = A->JA[k]-1;
            if (acol==col) {	/* If they match, put it in the global matrix */
              A->val[k] = A->val[k] + ALoc[i*dof_per_f+j];
            }
          }
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /*!
   * \fn eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global matrix:
   *        For each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and
   *        off-diagonals zero.  Then, make the corresponding column entries 0.
   *        For the RHS, if it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a CSR matrix and a single finite-element space.
   *       If bc=NULL and b=NULL, only eliminate the BC in the matrix.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             CSR stiffness matrix
   * \param time          Physical time if time-dependent
   *
   * \return A            Global CSR matrix with boundaries eliminated
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  INT i,j,cola,colb;
  INT ndof = FE->ndof;
  SHORT status;

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    status = ERROR_MAT_DOF;
    check_error(status, __FUNCTION__);
  }

  // Fix RHS First b_interior = (b - A(0;u_bdry)^T)_interior
  //               b_bdry = u_bdry
  if(b!=NULL)
    eliminate_DirichletBC_RHS(bc,FE,mesh,b,A,time);

  // Now fix A
  // Loop over rows of A
  for (i=0; i<ndof; i++) {
    cola = A->IA[i]-1;
    colb = A->IA[i+1]-1;
    if(FE->dirichlet[i]==1) { // Boundary Row
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) {
        if(A->JA[j]==i+1)
          A->val[j] = 1.0;
        else
          A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) {
        if(FE->dirichlet[A->JA[j]-1]==1) {
          // Zero out matrix entry
          A->val[j] = 0.0;
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /*!
   * \fn eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL),fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global RHS vector:
   *        If it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a CSR matrix and a single finite-element space.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             CSR stiffness matrix (without boundaries eliminated)
   * \param time          Physical time if time-dependent
   *
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  SHORT status;
  INT i;
  INT ndof = FE->ndof;
  REAL* ub = (REAL *) calloc(ndof,sizeof(REAL));

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    status = ERROR_MAT_DOF;
    check_error(status, __FUNCTION__);
  }

  // Get solution vector that's 0 on interior and boundary value on boundary
  for(i=0; i<ndof; i++) {
    if(FE->dirichlet[i]==1) {
      ub[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i);
    } else {
      ub[i] = 0.0;
    }
  }

  // b = b - Aub
  dcsr_aAxpy_1(-1.0,A,ub,b->val);

  // Fix boundary values
  for(i=0;i<ndof;i++) {
    if(FE->dirichlet[i]==1) {
      b->val[i] = ub[i];
    }
  }

  if(ub) free(ub);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void block_eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,
                                 void *A,INT Atype,REAL time)
{
  /*!
   * \fn block_eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,
                                 void *A,INT Atype,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global matrix:
   *        For each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and
   *        off-diagonals zero.  Then, make the corresponding column entries 0.
   *        For the RHS, if it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a BLOCK finite-element space and either block or CSR matrix.
   *       If bc=NULL and b=NULL, only eliminate the BC in the matrix.
   *       Assume all matrices are indexed at 1.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             Stiffness matrix
   * \param Atype         0 - CSR matrix; 1 - block CSR matrix
   * \param time          Physical time if time-dependent
   *
   * \return A            Global CSR matrix with boundaries eliminated
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  if(Atype==0) {
    dCSRmat *Atemp = (dCSRmat *) A;
    eliminate_DirichletBC_blockFE(bc,FE,mesh,b,Atemp,time);
  } else if(Atype==1) {
    block_dCSRmat *Atemp = (block_dCSRmat *) A;
    eliminate_DirichletBC_blockFE_blockA(bc,FE,mesh,b,Atemp,time);
  } else {
    printf("Wrong type of matrix.  Not eliminating anything...\n\n");
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void block_eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,
                                     dvector *b,void *A,INT Atype,REAL time)
{
  /*!
   * \fn block_eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,
                                     dvector *b,void *A,INT Atype,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global RHS vector:
   *        If it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a BLOCK finite-element space and either block or CSR matrix.
   *       Assume all matrices are indexed at 1.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             Stiffness matrix (without boundaries eliminated)
   * \param Atype         0 - CSR matrix; 1 - block CSR matrix
   * \param time          Physical time if time-dependent
   *
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  if(Atype==0) {
    dCSRmat *Atemp = (dCSRmat *) A;
    eliminate_DirichletBC_RHS_blockFE(bc,FE,mesh,b,Atemp,time);
  } else if(Atype==1) {
    block_dCSRmat *Atemp = (block_dCSRmat *) A;
    eliminate_DirichletBC_RHS_blockFE_blockA(bc,FE,mesh,b,Atemp,time);
  } else {
    printf("Wrong type of matrix.  Not eliminating anything...\n\n");
  }
  return;
}
/******************************************************************************************************/


/******************************************************************************************************/
void eliminate_DirichletBC_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /*!
   * \fn eliminate_DirichletBC_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global matrix:
   *        For each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and
   *        off-diagonals zero.  Then, make the corresponding column entries 0.
   *        For the RHS, if it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a CSR matrix and a BLOCK finite-element space.
   *       If bc=NULL and b=NULL, only eliminate the BC in the matrix.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             CSR stiffness matrix
   * \param time          Physical time if time-dependent
   *
   * \return A            Global CSR matrix with boundaries eliminated
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  SHORT status;
  INT i,j,cola,colb;
  INT ndof = FE->ndof;

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    status = ERROR_MAT_DOF;
    check_error(status, __FUNCTION__);
  }
  
  // Fix RHS First b_interior = (b - A(0;u_bdry)^T)_interior
  //               b_bdry = u_bdry
  if(b!=NULL)
    eliminate_DirichletBC_RHS_blockFE(bc,FE,mesh,b,A,time);

  // Now fix A
  // Loop over rows of A
  for (i=0; i<ndof; i++) {
    cola = A->IA[i]-1;
    colb = A->IA[i+1]-1;
    if(FE->dirichlet[i]==1) { // Boundary Row
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) {
        if(A->JA[j]==i+1)
          A->val[j] = 1.0;
        else
          A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) {
        if(FE->dirichlet[A->JA[j]-1]==1) {
          // Zero out matrix entry
          A->val[j] = 0.0;
        }
      }
    }
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_RHS_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time) 
{
  /*!
   * \fn eliminate_DirichletBC_RHS_blockFE(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,dCSRmat *A,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global RHS vector:
   *        If it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a CSR matrix and a BLOCK finite-element space.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             CSR stiffness matrix (without boundaries eliminated)
   * \param time          Physical time if time-dependent
   *
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  SHORT status;
  INT i,j;
  INT ndof = FE->ndof;
  INT ndof_local=0, entry = 0;
  REAL* ub = (REAL *) calloc(ndof,sizeof(REAL));

  // Error Check
  if(A->row!=ndof || A->col!=ndof) {
    status = ERROR_MAT_DOF;
    check_error(status, __FUNCTION__);
  }

  // Get solution vector that's 0 on interior and boundary value on boundary
  for(j=0;j<FE->nspaces;j++) {
    ndof_local = FE->var_spaces[j]->ndof;
    for(i=0; i<ndof_local; i++) {
      if(FE->dirichlet[entry+i]==1) {
        ub[entry + i] = blockFE_Evaluate_DOF(bc,FE,mesh,time,j,i);
      } else {
        ub[entry + i] = 0.0;
      }
    }
    entry += ndof_local;
  }

  // b = b - Aub
  dcsr_aAxpy_1(-1.0,A,ub,b->val);
  
  // Fix boundary values
  for(i=0;i<ndof;i++) {
    if(FE->dirichlet[i]==1) {
      b->val[i] = ub[i];
    }
  }

  if(ub) free(ub);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_blockFE_blockA(void (*bc)(REAL *, REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,
                                          block_dCSRmat *A,REAL time)
{
  /*!
   * \fn eliminate_DirichletBC_blockFE_blockA(void (*bc)(REAL *, REAL *,REAL),void *FE,trimesh *mesh,dvector *b,
                                          block_dCSRmat *A,REAL time)   *
   * \brief Eliminates the Dirichlet boundaries from the global matrix:
   *        For each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and
   *        off-diagonals zero.  Then, make the corresponding column entries 0.
   *        For the RHS, if it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a BLOCK CSR matrix and a BLOCK finite-element space.
   *       If bc=NULL and b=NULL, only eliminate the BC in the matrix.
   *       We also assume blocks of A are indexed at 1.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             block CSR stiffness matrix
   * \param time          Physical time if time-dependent
   *
   * \return A            Global block CSR matrix with boundaries eliminated
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  INT i,j,k,cola,colb;
  INT nsp = FE->nspaces;

  INT nrows;
  INT rowshift, colshift;
  INT l;

  // Fix RHS First b_interior = (b - A(0;u_bdry)^T)_interior
  //               b_bdry = u_bdry
  if(b!=NULL)
    eliminate_DirichletBC_RHS_blockFE_blockA(bc,FE,mesh,b,A,time);

  // Find dof shifts needed for each block
  INT *dofshift;
  dofshift = (INT *)calloc(nsp,sizeof(INT));
  for(i=0;i<nsp;i++){
    for(j=0;j<nsp;j++){
      if(dofshift[j]==0){
        if(A->blocks[i+j*nsp] != NULL){
          dofshift[j] = A->blocks[i+j*nsp]->row;
          //break;
        }
      }
    }
  }
  // Error Check
  for(i=0;i<nsp;i++) {
    if(dofshift[i]==0) {
      printf("ERROR HAZMATH DANGER: in function %s: NULL BLOCK ROW in A.\n",__FUNCTION__);
    }
  }

  colshift = 0;
  // Loop over blocks of A
  for(i=0;i<nsp;i++){ // Loop over block cols
    rowshift = 0;
    for(j=0;j<nsp;j++){ // Loop over block rows
      if(A->blocks[i+j*nsp]!=NULL){
        nrows = A->blocks[i+j*nsp]->row;
        for(k=0;k<nrows;k++){ // Loop over matrix rows
          cola = A->blocks[i+j*nsp]->IA[k]-1;
          colb = A->blocks[i+j*nsp]->IA[k+1]-1;
          if(FE->dirichlet[k+rowshift] == 1) { // Boundary Row
            for(l=cola;l<colb;l++){ // Loop over matrix columns
              if(A->blocks[i+j*nsp]->JA[l]==k+1 && i==j) // diagonal Entry
                A->blocks[i+j*nsp]->val[l] = 1.0;
              else
                A->blocks[i+j*nsp]->val[l] = 0.0;
            }//end for(l)
          } else { // Non-boundary Row
            for(l=cola;l<colb;l++){ // Loop over matrix columns
              if(FE->dirichlet[A->blocks[i+j*nsp]->JA[l]+colshift-1]==1) { // Boundary Column
                // If column is boundary column, set to zero
                A->blocks[i+j*nsp]->val[l] = 0.0;
              }//end if(boundary column)
            }//end for(l)
          }//end Non-boundary Row
        }//end for(k)
      }//end if(A->blocks[i+j*nsp]!=NULL)

      // Update dof shift for rows
      rowshift += dofshift[j];
    }//end for(j)

    // Update dof shift for columns
    colshift += dofshift[i];
  }//end for(i)

  free(dofshift);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void eliminate_DirichletBC_RHS_blockFE_blockA(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,block_dCSRmat *A,REAL time) 
{
  /*!
   * \fn eliminate_DirichletBC_RHS_blockFE_blockA(void (*bc)(REAL *,REAL *,REAL),block_fespace *FE,trimesh *mesh,dvector *b,block_dCSRmat *A,REAL time)
   *
   * \brief Eliminates the Dirichlet boundaries from the global RHS vector:
   *        If it's a boundary DOF, set to the actual boundary value.
   *        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
   *
   * \note We assume a BLOCK CSR matrix and a BLOCK finite-element space.
   *       We also assume blocks of A are indexed at 1.
   *
   * \param bc            Function to get boundary condition at given coordinates.
   * \param FE            block FE Space
   * \param mesh          Mesh struct
   * \param b             RHS vector
   * \param A             block CSR stiffness matrix (without boundaries eliminated)
   * \param time          Physical time if time-dependent
   *
   * \return b            Global RHS vector with boundaries eliminated
   *
   */

  INT i,j;
  INT ndof = FE->ndof;
  INT ndof_local=0, entry = 0;
  REAL* ub = (REAL *) calloc(ndof,sizeof(REAL));

  // Get solution vector that's 0 on interior and boundary value on boundary
  for(j=0;j<FE->nspaces;j++) {
    ndof_local = FE->var_spaces[j]->ndof;
    for(i=0; i<ndof_local; i++) {
      if(FE->dirichlet[entry+i]==1) {
        ub[entry + i] = blockFE_Evaluate_DOF(bc,FE,mesh,time,j,i);
      } else {
        ub[entry + i] = 0.0;
      }
    }
    entry += ndof_local;
  }

  // Shift indices for HAZMATH utilities
  for(i=0;i<(FE->nspaces)*(FE->nspaces);i++) {
    if(A->blocks[i] != NULL)
      dcsr_shift(A->blocks[i],-1);
  }

  // b = b - Aub
  bdcsr_aAxpy(-1.0,A,ub,b->val);
  
  // Shift back
  for(i=0;i<(FE->nspaces)*(FE->nspaces);i++) {
    if(A->blocks[i] != NULL)
      dcsr_shift(A->blocks[i],1);
  }

  // Fix boundary values
  for(i=0;i<ndof;i++) {
    if(FE->dirichlet[i]==1) {
      b->val[i] = ub[i];
    }
  }

  if(ub) free(ub);

  return;
}
/******************************************************************************************************/

