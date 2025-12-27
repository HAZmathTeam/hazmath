/*! \file src/assemble/assemble_utils.c
*
* \brief This code will contain all the tools needed to build stiffness matrices
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 4/22/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
* \note modified by James Adler 02/22/2019 for 0-1 fix
*/

#include "hazmath.h"

/******************************************************************************************************/
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
void create_CSR_rows(dCSRmat *A, fespace *FE)
{

  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=0;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Loop over all Elements connected to particular DOF
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el.JA[j];
      k_a = FE->el_dof->IA[if1];
      k_b = FE->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE->el_dof->JA[k];
        if (ix[mydof]!=i) { /* We haven't been here  */
        icp++;
        ix[mydof] = i;
      }
    }
  }
}
A->IA[ndof] = icp;
A->nnz = icp;

if(ix) free(ix);
icsr_free(&dof_el);

return;
}
/******************************************************************************************************/

/******************************************************************************************************/
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
void create_CSR_rows_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2;
  icsr_trans(FE2->el_dof,&dof_el_2);

  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;

  INT* ix = (INT *) calloc(ncols,sizeof(INT));
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
        ix[mydof] = i;
      }
    }
  }
}
A->IA[nrows] = icp;
A->nnz = icp;

if(ix) free(ix);
icsr_free(&dof_el_2);

return;
}
/******************************************************************************************************/

/******************************************************************************************************/
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
void create_CSR_rows_withBC(dCSRmat *A, fespace *FE)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=0;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Check if this is a boundary edge
    if (FE->dirichlet[i]==1) {	/* Only 1 nonzero this row */
      icp++;
    } else {			/* It's not a boundary and compute as normal */
    // Loop over all Elements connected to particular edge
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el.JA[j];
      k_a = FE->el_dof->IA[if1];
      k_b = FE->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE->el_dof->JA[k];
        if (ix[mydof]!=i && FE->dirichlet[mydof]==0) { /* We haven't been here AND it's not a boundary */
        icp++;
        ix[mydof] = i;
      }
    }
  }
}
}
A->IA[ndof] = icp;
A->nnz = icp;

if(ix) free(ix);
icsr_free(&dof_el);

return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn void create_CSR_rows_flag(dCSRmat *A,fespace *FE,INT flag0,INT flag1)
*
* \brief Computes the "possible" number of nonzeros for the global stiffness matrix
*        Also takes into account special boundary conditions.  Here if the boundary
*        of the DOF is equal to "flag" (whatever that corresponds to)
*        then we add something to the matrix.
*        Builds the IA array for A and assigns nnz.
*
* \param FE            FE Space
* \param A             dCSRmat Stiffness Matrix
* \param flag0,flag1   Indicates range of boundaries we care about
*
* \return A.IA         Row structure of CSR matrix
* \return A.nnz        Number of nonzeros of A
*/
void create_CSR_rows_flag(dCSRmat *A, fespace *FE,INT flag0,INT flag1)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and count possible nonzeros in A
  // Also build A->IA, while you're at it...
  icp=0;
  for (i=0; i<ndof; i++) {
    A->IA[i] = icp;
    // Check if this is a boundary dof we want
    if (FE->dof_flag[i]>=flag0 && FE->dof_flag[i]<=flag1) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1];
      for (j=j_a; j<j_b; j++) {
        if1 = dof_el.JA[j];
        k_a = FE->el_dof->IA[if1];
        k_b = FE->el_dof->IA[if1+1];
        for (k=k_a; k<k_b; k++) {
          mydof = FE->el_dof->JA[k];
          if (ix[mydof]!=i && (FE->dof_flag[mydof]>=flag0 && FE->dof_flag[mydof]<=flag1)) { /* We haven't been here AND it is a boundary */
          icp++;
          ix[mydof] = i;
        }
      }
    }
  }
}
A->IA[ndof] = icp;
A->nnz = icp;

if(ix) free(ix);
icsr_free(&dof_el);

return;
}
/******************************************************************************************************/

/******************************************************************************************************/
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
void create_CSR_cols(dCSRmat *A, fespace *FE)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and build A->JA
  icp=0;
  for (i=0; i<ndof; i++) {
    // Loop over all Elements connected to particular edge
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el.JA[j];
      k_a = FE->el_dof->IA[if1];
      k_b = FE->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE->el_dof->JA[k];
        if (ix[mydof]!=i) { /* We haven't been here */
        A->JA[icp] = mydof;
        icp++;
        ix[mydof] = i;
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
void create_CSR_cols_FE1FE2(dCSRmat *A, fespace *FE1, fespace *FE2)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map of the test space
  iCSRmat dof_el_2;
  icsr_trans(FE2->el_dof,&dof_el_2);

  INT nrows = FE2->ndof;
  INT ncols = FE1->ndof;

  INT* ix = (INT *) calloc(ncols,sizeof(INT));
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

return;
}
/******************************************************************************************************/

/******************************************************************************************************/
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
void create_CSR_cols_withBC(dCSRmat *A, fespace *FE)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and build A->JA
  icp=0;
  for (i=0; i<ndof; i++) {
    // Check if this is a boundary edge
    if (FE->dirichlet[i]==1) {	/* Only 1 nonzero this row */
      A->JA[icp]=i;
      icp++;
    } else {			/* It's not a boundary and compute as normal */
    // Loop over all Elements connected to particular edge
    j_a = dof_el.IA[i];
    j_b = dof_el.IA[i+1];
    for (j=j_a; j<j_b; j++) {
      if1 = dof_el.JA[j];
      k_a = FE->el_dof->IA[if1];
      k_b = FE->el_dof->IA[if1+1];
      for (k=k_a; k<k_b; k++) {
        mydof = FE->el_dof->JA[k];
        if (ix[mydof]!=i && FE->dirichlet[mydof]==0) { /* We haven't been here AND it's not a boundary */
        A->JA[icp] = mydof;
        icp++;
        ix[mydof] = i;
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
/*!
* \fn void create_CSR_cols_flag(dCSRmat *A, fespace *FE,INT flag0,INT flag1)
*
* \brief Finds the column sparsity structure of the Global Stiffness Matrix
*        Also takes into account special boundary conditions.  Here if the boundary
*        of the DOF is equal to "flag" (whatever that corresponds to)
*        then we add something to the matrix.
*        Builds the JA array for A.
*
* \param FE            FE Space
* \param A             dCSRmat Stiffness Matrix
* \param flag0,flag1   Indicates which range of boundaries DOF to grab
*
* \return A.JA         Columns of CSR matrix
*
*/
void create_CSR_cols_flag(dCSRmat *A, fespace *FE,INT flag0,INT flag1)
{
  INT i,j,k,j_a,j_b,k_a,k_b,mydof,if1,icp;

  // We will need the DOF to element map first
  iCSRmat dof_el;
  icsr_trans(FE->el_dof,&dof_el);

  INT ndof = FE->ndof;

  INT* ix = (INT *) calloc(ndof,sizeof(INT));
  for (i=0; i<ndof; i++) {
    ix[i] = -1;
  }

  // Loop over all DOF and build A->JA
  icp=0;
  for (i=0; i<ndof; i++) {
    // Check if this is a boundary edge
    if (FE->dof_flag[i]>=flag0 && FE->dof_flag[i]<=flag1) {
      // Loop over all Elements connected to particular edge
      j_a = dof_el.IA[i];
      j_b = dof_el.IA[i+1];
      for (j=j_a; j<j_b; j++) {
        if1 = dof_el.JA[j];
        k_a = FE->el_dof->IA[if1];
        k_b = FE->el_dof->IA[if1+1];
        for (k=k_a; k<k_b; k++) {
          mydof = FE->el_dof->JA[k];
          if (ix[mydof]!=i && (FE->dof_flag[mydof]>=flag0 && FE->dof_flag[mydof]<=flag1)) { /* We haven't been here AND it is a boundary */
          A->JA[icp] = mydof;
          icp++;
          ix[mydof] = i;
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
void LocaltoGlobal(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
{
  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE->dof_per_elm; i++) { /* Rows of Local Stiffness */
    row = dof_on_elm[i];
    // Adjust Right-hand side globally
    if(bLoc!=NULL)
    b->val[row] = b->val[row] + bLoc[i];

    for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
      col = dof_on_elm[j];

      col_a = A->IA[row];
      col_b = A->IA[row+1];
      for (k=col_a; k<col_b; k++) { /* Columns of A */
        acol = A->JA[k];
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
void LocaltoGlobal_FE1FE2(INT *dof_on_elm1,fespace *FE1,INT *dof_on_elm2,fespace *FE2,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
{
  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE2->dof_per_elm; i++) { /* Rows of Local Stiffness (test space)*/
    row = dof_on_elm2[i];
    // Adjust Right-hand side globally
    if(bLoc!=NULL)
    b->val[row] = b->val[row] + bLoc[i];

    for (j=0; j<FE1->dof_per_elm; j++) { /* Columns of Local Stiffness (trial space)*/
      col = dof_on_elm1[j];

      col_a = A->IA[row];
      col_b = A->IA[row+1];
      for (k=col_a; k<col_b; k++) { /* Columns of A */
        acol = A->JA[k];
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
void block_LocaltoGlobal(INT *dof_on_elm,block_fespace *FE,dvector *b,block_dCSRmat *A,REAL *ALoc,REAL *bLoc)
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

        /* Columns of Local Stiffness (trial space)*/
        for (j=0; j<dof_per_elm_trial; j++) {
          local_col = dof_on_elm[local_col_index + j];
          /* Columns of A */
          if(A->blocks[block_row*nblocks+block_col]) {
            col_a = A->blocks[block_row*nblocks+block_col]->IA[local_row];
            col_b = A->blocks[block_row*nblocks+block_col]->IA[local_row+1];
            for (k=col_a; k<col_b; k++) {
              acol = A->blocks[block_row*nblocks+block_col]->JA[k];
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
void LocaltoGlobal_withBC(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc)
{
  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<FE->dof_per_elm; i++) { /* Rows of Local Stiffness */
    row = dof_on_elm[i];
    if (FE->dirichlet[row]==0) { /* Only if not on a boundary */
      // Adjust Right-hand side globally
      if(bLoc!=NULL)
      b->val[row] = b->val[row] + bLoc[i];

      for (j=0; j<FE->dof_per_elm; j++) { /* Columns of Local Stiffness */
        col = dof_on_elm[j];
        if (FE->dirichlet[col]==0) { /* Only do stuff if hit a non-boundary edge */
          col_a = A->IA[row];
          col_b = A->IA[row+1];
          for (k=col_a; k<col_b; k++) { /* Columns of A */
            acol = A->JA[k];
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
/*!
* \fn LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag0,INT flag1)
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
* \param flag0,flag1   Indicates which range of boundaries DOF to grab
*
* \return A            Global CSR matrix
* \return b            Global RHS vector
*
*/
void LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag0,INT flag1)
{
  INT i,j,k,row,col,col_a,col_b,acol;

  for (i=0; i<dof_per_f; i++) { /* Rows of Local Stiffness */
    row = dof_on_f[i];
    if (FE->dof_flag[row]>=flag0 && FE->dof_flag[row]<=flag1) { /* Only if on special boundary */
      if(bLoc!=NULL)
      b->val[row] = b->val[row] + bLoc[i];

      for (j=0; j<dof_per_f; j++) { /* Columns of Local Stiffness */
        col = dof_on_f[j];
        if (FE->dof_flag[col]>=flag0 && FE->dof_flag[col]<=flag1) { /* Only do stuff if hit a special boundary */
          col_a = A->IA[row];
          col_b = A->IA[row+1];
          for (k=col_a; k<col_b; k++) { /* Columns of A */
            acol = A->JA[k];
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
/*!
* \fn block_LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,block_fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag0,INT flag1)
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
* \param flag0,flag1   Indicates which range of boundaries DOF to grab
*
* \return A            Global CSR matrix
* \return b            Global RHS vector
*
*/
void block_LocaltoGlobal_face(INT *dof_on_f,INT dof_per_f,INT* dof_per_face_blk,block_fespace *FE,dvector *b,block_dCSRmat *A,REAL *ALoc,REAL *bLoc,INT flag0,INT flag1)
{
  INT i,j,k,col_a,col_b,acol,block_row,block_col;
  INT local_row,local_col;

  // Loop over all the blocks
  INT nblocks = FE->nspaces;
  INT dof_per_face_test = 0;
  INT dof_per_face_trial = 0;
  INT local_row_index = 0;
  INT local_col_index = 0;
  INT global_row_index = 0;

  // Loop through all the blocks
  for(block_row=0;block_row<nblocks;block_row++) {
    dof_per_face_test = dof_per_face_blk[block_row];

    for(block_col=0;block_col<nblocks;block_col++) {
      dof_per_face_trial = dof_per_face_blk[block_col];

      for (i=0; i<dof_per_face_test; i++) { /* Rows of Local Stiffness */
        local_row = dof_on_f[local_row_index+i];
        // Update RHS
        if(bLoc!=NULL && block_col==0)  b->val[local_row+global_row_index] += bLoc[local_row_index+i];

        for (j=0; j<dof_per_face_trial; j++) { /* Columns of Local Stiffness */
          local_col = dof_on_f[local_col_index+j];
          if(A->blocks[block_row*nblocks+block_col]) {
            col_a = A->blocks[block_row*nblocks+block_col]->IA[local_row];
            col_b = A->blocks[block_row*nblocks+block_col]->IA[local_row+1];
            for (k=col_a; k<col_b; k++) { /* Columns of A */
              acol = A->blocks[block_row*nblocks+block_col]->JA[k];
              if (acol==local_col) {	/* If they match, put it in the global matrix */
                A->blocks[block_row*nblocks+block_col]->val[k] += ALoc[(local_row_index+i)*dof_per_f+(local_col_index+j)];
              }
            }
          }
        }
      }
      local_col_index += dof_per_face_trial;
    }
    local_col_index = 0;
    global_row_index += FE->var_spaces[block_row]->ndof;
    local_row_index += dof_per_face_test;
  }
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
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
void eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
{
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
    cola = A->IA[i];
    colb = A->IA[i+1];
    if(FE->dirichlet[i]==1) { // Boundary Row
      // Loop over columns and 0 out everything in a row except for diagonal
      for(j=cola; j<colb; j++) {
        if(A->JA[j]==i)
        A->val[j] = 1.0;
        else
        A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) {
        if(FE->dirichlet[A->JA[j]]==1) {
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
/*!
* \fn eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
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
void eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL,void *),fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
{
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
  dcsr_aAxpy(-1.0,A,ub,b->val);

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
/*!
* \fn block_eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,
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
*       Assume all matrices are indexed at 0.
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
void block_eliminate_DirichletBC(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,
void *A,INT Atype,REAL time)
{
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
/*!
* \fn block_eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,
dvector *b,void *A,INT Atype,REAL time)
*
* \brief Eliminates the Dirichlet boundaries from the global RHS vector:
*        If it's a boundary DOF, set to the actual boundary value.
*        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
*
* \note We assume a BLOCK finite-element space and either block or CSR matrix.
*       Assume all matrices are indexed at 0.
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
void block_eliminate_DirichletBC_RHS(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,
dvector *b,void *A,INT Atype,REAL time)
{
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
/*!
* \fn eliminate_DirichletBC_blockFE(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
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
void eliminate_DirichletBC_blockFE(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
{
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
    cola = A->IA[i];
    colb = A->IA[i+1];
    if(FE->dirichlet[i]==1) { // Boundary Row
      // Loop over columns and 0 out row except for diagonal
      for(j=cola; j<colb; j++) {
        if(A->JA[j]==i)
        A->val[j] = 1.0;
        else
        A->val[j] = 0.0;
      }
    } else { // Non-boundary-row
      // Loop over columns and 0 out if column is boundary
      for(j=cola; j<colb; j++) {
        if(FE->dirichlet[A->JA[j]]==1) {
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
/*!
* \fn eliminate_DirichletBC_RHS_blockFE(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
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
void eliminate_DirichletBC_RHS_blockFE(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,dCSRmat *A,REAL time)
{
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
  dcsr_aAxpy(-1.0,A,ub,b->val);

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
/*!
* \fn eliminate_DirichletBC_blockFE_blockA(void (*bc)(REAL *, REAL *,REAL,void *),void *FE,mesh_struct *mesh,dvector *b,
block_dCSRmat *A,REAL time)   *
* \brief Eliminates the Dirichlet boundaries from the global matrix:
*        For each row in A that corresponds to a Dirichlet boundary, make the diagonal 1 and
*        off-diagonals zero.  Then, make the corresponding column entries 0.
*        For the RHS, if it's a boundary DOF, set to the actual boundary value.
*        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
*
* \note We assume a BLOCK CSR matrix and a BLOCK finite-element space.
*       If bc=NULL and b=NULL, only eliminate the BC in the matrix.
*       We also assume blocks of A are indexed at 0.
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
void eliminate_DirichletBC_blockFE_blockA(void (*bc)(REAL *, REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,
block_dCSRmat *A,REAL time)
{
  INT i,j,k,cola,colb;
  //INT nsp = FE->nspaces;
  INT nsp = A->brow;

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
          cola = A->blocks[i+j*nsp]->IA[k];
          colb = A->blocks[i+j*nsp]->IA[k+1];
          if(FE->dirichlet[k+rowshift] == 1) { // Boundary Row
            for(l=cola;l<colb;l++){ // Loop over matrix columns
              if(A->blocks[i+j*nsp]->JA[l]==k && i==j) // diagonal Entry
              A->blocks[i+j*nsp]->val[l] = 1.0;
              else
              A->blocks[i+j*nsp]->val[l] = 0.0;
            }//end for(l)
          } else { // Non-boundary Row
            for(l=cola;l<colb;l++){ // Loop over matrix columns
              if(FE->dirichlet[A->blocks[i+j*nsp]->JA[l]+colshift]==1) { // Boundary Column
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
/*!
* \fn eliminate_DirichletBC_RHS_blockFE_blockA(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,block_dCSRmat *A,REAL time)
*
* \brief Eliminates the Dirichlet boundaries from the global RHS vector:
*        If it's a boundary DOF, set to the actual boundary value.
*        If it's a DOF connected to a boundary row, set f(DOF) = f(DOF) - A(DOF,DOF_BC)*u(DOF)
*
* \note We assume a BLOCK CSR matrix and a BLOCK finite-element space.
*       We also assume blocks of A are indexed at 0.
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
void eliminate_DirichletBC_RHS_blockFE_blockA(void (*bc)(REAL *,REAL *,REAL,void *),block_fespace *FE,mesh_struct *mesh,dvector *b,block_dCSRmat *A,REAL time)
{
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

  // b = b - Aub
  bdcsr_aAxpy(-1.0,A,ub,b->val);

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
/*!
* \fn generate_periodic_P(fespace* FE, dCSRmat* P_periodic)
*
* \brief generate (prolongation )matrix P_periodic for applying periodic boundary conditions later
*
* \param FE            FE Space
*
* \return P_periodic   P_periodic matrix in dCSRformat
*
*/
void generate_periodic_P(fespace* FE, dCSRmat* P_periodic)
{

  // local variables
  INT i,j;
  INT ndof = FE->ndof;
  INT nfreedof;
  SHORT stop;

  INT* flag = calloc(FE->ndof, sizeof(INT));

  // set flag array to -1
  iarray_set(ndof, flag, -1);

  // loop 1: initialize flag array
  for (i=0; i<ndof; i++)
  {
    j = FE->periodic[i]; // get the corresponding boundary index

    if (j != -1)
    {
      if ( (flag[i]==-1) & (flag[j]==-1) ) // both boundaries have not been touched
      {
        flag[i] = MIN(i,j);
        flag[j] = MIN(i,j);
      }
      else if ( (flag[i]==-1) & (flag[j]>-1) ) // boundary j has been touched
      {
        flag[i] = flag[j];
      }
      else if ( (flag[i]>-1) & (flag[j]==-1) ) // boudanry i has been touched
      {
        flag[j] = flag[i];
      }
      else // both boundaries have not been touched
      {
        if (flag[i] < flag[j])
        {
          flag[flag[j]] = flag[i];
          flag[j] = flag[i];
        }
        else
        {
          flag[flag[i]] = flag[j];
          flag[i] = flag[j];
        }
      }
    }
  }

  // loop 2: fix flag array
  stop = 0;
  while (stop == 0){
    stop = 1;
    for (i=0; i<ndof; i++)
    {
      if (flag[i]==-1)
      {
        flag[i]=i;
      }
      else if (flag[i] != i)
      {
        if (flag[flag[i]] != flag[i])
        {
          stop = 0;
          flag[i] = flag[flag[i]];
        }
      }
    }
  }

  // count how many dofs after apply periodic boundary condition and form JA for P_periodic
  INT* JA = calloc(ndof, sizeof(INT));
  // loop 3: counting
  nfreedof = 0;
  for (i=0; i<ndof; i++)
  {
    if (flag[i]==i)
    {
      JA[i] = nfreedof;
      nfreedof = nfreedof+1;
    }
  }

  // loop 4: generating JA
  for (i=0; i<ndof; i++)
  {
    if (flag[i] !=i )
    {
      JA[i] = JA[flag[i]];
    }
  }

  //free flag
  free(flag);

  /*
  for (i=0; i<ndof; i++) printf("flag[%d]=%d\n",i, flag[i]);
  for (i=0; i<ndof; i++) printf("JA[%d]=%d\n",i, JA[i]);
  */

  // generating IA
  INT* IA = calloc(ndof+1, sizeof(INT));
  for (i=0; i<ndof+1; i++) IA[i]=i;

  // generating val
  REAL* val = calloc(ndof, sizeof(REAL));
  for (i=0; i<ndof; i++) val[i] = 1.0;

  // form P_periodic
  P_periodic->row = ndof;
  P_periodic->col = nfreedof;
  P_periodic->nnz = ndof;

  P_periodic->IA = IA;
  P_periodic->JA = JA;
  P_periodic->val = val;


}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn generate_periodic_R_scaled(dCSRmat* P_periodic, dCSRmat* R_periodic_scaled)
*
* \brief generate (scaled restriction) matrix R_periodic_scaled (R_periodic_scaled u = u_periodic)
*
* \param P_periodic    P_periodic matrix in dCSRformat
*
* \return R_periodic_scaled   R_periodic_scaled matrix in dCSRformat
*
*/
void generate_periodic_R_scaled(dCSRmat* P_periodic, dCSRmat* R_periodic_scaled)
{

  // local variables
  INT i, j, nnz_row;

  // transpose to get unscaled R_periodic
  dcsr_trans(P_periodic, R_periodic_scaled);

  // scale R_periodic_scaled
  for (i=0; i<R_periodic_scaled->row; i++)
  {
    // number of ones in each row
    nnz_row = R_periodic_scaled->IA[i+1]-R_periodic_scaled->IA[i];

    // scaling each row
    for (j=R_periodic_scaled->IA[i]; j<R_periodic_scaled->IA[i+1]; j++)
    {
      R_periodic_scaled->val[j] = R_periodic_scaled->val[j]/nnz_row;
    }

  }

}

/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn eliminate_PeriodicBC(dCSRmat* P_periodic, dCSRmat* A, dvector* b)
*
* \brief Eliminate periodic boundary conditions. A = P_periodic'*A*P_periodic and b = P_periodic'*b;
*
* \param P_periodic    P_periodic matrix in dCSR format
* \param A             stiffness matrix in dCSR format
* \param b             right hand side
*
* \return A            stiffness matrix after elimination
* \return b            right hand side after elimination
*
*/
void eliminate_PeriodicBC(dCSRmat* P_periodic, dCSRmat* A, dvector* b)
{

  // local variables
  dCSRmat Atemp;
  dvector btemp;

  // copy stiffness matrix
  dcsr_alloc(A->row, A->col, A->nnz, &Atemp);
  dcsr_cp(A, &Atemp);

  // copy right hand side
  if( b ){
    dvec_alloc(b->row, &btemp);
    dvec_cp (b, &btemp);
  }

  // free original A and b
  if (A->IA) free(A->IA);
  if (A->JA) free(A->JA);
  if (A->val) free(A->val);

  // transpose P_periodic
  dCSRmat R_periodic;
  dcsr_trans(P_periodic, &R_periodic);

  // eliminate boundary
  dcsr_rap(&R_periodic, &Atemp, P_periodic, A);

  if(b){
    b->row = R_periodic.row;
    b->val = (REAL* )realloc(b->val, b->row*sizeof(REAL));
    dcsr_mxv(&R_periodic, btemp.val, b->val);
  }

  // free
  dcsr_free(&Atemp);
  dcsr_free(&R_periodic);
  if(b){dvec_free(&btemp);}


}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn eliminate_PeriodicBC_nonoverwrite(dCSRmat* P_periodic, dCSRmat* A, dvector* b, dCSRmat* PTAP, dvector* PTb)
*
* \brief Eliminate periodic boundary conditions. PTAP = P_periodic'*A*P_periodic and PTb = P_periodic'*b;
*
* \param P_periodic    P_periodic matrix in dCSR format
* \param A             stiffness matrix in dCSR format
* \param b             right hand side
*
* \return A            stiffness matrix after elimination
* \return b            right hand side after elimination
*
*/
void eliminate_PeriodicBC_nonoverwrite(dCSRmat* P_periodic, dCSRmat* A, dvector* b, dCSRmat* PTAP, dvector* PTb)
{

  // transpose P_periodic
  dCSRmat R_periodic;
  dcsr_trans(P_periodic, &R_periodic);

  // free PTAP
  if (PTAP->IA)  free(PTAP->IA);
  if (PTAP->JA)  free(PTAP->JA);
  if (PTAP->val) free(PTAP->val);

  // eliminate boundary
  dcsr_rap(&R_periodic, A, P_periodic, PTAP);

  if(PTb){
    PTb->row = R_periodic.row;
    if (PTb->val == NULL) {
      PTb->val = (REAL* )calloc(PTb->row, sizeof(REAL));
    }
    dcsr_mxv(&R_periodic, b->val, PTb->val);
  }

  // free
  dcsr_free(&R_periodic);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn apply_PeriodicBC(dCSRmat* P_periodic, dvector* u)
*
* \brief Apply the boundary conditions so that the vector contains periodic BC dofs: u = P_periodic*u;
*
* \param P_periodic    P_periodic matrix in dCSR format
* \param u             dvector
*
* \return u            dvector after apply the periodic BC
*
*/
void apply_PeriodicBC(dCSRmat* P_periodic, dvector* u)
{
  // local variable
  dvector utemp;

  // copy
  dvec_alloc(u->row, &utemp);
  dvec_cp(u, &utemp);

  // apply peeriodic BC
  u->row = P_periodic->row;
  u->val = (REAL* )realloc(u->val, u->row*sizeof(REAL));
  dcsr_mxv_agg(P_periodic, utemp.val, u->val);

  // free
  dvec_free(&utemp);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn apply_PeriodicBC_nonoverwrite(dCSRmat* P_periodic, dvector* u_periodic, dvector* u)
*
* \brief Apply the boundary conditions so that the vector contains periodic BC dofs: u = P_periodic*u_periodic;
*
* \param P_periodic    P_periodic matrix in dCSR format
* \param u_periodic    dvector
* \param u             dvector
*
* \return u            dvector after apply the periodic BC
*
*/
void apply_PeriodicBC_nonoverwrite(dCSRmat* P_periodic, dvector* u_periodic, dvector* u)
{

  // apply peeriodic BC
  u->row = P_periodic->row;
  if (u->val == NULL){
    u->val = (REAL* )calloc(u->row, sizeof(REAL));
  }
  dcsr_mxv_agg(P_periodic, u_periodic->val, u->val);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn generate_periodic_P_blockFE(block_fespace* FE, block_dCSRmat* P_periodic)
*
* \brief generate (prolongation) matrix P_periodic for applying periodic boundary conditions later
*
* \param FE            block FE Space
*
* \return P_periodic   P_periodic matrix in dCSRformat
*
*/
void generate_periodic_P_blockFE(block_fespace* FE, block_dCSRmat* P_periodic)
{
  // local variables
  INT i, nspaces=FE->nspaces;

  // allocate block dCSRmat
  bdcsr_alloc_minimal(nspaces, nspaces, P_periodic);

  // allocate diagonal blocks (other blocks are set to be NULL)
  for (i=0; i<(nspaces*nspaces); i++) P_periodic->blocks[i] = NULL;
  for (i=0; i<nspaces; i++) P_periodic->blocks[i*nspaces+i] = (dCSRmat *)calloc(1,sizeof(dCSRmat));

  // generate P_periodic
  for (i=0; i<nspaces; i++) generate_periodic_P(FE->var_spaces[i], P_periodic->blocks[i*nspaces+i]);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn generate_periodic_R_scaled_blockFE(dCSRmat* P_periodic, dCSRmat* R_periodic_scaled)
*
* \brief generate (scaled restriction) matrix R_periodic_scaled (R_periodic_scaled u = u_periodic)
*
* \param P_periodic    P_periodic matrix in dCSRformat
*
* \return R_periodic_scaled   R_periodic_scaled matrix in dCSRformat
*
*/
void generate_periodic_R_scaled_blockFE(block_dCSRmat* P_periodic, block_dCSRmat* R_periodic_scaled)
{
  // local variables
  INT i, nspaces = P_periodic->brow;

  // allocate block dCSRmat
  bdcsr_alloc_minimal(nspaces, nspaces, R_periodic_scaled);

  // allocate diagonal blocks (other blocks are set to be NULL)
  for (i=0; i<(nspaces*nspaces); i++) R_periodic_scaled->blocks[i] = NULL;
  for (i=0; i<nspaces; i++) R_periodic_scaled->blocks[i*nspaces+i] = (dCSRmat *)calloc(1,sizeof(dCSRmat));

  // generate P_periodic
  for (i=0; i<nspaces; i++) generate_periodic_R_scaled(P_periodic->blocks[i*nspaces+i], R_periodic_scaled->blocks[i*nspaces+i]);

}

/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn eliminate_PeriodicBC_blockFE(block_dCSRmat* P_periodic, block_dCSRmat* A, dvector* b)
*
* \brief Eliminate periodic boundary conditions (block version). A = P_periodic'*A*P_periodic and b = P_periodic'*b;
*
* \param P_periodic    P_periodic matrix in block_dCSR format
* \param A             stiffness matrix in block_dCSR format
* \param b             right hand side
*
* \return A            stiffness matrix after elimination
* \return b            right hand side after elimination
*
*/
void eliminate_PeriodicBC_blockFE(block_dCSRmat* P_periodic, block_dCSRmat* A, dvector* b)
{

  // local variables
  INT i, j;
  block_dCSRmat Atemp;
  dCSRmat RAtemp;
  dvector btemp;

  // copy stiffness matrix
  bdcsr_alloc(A->brow, A->bcol, &Atemp);
  bdcsr_cp(A, &Atemp);

  // copy right hand side
  if( b ){
    dvec_alloc(b->row, &btemp);
    dvec_cp (b, &btemp);
  }

  // free original A and b
  for (i=0; i<(A->brow*A->bcol); i++)
  {
    if (A->blocks[i]){
      if (A->blocks[i]->IA) free(A->blocks[i]->IA);
      if (A->blocks[i]->JA) free(A->blocks[i]->JA);
      if (A->blocks[i]->val) free(A->blocks[i]->val);
    }
  }

  // transpose P_periodic
  block_dCSRmat R_periodic;
  bdcsr_trans(P_periodic, &R_periodic);

  // eliminate boundary (use the fact that both P and R are block diagonal)
  for (i=0; i<A->brow; i++)
  {

    for (j=0; j<A->bcol; j++)
    {

      if ( Atemp.blocks[i*A->brow+j] == NULL  )
      {
        A->blocks[i*A->brow+j] = NULL;
      }
      else
      {
        dcsr_mxm(R_periodic.blocks[i*R_periodic.brow+i], Atemp.blocks[i*A->brow+j], &RAtemp);
        dcsr_mxm(&RAtemp, P_periodic->blocks[j*P_periodic->brow+j], A->blocks[i*A->brow+j]);
        dcsr_free(&RAtemp);
      }

    }

  }

  if( b ){
    b->row = 0;
    for (i=0; i<R_periodic.brow; i++) b->row = b->row + R_periodic.blocks[i*R_periodic.brow+i]->row;
    b->val = (REAL* )realloc(b->val, b->row*sizeof(REAL));
    bdcsr_mxv(&R_periodic, btemp.val, b->val);
  }

  // free
  bdcsr_free(&Atemp);
  bdcsr_free(&R_periodic);
  if(b){dvec_free(&btemp);}

}
/******************************************************************************************************/

/*!
* \fn eliminate_PeriodicBC_blockFE_nonoverwrite(block_dCSRmat* P_periodic, block_dCSRmat* A, dvector* b)
*
* \brief Eliminate periodic boundary conditions (block version). A = P_periodic'*A*P_periodic and b = P_periodic'*b;
*
* \param P_periodic    P_periodic matrix in block_dCSR format
* \param A             stiffness matrix in block_dCSR format
* \param b             right hand side
*
* \return A            stiffness matrix after elimination
* \return b            right hand side after elimination
*
*/
void eliminate_PeriodicBC_blockFE_nonoverwrite(block_dCSRmat* P_periodic, block_dCSRmat* A, dvector* b, block_dCSRmat* PTAP, dvector* PTb)
{

  // local variables
  INT i, j;
  dCSRmat RAtemp;

  // free PTAP
  for (i=0; i<(A->brow*A->bcol); i++)
  {
    if (PTAP->blocks[i]){
      if (PTAP->blocks[i]->IA)  free(PTAP->blocks[i]->IA);
      if (PTAP->blocks[i]->JA)  free(PTAP->blocks[i]->JA);
      if (PTAP->blocks[i]->val) free(PTAP->blocks[i]->val);
    }
  }

  // transpose P_periodic
  block_dCSRmat R_periodic;
  bdcsr_trans(P_periodic, &R_periodic);

  // eliminate boundary (use the fact that both P and R are block diagonal)
  for (i=0; i<A->brow; i++)
  {

    for (j=0; j<A->bcol; j++)
    {

      if ( A->blocks[i*A->brow+j] == NULL  )
      {
        PTAP->blocks[i*A->brow+j] = NULL;
      }
      else
      {
        dcsr_mxm(R_periodic.blocks[i*R_periodic.brow+i], A->blocks[i*A->brow+j], &RAtemp);
        dcsr_mxm(&RAtemp, P_periodic->blocks[j*P_periodic->brow+j], PTAP->blocks[i*A->brow+j]);
        dcsr_free(&RAtemp);
      }

    }

  }

  if( PTb ){
    PTb->row = 0;
    for (i=0; i<R_periodic.brow; i++) PTb->row = PTb->row + R_periodic.blocks[i*R_periodic.brow+i]->row;
    if (PTb->val == NULL) PTb->val = (REAL* )calloc(PTb->row, sizeof(REAL));
    bdcsr_mxv(&R_periodic, b->val, PTb->val);
  }

  // free
  bdcsr_free(&R_periodic);

}
/******************************************************************************************************/

/*!
* \fn eliminate_PeriodicBC_blockFE_RHS_nonoverwrite(block_dCSRmat* P_periodic,dvector* b)
*
* \brief Eliminate periodic boundary conditions (block version) for RHS ONLY
*        b = P_periodic'*b;
*
* \param P_periodic    P_periodic matrix in block_dCSR format
* \param b             right hand side
*
* \return b            right hand side after elimination
*
*/
void eliminate_PeriodicBC_blockFE_RHS_nonoverwrite(block_dCSRmat* P_periodic, dvector* b, dvector* PTb)
{

  // local variables
  INT i;

  // transpose P_periodic
  block_dCSRmat R_periodic;
  bdcsr_trans(P_periodic, &R_periodic);

  PTb->row = 0;
  for (i=0; i<R_periodic.brow; i++) PTb->row = PTb->row + R_periodic.blocks[i*R_periodic.brow+i]->row;
  if (PTb->val == NULL) PTb->val = (REAL* )calloc(PTb->row, sizeof(REAL));
  bdcsr_mxv(&R_periodic, b->val, PTb->val);

  // free
  bdcsr_free(&R_periodic);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn apply_PeriodicBC_blockFE(block_dCSRmat* P_periodic, dvector* u)
*
* \brief Apply the boundary conditions so that the vector contains periodic BC dofs: u = P_periodic*u (block version)
*
* \param P_periodic    P_periodic matrix in block_dCSR format
* \param u             dvector
*
* \return u            dvector after apply the periodic BC
*
*/
void apply_PeriodicBC_blockFE(block_dCSRmat* P_periodic, dvector* u)
{
  // local variable
  INT i;
  dvector utemp;

  // copy
  dvec_alloc(u->row, &utemp);
  dvec_cp(u, &utemp);

  // apply peeriodic BC
  u->row = 0;
  for (i=0; i<P_periodic->brow;i++) u->row = u->row + P_periodic->blocks[i*P_periodic->brow+i]->row;
  u->val = (REAL* )realloc(u->val, u->row*sizeof(REAL));
  bdcsr_mxv(P_periodic, utemp.val, u->val);

  // free
  dvec_free(&utemp);

}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn apply_PeriodicBC_blockFE_nonoverwrite(block_dCSRmat* P_periodic, dvector* u_periodic, dvector* u)
*
* \brief Apply the boundary conditions so that the vector contains periodic BC dofs: u = P_periodic*u_periodic (block version)
*
* \param P_periodic    P_periodic matrix in block_dCSR format
* \param u_periodic    dvector
* \param u             dvector
*
* \return u            dvector after apply the periodic BC
*
*/
void apply_PeriodicBC_blockFE_nonoverwrite(block_dCSRmat* P_periodic, dvector* u_periodic, dvector* u)
{

  // local variable
  INT i;

  // apply peeriodic BC
  u->row = 0;
  for (i=0; i<P_periodic->brow;i++) u->row = u->row + P_periodic->blocks[i*P_periodic->brow+i]->row;
  if (u->val == NULL)  u->val = (REAL* )calloc(u->row, sizeof(REAL));
  bdcsr_mxv(P_periodic, u_periodic->val, u->val);

}
/******************************************************************************************************/
