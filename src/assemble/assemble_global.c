/*! \file src/assemble/assemble_global.c
 *
 * \brief This code will build global stiffness matrices for various PDE systems
 *
 * \note It is set up to be generic so that a user could input their own local assembly routines.
 * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
 *
 *  Created by James Adler and Xiaozhe Hu on 4/22/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 * \note modified by James Adler 11/11/2016
 */

#include "hazmat.h"

// Full Assembly Routines
/******************************************************************************************************/
void assemble_global(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  /*!
   * \fn assemble_global(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element
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
   *        A_ij = a( phi_j, phi_i)
   *        b_i  = <f,phi_i>
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param local_assembly Routine to get local matrices
   * \param FE             FE Space
   * \param mesh           Mesh Data
   * \param cq             Quadrature Nodes
   * \param rhs            Routine to get RHS function (NULL if only assembling matrix)
   * \param coeff          Function that gives coefficient (for now assume constant)
   * \param time           Physical Time if time dependent
   *
   * \return A              Global stiffness CSR matrix
   * \return b              Global RHS vector
   *
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  create_CSR_rows(A,FE);

  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  create_CSR_cols(A,FE);

  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
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
  for (i=0; i<FE->nelm; i++) {
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
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,time);
    if(rhs!=NULL)
      FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobal(dof_on_elm,FE,b,A,ALoc,bLoc);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_withBC(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*bc)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
  /*!
   * \fn assemble_global_withBC(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*bc)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element
   *        types (eg. P0, P1, P2, Nedelec, and Raviart-Thomas).
   *        Also takes care of Dirichlet boundary conditions.  If the node is a boundary the row will be
   *        zeroed out except for the diagonal entry being 1.  The corresponding column will also be 0 and
   *        the right-hand side adjusted.
   *
   *        For this problem we compute:
   *
   *        Lu = f  ---->   a(u,v) = <f,v>
   *
   *        which gives Ax = b,
   *
   *        A_ij = a( phi_j, phi_i)
   *        b_i  = <f,phi_i>
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param local_assembly Routine to get local matrices
   * \param FE             FE Space
   * \param mesh           Mesh Data
   * \param cq             Quadrature Nodes
   * \param rhs            Routine to get RHS function (NULL if only assembling matrix)
   * \param bc             Routine to get boundary condition function (NULL if only assembling matrix)
   * \param coeff          Function that gives coefficient (for now assume constant)
   * \param time           Physical Time if time dependent
   *
   * \return A              Global stiffness CSR matrix
   * \return b              Global RHS vector
   *
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  create_CSR_rows_withBC(A,FE);

  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  create_CSR_cols_withBC(A,FE);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }

  // Now Build Global Matrix entries
  // First deal with boundary rows
  for (i=0; i<FE->ndof; i++) {
    if (FE->dof_bdry[i]==1) { /* This is a boundary row.  Just make identity and fix right hand side */
      A->val[A->IA[i]-1] = 1.0;
      if(rhs!=NULL)
        b->val[i] = FE_Evaluate_DOF(bc,FE,mesh,time,i);
    }
  }

  // Now adjust other rows
  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL)
    bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {
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
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,coeff,time);
    if(rhs!=NULL)
      FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobal_withBC(dof_on_elm,FE,b,A,ALoc,bLoc);
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_FE1FE2(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE1, fespace *FE2, trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) 
{
  /*!
   * \fn assemble_global_FE1FE2(dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,fespace *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE1, fespace *FE2, trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global stiffness matrix and rhs for any a(u,v) = <f,v> bilinear form using various element
   *        types (eg. P0, P1, P2, Nedelec, and Raviart-Thomas).
   *        Here we can assume u and v come from different FE spaces
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
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param local_assembly Routine to get local matrices
   * \param FE1 	       Finite-Element Space Struct for trial functions (u)
   * \param FE2            Finite-Element Space Struct for test functions (v)
   * \param mesh           Mesh Data
   * \param cq             Quadrature Nodes
   * \param rhs            Routine to get RHS function (NULL if only assembling matrix)
   * \param coeff          Function that gives coefficient (for now assume constant)
   * \param time           Physical Time if time dependent
   *
   * \return A              Global stiffness CSR matrix
   * \return b              Global RHS vector
   *
   */

  INT dof_per_elm1 = FE1->dof_per_elm;
  INT dof_per_elm2 = FE2->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j;

  // Allocate Row Array
  A->row = FE2->ndof; // test functions
  A->col = FE1->ndof; // trial functions
  A->IA = (INT *) calloc(A->row+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE2->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  create_CSR_rows_FE1FE2(A,FE1,FE2);

  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  create_CSR_cols_FE1FE2(A,FE1,FE2);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }

  // Now Build Global Matrix entries

  /* Loop over all Elements and build local matrix and rhs */
  INT local_size = dof_per_elm2*dof_per_elm1;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL)
    bLoc = (REAL *) calloc(dof_per_elm2,sizeof(REAL));

  INT* dof_on_elm1 = (INT *) calloc(dof_per_elm1,sizeof(INT));
  INT* dof_on_elm2 = (INT *) calloc(dof_per_elm2,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  // Loop over elements
  for (i=0; i<FE1->nelm; i++) {
    // Zero out local matrices
    for (j=0; j<local_size; j++) {
      ALoc[j]=0;
    }
    if(rhs!=NULL) {
      for (j=0; j<dof_per_elm2; j++) {
        bLoc[j]=0;
      }
    }

    // Find DOF for given Element
    rowa = FE1->el_dof->IA[i]-1;
    rowb = FE1->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm1[jcntr] = FE1->el_dof->JA[j];
      jcntr++;
    }
    rowa = FE2->el_dof->IA[i]-1;
    rowb = FE2->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm2[jcntr] = FE2->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE1,FE2,mesh,cq,dof_on_elm1,dof_on_elm2,v_on_elm,i,coeff,time);
    if(rhs!=NULL)
      FEM_RHS_Local(bLoc,FE2,mesh,cq,dof_on_elm2,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    LocaltoGlobal_FE1FE2(dof_on_elm1,FE1,dof_on_elm2,FE2,b,A,ALoc,bLoc);
  }

  if(dof_on_elm1) free(dof_on_elm1);
  if(dof_on_elm2) free(dof_on_elm2);
  if(v_on_elm) free(v_on_elm);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_block(block_dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,REAL),void (*local_rhs_assembly)(REAL *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
  /*!
   * \fn assemble_global_block(block_dCSRmat* A,dvector *b,void (*local_assembly)(REAL *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,REAL),void (*local_rhs_assembly)(REAL *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
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
   * \note All matrices are assumed to be blocks and indexed at 1 in the CSR formatting.
   *
   * \param local_assembly     Routine to get local matrices
   * \param local_rhs_assembly Routine to get local rhs vectors
   * \param FE                 block FE Space
   * \param mesh               Mesh Data
   * \param cq                 Quadrature Nodes
   * \param rhs                Routine to get RHS function (NULL if only assembling matrix)
   * \param bc                 Routine to get boundary condition function (NULL if only assembling matrix)
   * \param coeff              Function that gives coefficient (for now assume constant)
   * \param time               Physical Time if time dependent
   *
   * \return A                 Global stiffness BLOCK CSR matrix
   * \return b                 Global RHS vector
   *
   */

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
      rowa = FE->var_spaces[k]->el_dof->IA[i]-1;
      rowb = FE->var_spaces[k]->el_dof->IA[i+1]-1;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,time);
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

/******************************************************************************************************/
void assemble_global_Jacobian(block_dCSRmat* A,dvector *b,dvector *old_sol,void (*local_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,REAL),void (*local_rhs_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
  /*!
   * \fn assemble_global_Jacobian(block_dCSRmat* A,dvector *b,dvector *old_sol,void (*local_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,REAL),void (*local_rhs_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global stiffness BLOCK matrix and rhs for any a(u,v) = <f,v> bilinear form using various element
   *        types (eg. P0, P1, P2, Nedelec, and Raviart-Thomas).
   *        Here we assume a system and thus a block FE space and that this is from
   *        the assembly of a nonlinear problem (computing the Jacobian).
   *        If it is a linear system, just add NULL for old_sol.
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
   * \note All matrices are assumed to be blocks and indexed at 1 in the CSR formatting.
   *
   * \param old_sol            FE approximation of previous nonlinear solution
   * \param local_assembly     Routine to get local matrices
   * \param local_rhs_assembly Routine to get local rhs vectors
   * \param FE                 block FE Space
   * \param mesh               Mesh Data
   * \param cq                 Quadrature Nodes
   * \param rhs                Routine to get RHS function (NULL if only assembling matrix)
   * \param coeff              Function that gives coefficient (for now assume constant)
   * \param time               Physical Time if time dependent
   *
   * \return A                 Global stiffness BLOCK CSR matrix (Jacobian)
   * \return b                 Global RHS vector (Nonlinear residual)
   *
   */

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
      rowa = FE->var_spaces[k]->el_dof->IA[i]-1;
      rowb = FE->var_spaces[k]->el_dof->IA[i+1]-1;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local Stiffness Matrix for given Element
    (*local_assembly)(ALoc,old_sol,FE,mesh,cq,dof_on_elm,v_on_elm,i,time);
    if(rhs!=NULL)
      (*local_rhs_assembly)(bLoc,old_sol,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

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

// Assembles Global RHS vectors (if needed separately)
/******************************************************************************************************/
void assemble_global_RHS(dvector *b,fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
  /*!
   * \fn assemble_global_RHS(dvector *b,fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global rhs for any a(u,v) = <f,v> bilinear form using various element types
   *        (eg. P1, P2, Nedelec, and Raviart-Thomas).
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   *        For this problem we compute RHS of:
   *
   *        Lu = f  ---->   a(u,v) = <f,v>
   *
   *        which gives Ax = b,
   *
   *        b_i  = <f,phi_i>
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param FE             FE Space
   * \param mesh           Mesh Data
   * \param cq             Quadrature Nodes
   * \param rhs            Routine to get RHS function (NULL if only assembling matrix)
   * \param time           Physical Time if time dependent
   *
   * \return b             Global RHS vector
   *
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,row;

  // Allocate Arrays
  b->row = FE->ndof;
  if(b->val) {
    dvec_set(b->row,b,0.0);
  } else {
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  /* Loop over all Elements and build local rhs */
  REAL* bLoc= (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE->nelm; i++) {
    
    for (j=0; j<dof_per_elm; j++) {
      bLoc[j]=0;
    }

    // Find DOF for given Element
    rowa = FE->el_dof->IA[i]-1;
    rowb = FE->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      dof_on_elm[jcntr] = FE->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local RHS for given Element
    FEM_RHS_Local(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    for (j=0; j<dof_per_elm; j++) {
      row = dof_on_elm[j]-1;
      b->val[row] = b->val[row] + bLoc[j];
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_RHS_block(dvector *b,void (*local_rhs_assembly)(REAL *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
  /*!
   * \fn assemble_global_RHS(dvector *b,fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global rhs for any a(u,v) = <f,v> bilinear form using various element types
   *        (eg. P1, P2, Nedelec, and Raviart-Thomas).
   *        Here we assume a system and thus a block FE space.
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   *        For this problem we compute RHS of:
   *
   *        Lu = f  ---->   a(u,v) = <f,v>
   *
   *        which gives Ax = b,
   *
   *        b_i  = <f,phi_i>
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param local_rhs_assembly  Routine to assemble local RHS
   * \param FE                  block FE Space
   * \param mesh                Mesh Data
   * \param cq                  Quadrature Nodes
   * \param rhs                 Routine to get RHS function (NULL if only assembling matrix)
   * \param time                Physical Time if time dependent
   *
   * \return b             Global RHS vector
   *
   */

  INT dof_per_elm = 0;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,k,row;

  // Get block data first
  INT nblocks = FE->nspaces;
  
  // Allocate arrays
  b->row = FE->ndof;
  if(b->val) {
    dvec_set(b->row,b,0.0);
  } else {
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Loop over each block and get dof_per_elm
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  /* Loop over all Elements and build local rhs */
  REAL* bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  // Loop over elements
  for (i=0; i<mesh->nelm; i++) {
    // Zero out local matrices
    for (j=0; j<dof_per_elm; j++) {
      bLoc[j]=0;
    }

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(k=0;k<nblocks;k++) {
      rowa = FE->var_spaces[k]->el_dof->IA[i]-1;
      rowb = FE->var_spaces[k]->el_dof->IA[i+1]-1;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local RHS for given Element
    (*local_rhs_assembly)(bLoc,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    jcntr = 0;
    rowa = 0;
    for(k=0;k<nblocks;k++) {
      for(j=0;j<FE->var_spaces[k]->dof_per_elm;j++) {
        row = dof_on_elm[jcntr]-1;
        b->val[row+rowa]+=bLoc[jcntr];
      }
      rowa += FE->var_spaces[k]->ndof;
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_RHS_Jacobian(dvector *b,dvector *old_sol,void (*local_rhs_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time) 
{
  /*!
   * \fn assemble_global_RHS_Jacobian(dvector *b,dvector *old_sol,void (*local_rhs_assembly)(REAL *,dvector *,block_fespace *,trimesh *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL),REAL),block_fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),REAL time)
   *
   * \brief Computes the global rhs for any a(u,v) = <f,v> bilinear form using various element types
   *        (eg. P1, P2, Nedelec, and Raviart-Thomas).
   *        Here we assume a system and thus a block FE space and that this is from
   *        the assembly of a nonlinear problem (computing the Jacobian).
   *        If it is a linear system, just add NULL for old_sol.
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   *        For this problem we compute RHS of:
   *
   *        Lu = f  ---->   a(u,v) = <f,v>
   *
   *        which gives Ax = b,
   *
   *        b_i  = <f,phi_i>
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   * \param old_sol             FEM solution at previous Newton Step
   * \param local_rhs_assembly  Routine to assemble local RHS
   * \param FE                  block FE Space
   * \param mesh                Mesh Data
   * \param cq                  Quadrature Nodes
   * \param rhs                 Routine to get RHS function (NULL if only assembling matrix)
   * \param time                Physical Time if time dependent
   *
   * \return b                  Global RHS vector
   *
   */

  INT dof_per_elm = 0;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,k,row;

  // Get block data first
  INT nblocks = FE->nspaces;
  
  // Allocate arrays
  b->row = FE->ndof;
  if(b->val) {
    dvec_set(b->row,b,0.0);
  } else {
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Loop over each block and get dof_per_elm
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  /* Loop over all Elements and build local rhs */
  REAL* bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  // Loop over elements
  for (i=0; i<mesh->nelm; i++) {
    // Zero out local matrices
    for (j=0; j<dof_per_elm; j++) {
      bLoc[j]=0;
    }

    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(k=0;k<nblocks;k++) {
      rowa = FE->var_spaces[k]->el_dof->IA[i]-1;
      rowb = FE->var_spaces[k]->el_dof->IA[i+1]-1;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local RHS for given Element
    (*local_rhs_assembly)(bLoc,old_sol,FE,mesh,cq,dof_on_elm,v_on_elm,i,rhs,time);

    // Loop over DOF and place in appropriate slot globally
    jcntr = 0;
    rowa = 0;
    for(k=0;k<nblocks;k++) {
      for(j=0;j<FE->var_spaces[k]->dof_per_elm;j++) {
        row = dof_on_elm[jcntr]-1;
        b->val[row+rowa]+=bLoc[jcntr];
      }
      rowa += FE->var_spaces[k]->ndof;
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

// Assembly over Faces/Boundaries
/******************************************************************************************************/
void assemble_global_face(dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL),void (*rhs)(REAL *,REAL *, REAL),REAL time,INT flag) 
{
  /*!
   * \fn assemble_global_face(dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,trimesh *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL),void (*rhs)(REAL *,REAL *, REAL),REAL time,INT flag)
   *
   * \brief Computes the global stiffness matrix for any "boundary" bilinear form using various element types
   *        (eg. P1, P2, Nedelec, and Raviart-Thomas).
   *        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
   *        a(u,v)_i, where i denotes a set of faces (or edges) with in a boundary region marked with flag
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *        This allows for several matrices to be assembled then added or concatenated together.
   *
   *        For this problem we compute:
   *
   *        Lu = f  ---->   a(u,v)_bdry = <f,v>_bdry
   *
   *        which gives Ax = b,
   *
   *        A_ij = a( phi_j, phi_i)_bdry
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param old_sol                 FE approximation of previous solution if needed
   * \param local_assembly_face     Routine to get local matrices over each face
   * \param local_rhs_assembly_face Routine to get local rhs vectors over each face
   * \param FE                      FE Space
   * \param mesh                    Mesh Data
   * \param cq                      Quadrature Nodes
   * \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
   * \param coeff                   Function that gives coefficient (for now assume constant)
   * \param time                    Physical Time if time dependent
   * \param flag                    Marker for which faces are included in boundary integration
   *
   * \return A                      Global stiffness CSR matrix
   * \return b                      Global RHS vector
   *
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT dof_per_face = 0;
  INT dim = mesh->dim;
  if(FE->FEtype>=1 && FE->FEtype<10) { // PX Elements 2 -> 2 or 3   3 -> 3 or 6
    dof_per_face = dim + (FE->FEtype-1)*(2*dim -3);
  } else if(FE->FEtype==20) { // Nedelec Elements
    dof_per_face = 2*dim - 3; // 3 edges per face in 3D; face is edge in 2D
  } else if(FE->FEtype==30) { // Raviart-Thomas Elements
    dof_per_elm = 1;
  } else {
    printf("Face integration isn't set up for the FEM space you chose\n");
    exit(0);
  }

  INT i,j,elm;

  // Allocate Row Array
  A->row = FE->ndof;
  A->col = FE->ndof;
  A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
  if(rhs!=NULL) {
    b->row = FE->ndof;
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  // Get Sparsity Structure First
  // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
  create_CSR_rows_flag(A,FE,flag);

  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  create_CSR_cols_flag(A,FE,flag);
  
  // Set values
  A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
  for (i=0; i<A->nnz; i++) {
    A->val[i] = 0;
  }

  // Now Build Global Matrix entries

  /* Loop over all Faces and build local matrix */
  INT local_size = dof_per_face*dof_per_face;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(rhs!=NULL)
    bLoc = (REAL *) calloc(dof_per_face,sizeof(REAL));

  // Get mappings for given element and face
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* dof_on_f = (INT *) calloc(dof_per_face,sizeof(INT));
  INT rowa,rowb,jcntr;

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans_1(mesh->el_f,&f_el);

  // Loop over boundary faces
  for (i=0; i<mesh->nface; i++) {
    // Only grab the faces on the flagged boundary
    if(mesh->f_bdry[i]==flag) {
      // Zero out local matrices
      for (j=0; j<local_size; j++) {
        ALoc[j]=0;
      }
      if(rhs!=NULL) {
        for (j=0; j<dof_per_face; j++) {
          bLoc[j]=0;
        }
      }

      // Find DOF for given Face
      rowa = FE->f_dof->IA[i]-1;
      rowb = FE->f_dof->IA[i+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        dof_on_f[jcntr] = FE->f_dof->JA[j];
        jcntr++;
      }
      if(rowb-rowa!=dof_per_face) {
        printf("Face Integration has error.  Do you know how many DOF's per Face you have??\n");
        exit(0);
      }

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i]-1;
      elm = f_el.JA[rowa]-1;

      // Find DOF on that element
      rowa = FE->el_dof->IA[elm]-1;
      rowb = FE->el_dof->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->el_dof->JA[j];
        jcntr++;
      }
      
      // Find vertices for given Element
      rowa = mesh->el_v->IA[elm]-1;
      rowb = mesh->el_v->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        v_on_elm[jcntr] = mesh->el_v->JA[j];
        jcntr++;
      }

      // Compute Local Stiffness Matrix for given Element
      (*local_assembly_face)(ALoc,old_sol,FE,mesh,cq,dof_on_f,dof_on_elm,v_on_elm,i,elm,coeff,time);
      if(rhs!=NULL)
        (*local_rhs_assembly_face)(bLoc,old_sol,FE,mesh,cq,dof_on_f,dof_on_elm,v_on_elm,i,elm,rhs,time);
      
      // Loop over DOF and place in appropriate slot globally
      LocaltoGlobal_face(dof_on_f,dof_per_face,FE,b,A,ALoc,bLoc,flag);
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(dof_on_f) free(dof_on_f);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
  icsr_free(&f_el);
  
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_RHS_face(dvector* b,dvector *old_sol,void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,trimesh *,INT *,INT,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,void (*rhs)(REAL *,REAL *, REAL),REAL time,INT flag) 
{
  /*!
   * \fn assemble_global_RHS_face(dvector* b,dvector *old_sol,void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,trimesh *,INT *,INT,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL),REAL),fespace *FE,trimesh *mesh,void (*rhs)(REAL *,REAL *, REAL),REAL time,INT flag)
   *
   * \brief Computes the RHS for any "boundary" bilinear form using various element types
   *        (eg. P1, P2, Nedelec, and Raviart-Thomas).
   *        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
   *        a(u,v)_i, where i denotes a set of faces (or edges) with in a boundary region marked with flag
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *        This allows for several matrices to be assembled then added or concatenated together.
   *
   *        For this problem we compute RHS of:
   *
   *        Lu = f  ---->   a(u,v)_bdry = <f,v>_bdry
   *
   *        which gives Ax = b,
   *
   *        A_ij = a( phi_j, phi_i)_bdry
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param old_sol                 FE approximation of previous solution if needed
   * \param local_rhs_assembly_face Routine to get local rhs vectors over each face
   * \param FE                      FE Space
   * \param mesh                    Mesh Data
   * \param cq                      Quadrature Nodes
   * \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
   * \param time                    Physical Time if time dependent
   * \param flag                    Marker for which faces are included in boundary integration
   *
   * \return A                      Global stiffness CSR matrix
   * \return b                      Global RHS vector
   *
   */

  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT dof_per_face = 0;
  INT dim = mesh->dim;
  if(FE->FEtype>=1 && FE->FEtype<10) { // PX Elements 2 -> 2 or 3   3 -> 3 or 6
    dof_per_face = dim + (FE->FEtype-1)*(2*dim -3);
  } else if(FE->FEtype==20) { // Nedelec Elements
    dof_per_face = 2*dim - 3; // 3 edges per face in 3D; face is edge in 2D
  } else if(FE->FEtype==30) { // Raviart-Thomas Elements
    dof_per_face = 1;
  } else {
    printf("Face integration isn't set up for the FEM space you chose\n");
    exit(0);
  }

  INT i,j,elm,row;

  // Allocate Row Array
  b->row = FE->ndof;
  if(b->val) {
    dvec_set(b->row,b,0.0);
  } else {
    b->val = (REAL *) calloc(b->row,sizeof(REAL));
  }

  /* Loop over all Faces and build local matrix */
  REAL* bLoc=NULL;
  bLoc = (REAL *) calloc(dof_per_face,sizeof(REAL));

  // Get mappings for given element and face
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* dof_on_f = (INT *) calloc(dof_per_face,sizeof(INT));
  INT rowa,rowb,jcntr;

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans_1(mesh->el_f,&f_el);

  // Loop over boundary faces
  for (i=0; i<mesh->nface; i++) {
    // Only grab the faces on the flagged boundary
    if(mesh->f_bdry[i]==flag) {
      // Zero out local matrices
      for (j=0; j<dof_per_face; j++) {
        bLoc[j]=0;
      }

      // Find DOF for given Face
      rowa = FE->f_dof->IA[i]-1;
      rowb = FE->f_dof->IA[i+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        dof_on_f[jcntr] = FE->f_dof->JA[j];
        jcntr++;
      }
      if(rowb-rowa!=dof_per_face) {
        printf("Face Integration has error.  Do you know how many DOF's per Face you have??\n");
        exit(0);
      }

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i]-1;
      elm = f_el.JA[rowa]-1;

      // Find DOF on that element
      rowa = FE->el_dof->IA[elm]-1;
      rowb = FE->el_dof->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->el_dof->JA[j];
        jcntr++;
      }
      
      // Find vertices for given Element
      rowa = mesh->el_v->IA[elm]-1;
      rowb = mesh->el_v->IA[elm+1]-1;
      jcntr = 0;
      for (j=rowa; j<rowb; j++) {
        v_on_elm[jcntr] = mesh->el_v->JA[j];
        jcntr++;
      }

      // Compute Local Stiffness Matrix for given Element
      (*local_rhs_assembly_face)(bLoc,old_sol,FE,mesh,dof_on_f,dof_per_face,dof_on_elm,v_on_elm,i,elm,rhs,time);
      
      // Loop over DOF and place in appropriate slot globally
      for (j=0; j<dof_per_face; j++) { /* Rows of Local Stiffness */
        row = dof_on_f[j]-1;
        if (FE->dof_bdry[row]==flag) { /* Only if on special boundary */
          b->val[row] = b->val[row] + bLoc[j];
        }
      }
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(dof_on_f) free(dof_on_f);
  if(bLoc) free(bLoc);
  icsr_free(&f_el);
  
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void assemble_global_Ned_GradH1_RHS(dvector *b,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,dvector* u) 
{
  /*!
   * \fn assemble_global_Ned_GradH1_RHS(dvector *b,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,dvector* u)
   *
   * \brief Computes the global rhs for any <u,grad v>, where u is a Nedelec FE function and v is in P1
   *        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
   *
   * \note All matrices are assumed to be indexed at 1 in the CSR formatting.
   *
   * \param FE_H1                   H1 FE Space for P1
   * \param FE_Ned                  H(curl) FE Space for Nedelec
   * \param mesh                    Mesh Data
   * \param cq                      Quadrature Nodes
   * \param u                       Nedelec function for RHS
   *
   * \return b                      Global RHS vector
   *
   */

  INT ed_per_elm = FE_Ned->dof_per_elm;
  INT v_per_elm = mesh->v_per_elm;
  INT i,j,row;

  // Allocate Arrays
  b->row = FE_H1->ndof;
  b->val = (REAL *) calloc(b->row,sizeof(REAL));

  /* Loop over all Elements and build local rhs */
  REAL* bLoc= (REAL *) calloc(v_per_elm,sizeof(REAL));

  INT* ed_on_elm = (INT *) calloc(ed_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT rowa,rowb,jcntr;
  for (i=0; i<FE_H1->nelm; i++) {
    
    for (j=0; j<v_per_elm; j++) {
      bLoc[j]=0;
    }

    // Find Edges for given Element
    rowa = FE_Ned->el_dof->IA[i]-1;
    rowb = FE_Ned->el_dof->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      ed_on_elm[jcntr] = FE_Ned->el_dof->JA[j];
      jcntr++;
    }

    // Find vertices for given Element
    rowa = mesh->el_v->IA[i]-1;
    rowb = mesh->el_v->IA[i+1]-1;
    jcntr = 0;
    for (j=rowa; j<rowb; j++) {
      v_on_elm[jcntr] = mesh->el_v->JA[j];
      jcntr++;
    }

    // Compute Local RHS for given Element
    Ned_GradH1_RHS_local(bLoc,FE_H1,FE_Ned,mesh,cq,ed_on_elm,v_on_elm,i,u);

    // Loop over DOF and place in appropriate slot globally
    for (j=0; j<v_per_elm; j++) {
      row = v_on_elm[j]-1;
      b->val[row] += bLoc[j];
    }
  }

  if(ed_on_elm) free(ed_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

