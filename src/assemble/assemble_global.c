/*! \file src/assemble/assemble_global.c
*
* \brief This code will build global stiffness matrices for various PDE systems
*
* \note It is set up to be generic so that a user could input their own local assembly routines.
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 4/22/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
* \note modified by James Adler 02/21/2019 for 0-1 fix
*
* \note modified by James Adler 03/25/2021 to make assembly faster.  Will now only
*       pass local data to local assembly and not entire mesh.  We will leave old
*       assemblies here for now
*/

#include "hazmath.h"

/**********************************************************************/
/*!
 * \brief Built-in local assembly (new signature): DuDv (stiffness) + RHS for single-space.
 */
void local_assembly_DuDv(REAL *ALoc, REAL *bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *), REAL time)
{
  if(ALoc!=NULL)
    assemble_DuDv_local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
  if(bLoc!=NULL && rhs!=NULL)
    FEM_RHS_Local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
}
/**********************************************************************/
/*!
 * \brief Built-in local assembly (new signature): mass matrix + RHS for single-space.
 */
void local_assembly_mass(REAL *ALoc, REAL *bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *), REAL time)
{
  if(ALoc!=NULL)
    assemble_mass_local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
  if(bLoc!=NULL && rhs!=NULL)
    FEM_RHS_Local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
}
/**********************************************************************/
/*!
 * \brief Built-in local assembly (new signature): RHS-only for single-space.
 */
void local_assembly_rhs_only(REAL *ALoc, REAL *bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *), REAL time)
{
  if(bLoc!=NULL && rhs!=NULL)
    FEM_RHS_Local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
}

// New Global Assembly Routines (using local data structs)
/******************************************************************************************************/
/*!
* \fn void assemble_global_single(...)
*
* \brief Unified global assembly for single FE space problems using local data structs.
*        Computes the global stiffness matrix and/or rhs.
*        Uses simplex_local_data and fe_local_data instead of passing global structs
*        to the callback.
*
*        Does NOT apply Dirichlet boundary conditions. Use
*        eliminate_DirichletBC() after assembly.
*
* \param A                 Output CSR matrix (NULL if RHS-only)
* \param b                 Output RHS vector (NULL if matrix-only)
* \param FE                FE space
* \param sc                Simplicial complex (mesh)
* \param cq                Quadrature data
* \param local_assembly    Callback using local data structs
* \param old_sol           Previous solution for Newton (NULL for linear)
* \param rhs               RHS function (NULL if no source term)
* \param coeff             Coefficient function (NULL if not needed)
* \param time              Physical time
*
*/
void assemble_global_single(dCSRmat* A,dvector *b,fespace *FE,scomplex *sc,qcoordinates *cq,void (*local_assembly)(REAL *,REAL *,REAL *,simplex_local_data *,fe_local_data *,void (*)(REAL *,REAL *,REAL,void *),void (*)(REAL *,REAL *,REAL,void *),REAL),dvector *old_sol,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dim = sc->dim;
  INT dof_per_elm = FE->dof_per_elm;
  INT i,j;

  // Build sparsity structure (if assembling matrix)
  if(A!=NULL) {
    A->row = FE->ndof;
    A->col = FE->ndof;
    if(A->IA==NULL) {
      A->IA = (INT *) calloc(FE->ndof+1,sizeof(INT));
      create_CSR_rows(A,FE);
      A->JA = (INT *) calloc(A->nnz,sizeof(INT));
      create_CSR_cols(A,FE);
    }
    if(A->val==NULL)
      A->val = (REAL *) calloc(A->nnz,sizeof(REAL));
    for (i=0; i<A->nnz; i++) A->val[i] = 0;
  }

  // Allocate or zero RHS vector
  if(b!=NULL) {
    b->row = FE->ndof;
    if(b->val) {
      dvec_set(b->row,b,0.0);
    } else {
      b->val = (REAL *) calloc(b->row,sizeof(REAL));
    }
  }

  // Allocate local arrays
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = NULL;
  REAL* bLoc = NULL;
  if(A!=NULL) ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  if(b!=NULL) bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  // Create temporary 1-block block_fespace for initialize_localdata_elm
  block_fespace temp_bfe;
  temp_bfe.nspaces = 1;
  temp_bfe.nun = 1;
  temp_bfe.ndof = FE->ndof;
  temp_bfe.nelm = FE->nelm;
  temp_bfe.var_spaces = (fespace **) calloc(1,sizeof(fespace *));
  temp_bfe.var_spaces[0] = FE;
  temp_bfe.dirichlet = NULL;
  temp_bfe.dof_flag = NULL;
  temp_bfe.simplex_data = NULL;
  temp_bfe.fe_data = NULL;

  // Initialize local data structs
  simplex_local_data elm_data;
  fe_local_data fe_data;
  memset(&elm_data, 0, sizeof(simplex_local_data));
  memset(&fe_data, 0, sizeof(fe_local_data));
  initialize_localdata_elm(&elm_data, &fe_data, sc, &temp_bfe, cq->nq1d);

  // Loop over elements
  for (i=0; i<FE->nelm; i++) {
    if(ALoc!=NULL) memset(ALoc, 0, local_size*sizeof(REAL));
    if(bLoc!=NULL) memset(bLoc, 0, dof_per_elm*sizeof(REAL));

    // Gather local mesh and FE data
    get_elmlocaldata(&elm_data, sc, i);
    get_felocaldata_elm(&fe_data, &temp_bfe, old_sol, i);

    // Compute local matrix and/or RHS
    (*local_assembly)(ALoc, bLoc, fe_data.u_local, &elm_data, &fe_data, rhs, coeff, time);

    // Scatter local to global
    if(A!=NULL)
      LocaltoGlobal(fe_data.local_dof, FE, b, A, ALoc, bLoc);
    else if(bLoc!=NULL && b!=NULL) {
      for(j=0;j<dof_per_elm;j++)
        b->val[fe_data.local_dof[j]] += bLoc[j];
    }
  }

  // Free local data
  free_simplexlocaldata(&elm_data);
  free_felocaldata(&fe_data);
  if(temp_bfe.var_spaces) free(temp_bfe.var_spaces);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn void assemble_global_system(...)
*
* \brief Unified global assembly for block FE space problems using local data structs.
*        Computes the global stiffness BLOCK matrix and/or rhs.
*        Uses simplex_local_data and fe_local_data instead of passing global structs
*        to the callback.
*
*        Does NOT apply Dirichlet boundary conditions. Use
*        eliminate_DirichletBC() after assembly.
*
* \param A                 Output block CSR matrix (NULL if RHS-only)
* \param b                 Output RHS vector (NULL if matrix-only)
* \param FE                Block FE space
* \param sc                Simplicial complex (mesh)
* \param cq                Quadrature data
* \param local_assembly    Callback using local data structs
* \param old_sol           Previous solution for Newton (NULL for linear)
* \param rhs               RHS function (NULL if no source term)
* \param coeff             Coefficient function (NULL if not needed)
* \param time              Physical time
*
*/
void assemble_global_system(block_dCSRmat* A,dvector *b,block_fespace *FE,scomplex *sc,qcoordinates *cq,void (*local_assembly)(REAL *,REAL *,REAL *,simplex_local_data *,fe_local_data *,void (*)(REAL *,REAL *,REAL,void *),void (*)(REAL *,REAL *,REAL,void *),REAL),dvector *old_sol,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT dof_per_elm = 0;
  INT i,j,k,testdof,trialdof;

  INT nblocks = FE->nspaces;

  // Build sparsity structure of block matrix (if assembling matrix)
  if(A!=NULL) {
    if(nblocks!=A->brow || nblocks!=A->bcol) {
      fprintf(stderr,"HAZMATH ERROR in %s: block matrix is %lldx%lld but FE has %lld spaces.\n",
              __FUNCTION__,(long long)A->brow,(long long)A->bcol,(long long)nblocks);
      exit(255);
    }
    for(i=0;i<nblocks;i++) {
      for(j=0;j<nblocks;j++) {
        if(A->blocks[i*nblocks+j]) {
          testdof = FE->var_spaces[i]->ndof;
          trialdof = FE->var_spaces[j]->ndof;
          if(A->blocks[i*nblocks+j]->IA==NULL) {
            A->blocks[i*nblocks+j]->row = testdof;
            A->blocks[i*nblocks+j]->col = trialdof;
            A->blocks[i*nblocks+j]->IA = (INT *) calloc(testdof+1,sizeof(INT));
            create_CSR_rows_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);
            A->blocks[i*nblocks+j]->JA = (INT *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(INT));
            create_CSR_cols_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);
          }
          if(A->blocks[i*nblocks+j]->val==NULL)
            A->blocks[i*nblocks+j]->val = (REAL *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(REAL));
          for (k=0; k<A->blocks[i*nblocks+j]->nnz; k++) {
            A->blocks[i*nblocks+j]->val[k] = 0;
          }
        }
      }
    }
  }

  // Allocate or zero RHS vector
  if(b!=NULL) {
    b->row = FE->ndof;
    if(b->val) {
      dvec_set(b->row,b,0.0);
    } else {
      b->val = (REAL *) calloc(b->row,sizeof(REAL));
    }
  }

  // Compute total dof_per_elm
  for(i=0;i<nblocks;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  // Allocate local arrays
  INT local_size = dof_per_elm*dof_per_elm;
  REAL* ALoc = NULL;
  REAL* bLoc = NULL;
  if(A!=NULL) ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  if(b!=NULL) bLoc = (REAL *) calloc(dof_per_elm,sizeof(REAL));

  // Initialize local data structs
  simplex_local_data elm_data;
  fe_local_data fe_data;
  memset(&elm_data, 0, sizeof(simplex_local_data));
  memset(&fe_data, 0, sizeof(fe_local_data));
  initialize_localdata_elm(&elm_data, &fe_data, sc, FE, cq->nq1d);

  // Loop over elements
  for (i=0; i<fem->ns_leaf; i++) {
    // Zero out local arrays
    if(ALoc!=NULL) memset(ALoc, 0, local_size*sizeof(REAL));
    if(bLoc!=NULL) memset(bLoc, 0, dof_per_elm*sizeof(REAL));

    // Gather local mesh and FE data
    get_elmlocaldata(&elm_data, sc, i);
    get_felocaldata_elm(&fe_data, FE, old_sol, i);

    // Compute local matrix and/or RHS
    (*local_assembly)(ALoc, bLoc, fe_data.u_local, &elm_data, &fe_data, rhs, coeff, time);

    // Scatter local to global
    block_LocaltoGlobal(fe_data.local_dof, FE, b, A, ALoc, bLoc);
  }

  free_simplexlocaldata(&elm_data);
  free_felocaldata(&fe_data);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/

// Assembly over Faces/Boundaries
/******************************************************************************************************/
/*!
* \fn assemble_global_face(dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),fespace *FE,scomplex *sc,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL,void *),void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
*
* \brief Computes the global stiffness matrix for any "boundary" bilinear form using various element types
*        (eg. P1, P2, Nedelec, and Raviart-Thomas).
*        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
*        a(u,v)_i, where i denotes a set of faces (or edges) within a boundary region marked with flag
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
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.
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
* \param flag0,flag1             Marker for which faces are included in boundary integration (range of faces from flag0 to flag1)
*
* \return A                      Global stiffness CSR matrix
* \return b                      Global RHS vector
*
*/
void assemble_global_face(dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),fespace *FE,scomplex *sc,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL,void *),void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = (dim + 1);
  INT dof_per_face = 0;
  if(FE->FEtype>=1 && FE->FEtype<10) { // PX Elements 2 -> 2 or 3   3 -> 3 or 6
    dof_per_face = dim + (FE->FEtype-1)*(2*dim -3);
  } else if(FE->FEtype==20) { // Nedelec Elements
    dof_per_face = 2*dim - 3; // 3 edges per face in 3D; face is edge in 2D
  } else if(FE->FEtype==30) { // Raviart-Thomas Elements
    dof_per_elm = 1;
  } else if(FE->FEtype==103) { // Mini
    dof_per_elm = dim; // Just linear part
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
  create_CSR_rows_flag(A,FE,flag0,flag1);

  // Columns of A -> JA
  A->JA = (INT *) calloc(A->nnz,sizeof(INT));
  create_CSR_cols_flag(A,FE,flag0,flag1);

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
  INT rowa;

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans(fem->el_f,&f_el);

  // Loop over boundary faces
  for (i=0; i<fem->nface; i++) {
    // Only grab the faces on the flagged boundary
    if(fem->f_flag[i]>=flag0 && fem->f_flag[i]<=flag1) {
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
      get_incidence_row(i,FE->f_dof,dof_on_f);

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i];
      elm = f_el.JA[rowa];

      // Find DOF on that element
      get_incidence_row(elm,FE->el_dof,dof_on_elm);

      // Find vertices for given Element
      get_incidence_row(elm,fem->el_v,v_on_elm);

      // Compute Local Stiffness Matrix for given Element
      (*local_assembly_face)(ALoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,i,elm,coeff,time);
      if(rhs!=NULL)
      (*local_rhs_assembly_face)(bLoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,i,elm,rhs,time);

      // Loop over DOF and place in appropriate slot globally
      LocaltoGlobal_face(dof_on_f,dof_per_face,FE,b,A,ALoc,bLoc,flag0,flag1);
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
/*!
* \fn void assemble_global_RHS_face(dvector* b,dvector *old_sol,void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),fespace *FE,scomplex *sc,qcoordinates *cq,void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
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
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.
* \note Assumes different type of integral for different Element type:
*       PX -> <f,v>_bdry
*       Ned -> <f,nxv>_bdry
*       RT  -> <f,n*v>_bdry
*
* \param old_sol                 FE approximation of previous solution if needed
* \param local_rhs_assembly_face Routine to get local rhs vectors over each face
* \param FE                      FE Space
* \param mesh                    Mesh Data
* \param cq                      Quadrature Nodes
* \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
* \param time                    Physical Time if time dependent
* \param flag0,flag1             Marker for which faces are included in boundary integration (range of faces from flag0 to flag1)
*
* \return b                      Global RHS vector
*
*/
void assemble_global_RHS_face(dvector* b,dvector *old_sol,void (*local_rhs_assembly_face)(REAL *,dvector *,fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),fespace *FE,scomplex *sc,qcoordinates *cq,void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT dof_per_elm = FE->dof_per_elm;
  INT v_per_elm = (dim + 1);
  INT dof_per_face = 0;
  if(FE->FEtype>=1 && FE->FEtype<10) { // PX Elements 2 -> 2 or 3   3 -> 3 or 6
    dof_per_face = dim + (FE->FEtype-1)*(2*dim -3);
  } else if(FE->FEtype==20) { // Nedelec Elements
    dof_per_face = 2*dim - 3; // 3 edges per face in 3D; face is edge in 2D
  } else if(FE->FEtype==30) { // Raviart-Thomas Elements
    dof_per_face = 1;
  } else if(FE->FEtype==61) { // Bubbles
    dof_per_face = 1;
  } else if(FE->FEtype==103) { // MINI
    dof_per_face = dim; // only linears contribute to face
  } else {
    printf("Face integration isn't set up for the FEM space you chose\n");
    exit(0);
  }

  INT i,j,elm,row,rowa;

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

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans(fem->el_f,&f_el);

  // Loop over boundary faces
  for (i=0; i<fem->nface; i++) {
    // Only grab the faces on the flagged boundary
    if(fem->f_flag[i]>=flag0 && fem->f_flag[i]<=flag1) {
      // Zero out local matrices
      for (j=0; j<dof_per_face; j++) {
        bLoc[j]=0;
      }

      // Find DOF for given Face
      get_incidence_row(i,FE->f_dof,dof_on_f);

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i];
      elm = f_el.JA[rowa];

      // Find DOF on that element
      get_incidence_row(elm,FE->el_dof,dof_on_elm);

      // Find vertices for given Element
      get_incidence_row(elm,fem->el_v,v_on_elm);

      // Compute Local Stiffness Matrix for given Element
      (*local_rhs_assembly_face)(bLoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,dof_per_face,i,elm,rhs,time);

      // Loop over DOF and place in appropriate slot globally
      for (j=0; j<dof_per_face; j++) { /* Rows of Local Stiffness */
        row = dof_on_f[j];
        b->val[row] += bLoc[j];
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

/*!
* \fn assemble_global_face_block(block_dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,block_fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,block_fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,scomplex *sc,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL,void *),void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
*
* \brief Computes the global block stiffness matrix for any "boundary" bilinear form using various element types
*        (eg. P1, P2, Nedelec, and Raviart-Thomas).
*        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
*        a(u,v)_i, where i denotes a set of faces (or edges) within a boundary region marked with flag
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
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.
*
* \param old_sol                 FE approximation of previous solution if needed
* \param local_assembly_face     Routine to get local matrices over each face
* \param local_rhs_assembly_face Routine to get local rhs vectors over each face
* \param FE                      block FE Space
* \param mesh                    Mesh Data
* \param cq                      Quadrature Nodes
* \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
* \param coeff                   Function that gives coefficient (for now assume constant)
* \param time                    Physical Time if time dependent
* \param flag0,flag1             Marker for which faces are included in boundary integration (range of faces from flag0 to flag1)
*
* \return A                      Block Global stiffness CSR matrix
* \return b                      Global RHS vector
*
*/
void assemble_global_face_block(block_dCSRmat* A,dvector* b,dvector *old_sol,void (*local_assembly_face)(REAL *,dvector *,block_fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),void (*local_rhs_assembly_face)(REAL *,dvector *,block_fespace *,scomplex *,qcoordinates *,INT *,INT *,INT *,INT,INT,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,scomplex *sc,qcoordinates *cq,void (*coeff)(REAL *,REAL *,REAL,void *),void (*rhs)(REAL *,REAL *, REAL,void *),REAL time,INT flag0,INT flag1)
{

  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT dof_per_elm = 0;
  INT v_per_elm = (dim + 1);
  INT i,j,k,testdof,trialdof,jcntr,rowb,elm;
  INT nspaces = FE->nspaces;

  // Get block data first
  INT nblocks = A->brow;
  // Check for errors
  if(nblocks!=A->bcol) {
    printf("Your block matrix is not square.  It is an %lld x %lld matrix.\n\n",(long long )A->brow,(long long )A->bcol);
    exit(0);
  }
  if(nblocks!=nspaces) {
    printf("You have %lld FEM spaces, but only %lldx%lld blocks.  They must be consistent.\n\n",(long long )nspaces,(long long )A->brow,(long long )A->bcol);
    exit(0);
  }
  if(rhs!=NULL) {
    b->row = FE->ndof;
    if(b->val) {
      dvec_set(b->row,b,0.0);
    } else {
      b->val = (REAL *) calloc(b->row,sizeof(REAL));
    }
  }

  // Loop over each block and build sparsity structure of matrices
  for(i=0;i<nblocks;i++) {
    for(j=0;j<nblocks;j++) {
      testdof = FE->var_spaces[i]->ndof;
      trialdof = FE->var_spaces[j]->ndof;
      if(A->blocks[i*nblocks+j]) {
        if(A->blocks[i*nblocks+j]->IA==NULL){
          A->blocks[i*nblocks+j]->row = testdof; // test functions
          A->blocks[i*nblocks+j]->col = trialdof; // trial functions
          A->blocks[i*nblocks+j]->IA = (INT *) calloc(testdof+1,sizeof(INT));

          // Get Sparsity Structure First
          // Non-zeros of A and IA (ignores cancellations, so maybe more than necessary)
          create_CSR_rows_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);

          // Columns of A -> JA
          A->blocks[i*nblocks+j]->JA = (INT *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(INT));
          create_CSR_cols_FE1FE2(A->blocks[i*nblocks+j],FE->var_spaces[j],FE->var_spaces[i]);
        }

        // Set values
        if(A->blocks[i*nblocks+j]->val==NULL)
          A->blocks[i*nblocks+j]->val = (REAL *) calloc(A->blocks[i*nblocks+j]->nnz,sizeof(REAL));
        for (k=0; k<A->blocks[i*nblocks+j]->nnz; k++) {
          A->blocks[i*nblocks+j]->val[k] = 0;
        }
      }
    }
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }

  // Loop over each block and get dof_per_face
  INT dof_per_face = 0;
  INT* dof_per_face_blk = (INT *) calloc(nspaces,sizeof(INT));
  INT FEtype;
  for(i=0;i<nblocks;i++) {
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
    } else if (FEtype==99) { // Constraint Space (single DoF)
      dof_per_face_blk[i] = 1;
      dof_per_face += 1;
    } else if (FEtype==103) { // Mini
       dof_per_face_blk[i] = dim;
       dof_per_face += dim;
    } else {
      printf("Block face integration isn't set up for the FEM space you chose\n");
      check_error(ERROR_FE_TYPE,__FUNCTION__);
    }
  }

  // Now Build Global Matrix entries

  /* Loop over all Faces and build local matrix */
  INT local_size = dof_per_face*dof_per_face;
  REAL* ALoc = (REAL *) calloc(local_size,sizeof(REAL));
  REAL* bLoc=NULL;
  if(b!=NULL) bLoc = (REAL *) calloc(dof_per_face,sizeof(REAL));

  // Get mappings for given element and face
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* dof_on_f = (INT *) calloc(dof_per_face,sizeof(INT));
  INT rowa;

  // We will need the face to element map
  iCSRmat f_el;
  icsr_trans(fem->el_f,&f_el);

  // Loop over boundary faces
  for (i=0; i<fem->nface; i++) {
    // Only grab the faces on the flagged boundary
    if(fem->f_flag[i]>=flag0 && fem->f_flag[i]<=flag1) {
      // Zero out local matrices
      for (j=0; j<local_size; j++) {
        ALoc[j]=0;
      }
      if(rhs!=NULL) {
        for (j=0; j<dof_per_face; j++) {
          bLoc[j]=0;
        }
      }

      // Find DOF for given face
      jcntr = 0;
      for(k=0;k<nblocks;k++) {
        rowa = FE->var_spaces[k]->f_dof->IA[i];
        rowb = FE->var_spaces[k]->f_dof->IA[i+1];
        for (j=rowa; j<rowb; j++) {
          dof_on_f[jcntr] = FE->var_spaces[k]->f_dof->JA[j];
          jcntr++;
        }
      }

      // Find the corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i];
      elm = f_el.JA[rowa];

      // Find DOF for given Element
      // Note this is "local" ordering for the given FE space of the block
      // Not global ordering of all DOF
      jcntr = 0;
      for(k=0;k<nblocks;k++) {
        rowa = FE->var_spaces[k]->el_dof->IA[elm];
        rowb = FE->var_spaces[k]->el_dof->IA[elm+1];
        for (j=rowa; j<rowb; j++) {
          dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
          jcntr++;
        }
      }

      // Find vertices for given Element
      get_incidence_row(elm,fem->el_v,v_on_elm);

      // Compute Local Stiffness Matrix for given Element
      (*local_assembly_face)(ALoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,dof_per_face,i,elm,coeff,time);
      if(b!=NULL) (*local_rhs_assembly_face)(bLoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,dof_per_face,i,elm,rhs,time);

      // Loop over DOF and place in appropriate slot globally
      block_LocaltoGlobal_face(dof_on_f,dof_per_face,dof_per_face_blk,FE,b,A,ALoc,bLoc,flag0,flag1);
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(dof_per_face_blk) free (dof_per_face_blk);
  if(v_on_elm) free(v_on_elm);
  if(dof_on_f) free(dof_on_f);
  if(ALoc) free(ALoc);
  if(bLoc) free(bLoc);
  icsr_free(&f_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn void assemble_global_RHS_face_block(dvector *b,void (*local_rhs_assembly)(REAL *,block_fespace *,scomplex *,qcoordinates *,INT *,INT *,INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,scomplex *sc,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the global stiffness matrix for any "boundary" bilinear form using various element types
*        (eg. P1, P2, Nedelec, and Raviart-Thomas).
*        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
*        a(u,v)_i, where i denotes a set of faces (or edges) within a boundary region marked with flag
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
* \param local_rhs_assembly  Routine to assemble local RHS
* \param FE                  block FE Space
* \param mesh                Mesh Data
* \param cq                  Quadrature Nodes
* \param rhs                 Routine to get RHS function (NULL if only assembling matrix)
* \param time                Physical Time if time dependent
* \param flag0,flag1         Marker for which faces are included in boundary integration (range of faces from flag0 to flag1)
*
* \return b             Global RHS vector
*
*/
void assemble_global_RHS_face_block(dvector *b, dvector *old_sol, void (*local_rhs_assembly_face)(REAL *, dvector *, block_fespace *,scomplex *,qcoordinates *,INT *,INT *, INT *, INT, INT, INT,void (*)(REAL *,REAL *,REAL,void *),REAL),block_fespace *FE,scomplex *sc,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time, INT flag0, INT flag1)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT dof_per_elm = 0;
  INT v_per_elm = (dim + 1);
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

  // Loop over each block and get dof_per_face
  INT dof_per_face = 0;
  INT* dof_per_face_blk = (INT*)calloc(FE->nspaces,sizeof(INT));
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
    } else if (FEtype==99) { // Constraint Space (1 DoF)
      dof_per_face_blk[i] = 1;
      dof_per_face += 1;
    } else if (FEtype==99) { // Mini
      dof_per_face_blk[i] = dim;
      dof_per_face += dim;
    } else {
      printf("Block face integration isn't set up for the FEM space you chose\n");
      check_error(ERROR_FE_TYPE,__FUNCTION__);
    }
  }

  /* Loop over all Elements and build local rhs */
  REAL* bLoc = (REAL *) calloc(dof_per_face,sizeof(REAL));

  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* dof_on_f = (INT *) calloc(dof_per_face,sizeof(INT));

  // Need face to element map
  iCSRmat f_el;
  icsr_trans(fem->el_f,&f_el);

  INT jcntr,rowa,rowb;
  INT elm;
  //INT dof_on_f_shift;
  //INT dof_on_elm_shift;

  // Loop over boundary faces
  for(i=0;i<fem->nface;i++) {
    // Only grab faces on the flagged boundary
    if(fem->f_flag[i]>=flag0 && fem->f_flag[i]<=flag1) {
      // Zero out local vector
      for(j=0;j<dof_per_face;j++) {
        bLoc[j] = 0;
      }

      // Find DOF for given face
      jcntr = 0;
      for(k=0;k<nblocks;k++) {
        rowa = FE->var_spaces[k]->f_dof->IA[i];
        rowb = FE->var_spaces[k]->f_dof->IA[i+1];
        for (j=rowa; j<rowb; j++) {
          dof_on_f[jcntr] = FE->var_spaces[k]->f_dof->JA[j];
          jcntr++;
        }
      }

      // Find corresponding element associated with the face
      // Assume only 1 matters (if on boundary)
      rowa = f_el.IA[i];
      elm = f_el.JA[rowa];

      // Find DOF for given Element
      // Note this is "local" ordering for the given FE space of the block
      // Not global ordering of all DOF
      jcntr = 0;
      for(k=0;k<nblocks;k++) {
        rowa = FE->var_spaces[k]->el_dof->IA[elm];
        rowb = FE->var_spaces[k]->el_dof->IA[elm+1];
        for (j=rowa; j<rowb; j++) {
          dof_on_elm[jcntr] = FE->var_spaces[k]->el_dof->JA[j];
          jcntr++;
        }
      }
      //// DOF for given Face
      //// then Find DOF on that element
      //jcntr = 0;
      //dof_on_f_shift = 0;
      //dof_on_elm_shift = 0;
      //for(k=0;k<nblocks;k++){
      //  get_incidence_row(i,FE->var_spaces[k]->f_dof,dof_on_f+dof_on_f_shift);
      //  dof_on_f_shift += dof_per_face_blk[k];
      //
      //  get_incidence_row(elm,FE->var_spaces[k]->el_dof,dof_on_elm+dof_on_elm_shift);
      //  dof_on_elm_shift += FE->var_spaces[k]->dof_per_elm;
      //}

      // Find vertices for given element
      get_incidence_row(elm,fem->el_v,v_on_elm);

      // Compute Local RHS for given element
      (*local_rhs_assembly_face)(bLoc,old_sol,FE,sc,cq,dof_on_f,dof_on_elm,v_on_elm,dof_per_face,i,elm,rhs,time);

      // Put Local RHS in correct location
      jcntr = 0;
      rowa = 0;
      for(k=0;k<nblocks;k++) {
        for(j=0;j<dof_per_face_blk[k];j++) {
          row = dof_on_f[jcntr];
          b->val[row+rowa]+=bLoc[jcntr];
          jcntr++;
        }
        rowa += FE->var_spaces[k]->ndof;
      }
    }
  }

  if(dof_on_f) free(dof_on_f);
  if(dof_per_face_blk) free(dof_per_face_blk);
  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);
  icsr_free(&f_el);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn assemble_global_Ned_GradH1_RHS(dvector *b,fespace *FE_H1,fespace *FE_Ned,scomplex *sc,qcoordinates *cq,dvector* u)
*
* \brief Computes the global rhs for any <u,grad v>, where u is a Nedelec FE function and v is in P1
*        DOES NOT take care of Dirichlet boundary conditions.  A separate routine will eliminate them later
*
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.
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
void assemble_global_Ned_GradH1_RHS(dvector *b,fespace *FE_H1,fespace *FE_Ned,scomplex *sc,qcoordinates *cq,dvector* u)
{
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;
  INT ed_per_elm = FE_Ned->dof_per_elm;
  INT v_per_elm = (dim + 1);
  INT i,j,row;

  // Allocate Arrays
  b->row = FE_H1->ndof;
  b->val = (REAL *) calloc(b->row,sizeof(REAL));

  /* Loop over all Elements and build local rhs */
  REAL* bLoc= (REAL *) calloc(v_per_elm,sizeof(REAL));

  INT* ed_on_elm = (INT *) calloc(ed_per_elm,sizeof(INT));
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  for (i=0; i<FE_H1->nelm; i++) {

    for (j=0; j<v_per_elm; j++) {
      bLoc[j]=0;
    }

    // Find Edges for given Element
    get_incidence_row(i,FE_Ned->el_dof,ed_on_elm);

    // Find vertices for given Element
    get_incidence_row(i,fem->el_v,v_on_elm);

    // Compute Local RHS for given Element
    Ned_GradH1_RHS_local(bLoc,FE_H1,FE_Ned,sc,cq,ed_on_elm,v_on_elm,i,u);

    // Loop over DOF and place in appropriate slot globally
    for (j=0; j<v_per_elm; j++) {
      row = v_on_elm[j];
      b->val[row] += bLoc[j];
    }
  }

  if(ed_on_elm) free(ed_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(bLoc) free(bLoc);

  return;
}
/******************************************************************************************************/
