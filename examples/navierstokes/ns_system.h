/*! \file examples/stokes/stokes_system.h
*
*  Created by Peter Ohm on 1/5/17.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the local assembly routines
*        for the Stokes example.
*
* \note Updated by James Adler on 02/04/21
*/

/*!
* \fn void local_assembly_NS(REAL* ALoc, REAL *bLoc, dvector *old_sol, block_fespace *FE, scomplex *sc, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, void (*rhs)(...), void (*coeff)(...), REAL time)
*
* \brief Computes the local stiffness matrix for the linearized Navier-Stokes system.
*  
*            <2*eps(u), eps(v)> + Re<u0*grad(u) + u*grad(u0),v> - <p, div v> = <f, v> - <2*eps(u0), eps(v)> - Re<u0*grad(u0),v> + <p0, div v>
*              - <div u, q> = <div u0, q>
* where u0 is the previous Newton iterate.
*
*        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient.
*
* \param FE            Block FE Space
* \param mesh          Mesh Data
* \param cq            Quadrature Nodes
* \param dof_on_elm    Specific DOF on element
* \param v_on_elm      Specific vertices on element
* \param elm           Current element
* \param time          Physical Time if time dependent
*
* \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,u3,p)
*
* \note Assumes 2D or 3D only
*
* \note We provide the computation for the 0-0 block:
*       <2 eps(u), eps(v)> =
*                <2 dx(u1),dx(v1)> + <dy(u1),dy(v1)>  +     <dx(u2),dy(v1)>
*                <dy(u1),dx(v2)>                      +     <dx(u2),dx(v2)> + <2 dy(u2),dy(v2)>
*
*       for Dirichlet boundary conditions and when div u = 0 exactly, then
*       <2 eps(u), eps(v)> = <grad u, grad v>
*
*
*/
void local_assembly_NS(REAL *ALoc,REAL* bLoc, dvector *old_sol, block_fespace *FE, scomplex *sc, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop indices
  INT i,j,idim,quad,test,trial;

  // Mesh and FE data
  INT dim = sc->dim;
  INT nspaces = FE->nspaces;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;
  INT ip = dim;

  // Keep track of DoF per element
  INT dof_per_elm = 0;
  for (i=0; i<nspaces;i++) dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT u1dof = FE->var_spaces[iu1]->ndof;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dof = FE->var_spaces[iu2]->ndof;
  INT u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  INT u3dof,u3dofpelm;
  if(dim==3) {
    INT u3dof = FE->var_spaces[iu3]->ndof;
    INT u3dofpelm = FE->var_spaces[iu3]->dof_per_elm;
  }
  INT pdof = FE->var_spaces[ip]->ndof;
  INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;

  // Pointer to current unknowns DoF on element
  INT* local_dof_on_elm;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry and RHS entry
  REAL kij = 0.0;
  REAL rij = 0.0;

   // Stuff for previous solution
  REAL* local_uprev = NULL;
  REAL u1_prev;
  REAL u2_prev;
  REAL u3_prev;
  REAL p_prev;
  REAL* du1_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL* du2_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL* du3_prev = NULL;
  if(dim==3) du3_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL divuk,epsu11,epsu21,epsu31,epsu22,epsu32,epsu33;
  REAL ukgraduk1,ukgraduk2,ukgraduk3;

  // Test and Trial Functions
  REAL u1,u1x,u1y,u1z,u2,u2x,u2y,u2z,u3,u3x,u3y,u3z,p;
  REAL v1,v1x,v1y,v1z,v2,v2x,v2y,v2z,v3,v3x,v3y,v3z,q;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Get Constants
  REAL Re = 0.0;
  get_reynolds_number(&Re);

  // Source Term Values if any
  REAL rhs_val[dim+1];

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(sc->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    (*rhs)(rhs_val,qx,time,&(sc->flags[elm]));

    //  Get the Basis Functions at each quadrature node and previous solutions
    // u = (u1,u2,u3,p) and v = (v1,v2,v3,q)
    local_dof_on_elm = dof_on_elm;
    local_uprev = old_sol->val;
    FE_Interpolation(&u1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu1],sc);
    FE_DerivativeInterpolation(du1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu1],sc);
    get_FEM_basis(FE->var_spaces[iu1]->phi,FE->var_spaces[iu1]->dphi,qx,v_on_elm,local_dof_on_elm,sc,FE->var_spaces[iu1]);
    local_dof_on_elm += u1dofpelm;
    local_uprev += u1dof;

    FE_Interpolation(&u2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu2],sc);
    FE_DerivativeInterpolation(du2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu2],sc);
    get_FEM_basis(FE->var_spaces[iu2]->phi,FE->var_spaces[iu2]->dphi,qx,v_on_elm,local_dof_on_elm,sc,FE->var_spaces[iu2]);
    local_dof_on_elm += u2dofpelm;
    local_uprev += u2dof;

    if(dim==3) {
      FE_Interpolation(&u3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu3],sc);
      FE_DerivativeInterpolation(du3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[iu3],sc);
      get_FEM_basis(FE->var_spaces[iu3]->phi,FE->var_spaces[iu3]->dphi,qx,v_on_elm,local_dof_on_elm,sc,FE->var_spaces[iu3]);
      local_dof_on_elm += u3dofpelm;
      local_uprev += u3dof;
    }

    // p
    FE_Interpolation(&p_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[ip],sc);
    get_FEM_basis(FE->var_spaces[ip]->phi,FE->var_spaces[ip]->dphi,qx,v_on_elm,local_dof_on_elm,sc,FE->var_spaces[ip]);
    local_dof_on_elm += pdofpelm;
    local_uprev += pdof;

    // Some precomputations
    divuk = du1_prev[0] + du2_prev[1];
    epsu11 = du1_prev[0];
    epsu21 = 0.5*(du1_prev[1] + du2_prev[0]);
    epsu22 = du2_prev[1];
    ukgraduk1 = u1_prev*du1_prev[0] + u2_prev*du1_prev[1];
    ukgraduk2 = u1_prev*du2_prev[0] + u2_prev*du2_prev[1];
    if(dim==3) {
      divuk += du3_prev[2];
      epsu31 = 0.5*(du1_prev[2] + du3_prev[0]);
      epsu32 = 0.5*(du2_prev[2] + du3_prev[1]);
      ukgraduk1 += u3_prev*du1_prev[2];
      ukgraduk2 += u3_prev*du2_prev[2];
      ukgraduk3 = u1_prev*du3_prev[0] + u2_prev*du3_prev[1] + u3_prev*du3_prev[2];
    }

    // Loop over block rows of test functions
    local_row_index = 0;

    // v1 block row:
    for (test=0; test<u1dofpelm; test++) {
      v1 = FE->var_spaces[iu1]->phi[test];
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
      if(dim==3) v1z = FE->var_spaces[iu1]->dphi[test*dim+2];

      local_col_index = 0;

      // Loop over block columns of trial functions
      // u1-v1 block:
      // 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)> + <dz(u1),dz(v1)> + Re*<u1k*dx(u1)+u2k*dy(u1) + u1*dx(u1k),v>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u1dofpelm;trial++){
        u1=FE->var_spaces[iu1]->phi[trial];
        u1x=FE->var_spaces[iu1]->dphi[trial*dim];
        u1y=FE->var_spaces[iu1]->dphi[trial*dim+1];
        kij = 2*u1x*v1x + u1y*v1y + Re*(u1_prev*u1x+u2_prev*u1y+u1*du1_prev[0])*v1;
        if(dim==3) {
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          kij+=u1z*v1z + Re*(u3_prev*u1z)*v1;
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block:
      // <dx(u2),dy(v1)> + Re*<u2*dy(u1k),v1>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u2dofpelm;trial++){
        u2=FE->var_spaces[iu2]->phi[trial];
        u2x=FE->var_spaces[iu2]->dphi[trial*dim];
        kij = u2x*v1y + Re*u2*du1_prev[1]*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block <dx(u3),dz(v1)> + Re*<u3*dz(u1k),v1>
      if(dim==3) {
        for (trial=0; trial<u3dofpelm;trial++){
          u3=FE->var_spaces[iu3]->phi[trial];
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          kij = u3x*v1z + Re*u3*du1_prev[2]*v1;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v1 block: -<p, dx(v1)>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<pdofpelm;trial++){
        p = FE->var_spaces[ip]->phi[trial];
        kij = -p*v1x;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v1 part of RHS
      rij = rhs_val[0]*v1 - 2*du1_prev[0]*v1x - 2*epsu21*v1y - Re*ukgraduk1*v1 + p_prev*v1x;
      if(dim==3) rij+= -2*epsu31*v1z;
      bLoc[(local_row_index+test)] += w*rij;
    } // End v1 block
    local_row_index += u1dofpelm;

    // v2 block row
    for (test=0; test<u2dofpelm; test++) {
      v2 = FE->var_spaces[iu2]->phi[test];
      v2x = FE->var_spaces[iu2]->dphi[test*dim];
      v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
      if(dim==3) v2z = FE->var_spaces[iu2]->dphi[test*dim+2];

      local_col_index = 0;

      // u1-v2 block
      // <dy(u1),dx(v2)> + Re*<u1*dx(u2k),v2>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u1dofpelm;trial++){
        u1=FE->var_spaces[iu1]->phi[trial];
        u1y=FE->var_spaces[iu1]->dphi[trial*dim+1];
        kij = u1y*v2x + Re*u1*du2_prev[0]*v2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block:
      // <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)> + <dz(u2),dz(v2)> + Re*<u_prev*grad(u2) + u2*dy(u2k),v2>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u2dofpelm;trial++){
        u2=FE->var_spaces[iu2]->phi[trial];
        u2x=FE->var_spaces[iu2]->dphi[trial*dim];
        u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
        kij = u2x*v2x + 2*u2y*v2y + Re*(u1_prev*u2x + u2_prev*u2y + u2*du2_prev[1])*v2;
        if(dim==3) {
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          kij+=u2z*v2z + Re*(u3_prev*u2z)*v2;
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block
      // <dy(u3),dz(v2)> + Re*<u3*dz(u2k),v2>
      if(dim==3) {
        for (trial=0; trial<u3dofpelm;trial++){
          u3=FE->var_spaces[iu3]->phi[trial];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          kij = u3y*v2z + Re*u3*du2_prev[2]*v2;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v2 block: -<p, dy(v2)>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<pdofpelm;trial++){
        p = FE->var_spaces[ip]->phi[trial];
        kij = -p*v2y;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v2 part of RHS
      rij = rhs_val[1]*v2 - 2*du2_prev[1]*v2y - 2*epsu21*v2x - Re*ukgraduk2*v2 + p_prev*v2y;
      if(dim==3) rij+= -2*epsu32*v2z;
      bLoc[(local_row_index+test)] += w*rij;
    } // End v2 block
    local_row_index += u2dofpelm;

    // v3 block row
    if(dim==3) {
      for (test=0; test<u3dofpelm;test++){
        v3=FE->var_spaces[iu3]->phi[test];
        v3x=FE->var_spaces[iu3]->dphi[test*dim];
        v3y=FE->var_spaces[iu3]->dphi[test*dim+1];
        v3z=FE->var_spaces[iu3]->dphi[test*dim+2];

        local_col_index = 0;

        // u1-v3 block
        // <dz(u1),dx(v3)> + Re*<u1*dx(u3k),v3>
        local_col_index = 0;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<u1dofpelm;trial++){
          u1=FE->var_spaces[iu1]->phi[trial];
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          kij = u1z*v3x + Re*u1*du3_prev[0]*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // u2-v3 block:
        // <dz(u2),dy(v3)> + Re*<u2*dy(u3k),v3>
        local_col_index += u1dofpelm;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<u2dofpelm;trial++){
          u2=FE->var_spaces[iu2]->phi[trial];
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          kij = u2z*v3y + Re*u2*du3_prev[1]*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // u3-v3 block
        // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)> + Re*<u_prev*grad(u3) + u3*dz(u3k),v3>
        local_col_index += u2dofpelm;
        for (trial=0; trial<u3dofpelm;trial++){
          u3=FE->var_spaces[iu3]->phi[trial];
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          kij = u3x*v3x + u3y*v3y + 2*u3z*v3z + Re*(u1_prev*u3x + u2_prev*u3y + u3_prev*u3z + u3*du3_prev[2])*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // p-v3 block: -<p, dz(v3)>
        local_col_index += u3dofpelm;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<pdofpelm;trial++){
          p = FE->var_spaces[ip]->phi[trial];
          kij = -p*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // v3 part of RHS
        rij = rhs_val[2]*v3 - 2*du3_prev[2]*v3z - 2*epsu31*v3x - 2*epsu32*v3y - Re*ukgraduk3*v3 + p_prev*v3z;
        bLoc[(local_row_index+test)] += w*rij;
      } // End v3 block
      local_row_index += u3dofpelm;
    } // End dim if

    // q block row
    for (test=0; test<pdofpelm;test++){
      q=FE->var_spaces[ip]->phi[test];
      local_col_index = 0;

      // u1-q block: -<dx(u1), q>
      for (trial=0; trial<u1dofpelm;trial++){
        u1x=FE->var_spaces[iu1]->dphi[trial*dim];
        kij = -u1x*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-q block: -<dy(u2), q>
      for (trial=0; trial<u2dofpelm;trial++){
        u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
        kij = -u2y*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      if(dim==3) {
        // u3-q block: -<dz(u3), q>
        for (trial=0; trial<u3dofpelm;trial++){
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          kij = -u3z*q;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-q block: 0

      // q part of RHS
      rij = divuk*q;
      bLoc[(local_row_index+test)] += w*rij;
    } // End q block
  } // End quadrature loop

  if(du1_prev) free(du1_prev);
  if(du2_prev) free(du2_prev);
  if(du3_prev) free(du3_prev);

  return;
}

/*!
* \fn void meanzero_pressure(block_dCSR* A, dvector* b,mesh_struct* mesh,fespace* FEp,qcoordinates* cq)
*
* \brief Deals with pressure singularity by adding constraint that
*         int p = int p_true (0 if mean zero constraint)
*        Do this by adding a row and column to the p-p block that computes this
*        integral.  Note that this increases the DoF by 1 and the rhs by 1.
*
* \param A            Stokes system matrix
* \param b            System RHS
* \param mesh         Mesh Information
* \param FEp          FE space for pressure
* \param cq           Quadrature needed for anything but P0
*
* \return A, b        Modified system
*
*/
void meanzero_pressure(block_dCSRmat* A, dvector* b,scomplex* sc,fespace* FEp,qcoordinates* cq)
{

  INT ndof = FEp->ndof;
  INT i,newrows;

  // Reallocate rhs appropriately and add in value of mean if not zero
  newrows = b->row+1;

  // Copy current values into new longer vector
  REAL* btmp = (REAL *) calloc(newrows,sizeof(REAL));
  for(i=0;i<b->row;i++) btmp[i] = b->val[i];
  // Get mean value
  REAL pmeanval = 0.0;
  pmean(&pmeanval);
  btmp[b->row] = pmeanval;
  b->val=(REAL *)realloc(b->val,sizeof(REAL)*(newrows));

  // Copy btmp back into b
  b->row++;
  for(i=0;i<newrows;i++) b->val[i] = btmp[i];

  // Reallocate matrix
  // In P0 -> Add el_vol for each entry
  if(FEp->FEtype==0) {
    bdcsr_extend(A,sc->fem->el_vol,sc->fem->el_vol,sc->dim,1.0,1.0);
  } else {
    // In P1 or higher, the extra row is Mp*1, where Mp is the mass matrix for
    // the pressure space and 1 is the vector of ones.
    dCSRmat Mp = dcsr_create(0,0,0);
    assemble_global_single(&Mp, NULL, FEp, sc, cq,
                           local_assembly_mass, NULL, NULL, one_coeff_scal, 0.0);
    dvector* ones = dvec_create_p(ndof);
    dvec_set(ndof,ones,1.0);
    REAL* constraint = (REAL *) calloc(ndof,sizeof(REAL));
    dcsr_mxv(&Mp,ones->val,constraint);
    bdcsr_extend(A,constraint,constraint,sc->dim,1.0,1.0);

    dcsr_free(&Mp);
    if(ones) free(ones);
    if(constraint) free(constraint);
  }

  if(btmp) free(btmp);

  return;
}
