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
* \fn void local_assembly_Stokes(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
*
* \brief Computes the local stiffness matrix for the Stokes system.
*        For this problem we compute LHS of:
*
*        <2 eps(u), eps(v)> - <p, div v> = <f, v>
*                   - <div u, q> = 0
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
* \note For this example we assume viscosity is 1.
*
*/
void local_assembly_Stokes(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  // Loop indices
  INT i,j,idim,quad,test,trial;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;
  INT ip = dim;

  // Keep track of DoF per element
  INT dof_per_elm = 0;
  for (i=0; i<nspaces;i++) dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  INT u3dofpelm;
  if(dim==3) u3dofpelm = FE->var_spaces[iu3]->dof_per_elm;
  INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;

  // Pointer to current unknowns DoF on element
  INT* local_dof_on_elm;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Test and Trial Functions
  REAL u1,u1x,u1y,u1z,u2,u2x,u2y,u2z,u3,u3x,u3y,u3z,p;
  REAL v1,v1x,v1y,v1z,v2,v2x,v2y,v2z,v3,v3x,v3y,v3z,q;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u = (u1,u2,u3,p) and v = (v1,v2,v3,q)
    local_dof_on_elm = dof_on_elm;
    for(i=0;i<FE->nspaces;i++) {
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);
      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
    }

    // Loop over block rows of test functions
    local_row_index = 0;
    // u-v block: <2 eps(u),eps(v)>
    // p-v block: -<p,div v>
    // u-q block: -<div u, q>
    // p-q block: 0

    // v1 block row:
    for (test=0; test<u1dofpelm; test++) {
      v1 = FE->var_spaces[iu1]->phi[test];
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
      if(dim==3) v1z = FE->var_spaces[iu1]->dphi[test*dim+2];

      // Loop over block columns of trial functions
      // u1-v1 block: 2*<u1x,v1x> + <u1y,v1y>
      local_col_index = 0;
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        kij = 2.0*u1x*v1x + u1y*v1y;
        if(dim==3) {
          u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
          kij += u1z*v1z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block: <u2x,v1y>
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        kij = u2x*v1y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block: <u3x,v1z>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3x = FE->var_spaces[iu3]->dphi[trial*dim];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          kij = u3x*v1z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v1 block: -<p,v1x>
      for (trial=0; trial<pdofpelm; trial++) {
        p = FE->var_spaces[ip]->phi[trial];
        kij = -p*v1x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
    } // End v1 block
    local_row_index += u1dofpelm;

    // v2 block row
    for (test=0; test<u2dofpelm; test++) {
      v2 = FE->var_spaces[iu2]->phi[test];
      v2x = FE->var_spaces[iu2]->dphi[test*dim];
      v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
      if(dim==3) v2z = FE->var_spaces[iu2]->dphi[test*dim+2];

      local_col_index = 0;

      // u1-v2 block: <u1y,v2x>
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        kij = u1y*v2x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block: <u2x,v2x> + 2*<u2y,v2y>
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        kij = u2x*v2x + 2.0*u2y*v2y;
        if(dim==3) {
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          kij+=u2z*v2z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block: <u3y,v2z>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3y = FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          kij = u3y*v2z;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v2 block: -<p,v2y>
      for (trial=0; trial<pdofpelm; trial++) {
        p = FE->var_spaces[ip]->phi[trial];
        kij = -p*v2y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*kij;
      }
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
        // <dz(u1),dx(v3)>
        for (trial=0; trial<u1dofpelm;trial++){
          u1x=FE->var_spaces[iu1]->dphi[trial*dim];
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          kij = u1z*v3x;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u1dofpelm;

        // u2-v3 block:
        // <dz(u2),dy(v3)>
        for (trial=0; trial<u2dofpelm;trial++){
          u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          kij = u2z*v3y;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u2dofpelm;

        // u3-v3 block
        // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)>
        for (trial=0; trial<u3dofpelm;trial++){
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          kij = u3x*v3x + u3y*v3y + 2.0*u3z*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;

        // p-v3 block: -<p, dz(v3)>
        for (trial=0; trial<pdofpelm;trial++){
          p = FE->var_spaces[ip]->phi[trial];
          kij = -p*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
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
    } // End q block
  }

  return;
}
