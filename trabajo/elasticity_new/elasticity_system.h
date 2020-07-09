/*! \file examples/elasticity/elasticity_system.h
 *
 *  Created by James Adler on 7/4/20.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routines
 *        for the elasticity example.
 *
 */

/*!
 * \fn void elasticity_mixed_system(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
 *
 * \brief Computes the local stiffness matrix for the Elasticity system in mixed formulation.
 *        For this problem we compute LHS system:
 *
 *        2*mu*<eps(u), eps(v)> - <p, div v>
 *                   -<div u, q> - (1/lam)*<p,q>
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
 * \note We provide the computation for the 0-0 block in 2D:
 *       <2 eps(u), eps(v)> =
 *                <2 dx(u1),dx(v1)> + <dy(u1),dy(v1)>  +     <dx(u2),dy(v1)>
 *                <dy(u1),dx(v2)>                      +     <dx(u2),dx(v2)> + <2 dy(u2),dy(v2)>
 *
 *       for Dirichlet boundary conditions and when div u = 0 exactly, then
 *       <2 eps(u), eps(v)> = <grad u, grad v>
 */
void elasticity_mixed_system(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  // Loop indices
  INT i,quad,test,trial;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;
  INT ip = dim;

  // DoF per element
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  INT u3dofpelm;
  if(dim==3) u3dofpelm = FE->var_spaces[iu3]->dof_per_elm;
  INT pdofpelm = FE->var_spaces[ip]->dof_per_elm;

  INT* local_dof_on_elm;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL aij = 0.0;

  // Test and Trial Functions
  REAL u1,u1x,u1y,u1z,u2,u2x,u2y,u2z,u3,u3x,u3y,u3z,p;
  REAL v1,v1x,v1y,v1z,v2,v2x,v2y,v2z,v3,v3x,v3y,v3z,q;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Get Lame Parameters (assume constant for now)
  REAL mu;
  get_mu(&mu,NULL,time,NULL);
  REAL lam;
  get_lam(&lam,NULL,time,NULL);
  // Watch out for case of lam = infinity
  REAL laminv = 1.0/lam;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    local_dof_on_elm = dof_on_elm;
    // u1
    get_FEM_basis(FE->var_spaces[iu1]->phi,FE->var_spaces[iu1]->dphi,qx,v_on_elm,dof_on_elm,mesh,FE->var_spaces[iu1]);
    local_dof_on_elm += u1dofpelm;
    // u2
    get_FEM_basis(FE->var_spaces[iu2]->phi,FE->var_spaces[iu2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[iu2]);
    local_dof_on_elm += u2dofpelm;
    // u3
    if(dim==3){
      get_FEM_basis(FE->var_spaces[iu3]->phi,FE->var_spaces[iu3]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[iu3]);
      local_dof_on_elm += u3dofpelm;
    }
    // p
    get_FEM_basis(FE->var_spaces[ip]->phi,FE->var_spaces[ip]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[ip]);

    // Loop over block rows of test functions
    local_row_index = 0;
    // u-v block: 2*mu*<eps(u),eps(v)> + lam*<div u,div v>
    // p-v block: -<p,div v>
    // u-q block: -<div u, q>
    // p-q block: -(1/lam)*<p,q>


    // v1 block row:
    for (test=0; test<u1dofpelm; test++) {
      v1 = FE->var_spaces[iu1]->phi[test];
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];

      // Loop over block columns of trial functions
      local_col_index = 0;

      // u1-v1 block: mu*(2*<u1x,v1x> + <u1y,v1y>)
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = mu*(2.0*u1x*v1x + u1y*v1y);
        if(dim==3) {
          u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij += mu*(u1z*v1z);
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block: mu*<u2x,v1y>
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = mu*(u2x*v1y);
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block: mu*<u3x,v1z>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3x = FE->var_spaces[iu3]->dphi[trial*dim];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3x*v1z);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v1 block: -<p,v1x>
      for (trial=0; trial<pdofpelm; trial++) {
        p = FE->var_spaces[ip]->phi[trial];
        aij = -p*v1x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
    } // End v1 block
    local_row_index += u1dofpelm;

    for (test=0; test<u2dofpelm; test++) {
      v2 = FE->var_spaces[iu2]->phi[test];
      v2x = FE->var_spaces[iu2]->dphi[test*dim];
      v2y = FE->var_spaces[iu2]->dphi[test*dim+1];

      local_col_index = 0;

      // u1-v2 block: mu*<u1y,v2x>
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = mu*u1y*v2x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block: mu*(<u2x,v2x> + 2*<u2y,v2y>)
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = mu*(u2x*v2x + 2.0*u2y*v2y);
        if(dim==3) {
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij+=mu*u2z*v2z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block: mu*<u3y,v2z>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3y = FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3y*v2z);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v2 block: -<p,v2y>
      for (trial=0; trial<pdofpelm; trial++) {
        p = FE->var_spaces[ip]->phi[trial];
        aij = -p*v2y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
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
        // mu*<dz(u1),dx(v3)>
        for (trial=0; trial<u1dofpelm;trial++){
          u1x=FE->var_spaces[iu1]->dphi[trial*dim];
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij = mu*u1z*v3x;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u1dofpelm;

        // u2-v3 block:
        // mu*<dz(u2),dy(v3)>
        for (trial=0; trial<u2dofpelm;trial++){
          u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij = mu*u2z*v3y;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u2dofpelm;

        // u3-v3 block
        // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)>
        for (trial=0; trial<u3dofpelm;trial++){
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3x*v3x + u3y*v3y + 2.0*u3z*v3z);
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;

        // p-v3 block: -<p, dz(v3)>
        for (trial=0; trial<pdofpelm;trial++){
          p = FE->var_spaces[ip]->phi[trial];
          aij = -p*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
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
        aij = -u1x*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
      }
      local_col_index += u1dofpelm;

      // u2-q block: -<dy(u2), q>
      for (trial=0; trial<u2dofpelm;trial++){
        u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = -u2y*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
      }
      local_col_index += u2dofpelm;

      if(dim==3) {
        // u3-q block: -<dz(u3), q>
        for (trial=0; trial<u3dofpelm;trial++){
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = -u3z*q;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
      }

      // p-q block: -1/lam <p,q>
      for (trial=0; trial<pdofpelm;trial++){
        p = FE->var_spaces[ip]->phi[trial];
        aij = -laminv*p*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
      }

    } // End q block
  } // End quadrature

  return;
}
/*****************************************************************************************************/

/*!
 * \fn void elasticity_primal_system(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
 *
 * \brief Computes the local stiffness matrix for the Elasticity system in primal form
 *        For this problem we compute LHS system:
 *
 *        2*mu*<eps(u), eps(v)> + lam*<div u, div v>
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
 * \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,u3)
*
 * \note We provide the computation for the 0-0 block in 2D:
 *       <2 eps(u), eps(v)> =
 *                <2 dx(u1),dx(v1)> + <dy(u1),dy(v1)>  +     <dx(u2),dy(v1)>
 *                <dy(u1),dx(v2)>                      +     <dx(u2),dx(v2)> + <2 dy(u2),dy(v2)>
 *
 *       for Dirichlet boundary conditions and when div u = 0 exactly, then
 *       <2 eps(u), eps(v)> = <grad u, grad v>
 */
void elasticity_primal_system(REAL* ALoc, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  // Loop indices
  INT i,quad,test,trial;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT nspaces = FE->nspaces;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;

  // DoF per element
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT u1dofpelm = FE->var_spaces[iu1]->dof_per_elm;
  INT u2dofpelm = FE->var_spaces[iu2]->dof_per_elm;
  INT u3dofpelm;
  if(dim==3) u3dofpelm = FE->var_spaces[iu3]->dof_per_elm;

  INT* local_dof_on_elm;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL aij = 0.0;

  // Test and Trial Functions
  REAL u1,u1x,u1y,u1z,u2,u2x,u2y,u2z,u3,u3x,u3y,u3z;
  REAL v1,v1x,v1y,v1z,v2,v2x,v2y,v2z,v3,v3x,v3y,v3z;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Get Lame Parameters (assume constant for now)
  REAL mu;
  get_mu(&mu,NULL,time,NULL);
  REAL lam;
  get_lam(&lam,NULL,time,NULL);

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    local_dof_on_elm = dof_on_elm;
    // u1
    get_FEM_basis(FE->var_spaces[iu1]->phi,FE->var_spaces[iu1]->dphi,qx,v_on_elm,dof_on_elm,mesh,FE->var_spaces[iu1]);
    local_dof_on_elm += u1dofpelm;
    // u2
    get_FEM_basis(FE->var_spaces[iu2]->phi,FE->var_spaces[iu2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[iu2]);
    local_dof_on_elm += u2dofpelm;
    // u3
    if(dim==3){
      get_FEM_basis(FE->var_spaces[iu3]->phi,FE->var_spaces[iu3]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[iu3]);
      local_dof_on_elm += u3dofpelm;
    }

    // Loop over block rows of test functions
    // u-v block: 2*mu*<eps(u),eps(v)> + lam*<div u,div v>
    local_row_index = 0;

    // v1 block row:
    for (test=0; test<u1dofpelm; test++) {
      v1 = FE->var_spaces[iu1]->phi[test];
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];

      // Loop over block columns of trial functions
      local_col_index = 0;

      // u1-v1 block: mu*(2*<u1x,v1x> + <u1y,v1y>) + lam*<u1x,v1x>
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = mu*(2.0*u1x*v1x + u1y*v1y) + lam*u1x*v1x;
        if(dim==3) {
          u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij += mu*(u1z*v1z);
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block: mu*<u2x,v1y> + lam*<u2y,v1x>
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = mu*(u2x*v1y) + lam*u2y*v1x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block: mu*<u3x,v1z> + lam*<u3z,v1x>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3x = FE->var_spaces[iu3]->dphi[trial*dim];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3x*v1z) + lam*u3z*v1x;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
      }
    } // End v1 block
    local_row_index += u1dofpelm;

    for (test=0; test<u2dofpelm; test++) {
      v2 = FE->var_spaces[iu2]->phi[test];
      v2x = FE->var_spaces[iu2]->dphi[test*dim];
      v2y = FE->var_spaces[iu2]->dphi[test*dim+1];

      local_col_index = 0;

      // u1-v2 block: mu*<u1y,v2x> + lam*<u1x,v2y>
      for (trial=0; trial<u1dofpelm; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = mu*u1y*v2x + lam*u1x*v2y;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block: mu*(<u2x,v2x> + 2*<u2y,v2y>) + lam*<u2y,v2y>
      for (trial=0; trial<u2dofpelm; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = mu*(u2x*v2x + 2.0*u2y*v2y) + lam*u2y*v2y;
        if(dim==3) {
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij+=mu*u2z*v2z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block: mu*<u3y,v2z> + lam*<u3z,v2y>
      if(dim==3) {
        for(trial=0;trial<u3dofpelm;trial++) {
          u3y = FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3y*v2z) + lam*u3z*v2y;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
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
        // mu*<dz(u1),dx(v3)> + lam*<u1x,v3z>
        for (trial=0; trial<u1dofpelm;trial++){
          u1x=FE->var_spaces[iu1]->dphi[trial*dim];
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij = mu*u1z*v3x + lam*u1x*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u1dofpelm;

        // u2-v3 block:
        // mu*<dz(u2),dy(v3)> + lam*<u2y,v3z>
        for (trial=0; trial<u2dofpelm;trial++){
          u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij = mu*u2z*v3y + lam*u2y*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u2dofpelm;

        // u3-v3 block
        // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)> + lam*<u3z,v3z>
        for (trial=0; trial<u3dofpelm;trial++){
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = mu*(u3x*v3x + u3y*v3y + 2.0*u3z*v3z) + lam*u3z*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += u3dofpelm;
      } // End v3 block
      local_row_index += u3dofpelm;
    } // End dim if
  } // End quadrature

  return;
}
/*****************************************************************************************************/


/*!
* \fn local_pressure_bdry(REAL* bLoc,dvector* old_sol,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_f,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the local boundary terms for Dirichlet conditions on p due to
*        integration by parts: <grad p, v> = -<p, div v> + <p,v*n>_boundary
*
*        Here, we compute as a rhs term, on any boundary for which v*n isn't zero:
*
*       -<p,v*n>_boundary
*
* \param old_sol       Solution at previous Newton step (NULL if linear)
* \param FE            FE Space
* \param mesh          Mesh Data
* \param cq            Quadrature Nodes
* \param dof_on_f      Specific DOF on face
* \param dof_on_elm    Specific DOF on element
* \param v_on_elm      Specific vertices on element
* \param dof_per_f     # DOF per face
* \param face          Current face
* \param elm           Current element
* \param rhs           Function for Boundary RHS (steady-state part)
* \param time          Physical Time if time dependent
*
* \return bLoc         Local RHS
*
* \note This assumes that we are using the mixed formulation of elasticity, and
*       we assemble for the entire block finite-element space.
*
*/
void local_pressure_bdry(REAL* bLoc,dvector* old_sol,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_f,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Loop Indices
  INT i,j,quad,test,doft;

  // Mesh and FE data
  INT dim = mesh->dim;

  INT* local_dof_on_elm;
  INT u1dof = FE->var_spaces[0]->ndof;
  INT u1dofpelm = FE->var_spaces[0]->dof_per_elm;
  INT u1dofpface = FE->var_spaces[0]->dof_per_face;
  INT u2dof = FE->var_spaces[1]->ndof;
  INT u2dofpelm = FE->var_spaces[1]->dof_per_elm;
  INT u2dofpface = FE->var_spaces[1]->dof_per_face;
  INT u3dof,u3dofpelm,u3dofpface;
  if(dim==3) {
    INT u3dof = FE->var_spaces[2]->ndof;
    INT u3dofpelm = FE->var_spaces[2]->dof_per_elm;
    INT u3dofpface = FE->var_spaces[2]->dof_per_face;
  }
  INT pdof = FE->var_spaces[dim]->ndof;
  INT pdofpelm = FE->var_spaces[dim]->dof_per_elm;
  INT pdofpface = FE->var_spaces[dim]->dof_per_face;

  // Quadrature Weights and Nodes
  qcoordinates *cq_face = allocateqcoords_bdry(cq->nq1d,1,dim,2);
  quad_face(cq_face,mesh,cq->nq1d,face);
  INT maxdim=4;
  REAL qx[maxdim];
  REAL w;

  // RHS Entry
  REAL rij = 0.0;

  // Solution on boundary
  REAL rhs_func[dim+1];
  REAL pres_val;

  // Normal vector on face
  REAL n1,n2,n3;
  n1 = mesh->f_norm[face*dim];
  n2 = mesh->f_norm[face*dim+1];
  if(dim==3) n3 = mesh->f_norm[face*dim+2];

  // Keep track of local indexing
  INT local_row_index, local_col_index, local_dofonelm_index;

  // Test Functions
  REAL v1,v2,v3;

  //  Sum over quadrature points
  for (quad=0;quad<cq_face->n;quad++) {
    qx[0] = cq_face->x[quad];
    qx[1] = cq_face->y[quad];
    if(dim==3)
    qx[2] = cq_face->z[quad];
    w = cq_face->w[quad];

    // Get pressure at Quadrature
    (*rhs)(rhs_func,qx,time,&(mesh->el_flag[elm]));
    pres_val = rhs_func[dim];

    /// Get basis functions for v1 and v2
    //  Get the Basis Functions and previous solutions at each quadrature node
    // u = (u1,u2) and v = (v1,v2)
    local_dof_on_elm = dof_on_elm;
    get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);
    local_dof_on_elm += u1dofpelm;
    get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
    local_dof_on_elm += u2dofpelm;
    if(dim==3) {
      get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
    }

    // Loop over Test Functions (dof on faces)
    // v1 block row
    local_row_index = 0;
    local_dofonelm_index=0;
    for (test=0; test<u1dofpface;test++){
      // Make sure ordering for global matrix is right
      for(j=0;j<u1dofpelm;j++) {
        if(dof_on_f[local_row_index+test]==dof_on_elm[local_dofonelm_index+j]) {
          doft = j;
        }
      }
      v1=FE->var_spaces[0]->phi[doft];
      rij = -pres_val*v1*n1;
      bLoc[(local_row_index+test)] += w*rij;
    }
    local_row_index += u1dofpface;
    local_dofonelm_index += u1dofpelm;

    // v2 block row
    for (test=0; test<u2dofpface;test++){
      // Make sure ordering for global matrix is right
      for(j=0;j<u2dofpelm;j++) {
        if(dof_on_f[local_row_index+test]==dof_on_elm[local_dofonelm_index+j]) {
          doft = j;
        }
      }
      v2=FE->var_spaces[1]->phi[doft];
      rij = -pres_val*v2*n2;
      bLoc[(local_row_index+test)] += w*rij;
    }
    local_row_index += u2dofpface;
    local_dofonelm_index += u2dofpelm;

    // v3 block row
    if(dim==3) {
      for (test=0; test<u3dofpface;test++){
        // Make sure ordering for global matrix is right
        for(j=0;j<u3dofpelm;j++) {
          if(dof_on_f[local_row_index+test]==dof_on_elm[local_dofonelm_index+j]) {
            doft = j;
          }
        }
        v3=FE->var_spaces[2]->phi[doft];
        rij = -pres_val*v3*n3;
        bLoc[(local_row_index+test)] += w*rij;
      }
    }
  } // End quadrature loop

  if(cq_face){
    free_qcoords(cq_face);
    free(cq_face);
    cq_face=NULL;
  }
  return;
}
/*****************************************************************************************************/
