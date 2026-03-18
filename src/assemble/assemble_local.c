/*! \file src/assemble/assemble_local.c
 *
 * \brief This code will build local mass and stiffness matrices for various PDE systems
 *        Set up for just a few generic systems
 *        In general, the user would write their own specialized ones.
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 4/22/15.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note modified by James Adler 11/2/2016
 *  \note updated  by James Adler 02/21/2019 for 0-1 fix
 *
 */

#include "hazmath.h"

/******************************************************************************************************/
/*!
 * \fn void assemble_DuDv_local(REAL* ALoc,REAL* bLoc,REAL* u_local,simplex_local_data *elm_data,fe_local_data *fe_data,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local stiffness matrix for coeff*<Du,Dv> bilinear form using local data structs.
 *
 * \param ALoc         Local Stiffness Matrix (output)
 * \param bLoc         unused (kept for uniform callback signature)
 * \param u_local      unused (kept for uniform callback signature)
 * \param elm_data     Local simplex/mesh data
 * \param fe_data      Local FE data
 * \param rhs          unused (kept for uniform callback signature)
 * \param coeff        Coefficient function
 * \param time         Physical time
 *
 */
void assemble_DuDv_local(REAL* ALoc,REAL* bLoc,REAL* u_local,
    simplex_local_data *elm_data,fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;
  INT dof_per_elm = fe_data->n_dof_per_space[0];
  INT FEtype = fe_data->fe_types[0];
  INT flag = elm_data->flag;

  // Loop Indices
  INT quad,test,trial,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim+1];

  // Stiffness Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  // Vector Derivatives: Gradients (PX) and 3D Curls (3D Ned)
  if(FEtype<20 || (FEtype>=20 && FEtype<30 && dim==3)) {

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time,&flag);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_elm; trial++) {
          kij=0.0;
          for(idim=0;idim<dim;idim++)
            kij += coeff_val*(fe_data->dphi[0][test*dim+idim]*fe_data->dphi[0][trial*dim+idim]);
          ALoc[test*dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else { // Scalar Derivatives: Divergence (RT) and 2D Curls (2D Ned)

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time,&flag);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_elm; trial++) {
          kij = coeff_val*(fe_data->dphi[0][test]*fe_data->dphi[0][trial]);
          ALoc[test*dof_per_elm+trial] += w*kij;
        }
      }
    }
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void assemble_mass_local(REAL* ALoc,REAL* bLoc,REAL* u_local,simplex_local_data *elm_data,fe_local_data *fe_data,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local mass matrix for coeff*<u,v> bilinear form using local data structs.
 *
 * \param ALoc         Local Mass Matrix (output)
 * \param bLoc         unused (kept for uniform callback signature)
 * \param u_local      unused (kept for uniform callback signature)
 * \param elm_data     Local simplex/mesh data
 * \param fe_data      Local FE data
 * \param rhs          unused (kept for uniform callback signature)
 * \param coeff        Coefficient function
 * \param time         Physical time
 *
 */
void assemble_mass_local(REAL* ALoc,REAL* bLoc,REAL* u_local,
    simplex_local_data *elm_data,fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;
  INT dof_per_elm = fe_data->n_dof_per_space[0];
  INT scal_or_vec = fe_data->scal_or_vec[0];
  INT flag = elm_data->flag;

  // Loop Indices
  INT quad,test,trial,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim+1];

  // Mass Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  // Vector Functions
  if(scal_or_vec) {

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time,&flag);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_elm; trial++) {
          kij=0.0;
          for(idim=0;idim<dim;idim++)
            kij += coeff_val*(fe_data->phi[0][test*dim+idim]*fe_data->phi[0][trial*dim+idim]);
          ALoc[test*dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time,&flag);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_elm; trial++) {
          kij = coeff_val*(fe_data->phi[0][test]*fe_data->phi[0][trial]);
          ALoc[test*dof_per_elm+trial] += w*kij;
        }
      }
    }
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void FEM_RHS_Local(REAL* ALoc,REAL* bLoc,REAL* u_local,simplex_local_data *elm_data,fe_local_data *fe_data,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local RHS vector b_i = <f, phi_i> using local data structs.
 *
 * \param ALoc         unused (kept for uniform callback signature)
 * \param bLoc         Local RHS Vector (output)
 * \param u_local      unused (kept for uniform callback signature)
 * \param elm_data     Local simplex/mesh data
 * \param fe_data      Local FE data
 * \param rhs          RHS function
 * \param coeff        unused (kept for uniform callback signature)
 * \param time         Physical time
 *
 */
void FEM_RHS_Local(REAL* ALoc,REAL* bLoc,REAL* u_local,
    simplex_local_data *elm_data,fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;
  INT dof_per_elm = fe_data->n_dof_per_space[0];
  INT scal_or_vec = fe_data->scal_or_vec[0];
  INT flag = elm_data->flag;

  // Loop Indices
  INT quad,test,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim+1];

  // Right-hand side function at Quadrature Nodes
  REAL rhs_val_scalar;
  REAL rhs_val_vector[dim];

  if(scal_or_vec==0) { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      (*rhs)(&rhs_val_scalar,qx,time,&flag);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over test functions and integrate rhs
      for (test=0; test<dof_per_elm;test++) {
        bLoc[test] += w*rhs_val_scalar*fe_data->phi[0][test];
      }
    }
  } else { // Vector Functions

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      (*rhs)(rhs_val_vector,qx,time,&flag);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis_at_quadpt(elm_data, fe_data, 0, quad);

      // Loop over test functions and integrate rhs
      for (test=0; test<dof_per_elm;test++) {
        for(idim=0;idim<dim;idim++) {
          bLoc[test] += w*(rhs_val_vector[idim]*fe_data->phi[0][test*dim+idim]);
        }
      }
    }
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void FEM_Block_RHS_Local(REAL* ALoc,REAL* bLoc,REAL* u_local,simplex_local_data *elm_data,fe_local_data *fe_data,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local RHS vector for a block FEM system using local data structs.
 *
 * \param ALoc         unused (kept for uniform callback signature)
 * \param bLoc         Local RHS Vector (output)
 * \param u_local      unused (kept for uniform callback signature)
 * \param elm_data     Local simplex/mesh data
 * \param fe_data      Local FE data
 * \param rhs          RHS function (in FEM block ordering)
 * \param coeff        unused (kept for uniform callback signature)
 * \param time         Physical time
 *
 */
void FEM_Block_RHS_Local(REAL* ALoc,REAL* bLoc,REAL* u_local,
    simplex_local_data *elm_data,fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Loop Indices
  INT i,quad,test;

  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;
  INT nspaces = fe_data->nspaces;
  INT flag = elm_data->flag;
  INT nun=0;

  for(i=0;i<nspaces;i++) {
    if(fe_data->scal_or_vec[i]==0) /* Scalar Element */
      nun++;
    else /* Vector Element */
      nun += dim;
  }

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Right-hand side function at Quadrature Nodes
  REAL rhs_val[nun];

  // Offsets for block indexing
  INT dof_offset = 0;
  INT rhs_offset = 0;

  for(i=0;i<nspaces;i++) {
    INT dof_per_elm_i = fe_data->n_dof_per_space[i];
    INT scal_or_vec_i = fe_data->scal_or_vec[i];

    //  Sum over quadrature points
    for (quad=0;quad<nq;quad++) {
      qx[0] = elm_data->quad_local->x[quad*dim + 0];
      if(dim==2 || dim==3)
        qx[1] = elm_data->quad_local->x[quad*dim + 1];
      if(dim==3)
        qx[2] = elm_data->quad_local->x[quad*dim + 2];
      w = elm_data->quad_local->w[quad];
      (*rhs)(rhs_val,qx,time,&flag);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis_at_quadpt(elm_data, fe_data, i, quad);

      if(scal_or_vec_i==0) { // Scalar
        // Loop over test functions and integrate rhs
        for (test=0; test<dof_per_elm_i;test++) {
          bLoc[dof_offset+test] += w*rhs_val[rhs_offset]*fe_data->phi[i][test];
        }
      } else { // Vector
        INT idim;
        for (test=0; test<dof_per_elm_i;test++) {
          for(idim=0;idim<dim;idim++) {
            bLoc[dof_offset+test] += w*(rhs_val[rhs_offset+idim]*fe_data->phi[i][test*dim+idim]);
          }
        }
      }
    }

    dof_offset += dof_per_elm_i;
    if(scal_or_vec_i==0)
      rhs_offset++;
    else
      rhs_offset += dim;
  }

  return;
}
/******************************************************************************************************/
/******************************************************************************************************/
/*!
* \fn void assemble_symmetricDuDv_local(REAL* ALoc, block_fespace *FE, scomplex *sc, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
*
* \brief Computes the local stiffness matrix for a symmetric gradient.
*        For this problem we compute LHS of:
*
*        <eps(u), eps(v)> = <f, v>
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
* \note Assumes 2D or 3D only. Only works for Scalar Elements in block form
* \note We provide the computation for the 0-0 block in 2D:
*       <eps(u), eps(v)> =
*                <dx(u1),dx(v1)> + 0.5*<dy(u1),dy(v1)>  +  0.5*<dx(u2),dy(v1)>
*                0.5*<dy(u1),dx(v2)> + 0.5*<dx(u2),dx(v2)> + <dy(u2),dy(v2)>
*
*       for Dirichlet boundary conditions and when div u = 0 exactly, then
*       <2 eps(u), eps(v)> = <grad u, grad v>
*
*/
void assemble_symmetricDuDv_local(REAL* ALoc, block_fespace *FE, scomplex *sc, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  // Loop indices
  INT i,quad,test,trial;

  // Mesh and FE data
  INT dim = sc->dim;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;

  // DoF per element
  INT dof_per_elm = 0;
  INT dof_per_elm_arr[dim];
  for (i=0; i<dim;i++) {
    dof_per_elm_arr[i] = FE->var_spaces[i]->dof_per_elm;
    dof_per_elm += dof_per_elm_arr[i];
  }
  INT* local_dof_on_elm;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL aij = 0.0;

  // Test and Trial Functions
  REAL u1x,u1y,u1z,u2x,u2y,u2z,u3x,u3y,u3z;
  REAL v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(sc->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    local_dof_on_elm = dof_on_elm;
    for(i=0;i<dim;i++){
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,sc,FE->var_spaces[i]);
      // Shift local dof to next finite element space
      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
    }

    // Loop over block rows of test functions
    local_row_index = 0;

    // v1 block row:
    for (test=0; test<dof_per_elm_arr[iu1]; test++) {
      v1x = FE->var_spaces[iu1]->dphi[test*dim];
      v1y = FE->var_spaces[iu1]->dphi[test*dim+1];
      if(dim==3) v1z = FE->var_spaces[iu1]->dphi[test*dim+2];

      // Loop over block columns of trial functions
      local_col_index = 0;

      // u1-v1 block: <u1x,v1x> + 0.5*<u1y,v1y>
      for (trial=0; trial<dof_per_elm_arr[iu1]; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = u1x*v1x + 0.5*u1y*v1y;
        if(dim==3) {
          u1z = FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij += 0.5*u1z*v1z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += dof_per_elm_arr[iu1];

      // u2-v1 block: 0.5*<u2x,v1y>
      for (trial=0; trial<dof_per_elm_arr[iu2]; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = 0.5*(u2x*v1y);
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += dof_per_elm_arr[iu2];

      // u3-v1 block: 0.5*<u3x,v1z>
      if(dim==3) {
        for(trial=0;trial<dof_per_elm_arr[iu3];trial++) {
          u3x = FE->var_spaces[iu3]->dphi[trial*dim];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = 0.5*(u3x*v1z);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += dof_per_elm_arr[iu3];
      }
    } // End v1 block
    local_row_index += dof_per_elm_arr[iu1];

    for (test=0; test<dof_per_elm_arr[iu2]; test++) {
      v2x = FE->var_spaces[iu2]->dphi[test*dim];
      v2y = FE->var_spaces[iu2]->dphi[test*dim+1];
      if(dim==3) v2z = FE->var_spaces[iu2]->dphi[test*dim+2];


      local_col_index = 0;

      // u1-v2 block: 0.5*<u1y,v2x>
      for (trial=0; trial<dof_per_elm_arr[iu1]; trial++) {
        u1x = FE->var_spaces[iu1]->dphi[trial*dim];
        u1y = FE->var_spaces[iu1]->dphi[trial*dim+1];
        aij = 0.5*u1y*v2x;
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += dof_per_elm_arr[iu1];

      // u2-v2 block: 0.5<u2x,v2x> + <u2y,v2y>
      for (trial=0; trial<dof_per_elm_arr[iu2]; trial++) {
        u2x = FE->var_spaces[iu2]->dphi[trial*dim];
        u2y = FE->var_spaces[iu2]->dphi[trial*dim+1];
        aij = 0.5*u2x*v2x + u2y*v2y;
        if(dim==3) {
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij+=0.5*u2z*v2z;
        }
        ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
      }
      local_col_index += dof_per_elm_arr[iu2];

      // u3-v2 block: 0.5*<u3y,v2z>
      if(dim==3) {
        for(trial=0;trial<dof_per_elm_arr[iu3];trial++) {
          u3y = FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z = FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = 0.5*(u3y*v2z);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
        local_col_index += dof_per_elm_arr[iu3];
      }
    } // End v2 block
    local_row_index += dof_per_elm_arr[iu2];

    // v3 block row
    if(dim==3) {
      for (test=0; test<dof_per_elm_arr[iu3];test++){
        v3x=FE->var_spaces[iu3]->dphi[test*dim];
        v3y=FE->var_spaces[iu3]->dphi[test*dim+1];
        v3z=FE->var_spaces[iu3]->dphi[test*dim+2];

        local_col_index = 0;

        // u1-v3 block
        // 0.5*<dz(u1),dx(v3)>
        for (trial=0; trial<dof_per_elm_arr[iu1];trial++){
          u1x=FE->var_spaces[iu1]->dphi[trial*dim];
          u1z=FE->var_spaces[iu1]->dphi[trial*dim+2];
          aij = 0.5*u1z*v3x;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += dof_per_elm_arr[iu1];

        // u2-v3 block:
        // 0.5*<dz(u2),dy(v3)>
        for (trial=0; trial<dof_per_elm_arr[iu2];trial++){
          u2y=FE->var_spaces[iu2]->dphi[trial*dim+1];
          u2z=FE->var_spaces[iu2]->dphi[trial*dim+2];
          aij = 0.5*u2z*v3y;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += dof_per_elm_arr[iu2];

        // u3-v3 block
        // 0.5*<dx(u3),dx(v3)> + 0.5*<dy(u3),dy(v3)> + <dz(u3),dz(v3)>
        for (trial=0; trial<dof_per_elm_arr[iu3];trial++){
          u3x=FE->var_spaces[iu3]->dphi[trial*dim];
          u3y=FE->var_spaces[iu3]->dphi[trial*dim+1];
          u3z=FE->var_spaces[iu3]->dphi[trial*dim+2];
          aij = 0.5*u3x*v3x + 0.5*u3y*v3y + u3z*v3z;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*aij;
        }
        local_col_index += dof_per_elm_arr[iu3];
      } // End v3 block
    } // End dim if
  } // End quadrature

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void Ned_GradH1_RHS_local(REAL* bLoc,fespace *FE_H1,fespace *FE_Ned,scomplex *sc,qcoordinates *cq,INT *ed_on_elm,INT *v_on_elm,INT elm,dvector* u)
 *
 * \brief Computes the local weak formulation of <E,grad(q)> where E is a given
 *         Nedelec approximation and q in H_0^1 (linears)
 *
 * \param FE_H1         FE Space for H1 elements
 * \param FE_Ned        FE Space for Nedelec elements
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param ed_on_elm     Specific edges on the given element
 * \param v_on_elm      Specific vertices on the given element
 * \param elm           Current element
 * \param u             FEM Function that gives coefficient
 *
 * \return bLoc         Local RHS vector (Full Matrix)
 *
 * \note                Assuming 2D and 3D only
 *
 */
void Ned_GradH1_RHS_local(REAL* bLoc,fespace *FE_H1,fespace *FE_Ned,scomplex *sc,qcoordinates *cq,INT *ed_on_elm,INT *v_on_elm,INT elm,dvector* u)
{
  // Mesh and FE data
  INT dim = sc->dim;

  // Loop Indices
  INT quad,test;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim+1];

  // Right-hand side function at Quadrature Nodes
  REAL ucoeff[3];

  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    // Get FEM function at quadrature nodes
    FE_Interpolation(ucoeff,u->val,qx,ed_on_elm,v_on_elm,FE_Ned,sc);

    //  Get the Basis Functions at each quadrature node
    PX_H1_basis(FE_H1->phi,FE_H1->dphi,qx,v_on_elm,FE_H1->FEtype,sc);

    // Loop over test functions and integrate rhs
    for (test=0; test<FE_H1->dof_per_elm;test++) {
      bLoc[test] += w*(ucoeff[0]*FE_H1->dphi[test]+ucoeff[1]*FE_H1->dphi[test]);
      if(dim==3) bLoc[test] += w*ucoeff[2]*FE_H1->dphi[test];
    }
  }

  return;
}
/******************************************************************************************************/

/****** Boundary Assemblies *********************/
/******************************************************************************************************/
/*!
* \fn void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,scomplex *sc,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the local assembly of a RHS for any "boundary" bilinear form using various element types
*        (eg. P1, P2, Nedelec, and Raviart-Thomas).
*        This does integration over a surface or boundary (i.e., faces of your domain: faces in 3D, edges in 2D)
*        a(u,v)_i, where i denotes a set of faces (or edges) with in a boundary region marked with flag
*
*
*        For this problem we compute local RHS of:
*
*        Lu = f  ---->   a(u,v)_bdry = <f,v>_bdry
*
*        which gives Ax = b,
*
*        A_ij = a( phi_j, phi_i)_bdry
*
* \note All matrices are assumed to be indexed at 0 in the CSR formatting.

* \note Assumes different type of integral for different Element type:
*       Scalar -> <f,v>_bdry
*       Vector -> <f,n*v>_bdry
*
* \note This reallocates the quadrature on each face for every call to this function.
*       This is not optimal, and the quadrature should be set outside if you are
*       building your own local assembly.
*
* \param old_sol                 FE approximation of previous solution if needed
* \param FE                      FE Space
* \param mesh                    Mesh Data
* \param dof_on_f                DOF on the given face
* \param dof_on_elm              DOF on given element
* \param v_on_elm                Vertices on given element
* \param dof_per_f               # of DOF per face
* \param face                    Given Face
* \param elm                     Given Element
* \param rhs                     Routine to get RHS function (NULL if only assembling matrix)
* \param time                    Physical Time if time dependent
*
* \return bLoc                   Local RHS vector
*
*/
void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,scomplex *sc,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{
  // Mesh and FE data
  sc_fem *fem = sc->fem;
  INT dim = sc->dim;

  // Loop Indices
  INT j,quad,test,doft;

  // Quadrature Weights and Nodes
  qcoordinates *cq_face = allocateqcoords_bdry(cq->nq1d,1,dim,2);
  quad_face(cq_face,sc,cq->nq1d,face);
  REAL qx[dim+1];
  REAL w;

  // Get normal vector components on face if needed
  REAL nx=0.0,ny=0.0,nz=0.0;
  if(FE->FEtype>=20) {
    nx = fem->f_norm[face*dim];
    ny = fem->f_norm[face*dim+1];
    if(dim==3) nz = fem->f_norm[face*dim+2];
  }

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Value of f
  REAL rhs_val;

  if(FE->scal_or_vec==0) { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq_face->n;quad++) {
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time,&(fem->f_flag[face]));

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,sc,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_face;test++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<FE->dof_per_elm;j++) {
          if(dof_on_f[test]==dof_on_elm[j]) {
            doft = j;
          }
        }
        kij = rhs_val*FE->phi[doft];
        bLoc[test] += w*kij;
      }
    }
  } else { // Vector Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq_face->n;quad++) {
    //for (quad=0;quad<1;quad++) {
    //  qx[0] = fem->f_mid[face*dim];
    //  qx[1] = fem->f_mid[face*dim+1];
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
//        qx[2] = fem->f_mid[face*dim+2];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time,&(fem->f_flag[face]));

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,sc,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<dof_per_face;test++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<FE->dof_per_elm;j++) {
          if(dof_on_f[test]==dof_on_elm[j]) {
            doft = j;
          }
        }
        kij = rhs_val*(nx*FE->phi[doft*dim] + ny*FE->phi[doft*dim+1]);
        if(dim==3) kij +=rhs_val*nz*FE->phi[doft*dim+2];
        bLoc[test] += w*kij;
      }
    }
  }

  return;
}
/******************************************************************************************************/
