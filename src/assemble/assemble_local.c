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
* \fn void assemble_symmetricDuDv_local(REAL* ALoc,REAL* bLoc,REAL* u_local,simplex_local_data *elm_data,fe_local_data *fe_data,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
*
* \brief Computes the local stiffness matrix for a symmetric gradient.
*        For this problem we compute LHS of:
*
*        <eps(u), eps(v)> = <f, v>
*
*        where eps(u) = (grad u + (grad u)^T)/2 is the symmetric gradient.
*
* \param ALoc         Local Stiffness Matrix (output)
* \param bLoc         unused (kept for uniform callback signature)
* \param u_local      unused (kept for uniform callback signature)
* \param elm_data     Local simplex/mesh data
* \param fe_data      Local FE data
* \param rhs          unused (kept for uniform callback signature)
* \param coeff        Coefficient function
* \param time         Physical Time if time dependent
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
void assemble_symmetricDuDv_local(REAL* ALoc,REAL* bLoc,REAL* u_local,
    simplex_local_data *elm_data,fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop indices
  INT i,quad,test,trial;

  // Mesh and FE data
  INT dim = elm_data->dim;
  INT nq = elm_data->quad_local->nq;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;

  // DoF per element
  INT dof_per_elm = 0;
  INT dof_per_elm_arr[dim];
  INT block_offset[dim];
  for (i=0; i<dim;i++) {
    dof_per_elm_arr[i] = fe_data->n_dof_per_space[i];
    block_offset[i] = dof_per_elm;
    dof_per_elm += dof_per_elm_arr[i];
  }

  // Quadrature Weights and Nodes
  REAL w;

  // Stiffness Matrix Entry
  REAL aij;

  INT a,b,c;

  /* The symmetric gradient bilinear form <eps(u),eps(v)> for block (a,b) is:
   * aij = 0.5*(delta_{ab} * sum_c dphi[a][test,c]*dphi[b][trial,c]
   *           + dphi[a][test,b]*dphi[b][trial,a])
   */

  // Sum over quadrature points
  for (quad=0;quad<nq;quad++) {
    w = elm_data->quad_local->w[quad];

    //  Get the Basis Functions at each quadrature node
    for(i=0;i<dim;i++){
      get_FEM_basis_at_quadpt(elm_data, fe_data, i, quad);
    }

    // Loop over block rows (a) and block columns (b)
    for (a=0; a<dim; a++) {
      for (test=0; test<dof_per_elm_arr[a]; test++) {
        for (b=0; b<dim; b++) {
          for (trial=0; trial<dof_per_elm_arr[b]; trial++) {
            // Off-diagonal block term: 0.5 * dphi[a][test,b] * dphi[b][trial,a]
            aij = 0.5 * fe_data->dphi[a][test*dim+b] * fe_data->dphi[b][trial*dim+a];
            // Diagonal block term: 0.5 * sum_c dphi[a][test,c] * dphi[b][trial,c]
            if (a == b) {
              for (c=0; c<dim; c++)
                aij += 0.5 * fe_data->dphi[a][test*dim+c] * fe_data->dphi[b][trial*dim+c];
            }
            ALoc[(block_offset[a]+test)*dof_per_elm + (block_offset[b]+trial)] += w*aij;
          }
        }
      }
    }
  } // End quadrature

  return;
}
/******************************************************************************************************/

/* Ned_GradH1_RHS_local removed — was unused (called only from assemble_global_Ned_GradH1_RHS) */

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
