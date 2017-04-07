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
 *
 */

#include "hazmath.h"

/******************************************************************************************************/
/*!
 * \fn void assemble_DuDv_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local stiffness matrix for coeff*<Du,Dv> = <f,v> bilinear form using various element types
 *        (eg. P1, P2 -> (grad u, grad v), Nedelec <curl u, curl v>, and Raviart-Thomas <div u, div v>).
 *
 *        For this problem we compute:
 *
 *        coeff*D*D u = f  ---->   coeff*<D u, D v> = <f,v>
 *
 *        which gives Ax = b,
 *
 *        A_ij = coeff*<D phi_j,D phi_i>
 *        b_i  = <f,phi_i>
 *
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return ALoc         Local Stiffness Matrix (Full Matrix)
 *
 * \note Assumes 2D or 3D only for Nedelec and Raviart-Thomas Elements
 *
 */
void assemble_DuDv_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  INT dim = mesh->dim;

  // Loop Indices
  INT quad,test,trial,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  // Vector Derivatives: Gradients (PX) and 3D Curls (3D Ned)
  if(FE->FEtype<20 || (FE->FEtype>=20 && FE->FEtype<30 && dim==3)) {

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij=0.0;
          for(idim=0;idim<dim;idim++)
            kij += coeff_val*(FE->dphi[test*dim+idim]*FE->dphi[trial*dim+idim]);
          ALoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else { // Scalar Derivatives: Divergence (RT) and 2D Curls (2D Ned)

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij = coeff_val*(FE->dphi[test]*FE->dphi[trial]);
          ALoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  }

  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void assemble_mass_local(REAL* MLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local mass matrix for coeff*<u,v> = <f,v> bilinear form using various element types
 *        (eg. P0, P1, P2, Nedelec, and Raviart-Thomas).
 *
 *        For this problem we compute:
 *
 *        coeff*u = f  ---->   coeff*<u,v> = <f,v>
 *
 *        which gives Mx = b,
 *
 *        M_ij = coeff*<phi_j,phi_i>
 *        b_i  = <f,phi_i>
 *
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return MLoc         Local Mass Matrix (Full Matrix)
 *
 * \note Assumes 2D or 3D only for Nedelec and Raviart-Thomas Elements
 *
 */
void assemble_mass_local(REAL* MLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  INT dim = mesh->dim;

  // Loop Indices
  INT quad,test,trial,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  // Vector Functions
  if(FE->FEtype>=20) {

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij=0.0;
          for(idim=0;idim<dim;idim++)
            kij += coeff_val*(FE->phi[test*dim+idim]*FE->phi[trial*dim+idim]);
          MLoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij = coeff_val*(FE->phi[test]*FE->phi[trial]);
          MLoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  }

  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void assemble_DuDvplusmass_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local stiffness matrix for coeff1*<Du,Dv> + coeff2*<u,v> = <f,v> bilinear form
 *        using various element types
 *        (eg. P1, P2 -> (grad u, grad v) + (u,v), Nedelec <curl u, curl v> + (u,v),
 *        and Raviart-Thomas <div u, div v>) + (u,v).
 *
 *        For this problem we compute:
 *
 *        coeff1*D*D u +coeff2 u = f ----> coeff1*<D u, D v> + coeff2*<u,v>= <f,v>
 *
 *        which gives Ax = b,
 *
 *        A_ij = coeff[0]*<D phi_j,D phi_i> + coeff[1]*<phi_j,phi_i>
 *        b_i  = <f,phi_i>
 *
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return ALoc         Local Stiffness Matrix (Full Matrix)
 *
 * \note Assumes 2D or 3D only for Nedelec and Raviart-Thomas Elements
 *
 */
void assemble_DuDvplusmass_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  INT dim = mesh->dim;

  // Error Check
  SHORT status=0;

  // Loop Indices
  INT quad,test,trial,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));
  // Stiffness Matrix Entry
  REAL kij;
  // Coefficient Value at Quadrature Nodes
  REAL coeff_val[2];

  // Scalar Functions and Vector Derivatives: PX
  if(FE->FEtype<20) {

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(coeff_val,qx,time);
      } else {
        coeff_val[0] = 1.0;
        coeff_val[1] = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij=coeff_val[1]*(FE->phi[test]*FE->phi[trial]);
          for(idim=0;idim<dim;idim++)
            kij += coeff_val[0]*(FE->dphi[test*dim+idim]*FE->dphi[trial*dim+idim]);
          ALoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else if (FE->FEtype>=20 && FE->FEtype<30) { // Nedelec (diff in 2D and 3D)

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(coeff_val,qx,time);
      } else {
        coeff_val[0] = 1.0;
        coeff_val[1] = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          if(dim==2) {
            kij = coeff_val[0]*(FE->dphi[test]*FE->dphi[trial]) +
                coeff_val[1]*(FE->phi[test*dim]*FE->phi[trial*dim] + FE->phi[test*dim+1]*FE->phi[trial*dim+1]);
          } else if (dim==3) {
            kij = coeff_val[0]*(FE->dphi[test*dim]*FE->dphi[trial*dim] +
                FE->dphi[test*dim+1]*FE->dphi[trial*dim+1] + FE->dphi[test*dim+2]*FE->dphi[trial*dim+2]) +
                coeff_val[1]*(FE->phi[test*dim]*FE->phi[trial*dim] + FE->phi[test*dim+1]*FE->phi[trial*dim+1] +
                FE->phi[test*dim+2]*FE->phi[trial*dim+2]);
          } else {
            status = ERROR_DIM;
            check_error(status, __FUNCTION__);
          }
          ALoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else if(FE->FEtype==30) { // RT elements

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      if(coeff!=NULL) {
        (*coeff)(coeff_val,qx,time);
      } else {
        coeff_val[0] = 1.0;
        coeff_val[1] = 1.0;
      }

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over Test Functions (Rows)
      for (test=0; test<FE->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->dof_per_elm; trial++) {
          kij = coeff_val[0]*(FE->dphi[test]*FE->dphi[trial]) +
              coeff_val[1]*(FE->phi[test*dim]*FE->phi[trial*dim] + FE->phi[test*dim+1]*FE->phi[trial*dim+1]);
          if(dim==3) kij += coeff_val[1]*FE->phi[test*dim+2]*FE->phi[trial*dim+2];
          ALoc[test*FE->dof_per_elm+trial] += w*kij;
        }
      }
    }
  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void assemble_symmetricDuDv_local(REAL* ALoc, block_fespace *FE, trimesh *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
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
 * \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,u3,p)
 *
 * \note Assumes 2D or 3D only. Only works for Scalar Elements
 *
 */
void assemble_symmetricDuDv_local(REAL* ALoc, block_fespace *FE, trimesh *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{
  
  // Loop indices
  INT i,j,idim,quad,test,trial;

  //Mesh and FE data
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm;
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL * qx = (REAL *) calloc(dim,sizeof(REAL));

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for(quad=0;quad<cq->nq_per_elm;quad++){
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    // Get the Basis Functions at each quadrature node
    local_dof_on_elm = dof_on_elm;
    for(i=0;i<FE->nspaces;i++){
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,dof_on_elm,mesh,FE->var_spaces[i]);
      // Shift local dof to next finite element space
      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
    }

    // Assemble blocks
    local_row_index = 0;
    // Loop over block rows
    for(i=0;i<FE->nspaces;i++){
      // Reset local column indexing
      local_col_index = 0;
      // Loop over block columns
      for(j=0;j<FE->nspaces;j++){
        if(i==j){
          // Loop over Test Functions (Rows)
          for(test=0;test<FE->var_spaces[i]->dof_per_elm;test++){
            // Loop over Trial Functions (Columns)
            for(trial=0;trial<FE->var_spaces[j]->dof_per_elm;trial++){
               kij = 2*0.5*(FE->var_spaces[i]->dphi[test*dim+i]*FE->var_spaces[j]->dphi[trial*dim+j]);
               for(idim=0;idim<dim;idim++){
                 kij += 2*0.5*(FE->var_spaces[i]->dphi[test*dim+idim]*FE->var_spaces[j]->dphi[trial*dim+idim]);
               }
               ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
            }
          }
        } else {
          // Loop over Test Functions (Rows)
          for(test=0;test<FE->var_spaces[i]->dof_per_elm;test++){
            // Loop over Trial Functions (Columns)
            for(trial=0;trial<FE->var_spaces[j]->dof_per_elm;trial++){
               kij = 2*0.5*(FE->var_spaces[i]->dphi[test*dim+j]*FE->var_spaces[j]->dphi[trial*dim+i]);
               ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
            }
          }
        }
        // Update local column indexing
        local_col_index += FE->var_spaces[i]->dof_per_elm;
      }
    // Update local row indexing
    local_row_index += FE->var_spaces[i]->dof_per_elm;
    }

  }

  // Free
  if(qx) free(qx);

  return;
}
/******************************************************************************************************/

/****** Boundary Assemblies *******************/
/******************************************************************************************************/
/*!
 * \fn void boundary_mass_local(REAL* MLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq, \
                       INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local weak formulation of the mass matrix on a boundary face (3D -> 2D surface; 2D -> 1D curve)
 *         For this problem we compute the left-hand side of:
 *
 *         <u,v>_bdryobstacle    for all v
 *
 * \param FE            FE Space
 * \param old_sol       Old solution on FE space (not used here).
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_f      Specific DOF on the given face
 * \param dof_on_elm    Specific DOF on the given element
 * \param v_on_elm      Specific vertices on the given element
 * \param face          Current face
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return MLoc         Local Boundary Matrix (Full Matrix)
 *
 * \note                Assuming 2D and 3D only
 *
 */
void boundary_mass_local(REAL* MLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq, \
                         INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_f = 0;

  // flag for errors
  SHORT status;

  // Loop Indices
  INT j,quad,test,trial,doft,dofb;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));

  // Stiffness Matrix Entry
  REAL kij;


  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  if(FE->FEtype>=0 && FE->FEtype<10) { // PX elements

    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = FE->FEtype+1;
    } else if(dim==3) {
      dof_per_f = 3*FE->FEtype;
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }

    //  Sum over midpoints of edges
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      PX_H1_basis(FE->phi,FE->dphi,qx,dof_on_elm,FE->FEtype,mesh);

      // Loop over Test Functions (Rows - vertices)
      for (test=0; test<dof_per_f;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_f; trial++) {
          // Make sure ordering for global matrix is right
          for(j=0;j<FE->dof_per_elm;j++) {
            if(dof_on_f[test]==dof_on_elm[j]) {
              doft = j;
            }
            if(dof_on_f[trial]==dof_on_elm[j]) {
              dofb = j;
            }
          }
          kij = coeff_val*(FE->phi[dofb]*FE->phi[doft]);
          MLoc[test*dof_per_f+trial]+=w*kij;
        }
      }
    }
  } else if(FE->FEtype==20) {

    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = 1;
    } else if(dim==3) {
      dof_per_f = 3;
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      ned_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh);

      // Loop over Test Functions (Rows - edges)
      for (test=0; test<dof_per_f;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_f; trial++) {
          // Make sure ordering for global matrix is right
          for(j=0;j<FE->dof_per_elm;j++) {
            if(dof_on_f[test]==dof_on_elm[j]) {
              doft = j;
            }
            if(dof_on_f[trial]==dof_on_elm[j]) {
              dofb = j;
            }
          }
          kij = coeff_val*(FE->phi[dofb*dim]*FE->phi[doft*dim] + FE->phi[dofb*dim+1]*FE->phi[doft*dim+1]);
          if(dim==3) kij += coeff_val*(FE->phi[dofb*dim+2]*FE->phi[doft*dim+2]);
          MLoc[test*dof_per_f+trial]+=w*kij;
        }
      }
    }
  } else if(FE->FEtype==30) { // Raviart-Thomas elements
    // Get DOF Per Face
    if(dim==2) {
      dof_per_f = 1;
    } else if(dim==3) {
      dof_per_f = 1;
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[face*cq->nq_per_elm+quad];
      qx[1] = cq->y[face*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[face*cq->nq_per_elm+quad];
      w = cq->w[face*cq->nq_per_elm+quad];

      if(coeff!=NULL) {
        (*coeff)(&coeff_val,qx,time);
      } else {
        coeff_val = 1.0;
      }

      //  Get the Basis Functions at each quadrature node
      rt_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh);

      /// Loop over Test Functions (Rows - edges)
      for (test=0; test<dof_per_f;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<dof_per_f; trial++) {
          // Make sure ordering for global matrix is right
          for(j=0;j<FE->dof_per_elm;j++) {
            if(dof_on_f[test]==dof_on_elm[j]) {
              doft = j;
            }
            if(dof_on_f[trial]==dof_on_elm[j]) {
              dofb = j;
            }
          }
          kij = coeff_val*(FE->phi[dofb*dim]*FE->phi[doft*dim] + FE->phi[dofb*dim+1]*FE->phi[doft*dim+1]);
          kij += coeff_val*(FE->phi[dofb*dim+2]*FE->phi[doft*dim+2]);
          MLoc[test*dof_per_f+trial]+=w*kij;
        }
      }
    }
  } else {
    status = ERROR_FE_TYPE;
    check_error(status, __FUNCTION__);
  }

  if(qx) free(qx);

  return;
}
/******************************************************************************************************/

/***** RHS Routines *********************************/

/******************************************************************************************************/
/*!
 * \fn void FEM_RHS_Local(REAL* bLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local Right hand side vector for Galerkin Finite Elements
 *        b_i  = <f,phi_i>
 *
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific Vertices on element
 * \param elm           Current element
 * \param rhs           Function that gives RHS
 * \param time          Physical Time if time dependent
 *
 * \return bLoc         Local RHS Vector
 *
 * \note Assumes 2D or 3D only for Nedelec and Raviart-Thomas Elements
 *
 */
void FEM_RHS_Local(REAL* bLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
  // Mesh and FE data
  INT dim = mesh->dim;

  // Loop Indices
  INT quad,test,idim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Right-hand side function at Quadrature Nodes
  REAL* rhs_val=NULL;

  if(FE->FEtype<20) { // Scalar Functions

    rhs_val = (REAL *) calloc(1,sizeof(REAL));

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(dim==2 || dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*rhs)(rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over test functions and integrate rhs
      for (test=0; test<FE->dof_per_elm;test++) {
        bLoc[test] += w*rhs_val[0]*FE->phi[test];
      }
    }
  } else { // Vector Functions

    rhs_val = (REAL *) calloc(dim,sizeof(REAL));

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];
      (*rhs)(rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

      // Loop over test functions and integrate rhs
      for (test=0; test<FE->dof_per_elm;test++) {
        for(idim=0;idim<dim;idim++) {
          bLoc[test] += w*(rhs_val[idim]*FE->phi[test*dim+idim]);
        }
      }
    }
  }

  if(qx) free(qx);
  if(rhs_val) free(rhs_val);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void FEM_RHS_Local(REAL* bLoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
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
 * \return bLoc         Local RHS Vector
 *
 *
 */
void FEM_Block_RHS_Local(REAL* bLoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
  // Loop Indices
  INT i,quad,test;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_elm = 0;
  INT nun=0;

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
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));
  // Stiffness Matrix Entry

  // Right-hand side function at Quadrature Nodes
  REAL* rhs_val= (REAL *) calloc(nun,sizeof(REAL));

  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    (*rhs)(rhs_val,qx,time);

    local_row_index=0;
    unknown_index=0;
    local_dof_on_elm=dof_on_elm;
    for(i=0;i<FE->nspaces;i++) {

      // Basis Functions and its derivatives if necessary
      get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);

      // Loop over test functions and integrate rhs
      if(FE->var_spaces[i]->FEtype<20) { // Scalar Element
        for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
          bLoc[(local_row_index+test)] += w*rhs_val[unknown_index]*FE->var_spaces[i]->phi[test];
        }
        unknown_index++;

      } else { // Vector Element
        for (test=0; test<FE->var_spaces[i]->dof_per_elm;test++) {
          bLoc[(local_row_index+test)] += w*(rhs_val[unknown_index]*FE->var_spaces[i]->phi[test*dim] +
              rhs_val[unknown_index+1]*FE->var_spaces[i]->phi[test*dim+1]);
          if(dim==3) bLoc[(local_row_index+test)] += w*rhs_val[unknown_index+2]*FE->var_spaces[i]->phi[test*dim+2];
        }
        unknown_index += dim;
      }

      local_dof_on_elm += FE->var_spaces[i]->dof_per_elm;
      local_row_index += FE->var_spaces[i]->dof_per_elm;
    }
  }

  if(qx) free(qx);
  if(rhs_val) free(rhs_val);
  return;
}

/******************************************************************************************************/
/*!
 * \fn void Ned_GradH1_RHS_local(REAL* bLoc,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,INT *ed_on_elm,INT *v_on_elm,INT elm,dvector* u)
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
void Ned_GradH1_RHS_local(REAL* bLoc,fespace *FE_H1,fespace *FE_Ned,trimesh *mesh,qcoordinates *cq,INT *ed_on_elm,INT *v_on_elm,INT elm,dvector* u)  
{
  // Mesh and FE data
  INT dim = mesh->dim;
  
  // Loop Indices
  INT quad,test;
  
  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));
  
  // Right-hand side function at Quadrature Nodes
  REAL* ucoeff = (REAL *) calloc(dim,sizeof(REAL));
  
  //  Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];
    
    // Get FEM function at quadrature nodes
    FE_Interpolation(ucoeff,u->val,qx,ed_on_elm,v_on_elm,FE_Ned,mesh);
    
    //  Get the Basis Functions at each quadrature node
    PX_H1_basis(FE_H1->phi,FE_H1->dphi,qx,v_on_elm,FE_H1->FEtype,mesh);

    // Loop over test functions and integrate rhs
    for (test=0; test<FE_H1->dof_per_elm;test++) {
      bLoc[test] += w*(ucoeff[0]*FE_H1->dphi[test]+ucoeff[1]*FE_H1->dphi[test]);
      if(dim==3) bLoc[test] += w*ucoeff[2]*FE_H1->dphi[test];
    }
  }

  if(qx) free(qx);
  if(ucoeff) free(ucoeff);

  return;
}
/******************************************************************************************************/

/****** Boundary Assemblies *********************/
/******************************************************************************************************/
/*!
* \fn void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
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
* \note All matrices are assumed to be indexed at 1 in the CSR formatting.
* \note Assumes different type of integral for different Element type:
*       Scalar -> <f,v>_bdry
*       Vector -> <f,n*v>_bdry
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
void FEM_RHS_Local_face(REAL* bLoc,dvector* old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_f,INT *dof_on_elm,INT *v_on_elm,INT dof_per_face,INT face,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
{
  // Mesh and FE data
  INT dim = mesh->dim;

  // Loop Indices
  INT j,quad,test,doft;

  // Quadrature Weights and Nodes
  qcoordinates *cq_face = allocateqcoords(cq->nq1d,1,dim);
  quad_edgeface(cq_face,mesh,cq->nq1d,face,dim-1);
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));
  REAL w;

  // Get normal vector components on face if needed
  REAL nx=0.0,ny=0.0,nz=0.0;
  if(FE->FEtype>=20) {
    nx = mesh->f_norm[face*dim];
    ny = mesh->f_norm[face*dim+1];
    if(dim==3) nz = mesh->f_norm[face*dim+2];
  }

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Value of f
  REAL rhs_val;

  if(FE->FEtype<20) { // Scalar Functions

    //  Sum over quadrature points
    for (quad=0;quad<cq_face->n;quad++) {
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

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
      qx[0] = cq_face->x[quad];
      qx[1] = cq_face->y[quad];
      if(dim==3)
        qx[2] = cq_face->z[quad];
      w = cq_face->w[quad];
      (*rhs)(&rhs_val,qx,time);

      //  Get the Basis Functions at each quadrature node
      get_FEM_basis(FE->phi,FE->dphi,qx,v_on_elm,dof_on_elm,mesh,FE);

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

  if(qx) free(qx);
  return;
}
/******************************************************************************************************/

/****** Special Assemblies *********************/

/******************************************************************************************************/
/*!
 * \fn void impedancebdry_local(REAL* ZLoc,dvector *old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *ed_on_f, \
                       INT *ed_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
 *
 * \brief Computes the local weak formulation of the Impedance boundary condition for Maxwell's Equations
 *         Uses midpoint rule to integrate on edges of boundary face
 *         For this problem we compute the left-hand side of:
 *
 *         <n x E,n x F>_bdryobstacle    for all F in H_imp(curl) (Nedelec)
 *
 * \param old_sol       Solution at previous step (not needed in this function)
 * \param FE            FE Space
 * \param mesh          Mesh Data
 * \param cq            Quadrature Nodes
 * \param ed_on_f       Specific edges on the given face
 * \param ed_on_elm     Specific edges on the given element
 * \param v_on_elm      Specific vertices on the given element
 * \param face          Current face
 * \param elm           Current element
 * \param coeff         Function that gives coefficient (for now assume constant)
 * \param time          Physical Time if time dependent
 *
 * \return ZLoc         Local Boundary Matrix (Full Matrix)
 *
 * \note                ASSUMING 3D ONLY
 *
 */
void impedancebdry_local(REAL* ZLoc,dvector *old_sol,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *ed_on_f, \
                         INT *ed_on_elm,INT *v_on_elm,INT face,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time)
{
  // Mesh and FE data
  INT ed_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;

  // Loop Indices
  INT j,quad,test,trial,ed,edt,edb;

  // Quadrature Weights and Nodes
  INT nq = 2*dim-3; // = ed_per_face
  REAL* qx = (REAL *) calloc(nq,sizeof(REAL));
  // 3D: Using triangle midpoint rule, so qx is midpoint of edges and w is |F|/3
  REAL w = mesh->f_area[face]/3.0;

  // Get normal vector components on face
  REAL nx = mesh->f_norm[face*dim];
  REAL ny = mesh->f_norm[face*dim+1];
  REAL nz = mesh->f_norm[face*dim+2];

  // Stiffness Matrix Entry
  REAL kij,kij1,kij2,kij3,kij4,kij5,kij6;

  // Basis Functions and its curl
  REAL* phi= (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));
  REAL* cphi = (REAL *) calloc(ed_per_elm*dim,sizeof(REAL));

  // Coefficient Value at Quadrature Nodes
  REAL coeff_val=0.0;

  //  Sum over midpoints of edges
  for (quad=0;quad<nq;quad++) {
    ed = ed_on_f[quad]-1;
    qx[0] = mesh->ed_mid[ed*dim];
    qx[1] = mesh->ed_mid[ed*dim+1];
    qx[2] = mesh->ed_mid[ed*dim+2];

    if(coeff!=NULL) {
      (*coeff)(&coeff_val,qx,time);
    } else {
      coeff_val = 1.0;
    }

    //  Get the Basis Functions at each quadrature node
    ned_basis(phi,cphi,qx,v_on_elm,ed_on_elm,mesh);

    // Loop over Test Functions (Rows - edges)
    for (test=0; test<nq;test++) {
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<nq; trial++) {
        // Make sure ordering for global matrix is right
        for(j=0;j<mesh->ed_per_elm;j++) {
          if(ed_on_f[test]==ed_on_elm[j]) {
            edt = j;
          }
          if(ed_on_f[trial]==ed_on_elm[j]) {
            edb = j;
          }
        }
        kij1 = phi[edb*dim+1]*nz - phi[edb*dim+2]*ny;
        kij2 = phi[edb*dim+2]*nx - phi[edb*dim]*nz;
        kij3 = phi[edb*dim]*ny - phi[edb*dim+1]*nx;
        kij4 = phi[edt*dim+1]*nz - phi[edt*dim+2]*ny;
        kij5 = phi[edt*dim+2]*nx - phi[edt*dim]*nz;
        kij6 = phi[edt*dim]*ny - phi[edt*dim+1]*nx;
        kij = coeff_val*(kij1*kij4+kij2*kij5+kij3*kij6);
        ZLoc[test*nq+trial]+=w*kij;
      }
    }
  }

  if (phi) free(phi);
  if(cphi) free(cphi);
  if(qx) free(qx);
  return;
}
/******************************************************************************************************/


