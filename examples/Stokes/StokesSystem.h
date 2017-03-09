/*! \file StokesSystem.h
 *
 *  Created by Peter Ohm on 1/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program solves Stokes PDE using finite elements
 *
 *      -laplace(u) + div(p) = f
 *                   grad(u) = 0
 *
 *
 *        in 2D or 3D
 *
 *        Along the boundary of the region, Dirichlet conditions are
 *        imposed for u and Neumann for p.  P2-P1 or P2-P0 can be used.
 *
 */

/********************** Local Assembly **************************************************/
void local_assembly_Stokes(REAL* ALoc, block_fespace *FE, trimesh *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
{

  /*!
   * \fn void local_assembly_Stokes(REAL* ALoc, block_fespace *FE, trimesh *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm, REAL time)
*
   * \brief Computes the local stiffness matrix for the Stokes system.
   *        For this problem we compute LHS of:
   *
   *        <grad u, grad v> - <p, div v> = <f, v>
   *                   - <div u, q> = 0
   *
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
   */

  // Loop indices
  INT i,j,idim;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm;

  INT dim = mesh->dim;

  // Loop Indices
  INT quad,test,trial;
  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(3,sizeof(REAL));

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u
    get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,dof_on_elm,mesh,FE->var_spaces[0]);
    local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
    get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);
    if(dim==3){
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);
    }
    // p
    local_dof_on_elm += FE->var_spaces[dim-1]->dof_per_elm;
    get_FEM_basis(FE->var_spaces[dim]->phi,FE->var_spaces[dim]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[dim]);
    
    // ux-vx block: <grad u, grad v>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij=0.0;
        for(idim=0;idim<dim;idim++){
          kij += FE->var_spaces[0]->dphi[test*dim+idim]*FE->var_spaces[0]->dphi[trial*dim+idim];
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // p-vx block: -<p, dx(ux)>
    local_row_index = 0;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
        kij = -FE->var_spaces[dim]->phi[trial]*(FE->var_spaces[0]->dphi[test*dim+0]);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // ux-q block: -<div u, q>
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij = -(FE->var_spaces[0]->dphi[trial*dim+0])*FE->var_spaces[dim]->phi[test];
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // uy-vy block <grad u, grad v>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij=0.0;
        for(idim=0;idim<dim;idim++){
          kij += FE->var_spaces[1]->dphi[test*dim+idim]*FE->var_spaces[1]->dphi[trial*dim+idim];
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // p-vy block: -<p, dy(uy)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_col_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
        kij = -FE->var_spaces[dim]->phi[trial]*(FE->var_spaces[1]->dphi[test*dim+1]);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // uy-q block: -<div u, q>
    local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
    if(dim==3) local_row_index += FE->var_spaces[2]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij = -(FE->var_spaces[1]->dphi[trial*dim+1])*FE->var_spaces[dim]->phi[test];
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    if(dim==3){
      // uz-vz block <grad u, grad v>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij=0.0;
          for(idim=0;idim<dim;idim++){
            kij += FE->var_spaces[2]->dphi[test*dim+idim]*FE->var_spaces[2]->dphi[trial*dim+idim];
          }
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // p-vz block: -<p, dz(uz)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[dim]->dof_per_elm;trial++){
          kij = -FE->var_spaces[dim]->phi[trial]*(FE->var_spaces[2]->dphi[test*dim+2]);
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // uz-q block: -<div u, q>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm + FE->var_spaces[2]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[dim]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = -(FE->var_spaces[2]->dphi[trial*dim+2])*FE->var_spaces[dim]->phi[test];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }
    }
  }

  if (qx) free(qx);

  return;
}
/*****************************************************************************************************/
