/*! \file examples/stokes/stokes_system.h
 *
 *  Created by Peter Ohm on 1/5/17.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routines
 *        for the Stokes example.
 *
 * \note Updated by James Adler on 9/16/19
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
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm;
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions at each quadrature node
    // u = (u1,u2,u3) and v = (v1,v2,v3)
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

    // u1-v1 block: 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)>
    local_row_index = 0;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[0]->dphi[test*dim]*FE->var_spaces[0]->dphi[trial*dim];
        for(idim=0;idim<dim;idim++){
          kij += 1.0*FE->var_spaces[0]->dphi[test*dim+idim]*FE->var_spaces[0]->dphi[trial*dim+idim];
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // u2-v1 block <dx(u2),dy(v1)>>
    local_row_index = 0;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[0]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+0];
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // p-v1 block: -<p, dx(u1)>
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

    // u1-v2 block <dy(u1),dx(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[1]->dphi[test*dim+0]*FE->var_spaces[0]->dphi[trial*dim+1];
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // u2-v2 block <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)>
    local_row_index = FE->var_spaces[0]->dof_per_elm;
    local_col_index = FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        kij = 1.0*FE->var_spaces[1]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+1];
        for(idim=0;idim<dim;idim++){
          kij += 1.0*FE->var_spaces[1]->dphi[test*dim+idim]*FE->var_spaces[1]->dphi[trial*dim+idim];
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
    }

    // p-v2 block: -<p, dy(uy)>
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

    // u1-q block: -<dx(u1), q>
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

    // u2-q block: -<dy(u2), q>
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
      // u3-v3 block <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij=FE->var_spaces[2]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+2];
          for(idim=0;idim<dim;idim++){
            kij += FE->var_spaces[2]->dphi[test*dim+idim]*FE->var_spaces[2]->dphi[trial*dim+idim];
          }
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // p-v3 block: -<p, dz(u3)>
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

      // u3-q block: -<dz(u3), q>
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

      // u1-v3 block <dz(u1),dx(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
          kij = FE->var_spaces[2]->dphi[test*dim+0]*FE->var_spaces[0]->dphi[trial*dim+2];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u3-v1 block <dx(u3),dz(v1)>
      local_row_index = 0;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = FE->var_spaces[0]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+0];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u2-v3 block <dz(u2),dy(v3)>
      local_row_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
          kij = FE->var_spaces[2]->dphi[test*dim+1]*FE->var_spaces[1]->dphi[trial*dim+2];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }

      // u3-v2 block <dy(u3),dz(v2)>
      local_row_index = FE->var_spaces[0]->dof_per_elm;
      local_col_index = FE->var_spaces[0]->dof_per_elm + FE->var_spaces[1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++){
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
          kij = FE->var_spaces[1]->dphi[test*dim+2]*FE->var_spaces[2]->dphi[trial*dim+1];
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
      }
    }
  }

  return;
}
/*****************************************************************************************************/
