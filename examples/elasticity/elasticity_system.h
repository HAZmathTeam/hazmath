/*! \file BiotSystem.h
 *
 *  Created by Peter Ohm on 02/04/20.
 *  Copyright 2015_HAZMAT__. All rights reserved.
 *
 * \brief This program solves the following PDE using finite elements
 *
 *        2*mu*<eps(u),eps(v)> + lam*<div u,div v> - alpha*<p,div v> = <f1,v>
 *        alpha*<div u,q>                                            = 0
 *
 */

/********************** Local Assembly **************************************************/
void Elasticity_system(REAL* ALoc,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time)
{
  /*!
   * \fn void Elasticity_system(REAL* ALoc,block_fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT *v_on_elm,INT elm,REAL time)
   *
   * \brief Computes the local stiffness matrix for Elasticity
   *        For this problem we compute LHS of:
   *
   *        2*mu*<eps(u),eps(v)> + lam*<div u,div v> - alpha*<p,div v> = <f1,v>
   *        alpha*<div u,q>                                            = 0
   *
   *        where eps(u) = 0.5*(grad u + (grad u)^T) =
   *                       u1x            0.5*(u2x+u1y)  0.5*(u3x+u1z)
   *                       0.5*(u1y+u2x)  u2y            0.5*(u3y+u2z)
   *                       0.5*(u1z+u3x)  0.5*(u2z+u3y)  u3z
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
   * \note Assumes 2D only
   *
   */

  // Flag for erros
  SHORT status;

  // Loop Indices
  INT quad,test,trial,i,j;

  // Mesh and FE data
  INT dim = mesh->dim;
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;

  // FE Space names
  INT bbl = 0; // Bubble
  INT ds1 = 1; // Displacement x
  INT ds2 = 2; // Displacement y
  INT ds3 = 3; // Displacement z
  INT prs = dim+1; // Pressure

  // For each unknown, find the DOF
  INT* local_dof_on_elm = NULL;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Stiffness Matrix Entry
  REAL aij = 0.0;
  
  REAL divbbl = 0.0;
  REAL divbbl2= 0.0;

  // Data
  REAL* K = (REAL *) calloc(dim*dim,sizeof(REAL));
  REAL mu   =0;
  REAL mu_f =0;
  REAL lam  =0;
  REAL alpha=0;
  REAL tau  =1;

  INT local_row_index, local_col_index;

  if(dim==2) {

    //  Sum over quadrature points
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      w     = cq->w[elm*cq->nq_per_elm+quad];

      //  Get the Basis Functions at each quadrature node
      local_dof_on_elm = dof_on_elm;
      for( i=0; i<FE->nspaces; i++) {
        if(i>0) local_dof_on_elm = local_dof_on_elm + FE->var_spaces[i-1]->dof_per_elm;
        get_FEM_basis(FE->var_spaces[i]->phi,FE->var_spaces[i]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[i]);
      }

      // Data
      //conductivity2D(K,qx,time,NULL); // Actually the inverse of K// TODO
      get_mu(&mu,qx,time,NULL);
      mu_f = 1.0;//viscosity(&mu_f,qx,time,NULL);
      get_lam(&lam,qx,time,NULL);
      get_alpha(&alpha,qx,time,NULL);

      // b-b block: 2*mu*<eps(b),eps(b)> + lam*<div b,div b>
      // u-b block:
      // u-b block:
      // b-v block:
      // u-v block: 2*mu*<eps(u),eps(v)> + lam*<div u,div v>
      // p-v block: -alpha*<p,div v>
      // b-q block:
      // u-q block:
      // p-q block: 0

      // b-b block
      local_col_index = 0;
      local_row_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[bbl]->dof_per_elm; test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[bbl]->dof_per_elm; trial++) {
          divbbl = 0.0;
          divbbl2= 0.0;
          aij = 0.0;
          for (i=0; i<dim; i++) {
            divbbl  += FE->var_spaces[bbl]->dphi[trial*dim*dim + i*dim + i];
            divbbl2 += FE->var_spaces[bbl]->dphi[test*dim*dim + i*dim + i];
            for (j=0; j<dim; j++) {
              aij += 0.25*
                    (FE->var_spaces[bbl]->dphi[test*dim*dim+i*dim+j]+FE->var_spaces[bbl]->dphi[test*dim*dim+j*dim+i])*
                    (FE->var_spaces[bbl]->dphi[trial*dim*dim+i*dim+j]+FE->var_spaces[bbl]->dphi[trial*dim*dim+j*dim+i]);
            }
          }
          aij = 2*mu*aij + lam*(divbbl*divbbl2);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }
      
      // u1-b block: 2*mu*<eps(u1),eps(b)> + lam*<div1 b,div b>
      local_col_index += FE->var_spaces[bbl]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[bbl]->dof_per_elm; test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds1]->dof_per_elm; trial++) {
          divbbl = 0.0;
          aij = 0.0;
          j = 0; // x direction
          for (i=0; i<dim; i++) {
            divbbl += FE->var_spaces[bbl]->dphi[test*dim*dim + i*dim + i];
            aij += 0.5*(FE->var_spaces[bbl]->dphi[test*dim*dim+i*dim+j]*FE->var_spaces[ds1]->dphi[trial*dim+i]);
            aij += 0.5*(FE->var_spaces[bbl]->dphi[test*dim*dim+j*dim+i]*FE->var_spaces[ds1]->dphi[trial*dim+i]);
          }
          aij = 2*mu*aij + lam*(divbbl*FE->var_spaces[ds1]->dphi[trial*dim+j]);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u2-b block: 2*mu*<eps(u2),eps(b)> + lam*<div2 b,div b>
      local_col_index += FE->var_spaces[ds1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[bbl]->dof_per_elm; test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds2]->dof_per_elm; trial++) {
          divbbl = 0.0;
          aij = 0.0;
          j = 1; // y direction
          for (i=0; i<dim; i++) {
            divbbl += FE->var_spaces[bbl]->dphi[test*dim*dim + i*dim + i];
            aij += 0.5*(FE->var_spaces[bbl]->dphi[test*dim*dim+i*dim+j]*FE->var_spaces[ds2]->dphi[trial*dim+i]);
            aij += 0.5*(FE->var_spaces[bbl]->dphi[test*dim*dim+j*dim+i]*FE->var_spaces[ds2]->dphi[trial*dim+i]);
          }
          aij = 2*mu*aij + lam*(divbbl*FE->var_spaces[ds2]->dphi[trial*dim+j]);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // p-b block: -alpha*<p, div b>
      local_col_index += FE->var_spaces[ds2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[bbl]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[prs]->dof_per_elm; trial++) {
          aij = 0.0;
          for (i=0;i<dim;i++) {
            aij += (FE->var_spaces[bbl]->dphi[test*dim*dim+i*dim+i]*FE->var_spaces[prs]->phi[trial]);
          }
          aij = -alpha*aij;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }


      ////////////////////////////////////////////////////
      local_row_index += FE->var_spaces[bbl]->dof_per_elm;

      // b-v1 block: 2*mu*<eps(b),eps(v1)> + lam*<div b,div v1>
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds1]->dof_per_elm; test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[bbl]->dof_per_elm; trial++) {
          divbbl = 0.0;
          aij = 0.0;
          j = 0; // x direction
          for (i=0; i<dim; i++) {
            divbbl += FE->var_spaces[bbl]->dphi[trial*dim*dim + i*dim + i];
            aij += 0.5*(FE->var_spaces[bbl]->dphi[trial*dim*dim+i*dim+j]*FE->var_spaces[ds1]->dphi[test*dim+i]);
            aij += 0.5*(FE->var_spaces[bbl]->dphi[trial*dim*dim+j*dim+i]*FE->var_spaces[ds1]->dphi[test*dim+i]);
          }
          aij = 2*mu*aij + lam*(divbbl*FE->var_spaces[ds1]->dphi[test*dim+j]);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u1-v1 block: mu*(2*<u1x,v1x> + <u1y,v1y>) + lam*<u1x,v1x>
      local_col_index += FE->var_spaces[bbl]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds1]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds1]->dof_per_elm; trial++) {
          aij = mu*(2*FE->var_spaces[ds1]->dphi[trial*dim]*FE->var_spaces[ds1]->dphi[test*dim]
              + FE->var_spaces[ds1]->dphi[trial*dim+1]*FE->var_spaces[ds1]->dphi[test*dim+1])
              + lam*FE->var_spaces[ds1]->dphi[trial*dim]*FE->var_spaces[ds1]->dphi[test*dim];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u2-v1 block: mu*<u2x,v1y> + lam*<u2y,v1x>
      local_col_index += FE->var_spaces[ds1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds1]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds2]->dof_per_elm; trial++) {
          aij = mu*(FE->var_spaces[ds2]->dphi[trial*dim]*FE->var_spaces[ds1]->dphi[test*dim+1])
              + lam*FE->var_spaces[ds2]->dphi[trial*dim+1]*FE->var_spaces[ds1]->dphi[test*dim];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // p-v1 block: -alpha*<p,v1x>
      local_col_index += FE->var_spaces[ds2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds1]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[prs]->dof_per_elm; trial++) {
          aij = -alpha*FE->var_spaces[prs]->phi[trial]*FE->var_spaces[ds1]->dphi[test*dim];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }


      ////////////////////////////////////////////////////
      local_row_index += FE->var_spaces[ds1]->dof_per_elm;

      // b-v2 block: 2*mu*<eps(b),eps(v2)> + lam*<div b,div v2>
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds2]->dof_per_elm; test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[bbl]->dof_per_elm; trial++) {
          divbbl = 0.0;
          aij = 0.0;
          j = 1; // y direction
          for (i=0; i<dim; i++) {
            divbbl += FE->var_spaces[bbl]->dphi[trial*dim*dim + i*dim + i];
            aij += 0.5*(FE->var_spaces[bbl]->dphi[trial*dim*dim+i*dim+j]*FE->var_spaces[ds2]->dphi[test*dim+i]);
            aij += 0.5*(FE->var_spaces[bbl]->dphi[trial*dim*dim+j*dim+i]*FE->var_spaces[ds2]->dphi[test*dim+i]);
          }
          aij = 2*mu*aij + lam*(divbbl*FE->var_spaces[ds2]->dphi[test*dim+j]);
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u1-v2 block: mu*<u1y,v2x> + lam*<u1x,v2y>
      local_col_index += FE->var_spaces[bbl]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds2]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds1]->dof_per_elm; trial++) {
          aij = mu*(FE->var_spaces[ds1]->dphi[trial*dim+1]*FE->var_spaces[ds2]->dphi[test*dim])
              + lam*FE->var_spaces[ds1]->dphi[trial*dim]*FE->var_spaces[ds2]->dphi[test*dim+1];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u2-v2 block: mu*(<u2x,v2x> + 2*<u2y,v2y>) + lam*<u2y,v2y>
      local_col_index += FE->var_spaces[ds1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds2]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds2]->dof_per_elm; trial++) {
          aij = mu*(FE->var_spaces[ds2]->dphi[trial*dim]*FE->var_spaces[ds2]->dphi[test*dim]
              + 2*FE->var_spaces[ds2]->dphi[trial*dim+1]*FE->var_spaces[ds2]->dphi[test*dim+1])
              + lam*FE->var_spaces[ds2]->dphi[trial*dim+1]*FE->var_spaces[ds2]->dphi[test*dim+1];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // p-v2 block: -alpha*<p,v2y>
      local_col_index += FE->var_spaces[ds2]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[ds2]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[prs]->dof_per_elm; trial++) {
          aij = -alpha*FE->var_spaces[prs]->phi[trial]*FE->var_spaces[ds2]->dphi[test*dim+1];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }


      ////////////////////////////////////////////////////
      local_row_index += FE->var_spaces[ds2]->dof_per_elm;

      // b-q block: alpha*<div b,q>
      local_col_index = 0;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[prs]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[bbl]->dof_per_elm; trial++) {
          aij = 0.0;
          for (i=0;i<dim;i++) {
            aij += (FE->var_spaces[bbl]->dphi[trial*dim*dim+i*dim+i]*FE->var_spaces[prs]->phi[test]);
          }
          aij = alpha*aij;
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u1-q block: alpha*<u1x,q>
      local_col_index += FE->var_spaces[bbl]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[prs]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds1]->dof_per_elm; trial++) {
          aij = alpha*FE->var_spaces[ds1]->dphi[trial*dim]*FE->var_spaces[prs]->phi[test];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // u2-q block: alpha*<u2y,q>
      local_col_index += FE->var_spaces[ds1]->dof_per_elm;
      // Loop over Test Functions (Rows)
      for (test=0; test<FE->var_spaces[prs]->dof_per_elm;test++) {
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<FE->var_spaces[ds2]->dof_per_elm; trial++) {
          aij = alpha*FE->var_spaces[ds2]->dphi[trial*dim+1]*FE->var_spaces[prs]->phi[test];
          ALoc[(local_row_index+test)*dof_per_elm+(local_col_index+trial)] += w*aij;
        }
      }

      // p-q block: 0
      local_col_index += FE->var_spaces[ds2]->dof_per_elm;

    }
  } else if(dim==3) {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  if(qx) free(qx);
  if(K) free(K);
  return;
}
