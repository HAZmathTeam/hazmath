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
 * \fn void meanzero_pressure(block_dCSR* A, dvector* b, scomplex* sc, fespace* FEp, qcoordinates* cq)
 *
 * \brief Deals with pressure singularity by adding constraint that
 *         int p = int p_true (0 if mean zero constraint)
 *        Do this by adding a row and column to the p-p block that computes this
 *        integral.  Note that this increases the DoF by 1 and the rhs by 1.
 *
 * \param A            Stokes system matrix
 * \param b            System RHS
 * \param sc           Simplicial complex
 * \param FEp          FE space for pressure
 * \param cq           Quadrature needed for anything but P0
 *
 * \return A, b        Modified system
 *
 */
void meanzero_pressure(block_dCSRmat* A, dvector* b, scomplex* sc, fespace* FEp, qcoordinates* cq) {

  INT ndof = FEp->ndof;
  INT i, newrows;

  // Reallocate rhs appropriately and add in value of mean if not zero
  newrows = b->row + 1;

  // Copy current values into new longer vector
  REAL* btmp = (REAL*)calloc(newrows, sizeof(REAL));
  for (i = 0; i < b->row; i++) btmp[i] = b->val[i];
  // Get mean value
  REAL pmeanval = 0.0;
  pmean(&pmeanval);
  btmp[b->row] = pmeanval;
  b->val = (REAL*)realloc(b->val, sizeof(REAL) * (newrows));

  // Copy btmp back into b
  b->row++;
  for (i = 0; i < newrows; i++) b->val[i] = btmp[i];

  // Reallocate matrix
  // In P0 -> Add el_vol for each entry
  if (FEp->FEtype == 0) {
    bdcsr_extend(A, sc->fem->el_vol, sc->fem->el_vol, sc->dim, 1.0, 1.0);
  } else {
    // In P1 or higher, the extra row is Mp*1, where Mp is the mass matrix for
    // the pressure space and 1 is the vector of ones.
    dCSRmat Mp = dcsr_create(0,0,0);
    assemble_global_single(&Mp, NULL, FEp, sc, cq,
                           local_assembly_mass, NULL, NULL, one_coeff_scal, 0.0);
    dvector* ones = dvec_create_p(ndof);
    dvec_set(ndof, ones, 1.0);
    REAL* constraint = (REAL*)calloc(ndof, sizeof(REAL));
    dcsr_mxv(&Mp, ones->val, constraint);
    bdcsr_extend(A, constraint, constraint, sc->dim, 1.0, 1.0);

    dcsr_free(&Mp);
    if (ones) free(ones);
    if (constraint) free(constraint);
  }

  if (btmp) free(btmp);

  return;
}

/*!
 * \fn void local_assembly_Stokes(REAL* ALoc, simplex_local_data *elm_data, fe_local_data *fe_data, REAL time)
 *
 * \brief Computes the local stiffness matrix for the Stokes system using local data structs.
 *        Computes the local stiffness matrix for the Stokes system.
 *
 */
void local_assembly_Stokes(REAL* ALoc, simplex_local_data *elm_data, fe_local_data *fe_data, REAL time) {

  // Loop indices
  INT i, j, idim, quad, test, trial;

  // Mesh and FE data
  INT dim = elm_data->dim;
  INT nspaces = fe_data->nspaces;
  INT nq = elm_data->quad_local->nq;

  // Space indices
  INT iu1 = 0;
  INT iu2 = 1;
  INT iu3 = 2;
  INT ip = dim;

  // Keep track of DoF per element
  INT dof_per_elm = fe_data->n_dof;
  INT u1dofpelm = fe_data->n_dof_per_space[iu1];
  INT u2dofpelm = fe_data->n_dof_per_space[iu2];
  INT u3dofpelm;
  if (dim == 3) u3dofpelm = fe_data->n_dof_per_space[iu3];
  INT pdofpelm = fe_data->n_dof_per_space[ip];

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry
  REAL kij = 0.0;

  // Test and Trial Functions
  REAL u1, u1x, u1y, u1z, u2, u2x, u2y, u2z, u3, u3x, u3y, u3z, p;
  REAL v1, v1x, v1y, v1z, v2, v2x, v2y, v2z, v3, v3x, v3y, v3z, q;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Sum over quadrature points
  for (quad = 0; quad < nq; quad++) {
    qx[0] = elm_data->quad_local->x[quad * dim + 0];
    qx[1] = elm_data->quad_local->x[quad * dim + 1];
    if (dim == 3) qx[2] = elm_data->quad_local->x[quad * dim + 2];
    w = elm_data->quad_local->w[quad];

    //  Get the Basis Functions at each quadrature node
    for (i = 0; i < nspaces; i++) {
      get_FEM_basis_at_quadpt(elm_data, fe_data, i, quad);
    }

    // Loop over block rows of test functions
    local_row_index = 0;

    // v1 block row:
    for (test = 0; test < u1dofpelm; test++) {
      v1 = fe_data->phi[iu1][test];
      v1x = fe_data->dphi[iu1][test * dim];
      v1y = fe_data->dphi[iu1][test * dim + 1];
      if (dim == 3) v1z = fe_data->dphi[iu1][test * dim + 2];

      // u1-v1 block: 2*<u1x,v1x> + <u1y,v1y>
      local_col_index = 0;
      for (trial = 0; trial < u1dofpelm; trial++) {
        u1x = fe_data->dphi[iu1][trial * dim];
        u1y = fe_data->dphi[iu1][trial * dim + 1];
        kij = 2.0 * u1x * v1x + u1y * v1y;
        if (dim == 3) {
          u1z = fe_data->dphi[iu1][trial * dim + 2];
          kij += u1z * v1z;
        }
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block: <u2x,v1y>
      for (trial = 0; trial < u2dofpelm; trial++) {
        u2x = fe_data->dphi[iu2][trial * dim];
        u2y = fe_data->dphi[iu2][trial * dim + 1];
        kij = u2x * v1y;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block: <u3x,v1z>
      if (dim == 3) {
        for (trial = 0; trial < u3dofpelm; trial++) {
          u3x = fe_data->dphi[iu3][trial * dim];
          u3z = fe_data->dphi[iu3][trial * dim + 2];
          kij = u3x * v1z;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v1 block: -<p,v1x>
      for (trial = 0; trial < pdofpelm; trial++) {
        p = fe_data->phi[ip][trial];
        kij = -p * v1x;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
    } // End v1 block
    local_row_index += u1dofpelm;

    // v2 block row
    for (test = 0; test < u2dofpelm; test++) {
      v2 = fe_data->phi[iu2][test];
      v2x = fe_data->dphi[iu2][test * dim];
      v2y = fe_data->dphi[iu2][test * dim + 1];
      if (dim == 3) v2z = fe_data->dphi[iu2][test * dim + 2];

      local_col_index = 0;

      // u1-v2 block: <u1y,v2x>
      for (trial = 0; trial < u1dofpelm; trial++) {
        u1x = fe_data->dphi[iu1][trial * dim];
        u1y = fe_data->dphi[iu1][trial * dim + 1];
        kij = u1y * v2x;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block: <u2x,v2x> + 2*<u2y,v2y>
      for (trial = 0; trial < u2dofpelm; trial++) {
        u2x = fe_data->dphi[iu2][trial * dim];
        u2y = fe_data->dphi[iu2][trial * dim + 1];
        kij = u2x * v2x + 2.0 * u2y * v2y;
        if (dim == 3) {
          u2z = fe_data->dphi[iu2][trial * dim + 2];
          kij += u2z * v2z;
        }
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block: <u3y,v2z>
      if (dim == 3) {
        for (trial = 0; trial < u3dofpelm; trial++) {
          u3y = fe_data->dphi[iu3][trial * dim + 1];
          u3z = fe_data->dphi[iu3][trial * dim + 2];
          kij = u3y * v2z;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v2 block: -<p,v2y>
      for (trial = 0; trial < pdofpelm; trial++) {
        p = fe_data->phi[ip][trial];
        kij = -p * v2y;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
    } // End v2 block
    local_row_index += u2dofpelm;

    // v3 block row
    if (dim == 3) {
      for (test = 0; test < u3dofpelm; test++) {
        v3 = fe_data->phi[iu3][test];
        v3x = fe_data->dphi[iu3][test * dim];
        v3y = fe_data->dphi[iu3][test * dim + 1];
        v3z = fe_data->dphi[iu3][test * dim + 2];

        local_col_index = 0;

        // u1-v3 block: <dz(u1),dx(v3)>
        for (trial = 0; trial < u1dofpelm; trial++) {
          u1x = fe_data->dphi[iu1][trial * dim];
          u1z = fe_data->dphi[iu1][trial * dim + 2];
          kij = u1z * v3x;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u1dofpelm;

        // u2-v3 block: <dz(u2),dy(v3)>
        for (trial = 0; trial < u2dofpelm; trial++) {
          u2y = fe_data->dphi[iu2][trial * dim + 1];
          u2z = fe_data->dphi[iu2][trial * dim + 2];
          kij = u2z * v3y;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u2dofpelm;

        // u3-v3 block: <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)>
        for (trial = 0; trial < u3dofpelm; trial++) {
          u3x = fe_data->dphi[iu3][trial * dim];
          u3y = fe_data->dphi[iu3][trial * dim + 1];
          u3z = fe_data->dphi[iu3][trial * dim + 2];
          kij = u3x * v3x + u3y * v3y + 2.0 * u3z * v3z;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u3dofpelm;

        // p-v3 block: -<p, dz(v3)>
        for (trial = 0; trial < pdofpelm; trial++) {
          p = fe_data->phi[ip][trial];
          kij = -p * v3z;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
      } // End v3 block
      local_row_index += u3dofpelm;
    } // End dim if

    // q block row
    for (test = 0; test < pdofpelm; test++) {
      q = fe_data->phi[ip][test];
      local_col_index = 0;

      // u1-q block: -<dx(u1), q>
      for (trial = 0; trial < u1dofpelm; trial++) {
        u1x = fe_data->dphi[iu1][trial * dim];
        kij = -u1x * q;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u1dofpelm;

      // u2-q block: -<dy(u2), q>
      for (trial = 0; trial < u2dofpelm; trial++) {
        u2y = fe_data->dphi[iu2][trial * dim + 1];
        kij = -u2y * q;
        ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
      }
      local_col_index += u2dofpelm;

      if (dim == 3) {
        // u3-q block: -<dz(u3), q>
        for (trial = 0; trial < u3dofpelm; trial++) {
          u3z = fe_data->dphi[iu3][trial * dim + 2];
          kij = -u3z * q;
          ALoc[(local_row_index + test) * dof_per_elm + (local_col_index + trial)] += w * kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-q block: 0
    } // End q block
  }

  return;
}

/* Adapter for assemble_global_system: wraps local_assembly_Stokes + FEM_Block_RHS_Local */
void local_assembly_Stokes_unified(REAL *ALoc, REAL *bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *), REAL time)
{
  if(ALoc!=NULL)
    local_assembly_Stokes(ALoc, elm_data, fe_data, time);
  if(bLoc!=NULL && rhs!=NULL)
    FEM_Block_RHS_Local(ALoc, bLoc, u_local, elm_data, fe_data, rhs, coeff, time);
}
