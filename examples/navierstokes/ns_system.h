/*! \file examples/navierstokes/ns_system.h
*
*  Created by Peter Ohm on 1/5/17.
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the local assembly routines
*        for the Navier-Stokes example.
*
* \note Updated by James Adler on 02/04/21
*/

/*!
* \fn void local_assembly_NS(REAL *ALoc, REAL *bLoc, REAL *u_local, simplex_local_data *elm_data, fe_local_data *fe_data, void (*rhs)(...), void (*coeff)(...), REAL time)
*
* \brief Computes the local stiffness matrix for the linearized Navier-Stokes system
*        using local data structs.
*
*            <2*eps(u), eps(v)> + Re<u0*grad(u) + u*grad(u0),v> - <p, div v> = <f, v> - <2*eps(u0), eps(v)> - Re<u0*grad(u0),v> + <p0, div v>
*              - <div u, q> = <div u0, q>
* where u0 is the previous Newton iterate.
*
* \param ALoc         Local Stiffness Matrix (output)
* \param bLoc         Local RHS vector (output)
* \param u_local      Previous solution restricted to this element
* \param elm_data     Local simplex/mesh data
* \param fe_data      Local FE data
* \param rhs          RHS function
* \param coeff        Coefficient function (unused)
* \param time         Physical time
*
* \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,u3,p)
*
* \note Assumes 2D or 3D only
*
*/
void local_assembly_NS(REAL *ALoc,REAL* bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop indices
  INT i,j,d,idim,quad,test,trial;

  // Mesh and FE data
  INT dim = elm_data->dim;
  INT nspaces = fe_data->nspaces;

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
  if(dim==3) {
    u3dofpelm = fe_data->n_dof_per_space[iu3];
  }
  INT pdofpelm = fe_data->n_dof_per_space[ip];

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix Entry and RHS entry
  REAL kij = 0.0;
  REAL rij = 0.0;

  // Stuff for previous solution
  REAL u1_prev;
  REAL u2_prev;
  REAL u3_prev;
  REAL p_prev;
  REAL du1_prev[dim];
  REAL du2_prev[dim];
  REAL du3_prev[dim];
  REAL divuk,epsu11,epsu21,epsu31,epsu22,epsu32,epsu33;
  REAL ukgraduk1,ukgraduk2,ukgraduk3;

  // Compute offsets into u_local for each space
  INT offset_u1 = 0;
  INT offset_u2 = u1dofpelm;
  INT offset_u3 = u1dofpelm + u2dofpelm;
  INT offset_p;
  if(dim==3) offset_p = u1dofpelm + u2dofpelm + u3dofpelm;
  else offset_p = u1dofpelm + u2dofpelm;

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
  for (quad=0;quad<elm_data->quad_local->nq;quad++) {
    qx[0] = elm_data->quad_local->x[quad*dim];
    qx[1] = elm_data->quad_local->x[quad*dim+1];
    if(dim==3) qx[2] = elm_data->quad_local->x[quad*dim+2];
    w = elm_data->quad_local->w[quad];
    (*rhs)(rhs_val,qx,time,&elm_data->flag);

    // Get basis for all spaces
    for(i=0;i<nspaces;i++)
      get_FEM_basis_at_quadpt(elm_data, fe_data, i, quad);

    // Interpolate previous solution at this quadrature point
    u1_prev = 0.0;
    for(j=0;j<u1dofpelm;j++)
      u1_prev += u_local[offset_u1+j] * fe_data->phi[iu1][j];
    for(d=0;d<dim;d++) {
      du1_prev[d] = 0.0;
      for(j=0;j<u1dofpelm;j++)
        du1_prev[d] += u_local[offset_u1+j] * fe_data->dphi[iu1][j*dim+d];
    }

    u2_prev = 0.0;
    for(j=0;j<u2dofpelm;j++)
      u2_prev += u_local[offset_u2+j] * fe_data->phi[iu2][j];
    for(d=0;d<dim;d++) {
      du2_prev[d] = 0.0;
      for(j=0;j<u2dofpelm;j++)
        du2_prev[d] += u_local[offset_u2+j] * fe_data->dphi[iu2][j*dim+d];
    }

    if(dim==3) {
      u3_prev = 0.0;
      for(j=0;j<u3dofpelm;j++)
        u3_prev += u_local[offset_u3+j] * fe_data->phi[iu3][j];
      for(d=0;d<dim;d++) {
        du3_prev[d] = 0.0;
        for(j=0;j<u3dofpelm;j++)
          du3_prev[d] += u_local[offset_u3+j] * fe_data->dphi[iu3][j*dim+d];
      }
    }

    p_prev = 0.0;
    for(j=0;j<pdofpelm;j++)
      p_prev += u_local[offset_p+j] * fe_data->phi[ip][j];

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
      v1 = fe_data->phi[iu1][test];
      v1x = fe_data->dphi[iu1][test*dim];
      v1y = fe_data->dphi[iu1][test*dim+1];
      if(dim==3) v1z = fe_data->dphi[iu1][test*dim+2];

      local_col_index = 0;

      // Loop over block columns of trial functions
      // u1-v1 block:
      // 2*<dx(u1),dx(v1)> + <dy(u1),dy(v1)> + <dz(u1),dz(v1)> + Re*<u1k*dx(u1)+u2k*dy(u1) + u1*dx(u1k),v>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u1dofpelm;trial++){
        u1=fe_data->phi[iu1][trial];
        u1x=fe_data->dphi[iu1][trial*dim];
        u1y=fe_data->dphi[iu1][trial*dim+1];
        kij = 2*u1x*v1x + u1y*v1y + Re*(u1_prev*u1x+u2_prev*u1y+u1*du1_prev[0])*v1;
        if(dim==3) {
          u1z=fe_data->dphi[iu1][trial*dim+2];
          kij+=u1z*v1z + Re*(u3_prev*u1z)*v1;
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v1 block:
      // <dx(u2),dy(v1)> + Re*<u2*dy(u1k),v1>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u2dofpelm;trial++){
        u2=fe_data->phi[iu2][trial];
        u2x=fe_data->dphi[iu2][trial*dim];
        kij = u2x*v1y + Re*u2*du1_prev[1]*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v1 block <dx(u3),dz(v1)> + Re*<u3*dz(u1k),v1>
      if(dim==3) {
        for (trial=0; trial<u3dofpelm;trial++){
          u3=fe_data->phi[iu3][trial];
          u3x=fe_data->dphi[iu3][trial*dim];
          kij = u3x*v1z + Re*u3*du1_prev[2]*v1;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v1 block: -<p, dx(v1)>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<pdofpelm;trial++){
        p = fe_data->phi[ip][trial];
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
      v2 = fe_data->phi[iu2][test];
      v2x = fe_data->dphi[iu2][test*dim];
      v2y = fe_data->dphi[iu2][test*dim+1];
      if(dim==3) v2z = fe_data->dphi[iu2][test*dim+2];

      local_col_index = 0;

      // u1-v2 block
      // <dy(u1),dx(v2)> + Re*<u1*dx(u2k),v2>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u1dofpelm;trial++){
        u1=fe_data->phi[iu1][trial];
        u1y=fe_data->dphi[iu1][trial*dim+1];
        kij = u1y*v2x + Re*u1*du2_prev[0]*v2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-v2 block:
      // <dx(u2),dx(v2)> + 2*<dy(u2),dy(v2)> + <dz(u2),dz(v2)> + Re*<u_prev*grad(u2) + u2*dy(u2k),v2>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<u2dofpelm;trial++){
        u2=fe_data->phi[iu2][trial];
        u2x=fe_data->dphi[iu2][trial*dim];
        u2y=fe_data->dphi[iu2][trial*dim+1];
        kij = u2x*v2x + 2*u2y*v2y + Re*(u1_prev*u2x + u2_prev*u2y + u2*du2_prev[1])*v2;
        if(dim==3) {
          u2z=fe_data->dphi[iu2][trial*dim+2];
          kij+=u2z*v2z + Re*(u3_prev*u2z)*v2;
        }
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      // u3-v2 block
      // <dy(u3),dz(v2)> + Re*<u3*dz(u2k),v2>
      if(dim==3) {
        for (trial=0; trial<u3dofpelm;trial++){
          u3=fe_data->phi[iu3][trial];
          u3y=fe_data->dphi[iu3][trial*dim+1];
          kij = u3y*v2z + Re*u3*du2_prev[2]*v2;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }
        local_col_index += u3dofpelm;
      }

      // p-v2 block: -<p, dy(v2)>
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<pdofpelm;trial++){
        p = fe_data->phi[ip][trial];
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
        v3=fe_data->phi[iu3][test];
        v3x=fe_data->dphi[iu3][test*dim];
        v3y=fe_data->dphi[iu3][test*dim+1];
        v3z=fe_data->dphi[iu3][test*dim+2];

        local_col_index = 0;

        // u1-v3 block
        // <dz(u1),dx(v3)> + Re*<u1*dx(u3k),v3>
        local_col_index = 0;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<u1dofpelm;trial++){
          u1=fe_data->phi[iu1][trial];
          u1z=fe_data->dphi[iu1][trial*dim+2];
          kij = u1z*v3x + Re*u1*du3_prev[0]*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // u2-v3 block:
        // <dz(u2),dy(v3)> + Re*<u2*dy(u3k),v3>
        local_col_index += u1dofpelm;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<u2dofpelm;trial++){
          u2=fe_data->phi[iu2][trial];
          u2z=fe_data->dphi[iu2][trial*dim+2];
          kij = u2z*v3y + Re*u2*du3_prev[1]*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // u3-v3 block
        // <dx(u3),dx(v3)> + <dy(u3),dy(v3)> + 2*<dz(u3),dz(v3)> + Re*<u_prev*grad(u3) + u3*dz(u3k),v3>
        local_col_index += u2dofpelm;
        for (trial=0; trial<u3dofpelm;trial++){
          u3=fe_data->phi[iu3][trial];
          u3x=fe_data->dphi[iu3][trial*dim];
          u3y=fe_data->dphi[iu3][trial*dim+1];
          u3z=fe_data->dphi[iu3][trial*dim+2];
          kij = u3x*v3x + u3y*v3y + 2*u3z*v3z + Re*(u1_prev*u3x + u2_prev*u3y + u3_prev*u3z + u3*du3_prev[2])*v3;
          ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
        }

        // p-v3 block: -<p, dz(v3)>
        local_col_index += u3dofpelm;
        // Loop over Trial Functions (Columns)
        for (trial=0; trial<pdofpelm;trial++){
          p = fe_data->phi[ip][trial];
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
      q=fe_data->phi[ip][test];
      local_col_index = 0;

      // u1-q block: -<dx(u1), q>
      for (trial=0; trial<u1dofpelm;trial++){
        u1x=fe_data->dphi[iu1][trial*dim];
        kij = -u1x*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u1dofpelm;

      // u2-q block: -<dy(u2), q>
      for (trial=0; trial<u2dofpelm;trial++){
        u2y=fe_data->dphi[iu2][trial*dim+1];
        kij = -u2y*q;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }
      local_col_index += u2dofpelm;

      if(dim==3) {
        // u3-q block: -<dz(u3), q>
        for (trial=0; trial<u3dofpelm;trial++){
          u3z=fe_data->dphi[iu3][trial*dim+2];
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
