/*! \file LCSystem.h
 *  INCLUDE IN HAZMATH
 *
 *  Created by Anca Andrei on 04/01/20.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This contains all the local assembly routines
 *        for the LC Elastic example.
 *
 */

 // Compute Unit Length Constraint
 /**********************************/
 /*!
  * \fn void compute_LCelastic_unitlength(REAL* unitlength,REAL *u,block_fespace *FE,scomplex *sc,qcoordinates *cq)
  *
  * \brief Compute the unit length constraint for an LC simulation
  *        unitlength = ||n||_L2
  *
  * \param u 	         FE Approximation
  * \param FE          FE Space
  * \param sc          Simplicial Complex
  *
  * \return unitlength L2 norm of n
  *
  */
 void compute_LCelastic_unitlength(REAL* unitlength, REAL *u, block_fespace *FE, scomplex *sc, qcoordinates *cq) {
  INT i;

  // Get L2 norm of entire solution
  // 3 director components, 1 Lagrange multiplier
  REAL* solnormL2 = (REAL *) calloc(4, sizeof(REAL));
  L2norm_block(solnormL2,u,FE,sc,cq);

  // Grab just director portion
  REAL sum = 0;
  for(i=0;i<3;i++) sum += solnormL2[i]*solnormL2[i];

  // Grab each element and calculate the volume and add them up
  REAL vol = 0;
  INT elm;
  for (elm=0; elm<sc->fem->ns_leaf; elm++) {
    // add element volume to total
    vol += sc->fem->el_vol[elm];
  }
  *unitlength = sum/vol;

  if(solnormL2) free(solnormL2);
}


// Compute Energy
/**********************************/
/*!
 * \fn void compute_LCelastic_energy(REAL* energy,REAL *u,block_fespace *FE,scomplex *sc,qcoordinates *cq)
 *
 * \brief Compute the elastic energy for an LC simulation
 *        E(n) = 1/2 K1*||div n||^2 + 1/2 K3*<Z(n)*curl n,curl n>
 *
 * \param u 	      FE Approximation
 * \param FE          FE Space
 * \param sc          Simplicial Complex
 *
 * \return energy     Energy Value
 *
 */
void compute_LCelastic_energy(REAL* energy, REAL *u, block_fespace *FE, scomplex *sc, qcoordinates *cq) {

  // Counters
  INT elm, quad, i, j, rowa, rowb, jcntr;
  REAL splay = 0.0;
  REAL twist = 0.0;
  REAL bend = 0.0;

  // Mesh Stuff
  INT dim = sc->dim;
  INT* v_on_elm = (INT *) calloc(sc->dim+1, sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim, sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm, sizeof(INT));
  INT* local_dof_on_elm;

  // Solution Stuff
  REAL* local_uprev = NULL;
  REAL n1;
  REAL n2;
  REAL n3;
  REAL* dn1 = (REAL *) calloc(dim, sizeof(REAL));
  REAL* dn2 = (REAL *) calloc(dim, sizeof(REAL));
  REAL* dn3 = (REAL *) calloc(dim, sizeof(REAL));
  REAL divn;
  REAL* curln = (REAL *) calloc(3, sizeof(REAL));

  // Frank coefficients
  REAL* K = (REAL *) calloc(3, sizeof(REAL));
  // Coefficients
  get_frank_constants(K);
  REAL K1 = K[0];
  REAL K2 = K[1];
  REAL K3 = K[2];
  REAL kap = K2/K3;

  /* Loop over all Elements */
  for (elm=0; elm<sc->fem->ns_leaf; elm++) {
    // Find DOF for given Element
    // Note this is "local" ordering for the given FE space of the block
    // Not global ordering of all DOF
    jcntr = 0;
    for(i=0;i<nspaces;i++) {
      rowa = FE->var_spaces[i]->el_dof->IA[elm];
      rowb = FE->var_spaces[i]->el_dof->IA[elm+1];
      for (j=rowa; j<rowb; j++) {
        dof_on_elm[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
        jcntr++;
      }
    }

    // Find Vertices for given Element if not H1 elements
    get_incidence_row(elm,sc->fem->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(sc->dim==2 || sc->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(sc->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Interpolate FE solution to quadrature point
      // Get the Basis Functions and previous solutions at each quadrature node
      // n
      local_uprev = u;
      FE_Interpolation(&n1,local_uprev,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],sc);
      FE_DerivativeInterpolation(dn1,local_uprev,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],sc);
      local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
      local_uprev+=FE->var_spaces[0]->ndof;

      FE_Interpolation(&n2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],sc);
      FE_DerivativeInterpolation(dn2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],sc);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      local_uprev+=FE->var_spaces[1]->ndof;

      FE_Interpolation(&n3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],sc);
      FE_DerivativeInterpolation(dn3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],sc);

      divn     = dn1[0] + dn2[1];
      curln[0] = dn3[1];
      curln[1] = -dn3[0];
      curln[2] = dn2[0] - dn1[1];
      if(dim==3) {
        divn += dn3[2];
        curln[0] += -dn2[2];
        curln[1] += dn1[2];
      }

      // Computing the splay, twist, and bend.

      splay+=w*0.5*K1*divn*divn;
      twist+=w*0.5*K2*(n1*curln[0] + n2*curln[1] + n3*curln[2])*(n1*curln[0] + n2*curln[1] + n3*curln[2]);
      bend+=w*0.5*K3*((n2*curln[2]-n3*curln[1])*(n2*curln[2]-n3*curln[1]) +
          (n3*curln[0] - n1*curln[2])*(n3*curln[0] - n1*curln[2]) +
          (n1*curln[1] - n2*curln[0])*(n1*curln[1] - n2*curln[0]));
    }
  }

  energy[0] = splay + twist + bend;
  energy[1] = splay;
  energy[2] = twist;
  energy[3] = bend;

  if(dof_on_elm) free(dof_on_elm);
  if(v_on_elm) free(v_on_elm);
  if(qx) free(qx);
  if(dn1) free(dn1);
  if(dn2) free(dn2);
  if(dn3) free(dn3);
  if(curln) free(curln);
  if(K) free(K);

}

/*!
 * \fn void local_assembly_LCelastic(REAL *ALoc,REAL *bLoc, dvector *old_sol, block_fespace *FE, scomplex *sc, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
 *
 * \brief Computes the local stiffness matrix for the LC Elastic system.
 *        For this problem we compute LHS of:
 *
 *       K1<div(dn),div(v)> + K3<Z(nk)*curl(dn),curl(v)> +
 *       (K2-K3)*[<nk*curl(nk),dn*curl(v)+v*curl(dn)> +
 *                <dn*curl(nk),nk*curl(v) + v*curl(nk)> +
 *                <nk*curl(dn),v*curl(nk)>] + 2*int(lamk*dn*v)     +  2*int(dlam*nk*v)    = -L1(nk,lamk,v)
 *
 *       2*int(gam*nk*dn)                                          +  0                   = -L2(nk,lamk,gam)
 *
 * where
 *
 *       L1(n,lam,v)   = K1*<div n, div v> + K3<Z(n) curl n, curl v> + (K2-K3)*<n*curl n, v*curl n> + 2*int [(n*v)*lam]
 *       L2(n,lam,gam) = int [(n*n-1)*gam]
 *
 * \param old_sol       FE approximation of previous Newton step
 * \param FE            Block FE Space
 * \param sc            Simplicial Complex
 * \param cq            Quadrature Nodes
 * \param dof_on_elm    Specific DOF on element
 * \param v_on_elm      Specific vertices on element
 * \param elm           Current element
 * \param rhs           Function for rhs (NULL here)
 * \param time          Physical Time if time dependent
 *
 * \return ALoc         Local Stiffness Matrix (Full Matrix) ordered (u1,u2,u3,p)
 * \return bLoc         Local RHS ordered (u1,u2,u3,p)
 *
 *
 */
void local_assembly_LCelastic(REAL *ALoc,REAL* bLoc, REAL *u_local,
    simplex_local_data *elm_data, fe_local_data *fe_data,
    void (*rhs)(REAL *,REAL *,REAL,void *),
    void (*coeff)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop indices
  INT i, j, d, quad, test, trial;

  // Mesh and FE data
  INT dof_per_elm = fe_data->n_dof;
  INT dim = elm_data->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL qx[dim];

  // Stiffness Matrix and RHS Entry
  REAL kij = 0.0;
  REAL rij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Stuff for previous solution
  REAL nk1, nk2, nk3, lamk;
  REAL dnk1[dim];
  REAL dnk2[dim];
  REAL dnk3[dim];

  // Compute offsets into u_local for each space
  INT n1dofpelm = fe_data->n_dof_per_space[0];
  INT n2dofpelm = fe_data->n_dof_per_space[1];
  INT n3dofpelm = fe_data->n_dof_per_space[2];
  INT lamdofpelm = fe_data->n_dof_per_space[3];
  INT offset_n1 = 0;
  INT offset_n2 = n1dofpelm;
  INT offset_n3 = n1dofpelm + n2dofpelm;
  INT offset_lam = n1dofpelm + n2dofpelm + n3dofpelm;
  REAL divnk, curlnk1, curlnk2, curlnk3, nkdotcurlnk;
  REAL znk11, znk12, znk13, znk21, znk22, znk23, znk31, znk32, znk33;
  REAL znk_curl_nk1, znk_curl_nk2, znk_curl_nk3;
  REAL nk1x, nk1y, nk1z, nk2x, nk2y, nk2z, nk3x, nk3y, nk3z;

  // Trial Functions
  REAL n1, n2, n3, n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, lam;
  // Test Functions
  REAL v1, v2, v3, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, gam;

  // Frank coefficients
  REAL* K = (REAL *) calloc(3,sizeof(REAL));
  get_frank_constants(K);
  REAL K1=K[0];
  REAL K2=K[1];
  REAL K3=K[2];
  REAL K2_min_K3 = K2 - K3;
  REAL kappa = K2/K3;

  // Sum over quadrature points
  for (quad=0;quad<elm_data->quad_local->nq;quad++) {
    qx[0] = elm_data->quad_local->x[quad*dim];
    qx[1] = elm_data->quad_local->x[quad*dim+1];
    if(dim==3) qx[2] = elm_data->quad_local->x[quad*dim+2];
    w = elm_data->quad_local->w[quad];

    // Get basis functions for all spaces
    for(i=0;i<fe_data->nspaces;i++)
      get_FEM_basis_at_quadpt(elm_data, fe_data, i, quad);

    // Interpolate previous solution at quadrature point
    nk1 = 0.0;
    for(j=0;j<n1dofpelm;j++) nk1 += u_local[offset_n1+j] * fe_data->phi[0][j];
    for(d=0;d<dim;d++) {
      dnk1[d] = 0.0;
      for(j=0;j<n1dofpelm;j++) dnk1[d] += u_local[offset_n1+j] * fe_data->dphi[0][j*dim+d];
    }

    nk2 = 0.0;
    for(j=0;j<n2dofpelm;j++) nk2 += u_local[offset_n2+j] * fe_data->phi[1][j];
    for(d=0;d<dim;d++) {
      dnk2[d] = 0.0;
      for(j=0;j<n2dofpelm;j++) dnk2[d] += u_local[offset_n2+j] * fe_data->dphi[1][j*dim+d];
    }

    nk3 = 0.0;
    for(j=0;j<n3dofpelm;j++) nk3 += u_local[offset_n3+j] * fe_data->phi[2][j];
    for(d=0;d<dim;d++) {
      dnk3[d] = 0.0;
      for(j=0;j<n3dofpelm;j++) dnk3[d] += u_local[offset_n3+j] * fe_data->dphi[2][j*dim+d];
    }

    lamk = 0.0;
    for(j=0;j<lamdofpelm;j++) lamk += u_local[offset_lam+j] * fe_data->phi[3][j];

    // Some precomputations
    nk1x = dnk1[0];
    nk1y = dnk1[1];

    nk2x = dnk2[0];
    nk2y = dnk2[1];

    nk3x = dnk3[0];
    nk3y = dnk3[1];

    if(dim == 3){
      nk1z = dnk1[2];
      nk2z = dnk2[2];
      nk3z = dnk3[2];
    } else {
      nk1z = 0.0;
      nk2z = 0.0;
      nk3z = 0.0;
    }

    divnk = nk1x + nk2y + nk3z;
    curlnk1 = nk3y - nk2z;
    curlnk2 = nk1z - nk3x;
    curlnk3 = nk2x - nk1y;

    nkdotcurlnk = nk1*curlnk1 + nk2*curlnk2 + nk3*curlnk3;

    znk11 = 1 + (kappa - 1) * nk1 * nk1;
    znk12 = (kappa - 1) * nk1 * nk2;
    znk13 = (kappa - 1) * nk1 * nk3;
    znk21 = (kappa - 1) * nk2 * nk1;
    znk22 = 1 + (kappa - 1) * nk2 * nk2;
    znk23 = (kappa - 1) * nk2 * nk3;
    znk31 = (kappa - 1) * nk3 * nk1;
    znk32 = (kappa - 1) * nk3 * nk2;
    znk33 = 1 + (kappa - 1) * nk3 * nk3;

    znk_curl_nk1 = znk11*curlnk1 + znk12*curlnk2 + znk13*curlnk3;
    znk_curl_nk2 = znk21*curlnk1 + znk22*curlnk2 + znk23*curlnk3;
    znk_curl_nk3 = znk31*curlnk1 + znk32*curlnk2 + znk33*curlnk3;

    // v1 block row
    local_row_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<fe_data->n_dof_per_space[0];test++) {
      v1 = fe_data->phi[0][test];
      v1x = fe_data->dphi[0][test*dim];
      v1y = fe_data->dphi[0][test*dim+1];
      if(dim==3){
        v1z = fe_data->dphi[0][test*dim+2];
      } else {
        v1z = 0.0;
      }

      // n1 block column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[0];trial++) {
        n1 = fe_data->phi[0][trial];
        n1x = fe_data->dphi[0][trial*dim];
        n1y = fe_data->dphi[0][trial*dim+1];
        if(dim==3){
          n1z = fe_data->dphi[0][trial*dim+2];
        } else {
          n1z = 0.0;
        }

        kij = K1*n1x*v1x
            + K3*(v1z*znk22*n1z - v1y*znk32*n1z + v1y*znk33*n1y - v1z*znk23*n1y)
            + K2_min_K3*(nk2*v1z*n1*curlnk1 - nk3*v1y*n1*curlnk1 + nk2*n1z*v1*curlnk1 - nk3*n1y*v1*curlnk1)
            + K2_min_K3*(n1*curlnk1*v1*curlnk1)
            + lamk*n1*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += fe_data->n_dof_per_space[0];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[1];trial++){
        n2 = fe_data->phi[1][trial];
        n2x = fe_data->dphi[1][trial*dim];
        n2y = fe_data->dphi[1][trial*dim+1];
        if(dim==3){
          n2z = fe_data->dphi[1][trial*dim+2];
        } else {
          n2z = 0.0;
        }

        kij = K1*n2y*v1x
            + K3*(v1y*znk31*n2z - v1z*znk21*n2z + v1z*znk23*n2x - v1y*znk33*n2x)
            + K2_min_K3*(nkdotcurlnk*n2*v1z - nkdotcurlnk*v1*n2z + nk2*v1z*n2*curlnk2 - nk3*v1y*n2*curlnk2)
            + K2_min_K3*(nk3*n2x*v1*curlnk1 - nk1*n2z*v1*curlnk1 + n2*curlnk2*v1*curlnk1);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += fe_data->n_dof_per_space[1];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[2];trial++){
        n3 = fe_data->phi[2][trial];
        n3x = fe_data->dphi[2][trial*dim];
        n3y = fe_data->dphi[2][trial*dim+1];
        if(dim==3){
          n3z = fe_data->dphi[2][trial*dim+2];
        } else {
          n3z = 0.0;
        }

        kij = K1*n3z*v1x
            + K3*(v1z*znk21*n3y - v1y*znk31*n3y + v1y*znk32*n3x - v1z*znk22*n3x)
            + K2_min_K3*(nkdotcurlnk*v1*n3y - nkdotcurlnk*n3*v1y + nk2*v1z*n3*curlnk3 - nk3*v1y*n3*curlnk3)
            + K2_min_K3*(nk1*n3y*v1*curlnk1 - nk2*n3x*v1*curlnk1 + n3*curlnk3*v1*curlnk1);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column - done
      local_col_index += fe_data->n_dof_per_space[2];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[3];trial++){
        lam=fe_data->phi[3][trial];
        kij = lam*nk1*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v1 part of RHS
      rij = -(K1*divnk*v1x + K3*(znk_curl_nk2*v1z - znk_curl_nk3*v1y)
              + K2_min_K3*(nkdotcurlnk*v1*curlnk1) + lamk*nk1*v1);
      bLoc[(local_row_index+test)] += w*rij;
    }

    // v2 block row
    local_row_index += fe_data->n_dof_per_space[0];
    // Loop over Test Functions (Rows)
    for (test=0; test<fe_data->n_dof_per_space[1];test++) {
      v2 = fe_data->phi[1][test];
      v2x = fe_data->dphi[1][test*dim];
      v2y = fe_data->dphi[1][test*dim+1];
      if(dim==3){
        v2z = fe_data->dphi[1][test*dim+2];
      } else {
        v2z = 0.0;
      }

      // n1 block column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[0];trial++) {
        n1 = fe_data->phi[0][trial];
        n1x = fe_data->dphi[0][trial*dim];
        n1y = fe_data->dphi[0][trial*dim+1];
        if(dim==3){
          n1z = fe_data->dphi[0][trial*dim+2];
        } else {
          n1z = 0.0;
        }

        kij = K1*n1x*v2y
            + K3 * (v2x*znk32*n1z - v2z*znk12*n1z + v2z*znk13*n1y - v2x*znk33*n1y)
            + K2_min_K3*(nkdotcurlnk*v2*n1z - nkdotcurlnk*n1*v2z + nk3*v2x*n1*curlnk1 - nk1*v2z*n1*curlnk1)
            + K2_min_K3*(nk2*n1z*v2*curlnk2 - nk3*n1y*v2*curlnk2 + n1*curlnk1*v2*curlnk2);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += fe_data->n_dof_per_space[0];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[1];trial++){
        n2 = fe_data->phi[1][trial];
        n2x = fe_data->dphi[1][trial*dim];
        n2y = fe_data->dphi[1][trial*dim+1];
        if(dim==3){
          n2z = fe_data->dphi[1][trial*dim+2];
        } else {
          n2z = 0.0;
        }

        kij = K1*n2y*v2y
            + K3*(v2z*znk11*n2z - v2x*znk31*n2z + v2x*znk33*n2x - v2z*znk13*n2x)
            + K2_min_K3*(nk3*v2x*n2*curlnk2 - nk1*v2z*n2*curlnk2 + nk3*n2x*v2*curlnk2 - nk1*n2z*v2*curlnk2)
            + K2_min_K3*(n2*curlnk2*v2*curlnk2)
            + lamk*v2*n2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += fe_data->n_dof_per_space[1];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[2];trial++){
        n3 = fe_data->phi[2][trial];
        n3x = fe_data->dphi[2][trial*dim];
        n3y = fe_data->dphi[2][trial*dim+1];
        if(dim==3){
          n3z = fe_data->dphi[2][trial*dim+2];
        } else {
          n3z = 0.0;
        }

        kij = K1*n3z*v2y
            + K3*(v2x*znk31*n3y - v2z*znk11*n3y + v2z*znk12*n3x - v2x*znk32*n3x)
            + K2_min_K3*(nkdotcurlnk*n3*v2x - nkdotcurlnk*v2*n3x + nk3*v2x*n3*curlnk3 - nk1*v2z*n3*curlnk3)
            + K2_min_K3*(nk1*n3y*v2*curlnk2 - nk2*n3x*v2*curlnk2 + n3*curlnk3*v2*curlnk2);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column
      local_col_index += fe_data->n_dof_per_space[2];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[3];trial++){
        lam=fe_data->phi[3][trial];
        kij = lam*nk2*v2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v2 part of RHS block
      rij = -(K1*divnk*v2y + K3*(znk_curl_nk3*v2x - znk_curl_nk1*v2z)
              + K2_min_K3*(nkdotcurlnk*v2*curlnk2) + lamk*nk2*v2);
      bLoc[(local_row_index+test)] += w*rij;
    }

    // v3 block row
    local_row_index += fe_data->n_dof_per_space[1];
    // Loop over Test Functions (Rows)
    for (test=0; test<fe_data->n_dof_per_space[2];test++) {
      v3 = fe_data->phi[2][test];
      v3x = fe_data->dphi[2][test*dim];
      v3y = fe_data->dphi[2][test*dim+1];
      if(dim==3){
        v3z = fe_data->dphi[2][test*dim+2];
      } else {
        v3z = 0.0;
      }

      // n1 block Columns
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[0];trial++) {
        n1 = fe_data->phi[0][trial];
        n1x = fe_data->dphi[0][trial*dim];
        n1y = fe_data->dphi[0][trial*dim+1];
        if(dim==3){
          n1z = fe_data->dphi[0][trial*dim+2];
        } else {
          n1z = 0.0;
        }

        kij = K1*n1x*v3z
            + K3*(v3y*znk12*n1z - v3x*znk22*n1z + v3x*znk23*n1y - v3y*znk13*n1y)
            + K2_min_K3*(nkdotcurlnk*n1*v3y - nkdotcurlnk*v3*n1y + nk1*v3y*n1*curlnk1 - nk2*v3x*n1*curlnk1)
            + K2_min_K3*(nk2*n1z*v3*curlnk3 - nk3*n1y*v3*curlnk3 + n1*curlnk1*v3*curlnk3);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += fe_data->n_dof_per_space[0];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[1];trial++){
        n2 = fe_data->phi[1][trial];
        n2x = fe_data->dphi[1][trial*dim];
        n2y = fe_data->dphi[1][trial*dim+1];
        if(dim==3){
          n2z = fe_data->dphi[1][trial*dim+2];
        } else {
          n2z = 0.0;
        }

        kij = K1*n2y*v3z
            + K3*(v3x*znk21*n2z - v3y*znk11*n2z + v3y*znk13*n2x - v3x*znk23*n2x)
            + K2_min_K3*(nkdotcurlnk*v3*n2x - nkdotcurlnk*n2*v3x + nk1*v3y*n2*curlnk2 - nk2*v3x*n2*curlnk2)
            + K2_min_K3*(nk3*n2x*v3*curlnk3 - nk1*n2z*v3*curlnk3 + n2*curlnk2*v3*curlnk3);
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += fe_data->n_dof_per_space[1];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[2];trial++){
        n3 = fe_data->phi[2][trial];
        n3x = fe_data->dphi[2][trial*dim];
        n3y = fe_data->dphi[2][trial*dim+1];
        if(dim==3){
          n3z = fe_data->dphi[2][trial*dim+2];
        } else {
          n3z = 0.0;
        }

        kij = K1*n3z*v3z
            + K3*(v3y*znk11*n3y - v3x*znk21*n3y + v3x*znk22*n3x - v3y*znk12*n3x)
            + K2_min_K3*(nk1*v3y*n3*curlnk3 - nk2*v3x*n3*curlnk3 + nk1*n3y*v3*curlnk3 - nk2*n3x*v3*curlnk3)
            + K2_min_K3*(n3*curlnk3*v3*curlnk3)
            + lamk*v3*n3;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column
      local_col_index += fe_data->n_dof_per_space[2];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[3];trial++){
        lam=fe_data->phi[3][trial];
        kij = lam*nk3*v3;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v3 part of RHS block
      rij = -(K1*divnk*v3z + K3*(znk_curl_nk1*v3y - znk_curl_nk2*v3x)
              + K2_min_K3*(nkdotcurlnk*v3*curlnk3) + lamk*nk3*v3);
      bLoc[(local_row_index+test)] += w*rij;
    }

    // gam block row
    local_row_index += fe_data->n_dof_per_space[2];
    // Loop over Test Functions (Rows)
    for (test=0; test<fe_data->n_dof_per_space[3];test++){
      gam=fe_data->phi[3][test];

      // n1 block Column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[0];trial++){
        n1=fe_data->phi[0][trial];
        kij = gam*nk1*n1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += fe_data->n_dof_per_space[0];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[1];trial++){
        n2=fe_data->phi[1][trial];
        kij = gam*nk2*n2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += fe_data->n_dof_per_space[1];
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<fe_data->n_dof_per_space[2];trial++){
        n3=fe_data->phi[2][trial];
        kij = gam*nk3*n3;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // gamma part of RHS block
      rij = -0.5*gam*(nk1*nk1 + nk2*nk2 + nk3*nk3 - 1.0);
      bLoc[(local_row_index+test)] += w*rij;
    }
  }

  // Free stuff
  if (K) free(K);

  return;
}

// --------------------------------------------------------
// FE setup and boundary conditions setup

// Relabel boundary flags we'll have the following labels:
// 1: x=0, 2: x=1, 3: y=0, 4: y=1, 5: z=0, 6: z=1
void relabel_boundary(fespace* FE, INT dim) {

  INT i;
  for(i=0;i<FE->ndof;i++) {
    if(FE->dof_flag[i]==1) {
      if(FE->cdof->x[i]==0.0) {
        FE->dof_flag[i] = 1;
      }
      if(FE->cdof->x[i]==1.0) {
        FE->dof_flag[i] = 2;
      }
      if(FE->cdof->y[i]==0.0) {
        FE->dof_flag[i] = 3;
      }
      if(FE->cdof->y[i]==1.0) {
          FE->dof_flag[i] = 4;
      }
      if(dim==3){
        if(FE->cdof->z[i]==0.0) {
            FE->dof_flag[i] = 5;
        }
        if(FE->cdof->z[i]==1.0) {
            FE->dof_flag[i] = 6;
        }
      }
    }
  }
  return;
}

// Set up FE spaces for the LC Elasticity problem
void setup_FEspaces(INT order_n, INT order_lam,fespace* FE_nx, fespace* FE_ny, fespace* FE_nz, fespace* FE_lam, block_fespace *FE,block_dCSRmat *P_periodic,scomplex* sc,INT dim) {

  // Need Spaces for each component of the director plus lambda
  create_fespace(FE_nx,sc,order_n);
  create_fespace(FE_ny,sc,order_n);
  create_fespace(FE_nz,sc,order_n);
  create_fespace(FE_lam,sc,order_lam);

  // Relabel boundary flags for FE space to distinguish between different boundaries
  // Boundary codes
  // 1: x=0  2: x=1   3: y=0    4: y=1, 5: z=0, 6: z=1
  relabel_boundary(FE_nx, dim);
  relabel_boundary(FE_ny, dim);
  relabel_boundary(FE_nz, dim);
  relabel_boundary(FE_lam, dim);

  // Decide whether you want Dirichlet or Periodic boundary conditions.
  // periodic_flag = 1 means it's periodic bc
  // dirichlet_flag = 1 means it's dirichlet bc
  // set in data (*.h) file
  INT bc_flag = characterize_boundary_conditions();
  INT periodic_flag = bc_flag;
  INT dirichlet_flag = 0;
  if (bc_flag == 0){
    dirichlet_flag = 1;
  }

  // Lambda has Neumann boundaries regardless
  set_dirichlet_bdry(FE_lam,sc,-10,-10);

  // Dirichlet Boundaries (Dirichlet all around for n)
  if (dirichlet_flag == 1){
    set_dirichlet_bdry(FE_nx,sc,1,6);
    set_dirichlet_bdry(FE_ny,sc,1,6);
    set_dirichlet_bdry(FE_nz,sc,1,6);
  }

  // Periodic Boundaries on right-left; Dirichelt on top-bottom
  if (periodic_flag == 1){
    if(dim==2) {
      // For 2D, the y bounds are Dirichlet while x is periodic.
      set_periodic_bdry(FE_nx,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(FE_ny,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(FE_nz,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_periodic_bdry(FE_lam,sc,0.0,1.0,0.0,0.0,0.0,0.0);
      set_dirichlet_bdry(FE_nx,sc,3,4);
      set_dirichlet_bdry(FE_ny,sc,3,4);
      set_dirichlet_bdry(FE_nz,sc,3,4);
    } else if(dim==3) {
      // For 3D, the z bounds are Dirichlet while x & y are periodic
      set_periodic_bdry(FE_nx,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(FE_ny,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(FE_nz,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      set_periodic_bdry(FE_lam,sc,0.0,1.0,0.0,1.0,0.0,0.0);
      // Setting Dirichlet boundary conditions
      set_dirichlet_bdry(FE_nx,sc,5,6);
      set_dirichlet_bdry(FE_ny,sc,5,6);
      set_dirichlet_bdry(FE_nz,sc,5,6);
    } else {
      check_error(ERROR_DIM, __FUNCTION__);
    }
  }

  // Create Block System with ordering (n,lam)
  INT ndof = FE_nx->ndof + FE_ny->ndof + FE_nz->ndof +FE_lam->ndof;
  INT nspaces = 4;
  INT nun = 4;
  // Get Global FE Space
  initialize_fesystem(FE,nspaces,nun,ndof,sc->fem->ns_leaf);
  FE->var_spaces[0] = FE_nx;
  FE->var_spaces[1] = FE_ny;
  FE->var_spaces[2] = FE_nz;
  FE->var_spaces[3] = FE_lam;

  // Set Dirichlet Boundaries
  set_dirichlet_bdry_block(FE,sc);

  // Set Periodic Boundaries
  
  if (periodic_flag == 1){
    generate_periodic_P_blockFE(FE, P_periodic);
  }
  return;
}


