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
  * \fn void compute_LCelastic_unitlength(REAL* unitlength,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
  *
  * \brief Compute the unit length constraint for an LC simulation
  *        unitlength = ||n||_L2
  *
  * \param u 	         FE Approximation
  * \param FE          FE Space
  * \param mesh        Mesh Data
  *
  * \return unitlength L2 norm of n
  *
  */
 void compute_LCelastic_unitlength(REAL* unitlength,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq) {

   INT i;

   // Get L2 norm of entire solution
   REAL* solnormL2 = (REAL *) calloc(mesh->dim+2, sizeof(REAL));
   L2norm_block(solnormL2,u,FE,mesh,cq);

   // Grab just n portion
   REAL sum = 0;
   for(i=0;i<mesh->dim+1;i++) sum += solnormL2[i]*solnormL2[i];

   *unitlength = sqrt(sum);

   if(solnormL2) free(solnormL2);
 }


// Compute Energy
/**********************************/
/*!
 * \fn void compute_LCelastic_energy(REAL* energy,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
 *
 * \brief Compute the elastic energy for an LC simulation
 *        E(n) = 1/2 K1*||div n||^2 + 1/2 K3*<Z(n)*curl n,curl n>
 *
 * \param u 	      FE Approximation
 * \param FE          FE Space
 * \param mesh        Mesh Data
 *
 * \return energy     Energy Value
 *
 */
void compute_LCelastic_energy(REAL* energy,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq) {

  // Counters
  INT elm,quad,i,j,rowa,rowb,jcntr;
  REAL splay = 0.0;
  REAL twist = 0.0;
  REAL bend = 0.0;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* local_dof_on_elm;

  // Solution Stuff
  REAL* local_uprev = NULL;
  REAL n1;
  REAL n2;
  REAL n3;
  REAL* dn1 = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn2 = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn3 = (REAL *) calloc( dim, sizeof(REAL));
  REAL divn;
  REAL* curln = (REAL *) calloc( 3, sizeof(REAL));

  // Frank coefficients
  REAL* K = (REAL *) calloc(3,sizeof(REAL));
  // Coefficients
  get_frank_constants(K);
  REAL K1 = K[0];
  REAL K2 = K[1];
  REAL K3 = K[2];
  REAL kap = K2/K3;

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {

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
    get_incidence_row(elm,mesh->el_v,v_on_elm);

    // Loop over quadrature nodes on element
    for (quad=0;quad<cq->nq_per_elm;quad++) {
      qx[0] = cq->x[elm*cq->nq_per_elm+quad];
      if(mesh->dim==2 || mesh->dim==3)
        qx[1] = cq->y[elm*cq->nq_per_elm+quad];
      if(mesh->dim==3)
        qx[2] = cq->z[elm*cq->nq_per_elm+quad];
      w = cq->w[elm*cq->nq_per_elm+quad];

      // Interpolate FE solution to quadrature point
      //  Get the Basis Functions and previous solutions at each quadrature node
      // n
      local_uprev = u;
      FE_Interpolation(&n1,local_uprev,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
      FE_DerivativeInterpolation(dn1,local_uprev,qx,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
      local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
      local_uprev+=FE->var_spaces[0]->ndof;

      FE_Interpolation(&n2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      FE_DerivativeInterpolation(dn2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      local_uprev+=FE->var_spaces[1]->ndof;

      FE_Interpolation(&n3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
      FE_DerivativeInterpolation(dn3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);

      // Compute divs and curls
      divn     = dn1[0] + dn2[1];
      curln[0] = dn3[1];
      curln[1] = -dn3[0];
      curln[2] = dn2[0] - dn1[1];

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
 * \fn void local_assembly_LCelastic(REAL *ALoc,REAL *bLoc, dvector *old_sol, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time)
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
 * \param mesh          Mesh Data
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
void local_assembly_LCelastic(REAL *ALoc,REAL* bLoc, dvector *old_sol, block_fespace *FE, mesh_struct *mesh, qcoordinates *cq, INT *dof_on_elm, INT *v_on_elm, INT elm,void (*rhs)(REAL *,REAL *,REAL,void *),REAL time)
{

  // Loop indices
  INT i,j,quad,test,trial;

  // Mesh and FE data
  INT dof_per_elm = 0;
  for (i=0; i<FE->nspaces;i++)
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
  INT* local_dof_on_elm;
  INT dim = mesh->dim;

  // Quadrature Weights and Nodes
  REAL w;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));

  // Stiffness Matrix and RHS Entry
  REAL kij = 0.0;
  REAL rij = 0.0;

  // Keep track of local indexing
  INT local_row_index, local_col_index;

  // Stuff for previous solution
  REAL* local_uprev = NULL;
  REAL n1_prev;
  REAL n2_prev;
  REAL n3_prev;
  REAL lam_prev;
  REAL* dn1_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn2_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn3_prev = (REAL *) calloc( dim, sizeof(REAL));
  REAL divnk,curlnk1,curlnk2,curlnk3,nkdotcurlnk;
  REAL nkcrosscurlnk1,nkcrosscurlnk2,nkcrosscurlnk3;

  // Trial Functions
  REAL n1,n2,n3,lam,n1x,n1y,n2x,n2y,n3x,n3y;
  // Test Functions
  REAL v1,v2,v3,v1x,v1y,v2x,v2y,v3x,v3y,gam;

  // Frank coefficients
  REAL* K = (REAL *) calloc(3,sizeof(REAL));
  get_frank_constants(K);
  REAL K1=K[0];
  REAL K2=K[1];
  REAL K3=K[2];

  // Sum over quadrature points
  for (quad=0;quad<cq->nq_per_elm;quad++) {
    qx[0] = cq->x[elm*cq->nq_per_elm+quad];
    qx[1] = cq->y[elm*cq->nq_per_elm+quad];
    if(mesh->dim==3) qx[2] = cq->z[elm*cq->nq_per_elm+quad];
    w = cq->w[elm*cq->nq_per_elm+quad];

    //  Get the Basis Functions and previous solutions at each quadrature node
    // nk1, n1, and v1
    local_dof_on_elm = dof_on_elm;
    local_uprev = old_sol->val;
    FE_Interpolation(&n1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
    FE_DerivativeInterpolation(dn1_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
    get_FEM_basis(FE->var_spaces[0]->phi,FE->var_spaces[0]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[0]);

    // nk2, n2, and v2
    local_dof_on_elm += FE->var_spaces[0]->dof_per_elm;
    local_uprev += FE->var_spaces[0]->ndof;
    FE_Interpolation(&n2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
    FE_DerivativeInterpolation(dn2_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
    get_FEM_basis(FE->var_spaces[1]->phi,FE->var_spaces[1]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[1]);

    // nk3, n3, and v3
    local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
    local_uprev+=FE->var_spaces[1]->ndof;
    FE_Interpolation(&n3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
    FE_DerivativeInterpolation(dn3_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
    get_FEM_basis(FE->var_spaces[2]->phi,FE->var_spaces[2]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[2]);

    // lamk, lam, gam
    local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
    local_uprev += FE->var_spaces[2]->ndof;
    FE_Interpolation(&lam_prev,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh);
    get_FEM_basis(FE->var_spaces[3]->phi,FE->var_spaces[3]->dphi,qx,v_on_elm,local_dof_on_elm,mesh,FE->var_spaces[3]);

    // Some precomputations
    divnk = dn1_prev[0] + dn2_prev[1];
    curlnk1 = dn3_prev[1];
    curlnk2 = -dn3_prev[0];
    curlnk3 = dn2_prev[0]-dn1_prev[1];
    nkdotcurlnk = n1_prev*curlnk1 + n2_prev*curlnk2 + n3_prev*curlnk3;
    nkcrosscurlnk1 = n2_prev*curlnk3 - n3_prev*curlnk2;
    nkcrosscurlnk2 = -n1_prev*curlnk3 + n3_prev*curlnk1;
    nkcrosscurlnk3 = n1_prev*curlnk2 - n2_prev*curlnk1;

    // v1 block row
    local_row_index = 0;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[0]->dof_per_elm;test++) {
      v1=FE->var_spaces[0]->phi[test];
      v1x=FE->var_spaces[0]->dphi[test*dim];
      v1y=FE->var_spaces[0]->dphi[test*dim+1];

      // n1 block column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++) {
        n1=FE->var_spaces[0]->phi[trial];
        n1x=FE->var_spaces[0]->dphi[trial*dim];
        n1y=FE->var_spaces[0]->dphi[trial*dim+1];
        kij = K1*n1x*v1x
            + K2*(v1*curlnk1-v1y*n3_prev)*(-n3_prev*n1y+curlnk1*n1)
            + K3*(n1y*n2_prev*n2_prev*v1y + (n1y*n1_prev-n1*curlnk3)
                *(n1_prev*v1y - v1*curlnk3) + n1*curlnk2*v1*curlnk2
                  + nkcrosscurlnk2*(v1*n1y+n1*v1y))
                  + 2*lam_prev*n1*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += FE->var_spaces[0]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        n2=FE->var_spaces[1]->phi[trial];
        n2x=FE->var_spaces[1]->dphi[trial*dim];
        n2y=FE->var_spaces[1]->dphi[trial*dim+1];
        kij = K1*n2y*v1x
            + K2*(v1*curlnk1-v1y*n3_prev)*(n3_prev*n2x+curlnk2*n2)
            + K3*(- n2_prev * v1y * (n2_prev * n2x + n2 * curlnk3)
                  + (-n2x*n1_prev)*(n1_prev*v1y - v1*curlnk3)
                  - n2 *curlnk1*v1*curlnk2 - nkcrosscurlnk1*(n2*v1y)
                  - nkcrosscurlnk2*(v1*n2x));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += FE->var_spaces[1]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
        n3=FE->var_spaces[2]->phi[trial];
        n3x=FE->var_spaces[2]->dphi[trial*dim];
        n3y=FE->var_spaces[2]->dphi[trial*dim+1];
        kij = K2 *((v1*curlnk1-v1y*n3_prev)
                  *(n1_prev*n3y-n2_prev*n3x+n3*curlnk3)
                  + nkdotcurlnk*(v1*n3y-n3*v1y))
            + K3 *(-n2_prev*v1y*(n3_prev*n3x-n3*curlnk2)
                 +(n1_prev*v1y - v1*curlnk3)*(n3_prev*n3y+n3*curlnk1)
                 + v1*curlnk2*(-n1_prev*n3x-n2_prev*n3y)
                 + nkcrosscurlnk3*(-v1*n3x));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column - done
      local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
        lam=FE->var_spaces[3]->phi[trial];
        kij = 2*lam*n1_prev*v1;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v1 part of RHS
      rij = - K1*divnk*v1x
            - K2*nkdotcurlnk*(v1*curlnk1 - v1y*n3_prev)
            - K3*(nkcrosscurlnk1*(-v1y*n2_prev)
                  + nkcrosscurlnk2*(v1y*n1_prev-v1*curlnk3)
                  + nkcrosscurlnk3*v1*curlnk2)
            - 2*lam_prev*n1_prev*v1;
      bLoc[(local_row_index+test)] += w*rij;
    }

    // v2 block row
    local_row_index += FE->var_spaces[0]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[1]->dof_per_elm;test++) {
      v2=FE->var_spaces[1]->phi[test];
      v2x=FE->var_spaces[1]->dphi[test*dim];
      v2y=FE->var_spaces[1]->dphi[test*dim+1];

      // n1 block column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++) {
        n1=FE->var_spaces[0]->phi[trial];
        n1x=FE->var_spaces[0]->dphi[trial*dim];
        n1y=FE->var_spaces[0]->dphi[trial*dim+1];
        kij = K1*n1x*v2y
            + K2*((v2*curlnk2+v2x*n3_prev)*(-n1y*n3_prev+n1*curlnk1)) +
              K3*((v2x*n2_prev+v2*curlnk3)*(-n1y*n2_prev)
                - v2x*n1_prev*(n1y*n1_prev-n1*curlnk3)
                - v2*curlnk1*(n1*curlnk2)
                - nkcrosscurlnk1*(v2*n1y)
                - nkcrosscurlnk2*(v2x*n1));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += FE->var_spaces[0]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        n2=FE->var_spaces[1]->phi[trial];
        n2x=FE->var_spaces[1]->dphi[trial*dim];
        n2y=FE->var_spaces[1]->dphi[trial*dim+1];
        kij = K1*n2y*v2y
            + K2*((v2*curlnk2+v2x*n3_prev)*(n2x*n3_prev+n2*curlnk2))
            + K3*((v2x*n2_prev+v2*curlnk3)*(n2x*n2_prev+n2*curlnk3)
                  - v2x*n1_prev*(-n2x*n1_prev) - v2*curlnk1*(-n2*curlnk1)
                  + nkcrosscurlnk1*(v2*n2x+v2x*n2))
            + 2*v2*n2*lam_prev;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += FE->var_spaces[1]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
        n3=FE->var_spaces[2]->phi[trial];
        n3x=FE->var_spaces[2]->dphi[trial*dim];
        n3y=FE->var_spaces[2]->dphi[trial*dim+1];
        kij = K2 * ((v2*curlnk2+v2x*n3_prev)
                 * (n3y*n1_prev-n3x*n2_prev+n3*curlnk3)
                 + (v2x*n3-v2*n3x)*nkdotcurlnk)
            + K3 * ((v2x*n2_prev+v2*curlnk3)*(n3x*n3_prev-n3*curlnk2)
                    - v2x*n1_prev*(n3y*n3_prev+n3*curlnk1)
                    - v2*curlnk1*(n3x*n1_prev-n3y*n2_prev)
                    - nkcrosscurlnk3*(v2*n3y));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column
      local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
        lam=FE->var_spaces[3]->phi[trial];
        kij = 2*n2_prev*lam*v2;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v2 part of RHS block
      rij = - K1 * divnk*v2y
            - K2 * (v2*curlnk2+v2x*n3_prev)*nkdotcurlnk
            - K3 * (nkcrosscurlnk1*(v2x*n2_prev+v2*curlnk3)
                    - nkcrosscurlnk2*(v2x*n1_prev)
                    - nkcrosscurlnk3*(v2*curlnk1))
            - 2*lam_prev*n2_prev*v2;
      bLoc[(local_row_index+test)] += w*rij;
    }

    // v3 block row
    local_row_index += FE->var_spaces[1]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[2]->dof_per_elm;test++) {
      v3=FE->var_spaces[2]->phi[test];
      v3x=FE->var_spaces[2]->dphi[test*dim];
      v3y=FE->var_spaces[2]->dphi[test*dim+1];

      // n1 block Columns
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++) {
        n1=FE->var_spaces[0]->phi[trial];
        n1x=FE->var_spaces[0]->dphi[trial*dim];
        n1y=FE->var_spaces[0]->dphi[trial*dim+1];
        kij = K2 * ((v3*curlnk3+v3y*n1_prev - v3x*n2_prev)
                 * (-n1y*n3_prev+n1*curlnk1) + nkdotcurlnk*(v3y*n1-v3*n1y))
            + K3 * ((v3x*n3_prev-v3*curlnk2)*(-n1y*n2_prev)
                 + (v3y*n3_prev+v3*curlnk1)*(n1y*n1_prev-n1*curlnk3)
                 + (-v3x*n1_prev-v3y*n2_prev)*(n1*curlnk2)
                 - nkcrosscurlnk3*(v3x*n1));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += FE->var_spaces[0]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        n2=FE->var_spaces[1]->phi[trial];
        n2x=FE->var_spaces[1]->dphi[trial*dim];
        n2y=FE->var_spaces[1]->dphi[trial*dim+1];
        kij = K2 * ((v3*curlnk3+v3y*n1_prev-v3x*n2_prev)
                 * (n2x*n3_prev+n2*curlnk2) + nkdotcurlnk*(v3*n2x-v3x*n2))
            + K3 * ((v3x*n3_prev-v3*curlnk2)*(n2x*n2_prev+n2*curlnk3)
                 + (v3y*n3_prev+v3*curlnk1)*(-n2x*n1_prev)
                 + (-v3x*n1_prev-v3y*n2_prev)*(-n2*curlnk1)
                 - nkcrosscurlnk3*(v3y*n2));
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += FE->var_spaces[1]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
        n3=FE->var_spaces[2]->phi[trial];
        n3x=FE->var_spaces[2]->dphi[trial*dim];
        n3y=FE->var_spaces[2]->dphi[trial*dim+1];
        kij = K2 * ((v3*curlnk3+v3y*n1_prev-v3x*n2_prev)
                 * (n3y*n1_prev-n3x*n2_prev+n3*curlnk3))
            + K3 * ((v3x*n3_prev-v3*curlnk2)*(n3x*n3_prev-n3*curlnk2)
                 + (v3y*n3_prev+v3*curlnk1)*(n3y*n3_prev+n3*curlnk1)
                 + (-v3x*n1_prev-v3y*n2_prev)*(-n3x*n1_prev-n3y*n2_prev)
                 + nkcrosscurlnk1*(v3*n3x+v3x*n3)
                 + nkcrosscurlnk2*(v3*n3y+v3y*n3))
            + 2 * lam_prev*v3*n3;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // lam block column
      local_col_index += FE->var_spaces[2]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[3]->dof_per_elm;trial++){
        lam=FE->var_spaces[3]->phi[trial];
        kij = 2*n3_prev*lam*v3;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // v3 part of RHS block
      rij = - K2 * (v3y*n1_prev-v3x*n2_prev+v3*curlnk3)*nkdotcurlnk
            - K3 * (nkcrosscurlnk1*(v3x*n3_prev-v3*curlnk2)
                 + nkcrosscurlnk2*(v3y*n3_prev+v3*curlnk1)
                 + nkcrosscurlnk3*(-v3x*n1_prev-v3y*n2_prev))
            - 2 * lam_prev*n3_prev*v3;
      bLoc[(local_row_index+test)] += w*rij;
    }

    // gam block row
    local_row_index += FE->var_spaces[2]->dof_per_elm;
    // Loop over Test Functions (Rows)
    for (test=0; test<FE->var_spaces[3]->dof_per_elm;test++){
      gam=FE->var_spaces[3]->phi[test];

      // n1 block Column
      local_col_index = 0;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[0]->dof_per_elm;trial++){
        n1=FE->var_spaces[0]->phi[trial];
        kij = 2*n1_prev*n1*gam;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n2 block column
      local_col_index += FE->var_spaces[0]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[1]->dof_per_elm;trial++){
        n2=FE->var_spaces[1]->phi[trial];
        kij = 2*n2_prev*n2*gam;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // n3 block column
      local_col_index += FE->var_spaces[1]->dof_per_elm;
      // Loop over Trial Functions (Columns)
      for (trial=0; trial<FE->var_spaces[2]->dof_per_elm;trial++){
        n3=FE->var_spaces[2]->phi[trial];
        kij = 2*n3_prev*n3*gam;
        ALoc[(local_row_index+test)*dof_per_elm + (local_col_index+trial)] += w*kij;
      }

      // gamma part of RHS block
      rij = -gam * (n1_prev*n1_prev+n2_prev*n2_prev+n3_prev*n3_prev - 1.0);
      bLoc[(local_row_index+test)] += w*rij;
    }
  }

  // Free stuff
  if (qx) free(qx);
  if (K) free(K);
  if(dn1_prev) free(dn1_prev);
  if(dn2_prev) free(dn2_prev);
  if(dn3_prev) free(dn3_prev);

  return;
}

// --------------------------------------------------------
// RELABELING OUR BOUNDARIES FOR PERIODIC

// So for relabeling, I only need the finite element space.
// Relabel boundary flags we'll have the following labels:
// 1: x=0  2: x=1   3: y=0    4: y=1
void relabel_boundary(fespace* FE) {

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
    }
  }
  return;
}
