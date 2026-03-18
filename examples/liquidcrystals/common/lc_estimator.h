/*! \file lc_estimator.h
*  INCLUDE IN HAZMATH
*
*  Copyright 2015_HAZMATH__. All rights reserved.
*
* \brief This contains all the routines needed to compute the Error Estimator
*        for the LC Elastic examples, as well as some other temporary HAZMATH
*        routines that are needed for mesh generation.  These, once verified,
*        will go into hazmath directly.
*
* \note Works in 2D and 3D
*
*/

/*!
* \fn void LCerror_estimator(REAL* est,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq)
*
* \brief Compute the error estimator for each element
*        Th_T^2 = hT^2 ||−K1 grad(div(n)) + K3 curl(Z(n)curl(n)) + (K2 − K3)(n*curl(n))curl(n) + 2 K2 q0 curl(n)) + lam*n||_0,T^2
*               + ||n*n-1||_0,T^2 + sum_edgesonT ( hF || [[K1 div(n)N_f + K3 Z(n)curl(n) x N_f]] ||_0,f^2
*        Z(n) = I - (1-K2/K3)nn^T
*
* \param u 	        FE Approximation
* \param FE          FE Space
* \param mesh        Mesh Data
* \param cq          Quadrature
*
* \return est        error estimator on each element as array
*
*/
void LCerror_estimator(REAL* est,REAL *u,block_fespace *FE,mesh_struct *mesh,qcoordinates *cq) {

  // Counters
  INT elm,neighborelm,quad,i,j,rowa,rowb,jcntr,face,iface;
  INT haveneighbor = 0;
  REAL est_on_elm = 0.0;

  // Mesh Stuff
  INT dim = mesh->dim;
  INT twoders = 3*(dim-1); // uxx, uxy, uyy, uxz, uyz, uzz
  INT v_per_elm = mesh->v_per_elm;
  INT* v_on_elm = (INT *) calloc(v_per_elm,sizeof(INT));
  INT* v_on_elm_neighbor = (INT *) calloc(v_per_elm,sizeof(INT));
  INT f_per_elm = mesh->f_per_elm;
  INT* f_on_elm = (INT *) calloc(f_per_elm,sizeof(INT));
  REAL hT; // Element diamter
  REAL hF; // Face diameter
  iCSRmat* f_el = (iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  icsr_trans(mesh->el_f,f_el);

  // Quadrature Weights and Nodes
  REAL w, wface;
  REAL* qx = (REAL *) calloc(dim,sizeof(REAL));
  REAL* qxf = (REAL *) calloc(dim,sizeof(REAL));
  qcoordinates *cq_face = allocateqcoords_bdry(cq->nq1d,1,dim,2);

  // FEM Stuff
  INT nspaces = FE->nspaces;
  INT dof_per_elm = 0;
  INT dof_per_face = 0;
  for(i=0;i<FE->nspaces;i++) {
    dof_per_elm += FE->var_spaces[i]->dof_per_elm;
    dof_per_face += FE->var_spaces[i]->dof_per_face;
  }
  INT* dof_on_elm = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* dof_on_elm_neighbor = (INT *) calloc(dof_per_elm,sizeof(INT));
  INT* local_dof_on_elm;
  INT* local_dof_on_elm_neighbor;

  // Solution Stuff
  REAL* local_uprev = NULL;
  REAL n1, n1neigh;
  REAL n2, n2neigh;
  REAL n3, n3neigh;
  REAL lam;
  REAL* dn1 = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn2 = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn3 = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn1neigh = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn2neigh = (REAL *) calloc( dim, sizeof(REAL));
  REAL* dn3neigh = (REAL *) calloc( dim, sizeof(REAL));
  REAL* ddn1 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL* ddn2 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL* ddn3 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL divn,divnneigh,ndotcurln,ndotnm1;
  REAL* curln = (REAL *) calloc( 3, sizeof(REAL));
  REAL* curlnneigh = (REAL *) calloc( 3, sizeof(REAL));
  REAL* ndotcurlncurln = (REAL *) calloc( 3, sizeof(REAL));
  REAL* graddiv = (REAL *) calloc(3,sizeof(REAL));
  REAL* curlnx = (REAL *) calloc(3,sizeof(REAL));
  REAL* curlny = (REAL *) calloc(3,sizeof(REAL));
  REAL* curlnz = (REAL *) calloc(3,sizeof(REAL));
  REAL* cZc = (REAL *) calloc(3,sizeof(REAL));
  REAL* term1 = (REAL *) calloc(3,sizeof(REAL));
  REAL* term2 = (REAL *) calloc(3,sizeof(REAL));
  REAL Zn11,Zn12,Zn13,Zn22,Zn23,Zn33;
  REAL Zn11neigh,Zn12neigh,Zn13neigh,Zn22neigh,Zn23neigh,Zn33neigh;
  REAL Zn11x,Zn12x,Zn13x,Zn22x,Zn23x,Zn33x;
  REAL Zn11y,Zn12y,Zn13y,Zn22y,Zn23y,Zn33y;
  REAL Zn11z,Zn12z,Zn13z,Zn22z,Zn23z,Zn33z;
  REAL* Zcurln = (REAL *) calloc(3,sizeof(REAL));
  REAL* Zcurlnneigh = (REAL *) calloc(3,sizeof(REAL));

  // Frank coefficients
  REAL* K = (REAL *) calloc(3,sizeof(REAL));
  // Coefficients
  get_frank_constants(K);
  REAL K1 = K[0];
  REAL K2 = K[1];
  REAL K3 = K[2];
  REAL kap = K2/K3;
  REAL onemkap = 1.0 - kap;

  /* Loop over all Elements */
  for (elm=0; elm<mesh->nelm; elm++) {
    est_on_elm = 0.0;

    // Get cell diameter
    hT = pow(dim*(dim-1)*mesh->el_vol[elm],1.0/dim);

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

    // Interior calculations first
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
      P2_2ndDerivativeInterpolation(ddn1,local_uprev,qx,dof_on_elm,FE->var_spaces[0],mesh);
      local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
      local_uprev+=FE->var_spaces[0]->ndof;

      FE_Interpolation(&n2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      FE_DerivativeInterpolation(dn2,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
      P2_2ndDerivativeInterpolation(ddn2,local_uprev,qx,dof_on_elm,FE->var_spaces[1],mesh);
      local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
      local_uprev+=FE->var_spaces[1]->ndof;

      FE_Interpolation(&n3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
      FE_DerivativeInterpolation(dn3,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
      P2_2ndDerivativeInterpolation(ddn3,local_uprev,qx,dof_on_elm,FE->var_spaces[2],mesh);
      local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
      local_uprev+=FE->var_spaces[2]->ndof;

      // lambda
      FE_Interpolation(&lam,local_uprev,qx,local_dof_on_elm,v_on_elm,FE->var_spaces[3],mesh);

      // Compute divs and curls
      divn = dn1[0] + dn2[1];
      if(dim==3) divn += dn3[2];
      curln[0] = dn3[1];
      curln[1] = -dn3[0];
      curln[2] = dn2[0] - dn1[1];
      if(dim==3) {
        curln[0] += -dn2[2];
        curln[1] += dn1[2];
      }
      ndotcurln = n1*curln[0]+n2*curln[1]+n3*curln[2];
      for(i=0;i<dim;i++) ndotcurlncurln[i] = ndotcurln*curln[i];
      ndotnm1 = (n1*n1 + n2*n2 + n3*n3) - 1.0;
      // Z's
      Zn11 = 1.0 - onemkap*n1*n1;
      Zn12 = -onemkap*n1*n2;
      Zn13 = -onemkap*n1*n3;
      Zn22 = 1.0 - onemkap*n2*n2;
      Zn23 = -onemkap*n2*n3;
      Zn33 = 1.0 - onemkap*n3*n3;
      Zn11x = -onemkap*2.0*n1*dn1[0];
      Zn11y = -onemkap*2.0*n1*dn1[1];
      Zn12x = -onemkap*(n1*dn2[0]+dn1[0]*n2);
      Zn12y = -onemkap*(n1*dn2[1]+dn1[1]*n2);
      Zn13x = -onemkap*(n1*dn3[0]+dn1[0]*n3);
      Zn13y = -onemkap*(n1*dn3[1]+dn1[1]*n3);
      Zn22x = -onemkap*2.0*n2*dn2[0];
      Zn22y = -onemkap*2.0*n2*dn2[1];
      Zn23x = -onemkap*(n2*dn3[0]+dn2[0]*n3);
      Zn23y = -onemkap*(n2*dn3[1]+dn2[1]*n3);
      Zn33x = -onemkap*2.0*n3*dn3[0];
      Zn33y = -onemkap*2.0*n3*dn3[1];
      if(dim==3) {
        Zn11z = -onemkap*2.0*n1*dn1[2];
        Zn12z = -onemkap*(n1*dn2[2]+dn1[2]*n2);
        Zn13z = -onemkap*(n1*dn3[2]+dn1[2]*n3);
        Zn22z = -onemkap*2.0*n2*dn2[2];
        Zn23z = -onemkap*(n2*dn3[2]+dn2[2]*n3);
        Zn33z = -onemkap*2.0*n3*dn3[2];
      }
      // 2nd derivatives
      graddiv[0] = ddn1[0] + ddn2[1];
      graddiv[1] = ddn1[1] + ddn2[2];
      graddiv[2] = 0.0;
      curlnx[0] = ddn3[1];
      curlnx[1] = -ddn3[0];
      curlnx[2] = ddn2[0] - ddn1[1];
      curlny[0] = ddn3[2];
      curlny[1] = -ddn3[1];
      curlny[2] = ddn2[1] - ddn1[2];
      if(dim==3) {
        graddiv[0] += ddn3[3];
        graddiv[1] += ddn3[4];
        graddiv[2] += ddn1[3] + ddn2[4] + ddn3[5];
        curlnx[0] += -ddn2[3];
        curlnx[1] += ddn1[3];
        curlny[0] += -ddn2[4];
        curlny[1] += ddn1[4];
        curlnz[0] = ddn3[4] - ddn2[5];
        curlnz[1] = ddn1[5] - ddn3[3];
        curlnz[2] = ddn2[3] - ddn1[4];
      }
      cZc[0] = (Zn13y*curln[0] + Zn13*curlny[0] + Zn23y*curln[1] + Zn23*curlny[1] + Zn33y*curln[2] + Zn33*curlny[2]);
      cZc[1] = -(Zn13x*curln[0] + Zn13*curlnx[0] + Zn23x*curln[1] + Zn23*curlnx[1] + Zn33x*curln[2] + Zn33*curlnx[2]);
      cZc[2] = (Zn12x*curln[0] + Zn12*curlnx[0] + Zn22x*curln[1] + Zn22*curlnx[1] + Zn23x*curln[2] + Zn23*curlnx[2]) -
      (Zn11y*curln[0] + Zn11*curlny[0] + Zn12y*curln[1] + Zn12*curlny[1] + Zn13y*curln[2] + Zn13*curlny[2]);
      if(dim==3) {
        cZc[0] += -(Zn12z*curln[0] + Zn12*curlnz[0] + Zn22z*curln[1] + Zn22*curlnz[1] + Zn23z*curln[2] + Zn23*curlnz[2]);
        cZc[1] += (Zn11z*curln[0] + Zn11*curlnz[0] + Zn12z*curln[1] + Zn12*curlnz[1] + Zn13z*curln[2] + Zn13*curlnz[2]);
      }

      // Computing element estimator
      // hT^2 ||−K1 grad(div(n)) + K3 curl(Z(n)curl(n)) + (K2 − K3)(n*curl(n))curl(n) + 2 K2 q0 curl(n)) + lam*n||_0,T^2
      term1[0] = -K1*graddiv[0] + K3*cZc[0] + (K2-K3)*ndotcurlncurln[0] + lam*n1;
      term1[1] = -K1*graddiv[1] + K3*cZc[1] + (K2-K3)*ndotcurlncurln[1] + lam*n2;
      term1[2] = -K1*graddiv[2] + K3*cZc[2] + (K2-K3)*ndotcurlncurln[2] + lam*n3;

      est_on_elm+=w*(hT*hT*(term1[0]*term1[0]+term1[1]*term1[1]+term1[2]*term1[2]) + ndotnm1*ndotnm1);
    }

    // Face terms
    // Find Faces for given Element
    get_incidence_row(elm,mesh->el_f,f_on_elm);

    // Loop over each face and integrate
    for(iface=0;iface<f_per_elm;iface++) {
      face = f_on_elm[iface];

      // Get face diamter (length in 2D, sqrt(area/2) in 3D)
      hF = mesh->f_area[face];
      if(dim==3) hF = sqrt(2.0*hF);

      // Get quadrature on face
      quad_face(cq_face,mesh,cq->nq1d,face);

      // Need to find connecting elements on face if not on boundary
      rowa = f_el->IA[face];
      rowb = f_el->IA[face+1];
      if(rowb-rowa==1) {
        haveneighbor=0;
      } else if (rowb-rowa==2) {
        haveneighbor=1;
      } else {
        printf("\n!!!Something wrong with element neighbor count!!!\n");
      }

      if(haveneighbor) {
        for (j=rowa; j<rowb; j++) {
          if(f_el->JA[j]!=elm) {
            neighborelm = f_el->JA[j];
          }
        }

        // Get DoF on neighbor element
        jcntr = 0;
        for(i=0;i<nspaces;i++) {
          rowa = FE->var_spaces[i]->el_dof->IA[neighborelm];
          rowb = FE->var_spaces[i]->el_dof->IA[neighborelm+1];
          for (j=rowa; j<rowb; j++) {
            dof_on_elm_neighbor[jcntr] = FE->var_spaces[i]->el_dof->JA[j];
            jcntr++;
          }
        }

        // Find Vertices for neighbor element
        get_incidence_row(neighborelm,mesh->el_v,v_on_elm_neighbor);
      }

      // Loop over quadrature on face
      for (quad=0;quad<cq_face->nq_per_elm;quad++) {
        qxf[0] = cq_face->x[quad];
        if(mesh->dim==2 || mesh->dim==3)
        qxf[1] = cq_face->y[quad];
        if(mesh->dim==3)
        qxf[2] = cq_face->z[quad];
        wface = cq->w[quad];

        // Interpolate FE solution to quadrature point
        //  Get current solutions at each quadrature node
        // n
        local_uprev = u;
        FE_Interpolation(&n1,local_uprev,qxf,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
        FE_DerivativeInterpolation(dn1,local_uprev,qxf,dof_on_elm,v_on_elm,FE->var_spaces[0],mesh);
        local_dof_on_elm = dof_on_elm + FE->var_spaces[0]->dof_per_elm;
        if(haveneighbor) {
          FE_Interpolation(&n1neigh,local_uprev,qxf,dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[0],mesh);
          FE_DerivativeInterpolation(dn1neigh,local_uprev,qxf,dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[0],mesh);
          local_dof_on_elm_neighbor = dof_on_elm_neighbor + FE->var_spaces[0]->dof_per_elm;
        }
        local_uprev+=FE->var_spaces[0]->ndof;

        FE_Interpolation(&n2,local_uprev,qxf,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
        FE_DerivativeInterpolation(dn2,local_uprev,qxf,local_dof_on_elm,v_on_elm,FE->var_spaces[1],mesh);
        local_dof_on_elm += FE->var_spaces[1]->dof_per_elm;
        if(haveneighbor) {
          FE_Interpolation(&n2neigh,local_uprev,qxf,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[1],mesh);
          FE_DerivativeInterpolation(dn2neigh,local_uprev,qxf,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[1],mesh);
          local_dof_on_elm_neighbor += FE->var_spaces[1]->dof_per_elm;
        }
        local_uprev+=FE->var_spaces[1]->ndof;

        FE_Interpolation(&n3,local_uprev,qxf,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
        FE_DerivativeInterpolation(dn3,local_uprev,qxf,local_dof_on_elm,v_on_elm,FE->var_spaces[2],mesh);
        local_dof_on_elm += FE->var_spaces[2]->dof_per_elm;
        if(haveneighbor) {
          FE_Interpolation(&n3neigh,local_uprev,qxf,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[2],mesh);
          FE_DerivativeInterpolation(dn3neigh,local_uprev,qxf,local_dof_on_elm_neighbor,v_on_elm_neighbor,FE->var_spaces[2],mesh);
          local_dof_on_elm_neighbor += FE->var_spaces[2]->dof_per_elm;
        }
        local_uprev+=FE->var_spaces[2]->ndof;

        // Compute divs and curls
        divn = dn1[0] + dn2[1];
        if(dim==3) divn += dn3[2];
        curln[0] = dn3[1];
        curln[1] = -dn3[0];
        curln[2] = dn2[0] - dn1[1];
        if(dim==3) {
          curln[0] += -dn2[2];
          curln[1] += dn1[2];
        }
        // Z's
        Zn11 = 1.0 - onemkap*n1*n1;
        Zn12 = -onemkap*n1*n2;
        Zn13 = -onemkap*n1*n3;
        Zn22 = 1.0 - onemkap*n2*n2;
        Zn23 = -onemkap*n2*n3;
        Zn33 = 1.0 - onemkap*n3*n3;
        // Z*curln
        Zcurln[0] = Zn11*curln[0] + Zn12*curln[1] + Zn13*curln[2];
        Zcurln[1] = Zn12*curln[0] + Zn22*curln[1] + Zn23*curln[2];
        Zcurln[2] = Zn13*curln[0] + Zn23*curln[1] + Zn33*curln[2];

        // Compute divs and curls on neighbor
        if(haveneighbor) {
          divnneigh = dn1neigh[0] + dn2neigh[1];
          if(dim==3) divnneigh += dn3neigh[2];
          curlnneigh[0] = dn3neigh[1];
          curlnneigh[1] = -dn3neigh[0];
          curlnneigh[2] = dn2neigh[0] - dn1neigh[1];
          if(dim==3) {
            curlnneigh[0] += -dn2neigh[2];
            curlnneigh[1] += dn1neigh[2];
          }
          // Z's
          Zn11neigh = 1.0 - onemkap*n1neigh*n1neigh;
          Zn12neigh = -onemkap*n1neigh*n2neigh;
          Zn13neigh = -onemkap*n1neigh*n3neigh;
          Zn22neigh = 1.0 - onemkap*n2neigh*n2neigh;
          Zn23neigh = -onemkap*n2neigh*n3neigh;
          Zn33neigh = 1.0 - onemkap*n3neigh*n3neigh;
          // Z*curln
          Zcurlnneigh[0] = Zn11neigh*curlnneigh[0] + Zn12neigh*curlnneigh[1] + Zn13neigh*curlnneigh[2];
          Zcurlnneigh[1] = Zn12neigh*curlnneigh[0] + Zn22neigh*curlnneigh[1] + Zn23neigh*curlnneigh[2];
          Zcurlnneigh[2] = Zn13neigh*curlnneigh[0] + Zn23neigh*curlnneigh[1] + Zn33neigh*curlnneigh[2];
        }

        // Computing face estimator
        // hF || [[K1 div(n)norm(f) + K3 (Z*curl(n))xnorm(f)]]||_f^2
        term2[0] = K1*divn*mesh->f_norm[face*dim+0] + K3*(Zcurln[1]*mesh->f_norm[face*dim+2] - Zcurln[2]*mesh->f_norm[face*dim+1]);
        term2[1] = K1*divn*mesh->f_norm[face*dim+1] + K3*(Zcurln[2]*mesh->f_norm[face*dim+0] - Zcurln[0]*mesh->f_norm[face*dim+2]);
        term2[2] = K1*divn*mesh->f_norm[face*dim+2] + K3*(Zcurln[0]*mesh->f_norm[face*dim+1] - Zcurln[1]*mesh->f_norm[face*dim+0]);
        if(haveneighbor) {
          term2[0] += -(K1*divnneigh*mesh->f_norm[face*dim+0] + K3*(Zcurlnneigh[1]*mesh->f_norm[face*dim+2] - Zcurlnneigh[2]*mesh->f_norm[face*dim+1]));
          term2[1] += -(K1*divnneigh*mesh->f_norm[face*dim+1] + K3*(Zcurlnneigh[2]*mesh->f_norm[face*dim+0] - Zcurlnneigh[0]*mesh->f_norm[face*dim+2]));
          term2[2] += -(K1*divnneigh*mesh->f_norm[face*dim+2] + K3*(Zcurlnneigh[0]*mesh->f_norm[face*dim+1] - Zcurlnneigh[1]*mesh->f_norm[face*dim+0]));
        } else {
          term2[0] = 0.0;
          term2[1] = 0.0;
          term2[2] = 0.0;
        }

        est_on_elm+=wface*hF*(term2[0]*term2[0] + term2[1]*term2[1] + term2[2]*term2[2]);

      }
      est[elm] = sqrt(est_on_elm);
    }
  }

  if(dof_on_elm) free(dof_on_elm);
  if(dof_on_elm_neighbor) free(dof_on_elm_neighbor);
  if(v_on_elm) free(v_on_elm);
  if(v_on_elm_neighbor) free(v_on_elm_neighbor);
  if(f_on_elm) free(f_on_elm);
  icsr_free(f_el);
  if(qx) free(qx);
  if(qxf) free(qxf);
  if(dn1) free(dn1);
  if(dn2) free(dn2);
  if(dn3) free(dn3);
  if(dn1neigh) free(dn1neigh);
  if(dn2neigh) free(dn2neigh);
  if(dn3neigh) free(dn3neigh);
  if(ddn1) free(ddn1);
  if(ddn2) free(ddn2);
  if(ddn3) free(ddn3);
  if(curln) free(curln);
  if(curlnneigh) free(curlnneigh);
  if(cZc) free(cZc);
  if(Zcurln) free(Zcurln);
  if(Zcurlnneigh) free(Zcurlnneigh);
  if(curlnx) free(curlnx);
  if(curlny) free(curlny);
  if(curlnz) free(curlnz);
  if(graddiv) free(graddiv);
  if(ndotcurlncurln) free(ndotcurlncurln);
  if(term1) free(term1);
  if(term2) free(term2);
  if(K) free(K);
  if(cq_face){
    free_qcoords(cq_face);
    free(cq_face);
    cq_face=NULL;
  }

  return;
}

// MARKING Routines

// /*!
// * \fn ivector maximal_mark(scomplex *sc,dvector *estimator, REAL gamma)
// *
// * \brief mark elements for adaptive refinement using maximal marking strategy
// *           mark all elements with error estimator greater than gamma*max(estimator)
// *
// * \param sc 	       Pointer to the simplicial complex
// * \param estimator   Pointer to error estimator
// * \param gamma       Fraction of max error estimator to mark (eg. 0 -> mark all, 1 -> mark none)
// *
// * \return marked     ivector of elements to be marked, TRUE for marked
// *
// */
// ivector maximal_mark(scomplex *sc, REAL *estimator, REAL gamma) {
  
//   INT i,k=0;
//   INT ns = sc->ns;
//   REAL errmax;
//   ivector marked=ivec_create(sc->ns);
//   errmax=estimator[0];
//   for(i=1;i<ns;i++) {
//     if(fabs(estimator[i])>errmax) {
//       errmax=estimator[i];
//     }
//   }
//   for(i=0;i<ns;i++) {
//     if(fabs(estimator[i])>gamma*errmax){
//       marked.val[i]=TRUE;
//       k++;
//     } else {
//       marked.val[i]=FALSE;
//     }
//   }
  
//   printf("\nMarking elements for AMR using Maximal Mark strategy.\n\tMarking elements with %3.2f%% of the max error: --> %d of %d elements were marked to be refined.\n",gamma*100,k,ns);
//   return marked;
// }

// // Simplex complex routines that will eventually go into HAZMATH itself

// /*!
//  * \fn void scfinalize_temp(scomplex *sc,const INT set_bndry_codes)
//  *
//  * \brief Removes all hierachy and make sc to represent only the final
//  *        grid. computes connected components and connected components
//  *        on the boundary.
//  *
//  * \param sc: simplicial complex
//  * \param set_bndry_codes: if 0 then create the sparse matrix for all vertices;
//  *
//  * \return
//  *
//  * \note
//  *
//  * \author ludmil (20151010)
//  * \modified ludmil (20210831)
//  * \modified ludmil (20211121)
//  *
//  */
// void scfinalize_temp(scomplex *sc,const INT set_bndry_codes)
// {
//   // INT n=sc->n;
//   INT ns,j=-10,k=-10;
//   INT n1=sc->n+1;
//   /*
//       store the finest mesh in sc structure.
//       on input sc has all the hierarchy, on return sc only has the final mesh.
//   */
//   //  free(sc->parent_v->val);  sc->parent_v->val=NULL;
//   ns=0;
//   for (j=0;j<sc->ns;j++){
//     /*
//       On the last grid are all simplices that were not refined, so
//       these are the ones for which child0 and childn are not set.
//     */
//     if(sc->child0[j]<0 || sc->childn[j]<0){
//       for (k=0;k<n1;k++) {
// 	sc->nodes[ns*n1+k]=sc->nodes[j*n1+k];
//       }
//       sc->child0[ns]=-1;
//       sc->childn[ns]=-1;
//       sc->gen[ns]=sc->gen[j];
//       sc->flags[ns]=sc->flags[j];
//       ns++;
//     }
//   }
//   sc->ns=ns;
//   sc->nodes=realloc(sc->nodes,n1*sc->ns*sizeof(INT));
//   sc->nbr=realloc(sc->nbr,n1*sc->ns*sizeof(INT));
//   sc->vols=realloc(sc->vols,sc->ns*sizeof(REAL));
//   sc->child0=realloc(sc->child0,sc->ns*sizeof(INT));
//   sc->childn=realloc(sc->childn,sc->ns*sizeof(INT));
//   sc->gen=realloc(sc->gen,sc->ns*sizeof(INT));
//   sc->flags=realloc(sc->flags,sc->ns*sizeof(INT));
//   find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
//   // this also can be called separately
//   // set_bndry_codes should always be set to 1.
//   //  set_bndry_codes=1;
//   find_cc_bndry_cc(sc,(INT )1); //set_bndry_codes);
//   //
//   /* if(set_bndry_codes){ */
//   /*   for(j=0;j<sc->nv;++j){ */
//   /*     if(sc->bndry[j]>128) sc->bndry[j]-=128; */
//   /*   } */
//   /* } */
//   // clean up: TODO: DO NOT FREE ANYTHING UNTIL LUDMIL FIXES!
//   // icsr_free(sc->bndry_v);
//   // free(sc->bndry_v);
//   // sc->bndry_v=NULL;
//   // icsr_free(sc->parent_v);
//   // free(sc->parent_v);
//   // sc->parent_v=NULL;
//   return;
// }

// /*!
//  * \fn void get_initial_mesh_ni(nested_it* ni,INT dim,INT init_ref_levels)
//  *
//  * \brief Gets an initial uniform mesh to start the nested iteration process
//  *
//  * \param ni: nested iteration struct containing current mesh and simplicial complexes
//  * \param dim:  Dimension of problem
//  * \param init_ref_levels: Number of initial uniform refinements to get to first mesh
//  *
//  * \return ni: creates the initial simplicial complex and mesh for nested iteration
//  *
//  */
// void get_initial_mesh(nested_it* ni,INT dim,INT init_ref_levels)
// {

//   INT jlevel;
//   // Get the coarsest mesh on the cube in dimension dim and set the refinement type.
//   ni->sc_all=mesh_cube_init(dim,1,11);
//   scomplex* sc=ni->sc_all[0];
  
//   // uniform refine
//   if(dim==2){
//     for(jlevel=0;jlevel<init_ref_levels;++jlevel){
//       uniformrefine2d(sc);
//       sc_vols(sc);
//     }
//   } else if(dim==3){
//     for(jlevel=0;jlevel<init_ref_levels;++jlevel){
//       uniformrefine3d(sc);
//       sc_vols(sc);
//     }
//   } else {
//     check_error(ERROR_DIM, __FUNCTION__);
//   }
//   // Get boundary codes TODO: LTZ check on scfinalize version
//   scfinalize_temp(sc,(INT )1);
//   sc_vols(sc);

//   // Convert to mesh_struct for FEM assembly
//   ni->mesh=malloc(sizeof(mesh_struct));
//   ni->mesh[0]=sc2mesh(sc);
//   build_mesh_all(ni->mesh);

//   return;
// }
