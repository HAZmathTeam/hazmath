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
* \fn void LCerror_estimator(REAL* est,REAL *u,block_fespace *FE,scomplex *sc,qcoordinates *cq)
*
* \brief Compute the error estimator for each element
*        Th_T^2 = hT^2 ||−K1 grad(div(n)) + K3 curl(Z(n)curl(n)) + (K2 − K3)(n*curl(n))curl(n) + 2 K2 q0 curl(n)) + lam*n||_0,T^2
*               + ||n*n-1||_0,T^2 + sum_edgesonT ( hF || [[K1 div(n)N_f + K3 Z(n)curl(n) x N_f]] ||_0,f^2
*        Z(n) = I - (1-K2/K3)nn^T
*
* \param u 	        FE Approximation
* \param FE          FE Space
* \param sc          Simplicial Complex
* \param cq          Quadrature
*
* \return est        error estimator on each element as array
*
*/
void LCerror_estimator(REAL* est,REAL *u,block_fespace *FE,scomplex *sc,qcoordinates *cq) {

  // Counters
  INT elm,neighborelm,quad,i,j,d,rowa,rowb,jcntr,face,iface;
  INT haveneighbor = 0;
  REAL est_on_elm = 0.0;

  // Mesh and FE local data
  INT dim = sc->dim;
  INT twoders = 3*(dim-1); // uxx, uxy, uyy, uxz, uyz, uzz
  INT f_per_elm = dim+1;
  INT* f_on_elm = (INT *) calloc(f_per_elm,sizeof(INT));
  REAL hT; // Element diamter
  REAL hF; // Face diameter
  iCSRmat* f_el = (iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  icsr_trans(sc->fem->el_f,f_el);

  // Local data for element and neighbor
  simplex_local_data elm_data, nbr_data;
  fe_local_data fe_data, fe_nbr_data;
  memset(&elm_data, 0, sizeof(simplex_local_data));
  memset(&fe_data, 0, sizeof(fe_local_data));
  memset(&nbr_data, 0, sizeof(simplex_local_data));
  memset(&fe_nbr_data, 0, sizeof(fe_local_data));
  initialize_localdata_elm(&elm_data, &fe_data, sc, FE, cq->nq1d);
  initialize_localdata_elm(&nbr_data, &fe_nbr_data, sc, FE, cq->nq1d);

  // Wrap REAL* u in a temporary dvector for get_felocaldata_elm
  dvector sol_dvec;
  sol_dvec.row = 0;
  for (i = 0; i < FE->nspaces; i++) sol_dvec.row += FE->var_spaces[i]->ndof;
  sol_dvec.val = u;

  // Per-space DOF counts and offsets
  INT n1dofpelm = fe_data.n_dof_per_space[0];
  INT n2dofpelm = fe_data.n_dof_per_space[1];
  INT n3dofpelm = fe_data.n_dof_per_space[2];
  INT lamdofpelm = fe_data.n_dof_per_space[3];
  INT offset_n1 = 0;
  INT offset_n2 = n1dofpelm;
  INT offset_n3 = n1dofpelm + n2dofpelm;
  INT offset_lam = n1dofpelm + n2dofpelm + n3dofpelm;

  // Quadrature Weights and Nodes
  REAL w, wface;
  REAL qx[3], qxf[3];
  qcoordinates *cq_face = allocateqcoords_bdry(cq->nq1d,1,dim,2);

  // Solution Stuff
  REAL n1, n1neigh;
  REAL n2, n2neigh;
  REAL n3, n3neigh;
  REAL lam;
  REAL dn1[3], dn2[3], dn3[3];
  REAL dn1neigh[3], dn2neigh[3], dn3neigh[3];
  REAL* ddn1 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL* ddn2 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL* ddn3 = (REAL *) calloc( twoders, sizeof(REAL));
  REAL divn,divnneigh,ndotcurln,ndotnm1;
  REAL curln[3], curlnneigh[3], ndotcurlncurln[3];
  REAL graddiv[3], curlnx[3], curlny[3], curlnz[3];
  REAL cZc[3], term1[3], term2[3];
  REAL Zn11,Zn12,Zn13,Zn22,Zn23,Zn33;
  REAL Zn11neigh,Zn12neigh,Zn13neigh,Zn22neigh,Zn23neigh,Zn33neigh;
  REAL Zn11x,Zn12x,Zn13x,Zn22x,Zn23x,Zn33x;
  REAL Zn11y,Zn12y,Zn13y,Zn22y,Zn23y,Zn33y;
  REAL Zn11z,Zn12z,Zn13z,Zn22z,Zn23z,Zn33z;
  REAL Zcurln[3], Zcurlnneigh[3];

  // Frank coefficients
  REAL K[3];
  get_frank_constants(K);
  REAL K1 = K[0];
  REAL K2 = K[1];
  REAL K3 = K[2];
  REAL kap = K2/K3;
  REAL onemkap = 1.0 - kap;

  /* Loop over all Elements */
  for (elm=0; elm<sc->fem->ns_leaf; elm++) {
    est_on_elm = 0.0;

    // Get cell diameter
    hT = pow(dim*(dim-1)*sc->fem->el_vol[elm],1.0/dim);

    // Get element local data
    get_elmlocaldata(&elm_data, sc, elm);
    get_felocaldata_elm(&fe_data, FE, &sol_dvec, elm);

    // Interior calculations first
    // Loop over quadrature nodes on element
    for (quad=0;quad<elm_data.quad_local->nq;quad++) {
      for (d = 0; d < dim; d++) qx[d] = elm_data.quad_local->x[quad*dim+d];
      w = elm_data.quad_local->w[quad];

      // Get basis functions for all spaces at this quadrature point
      for (i = 0; i < fe_data.nspaces; i++)
        get_FEM_basis_at_quadpt(&elm_data, &fe_data, i, quad);

      // Interpolate n1, dn1
      n1 = 0.0;
      for (j = 0; j < n1dofpelm; j++) n1 += fe_data.u_local[offset_n1+j] * fe_data.phi[0][j];
      for (d = 0; d < dim; d++) {
        dn1[d] = 0.0;
        for (j = 0; j < n1dofpelm; j++) dn1[d] += fe_data.u_local[offset_n1+j] * fe_data.dphi[0][j*dim+d];
      }
      // P2 2nd derivatives (kept on old scomplex path for now)
      fe_2ndderiv_interp_to_x(ddn1, fe_data.u_local+offset_n1, qx, &fe_data, &elm_data, 0);

      // Interpolate n2, dn2
      n2 = 0.0;
      for (j = 0; j < n2dofpelm; j++) n2 += fe_data.u_local[offset_n2+j] * fe_data.phi[1][j];
      for (d = 0; d < dim; d++) {
        dn2[d] = 0.0;
        for (j = 0; j < n2dofpelm; j++) dn2[d] += fe_data.u_local[offset_n2+j] * fe_data.dphi[1][j*dim+d];
      }
      fe_2ndderiv_interp_to_x(ddn2, fe_data.u_local+offset_n2, qx, &fe_data, &elm_data, 1);

      // Interpolate n3, dn3
      n3 = 0.0;
      for (j = 0; j < n3dofpelm; j++) n3 += fe_data.u_local[offset_n3+j] * fe_data.phi[2][j];
      for (d = 0; d < dim; d++) {
        dn3[d] = 0.0;
        for (j = 0; j < n3dofpelm; j++) dn3[d] += fe_data.u_local[offset_n3+j] * fe_data.dphi[2][j*dim+d];
      }
      fe_2ndderiv_interp_to_x(ddn3, fe_data.u_local+offset_n3, qx, &fe_data, &elm_data, 2);

      // lambda
      lam = 0.0;
      for (j = 0; j < lamdofpelm; j++) lam += fe_data.u_local[offset_lam+j] * fe_data.phi[3][j];

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
    get_incidence_row(elm,sc->fem->el_f,f_on_elm);

    // Loop over each face and integrate
    for(iface=0;iface<f_per_elm;iface++) {
      face = f_on_elm[iface];

      // Get face diamter (length in 2D, sqrt(area/2) in 3D)
      hF = sc->fem->f_area[face];
      if(dim==3) hF = sqrt(2.0*hF);

      // Get quadrature on face
      quad_face(cq_face,sc,cq->nq1d,face);

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
        // Get neighbor local data
        get_elmlocaldata(&nbr_data, sc, neighborelm);
        get_felocaldata_elm(&fe_nbr_data, FE, &sol_dvec, neighborelm);
      }

      // Loop over quadrature on face
      for (quad=0;quad<cq_face->nq_per_elm;quad++) {
        qxf[0] = cq_face->x[quad];
        if(sc->dim==2 || sc->dim==3)
        qxf[1] = cq_face->y[quad];
        if(sc->dim==3)
        qxf[2] = cq_face->z[quad];
        wface = cq->w[quad];

        // Interpolate FE solution on element at face quad point
        fe_interpolation_to_x(&n1, fe_data.u_local+offset_n1, qxf, &fe_data, &elm_data, 0);
        fe_dinterp_to_x(dn1, fe_data.u_local+offset_n1, qxf, &fe_data, &elm_data, 0);
        if(haveneighbor) {
          fe_interpolation_to_x(&n1neigh, fe_nbr_data.u_local+offset_n1, qxf, &fe_nbr_data, &nbr_data, 0);
          fe_dinterp_to_x(dn1neigh, fe_nbr_data.u_local+offset_n1, qxf, &fe_nbr_data, &nbr_data, 0);
        }

        fe_interpolation_to_x(&n2, fe_data.u_local+offset_n2, qxf, &fe_data, &elm_data, 1);
        fe_dinterp_to_x(dn2, fe_data.u_local+offset_n2, qxf, &fe_data, &elm_data, 1);
        if(haveneighbor) {
          fe_interpolation_to_x(&n2neigh, fe_nbr_data.u_local+offset_n2, qxf, &fe_nbr_data, &nbr_data, 1);
          fe_dinterp_to_x(dn2neigh, fe_nbr_data.u_local+offset_n2, qxf, &fe_nbr_data, &nbr_data, 1);
        }

        fe_interpolation_to_x(&n3, fe_data.u_local+offset_n3, qxf, &fe_data, &elm_data, 2);
        fe_dinterp_to_x(dn3, fe_data.u_local+offset_n3, qxf, &fe_data, &elm_data, 2);
        if(haveneighbor) {
          fe_interpolation_to_x(&n3neigh, fe_nbr_data.u_local+offset_n3, qxf, &fe_nbr_data, &nbr_data, 2);
          fe_dinterp_to_x(dn3neigh, fe_nbr_data.u_local+offset_n3, qxf, &fe_nbr_data, &nbr_data, 2);
        }

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
        term2[0] = K1*divn*sc->fem->f_norm[face*dim+0] + K3*(Zcurln[1]*sc->fem->f_norm[face*dim+2] - Zcurln[2]*sc->fem->f_norm[face*dim+1]);
        term2[1] = K1*divn*sc->fem->f_norm[face*dim+1] + K3*(Zcurln[2]*sc->fem->f_norm[face*dim+0] - Zcurln[0]*sc->fem->f_norm[face*dim+2]);
        term2[2] = K1*divn*sc->fem->f_norm[face*dim+2] + K3*(Zcurln[0]*sc->fem->f_norm[face*dim+1] - Zcurln[1]*sc->fem->f_norm[face*dim+0]);
        if(haveneighbor) {
          term2[0] += -(K1*divnneigh*sc->fem->f_norm[face*dim+0] + K3*(Zcurlnneigh[1]*sc->fem->f_norm[face*dim+2] - Zcurlnneigh[2]*sc->fem->f_norm[face*dim+1]));
          term2[1] += -(K1*divnneigh*sc->fem->f_norm[face*dim+1] + K3*(Zcurlnneigh[2]*sc->fem->f_norm[face*dim+0] - Zcurlnneigh[0]*sc->fem->f_norm[face*dim+2]));
          term2[2] += -(K1*divnneigh*sc->fem->f_norm[face*dim+2] + K3*(Zcurlnneigh[0]*sc->fem->f_norm[face*dim+1] - Zcurlnneigh[1]*sc->fem->f_norm[face*dim+0]));
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

  free_simplexlocaldata(&elm_data);
  free_felocaldata(&fe_data);
  free_simplexlocaldata(&nbr_data);
  free_felocaldata(&fe_nbr_data);
  if(f_on_elm) free(f_on_elm);
  icsr_free(f_el);
  if(ddn1) free(ddn1);
  if(ddn2) free(ddn2);
  if(ddn3) free(ddn3);
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
