/*! \file src/solver/gmg_setup.c
 *
 *  Geometric Multigrid: SETUP phase
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 */

#include "hazmath.h"

/***********************************************************************************************/
/**
 * \fn static void build_linear_R (dCSRmat *tentp, INT nf1d, INT nc1d )
 *
 * \brief Build tentative R for piecewise linear elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param nf1d               Number of fine vertices in 1d
 * \param nc1d               Number of coarse vertices in 1d
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
void build_linear_R (dCSRmat *R,
                           INT      nf1d,
                           INT      nc1d)
{
    INT i;
    INT ci, cj;
    INT fdof, cdof;
    INT jstart=0;

    // allocate memory for R
    R->row = nc1d*nc1d;
    R->col = nf1d*nf1d;
    R->nnz = (nc1d-2)*(nc1d-2)*7 + 2*3 + 2*4 + 4*(nc1d-2)*5;
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;

    REAL stencil[] = { 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5 };

    IA[0] = 0;
    // Stencil for piece-wise linear
    // ia  = [ ptr ], generally 7
    // ja  = [ p-n, p-n+1, p-1,   p, p+1, p+n-1, p+n ]
    // val = [ 0.5,   0.5, 0.5, 1.0, 0.5,   0.5, 0.5 ]
    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d; ci++){
        cdof = ci + cj*nc1d;
        fdof = ci*2 + cj*nf1d*2;

        if(ci == 0){ // Left Edge
            if(cj == 0){ // bottom left corner
              JA[0+jstart] = fdof; // Coarse point
              JA[1+jstart] = fdof+1;
              JA[2+jstart] = fdof+nf1d;
              for(i=0; i<3; i++){
                  val[i+jstart] = stencil[i+3];
              }
              jstart = jstart+3; // Update where we start filling val from
            }else if(cj == nc1d-1){ // top left corner
              JA[0+jstart] = fdof-nf1d;
              JA[1+jstart] = fdof-nf1d+1;
              JA[2+jstart] = fdof; // Coarse point
              JA[3+jstart] = fdof+1;
              for(i=0; i<4; i++){
                  val[i+jstart] = stencil[i+1];
              }
              jstart = jstart+4; // Update where we start filling val from
            }else{
              JA[0+jstart] = fdof-nf1d;
              JA[1+jstart] = fdof-nf1d+1;
              JA[2+jstart] = fdof; // Coarse point
              JA[3+jstart] = fdof+1;
              JA[4+jstart] = fdof+nf1d;
              for(i=0; i<5; i++){
                  val[i+jstart] = stencil[i+1];
              }
              jstart = jstart+5; // Update where we start filling val from
            }
        }else if( ci == nc1d-1){ // Right Edge
            if(cj == 0){ // bottom right corner
              JA[0+jstart] = fdof-1;
              JA[1+jstart] = fdof; // Coarse point
              JA[2+jstart] = fdof+nf1d-1;
              JA[3+jstart] = fdof+nf1d;
              for(i=0; i<4; i++){
                  val[i+jstart] = stencil[i+2];
              }
              jstart = jstart+4; // Update where we start filling val from
            }else if(cj == nc1d-1){ // top right corner
              JA[0+jstart] = fdof-nf1d;
              JA[1+jstart] = fdof-1;
              JA[2+jstart] = fdof; // Coarse point
              for(i=0; i<3; i++){
                  val[i+jstart] = stencil[i+1];
              }
              jstart = jstart+3; // Update where we start filling val from
            }else{
              JA[0+jstart] = fdof-nf1d;
              JA[1+jstart] = fdof-1;
              JA[2+jstart] = fdof; // Coarse point
              JA[3+jstart] = fdof+nf1d-1;
              JA[4+jstart] = fdof+nf1d;
              for(i=0; i<5; i++){
                  val[i+jstart] = stencil[i+1];
              }
              jstart = jstart+5; // Update where we start filling val from
            }
        }else if( cj == 0 ){
            JA[0+jstart] = fdof-1;
            JA[1+jstart] = fdof; // Coarse point
            JA[2+jstart] = fdof+1;
            JA[3+jstart] = fdof+nf1d-1;
            JA[4+jstart] = fdof+nf1d;
            for(i=0; i<5; i++){
                val[i+jstart] = stencil[i+2];
            }
            jstart = jstart+5; // Update where we start filling val from
        }else if( cj == nc1d-1 ){
            JA[0+jstart] = fdof-nf1d;
            JA[1+jstart] = fdof-nf1d+1;
            JA[2+jstart] = fdof-1;
            JA[3+jstart] = fdof; // Coarse point
            JA[4+jstart] = fdof+1;
            for(i=0; i<5; i++){
                val[i+jstart] = stencil[i];
            }
            jstart = jstart+5; // Update where we start filling val from
        }else{
            JA[0+jstart] = fdof-nf1d;
            JA[1+jstart] = fdof-nf1d+1;
            JA[2+jstart] = fdof-1;
            JA[3+jstart] = fdof; // Coarse point
            JA[4+jstart] = fdof+1;
            JA[5+jstart] = fdof+nf1d-1;
            JA[6+jstart] = fdof+nf1d;
            for(i=0; i<7; i++){
                val[i+jstart] = stencil[i];
            }
            jstart = jstart+7; // Update where we start filling val from
        }// End of edge detection

        // Store IA
        IA[cdof+1] = jstart;
      }//ci
    }//cj
    return;
}

/***********************************************************************************************/
/**
 * \fn void build_constant_R (dCSRmat *tentp,
 *                              INT      nf1d,
 *                              INT      nc1d)
 *
 * \brief Build tentative P for piecewise constant elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param nf1d               Number of fine elements in 1d (number of vertices on an edge - 1)
 * \param nc1d               Number of coarse elements in 1d (number of vertices on an edge - 1)
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 11/29/2018
 *
 */
void build_constant_R (dCSRmat *R,
                             INT      nf1d,
                             INT      nc1d)
{
    int i;
    int ci, cj;
    int fdof, cdof;
    int jstart=0;
    // nf1d is number of elements we would have in 1d (number of vertices in 1d minus 1)

    // allocate memory for R
    R->row = nc1d*nc1d*2;
    R->col = nf1d*nf1d*2;
    R->nnz = nc1d*nc1d*2*4;
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;

    IA[0] = 0;
    // This loops over the coarse and fine elements. Don't question it.
    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d*2; ci=ci+2){
        // Lower Triangle and Upper Triangle that share Diagonal

        // Lower Triangle
        cdof = ci + cj*nc1d*2;
        fdof = ci*2 + cj*nf1d*4;
        //Fill
        JA[0+jstart] = fdof;
        JA[1+jstart] = fdof+1;
        JA[2+jstart] = fdof+2;
        JA[3+jstart] = fdof+nf1d*2;
        for(i=0; i<4; i++){
          val[i+jstart] = 1.0;
        }
        jstart = jstart+4;
        IA[cdof+1] = jstart;

        // Upper Triangle
        cdof = cdof+1;
        //Fill
        JA[0+jstart] = fdof+3;
        JA[1+jstart] = fdof+nf1d*2+1;
        JA[2+jstart] = fdof+nf1d*2+2;
        JA[3+jstart] = fdof+nf1d*2+3;
        for(i=0; i<4; i++){
          val[i+jstart] = 1.0;
        }
        jstart = jstart+4;
        IA[cdof+1] = jstart;

      }
    }
    return;
}

/***********************************************************************************************/
/**
 * \fn void build_face_R (dCSRmat *tentp,
 *                                  mesh_struct *fmesh,
 *                                  mesh_struct *cmesh,
 *                                  INT     nf1d,
 *                                  INT     nc1d)
 * \brief Build tentative P for RT0 elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param fmesh              Fine mesh_struct
 * \param cmesh              Coarse mesh_struct
 * \param nf1d               Number of fine elements in 1d (one less than the number of vertices in 1d)
 * \param nc1d               Number of coarse elements in 1d (one less than the number of vertices in 1d)
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
void build_face_R (dCSRmat *R,
                             mesh_struct  *fmesh,
                             mesh_struct  *cmesh,
                             INT      nf1d,
                             INT      nc1d)
{
    INT i,j,ii,jj,kk;
    INT ci, cj;
    INT cdof;
    INT felm, celm;
    INT fface, cface;
    INT felmList[8];
    INT jstart=0;
    INT index=0;
    INT locFaceId;
    // nf1d is number of elements we would have in 1d (number of vertices in 1d minus 1)

    INT rowa;
    INT rowb;
    INT rowaf;
    INT rowbf;

    // Basis stuff
    INT dim = fmesh->dim;
    REAL* x = (REAL *) calloc(dim,sizeof(REAL));
    REAL* phi = (REAL *) calloc(dim*cmesh->f_per_elm,sizeof(REAL));
    REAL* dphi = (REAL *) calloc(cmesh->f_per_elm,sizeof(REAL));
    REAL value;

    INT* v_on_elm = (INT*)calloc(cmesh->v_per_elm,sizeof(INT));
    INT* f_on_elm = (INT*)calloc(cmesh->f_per_elm,sizeof(INT));

    //Garbage
    INT* I;
    INT* J;
    REAL* V;
    INT* IJfilled;

    // allocate memory for R
    R->row = cmesh->nface;
    R->col = fmesh->nface;
    R->nnz = R->row*12; //TODO: FIX THIS
    //R->nnz = R->row*9; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    I = (INT *)calloc(R->nnz,sizeof(INT));
    J = (INT *)calloc(R->nnz,sizeof(INT));
    V = (REAL *)calloc(R->nnz,sizeof(REAL));
    IJfilled = (INT *)calloc(R->row*16,sizeof(INT));//TODO:BAD!!!!!
    for(i=0; i<R->row*16; i++) IJfilled[i] = -1; // Set all to -1 //TODO:BAD!!!!
    INT not_duplicate_entry = 0;

    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;
    // This loops over the coarse and fine elements. Don't question it.
    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d*2; ci=ci+2){
        // Lower Triangle and upper triangle box
        celm = ci + cj*nc1d*2;
        felm = ci*2 + cj*nf1d*4;
        // Make list of all fine elms to loop over
        felmList[0] = felm;         // Lower
        felmList[1] = felm+1;       // Lower
        felmList[2] = felm+2;       // Lower
        felmList[3] = felm+nf1d*2;    // Lower
        felmList[4] = felm+3;       // Upper
        felmList[5] = felm+nf1d*2+1;  // Upper
        felmList[6] = felm+nf1d*2+2;  // Upper
        felmList[7] = felm+nf1d*2+3;  // Upper
        for(ii=0;ii<8;ii++){
          felm = felmList[ii];
          if(ii==4) celm = celm+1; // Switch to upper triangle

          // Get Coarse DOF
          get_incidence_row(celm,cmesh->el_v,v_on_elm);
          get_incidence_row(celm,cmesh->el_f,f_on_elm);
          rowa = cmesh->el_f->IA[celm];
          rowb = cmesh->el_f->IA[celm+1];
          locFaceId = 0;
          for(j=rowa;j<rowb;j++){
            cface = cmesh->el_f->JA[j];
            // Get Fine DOF
            rowaf = fmesh->el_f->IA[felm];
            rowbf = fmesh->el_f->IA[felm+1];
            for(jj=rowaf;jj<rowbf;jj++){
              fface = fmesh->el_f->JA[jj];
              // Fill P
              value = 0.0;
              // Evaluate coarse mesh basis function on fine mesh dof.
              x[0] = fmesh->f_mid[fface*dim];
              x[1] = fmesh->f_mid[fface*dim+1];
              rt_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
              // Midpoint Rule: \int phi_{coarse} n_{fine}
              for(i=0;i<dim;i++) value += fmesh->f_norm[fface*dim+i]*phi[locFaceId*dim+i];
              // Scale midpoint rule based on face area and scale basis function by face area
              //value = value * fmesh->f_area[fface] / cmesh->f_area[cface];
              // TODO: It works without the scaling (why? I don't know)
              not_duplicate_entry = 1;
              for(kk=0;kk<16;kk++){
                if( IJfilled[cface*16 + kk] == fface ){
                  not_duplicate_entry = 0;
                  break;
                } else if(IJfilled[cface*16 + kk] == -1){
                  not_duplicate_entry = 1;
                  IJfilled[cface*16+kk] = fface;
                  break;
                }
              }
              if( kk == 16 ){ printf("\n\n\n--------------------------------------------------\nNODE NOT STORED!!!!\n");}
              //not_duplicate_entry = 1;
              //for(kk=0;kk<index;kk++){
              //  if(cface==I[kk] && fface==J[kk]) not_duplicate_entry = 0;
              //}
              if(not_duplicate_entry){
                if(ABS(value) > 1e-8){
                  I[index] = cface;
                  J[index] = fface;
                  V[index] = value;
                  index++;
                }
              }
            }
            locFaceId++;
          }//j
        }//ii
      }//ci
    }//cj

    // Fill the CSR matrix format
    IA[0] = jstart;
    for(cdof=0; cdof < R->row; cdof++){
      for(ii=0;ii<index;ii++){
        if(cdof == I[ii]){
          JA[jstart] = J[ii];
          val[jstart] = V[ii];
          jstart++;
        }
      }//ii
      IA[cdof+1] = jstart;
    }//i
    R->nnz = jstart;
    return;
}

/**
 * \fn void build_bubble_R (dCSRmat *R,
 *                          dCSRmat *Rblx,
 *                          dCSRmat *Rbly,
 *                          mesh_struct *fmesh,
 *                          mesh_struct *cmesh,
 *                          INT     nf1d,
 *                          INT     nc1d)
 * \brief Build tentative P for RT0 elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param fmesh              Fine mesh_struct
 * \param cmesh              Coarse mesh_struct
 * \param nf1d               Number of fine elements in 1d (one less than the number of vertices in 1d)
 * \param nc1d               Number of coarse elements in 1d (one less than the number of vertices in 1d)
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
void build_bubble_R (dCSRmat *R,
                     dCSRmat *Rblx,
                     dCSRmat *Rbly,
                     mesh_struct  *fmesh,
                     mesh_struct  *cmesh,
                     INT      nf1d,
                     INT      nc1d)
{
    // nf1d is number of elements we would have in 1d (number of vertices in 1d minus 1)
    INT i,j,ii,jj;
    INT kk;
    INT ed;
    INT ci, cj;
    INT felm, celm;
    INT fface, cface;
    INT felmList[8];
    INT jstart=0;
    INT index=0;
    INT locFaceId;

    INT rowa;
    INT rowb;
    INT rowaf;
    INT rowbf;

    // Quadrature  Stuff
    //INT nq1d = 4;
    //INT quad;
    //REAL qval;
    //qcoordinates *cqface = get_quadrature_boundary(fmesh,nq1d,2);//TODO:free

    // Basis stuff
    INT dim = fmesh->dim;
    REAL* x = (REAL *) calloc(dim,sizeof(REAL));
    REAL* phi = (REAL *) calloc(dim*cmesh->f_per_elm,sizeof(REAL));
    REAL* dphi = (REAL *) calloc(dim*dim*cmesh->f_per_elm,sizeof(REAL));
    REAL* Fphi = (REAL *) calloc(dim*cmesh->f_per_elm,sizeof(REAL));
    REAL* Fdphi = (REAL *) calloc(dim*dim*cmesh->f_per_elm,sizeof(REAL));
    REAL value;
    REAL lval;
    //REAL alpha;

    INT* v_on_elm = (INT*)calloc(cmesh->v_per_elm,sizeof(INT));
    INT* f_on_elm = (INT*)calloc(cmesh->f_per_elm,sizeof(INT));
    //INT* v_on_elmF = (INT*)calloc(cmesh->v_per_elm,sizeof(INT));
    //INT* f_on_elmF = (INT*)calloc(cmesh->f_per_elm,sizeof(INT));
    INT* v_on_f   = (INT*)calloc(dim,sizeof(INT));

    //Garbage
    INT* I;
    INT* J;
    REAL* V;
    INT* IJfilled;

    // allocate memory for R
    R->row = cmesh->nface;
    R->col = fmesh->nface;
    R->nnz = R->row*12; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));
    // Rename for convenience
    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;
    // allocate helper arrays.
    INT not_duplicate_entry = 0;
    I = (INT *)calloc(R->nnz,sizeof(INT));
    J = (INT *)calloc(R->nnz,sizeof(INT));
    V = (REAL *)calloc(R->nnz,sizeof(REAL));
    IJfilled = (INT *)calloc(R->row*16,sizeof(INT));//TODO:BAD!!!!!
    for(i=0; i<R->row*16; i++) IJfilled[i] = -1;

    // allocate memory for Rl1
    // allocate memory for Rl2
    INT vertex;
    INT* VertFilled = (INT *)calloc(R->row*9,sizeof(INT));//TODO:BAD!!!!!
    for(i=0; i<R->row*9; i++) VertFilled[i] = -1;
    INT indexl=0;
    //TODO: Fix allocation amount
    INT*  Il1 = (INT*) calloc(12*fmesh->nv,sizeof(INT));
    INT*  Jl1 = (INT*) calloc(12*fmesh->nv,sizeof(INT));
    REAL* Vl1 = (REAL*)calloc(12*fmesh->nv,sizeof(REAL));
    INT*  Il2 = (INT*) calloc(12*fmesh->nv,sizeof(INT));
    INT*  Jl2 = (INT*) calloc(12*fmesh->nv,sizeof(INT));
    REAL* Vl2 = (REAL*)calloc(12*fmesh->nv,sizeof(REAL));


    // This loops over the coarse and fine elements. Don't question it.
    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d*2; ci=ci+2){
        // Lower Triangle and upper triangle box
        celm = ci + cj*nc1d*2;
        felm = ci*2 + cj*nf1d*4;
        // Make list of all fine elms to loop over
        felmList[0] = felm;         // Lower
        felmList[1] = felm+1;       // Lower
        felmList[2] = felm+2;       // Lower
        felmList[3] = felm+nf1d*2;    // Lower
        felmList[4] = felm+3;       // Upper
        felmList[5] = felm+nf1d*2+1;  // Upper
        felmList[6] = felm+nf1d*2+2;  // Upper
        felmList[7] = felm+nf1d*2+3;  // Upper
        for(ii=0;ii<8;ii++){
          felm = felmList[ii];
          if(ii==4) celm = celm+1; // Switch to upper triangle
          // Get Coarse DOF
          get_incidence_row(celm,cmesh->el_v,v_on_elm);
          get_incidence_row(celm,cmesh->el_f,f_on_elm);
          rowa = cmesh->el_f->IA[celm];
          rowb = cmesh->el_f->IA[celm+1];
          locFaceId = 0;
          for(j=rowa;j<rowb;j++){
            cface = cmesh->el_f->JA[j];
            // Get Fine DOF
            rowaf = fmesh->el_f->IA[felm];
            rowbf = fmesh->el_f->IA[felm+1];
            for(jj=rowaf;jj<rowbf;jj++){
              fface = fmesh->el_f->JA[jj];
              // Fill R
              value = 0.0;

              // The following is a simpson's rule integration, with the endpoints of the function = 0, so only the midpoint is needed.
              x[0] = fmesh->f_mid[fface*dim];
              x[1] = fmesh->f_mid[fface*dim+1];
              bubble_face_basis(Fphi,Fdphi,x,v_on_elm,f_on_elm,cmesh);

              if( ABS(fmesh->f_norm[fface*dim+0]) > 1e-6){
                  i = 0;
              } else {
                  i = 1;
              }

              get_incidence_row(fface,fmesh->f_v,v_on_f);
              lval = 0.0;
              for(ed=0;ed<dim;ed++){
                x[0] = fmesh->cv->x[v_on_f[ed]];
                x[1] = fmesh->cv->y[v_on_f[ed]];
                // Evaluate Coarse Bubble at endpoints
                bubble_face_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
                // Save in R linear
                vertex = v_on_f[ed];

                not_duplicate_entry = 1;
                for(kk=0;kk<9;kk++){
                  if( VertFilled[cface*9 + kk] == vertex ){
                    not_duplicate_entry = 0;
                    break;
                  } else if(VertFilled[cface*9 + kk] == -1){
                    not_duplicate_entry = 1;
                    VertFilled[cface*9+kk] = vertex;
                    break;
                  }
                }
                if( kk == 9 ){ printf("\n\n\n--------------------------------------------------\nVERT NOT STORED!!!!\n");}
                //not_duplicate_entry = 1;
                //for(kk=0;kk<indexl;kk++){
                //  if(cface==Il1[kk] && vertex==Jl1[kk]) not_duplicate_entry = 0;
                //}
                if(not_duplicate_entry){
                  Il1[indexl] = cface;
                  Jl1[indexl] = vertex;
                  Vl1[indexl] = phi[locFaceId*dim+0];
                  Il2[indexl] = cface;
                  Jl2[indexl] = vertex;
                  Vl2[indexl] = phi[locFaceId*dim+1];
                  indexl++;
                }
                // Linear part on fine edge.
                lval += phi[locFaceId*dim+i]/2;
                //for(i=0;i<dim;i++) lval += phi[locFaceId*dim+i]*fmesh->f_norm[fface*dim+i];
              }

              value = ( Fphi[locFaceId*dim+i] - lval ) / ABS(fmesh->f_norm[fface*dim+i]);

              not_duplicate_entry = 1;
              for(kk=0;kk<16;kk++){
                if( IJfilled[cface*16 + kk] == fface ){
                  not_duplicate_entry = 0;
                  break;
                } else if(IJfilled[cface*16 + kk] == -1){
                  not_duplicate_entry = 1;
                  IJfilled[cface*16+kk] = fface;
                  break;
                }
              }
              if( kk == 16 ){ printf("\n\n\n--------------------------------------------------\nEDGE NOT STORED!!!!\n");}
              //not_duplicate_entry = 1;
              //for(kk=0;kk<index;kk++){
              //  if(cface==I[kk] && fface==J[kk]) not_duplicate_entry = 0;
              //}
              if(not_duplicate_entry){
                if(ABS(value) > 1e-8){
                  I[index] = cface;
                  J[index] = fface;
                  V[index] = value;
                  index++;
                }
              }

            }//jj (fface)
            locFaceId++;
          }//j (cface)
        }//ii
      }//ci
    }//cj

    //printf("index: %d\tnnz: %d\n",index,R->nnz);//TODO: Correct number of nnz in the matrix
    //printf("indexl: %d\tnnz: %d\n",indexl,R->row*12);//TODO: Correct number of nnz in the matrix

    // Fill the CSR matrix format
    INT cdof;
    IA[0] = jstart;
    for(cdof=0; cdof < R->row; cdof++){
      for(ii=0;ii<index;ii++){
        if(cdof == I[ii]){
          JA[jstart] = J[ii];
          val[jstart] = V[ii];
          jstart++;
        }
      }//ii
      IA[cdof+1] = jstart;
    }//cdof
    R->nnz = jstart;

    // allocate memory for R
    Rblx->row = cmesh->nface;
    Rblx->col = fmesh->nv;
    Rblx->nnz = R->row*12; //TODO: FIX THIS
    Rblx->IA  = (INT *)calloc(Rblx->row+1,sizeof(INT));
    Rblx->JA  = (INT *)calloc(Rblx->nnz, sizeof(INT));
    Rblx->val = (REAL *)calloc(Rblx->nnz, sizeof(REAL));
    // allocate memory for R
    Rbly->row = cmesh->nface;
    Rbly->col = fmesh->nv;
    Rbly->nnz = R->row*12; //TODO: FIX THIS
    Rbly->IA  = (INT *)calloc(Rblx->row+1,sizeof(INT));
    Rbly->JA  = (INT *)calloc(Rbly->nnz, sizeof(INT));
    Rbly->val = (REAL *)calloc(Rbly->nnz, sizeof(REAL));
    // Fill CSR matrix for Bubble Linear
    jstart = 0;
    Rblx->IA[0] = jstart;
    Rbly->IA[0] = jstart;
    for(cdof=0; cdof< Rblx->row; cdof++){
      for(ii=0;ii<indexl;ii++){
        if(cdof == Il1[ii]){
          Rblx->JA[jstart]  = Jl1[ii];
          Rblx->val[jstart] = Vl1[ii];
          Rbly->JA[jstart]  = Jl2[ii];
          Rbly->val[jstart] = Vl2[ii];
          jstart++;
        }
      }//ii
      Rblx->IA[cdof+1] = jstart;
      Rbly->IA[cdof+1] = jstart;
    }
    Rblx->nnz = jstart;
    Rbly->nnz = jstart;
//printf("___________________________________________\n");
//printf("%d  %d  %d\n",R->row,Rbly->row,Rblx->row);
//printf("%d  %d  %d\n",R->col,Rbly->col,Rblx->col);
//csr_print_matlab(stdout,R);
//printf("___________________________________________\n");

    //TODO: FREE!
}


/***********************************************************************************************/
/**
 * \fn SHORT gmg_setup_RT0 (dCSRmat *tentp,
 *                                  INT     ndof
 *                                  INT     levelNum)
 * \brief Build
 *
 * \param tentp              Pointer to the prolongation operators
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
INT gmg_setup_RT0(mesh_struct* fine_level_mesh)
{
//INT gmg_setup_RT0(GMG_data *mgl,
//                    GMG_param *param)
  // Things from param
  //SHORT max_levels = param->max_levels;
  SHORT max_levels = 2;
  //SHORT lvl = 0;

  //
  INT i;
  INT fSize, cSize;
  INT nf1d, nc1d;

  // To read from file
  FILE* cgfid;
  char cgridfile[500];

  mesh_struct** meshHeirarchy = (mesh_struct**)calloc(max_levels, sizeof(mesh_struct*));
  //meshHeirarchy[0] = mgl->fine_level_mesh;
  meshHeirarchy[0] = fine_level_mesh;

  dCSRmat R;

  // Pre-loop stuff
  for(i=0;i<max_levels-1;i++){
    fSize = sqrt(meshHeirarchy[i]->nv)-1;// Is there a type problem here?
    cSize = (fSize/2);
    sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",cSize+1);//Need to customize this to specific directory
    //sprintf(cgridfile,"/home/pohm01/HAZMAT/hazmath/examples/grids/2D/unitSQ_n%d.haz",cSize+1);//Need to customize this to specific directory
    cgfid = HAZ_fopen(cgridfile,"r");
    meshHeirarchy[i+1] = (mesh_struct*)calloc(1, sizeof(mesh_struct));
    initialize_mesh(meshHeirarchy[i+1]);
    creategrid_fread(cgfid,0,meshHeirarchy[i+1]);
    fclose(cgfid);

    // Let's build here
    nf1d = fSize;
    nc1d = cSize;
    printf("Let's build the Restriction\n");
    printf("Sizes are nf1d = %d, and nc1d = %d\n",nf1d, nc1d);

    //build_constant_R (&R,nf1d,nc1d);
    build_linear_R (&R,nf1d+1,nc1d+1);
    //build_face_R(&R, meshHeirarchy[i], meshHeirarchy[i+1], nf1d, nc1d);

    // Setup for next loop
  }

  // Free mesh
  for(i=1;i<max_levels;i++){
    free_mesh(meshHeirarchy[i]);
  }
  free(meshHeirarchy);

  return 0;

}



/***********************************************************************************************/
/**
 * \fn SHORT gmg_load_coarse_grids_from_file (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Load coarse grids from file
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Peter Ohm
 * \date   07/22/2019
 *
 */
SHORT gmg_load_coarse_grids_from_file( MG_blk_data *mgl,
                                       AMG_param *param)
{
  // local variables
  const SHORT prtlvl     = param->print_level;
  const SHORT cycle_type = param->cycle_type;
  const SHORT csolver    = param->coarse_solver;
  const SHORT max_levels = param->max_levels;

  INT lvl = 0;
  INT status = SUCCESS;
  INT dim = mgl[0].fine_level_mesh->dim;
  INT i,j,nf1d,nc1d,csize;

  // To Read from file
  FILE* cgfid;
  char cgridfile[500];

  for( lvl = 0; lvl < max_levels-1; lvl++){
    /*-- Build coarse level mesh --*/
    csize = (sqrt(mgl[lvl].fine_level_mesh->nv)-1)/2 + 1;
    mgl[lvl+1].fine_level_mesh = (mesh_struct*)calloc(1, sizeof(mesh_struct));
    initialize_mesh(mgl[lvl+1].fine_level_mesh);
    /*-- Need to customize this to specific directory --*/
    sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",csize);
    //sprintf(cgridfile,"/home/xiaozhehu/Work/Projects/HAZMATH/hazmath/examples/grids/2D/unitSQ_n%d.haz",csize);
    //sprintf(cgridfile,"/home/pohm01/HAZMAT/hazmath/examples/grids/2D/unitSQ_n%d.haz",csize);
    cgfid = HAZ_fopen(cgridfile,"r");
    creategrid_fread(cgfid,0,mgl[lvl+1].fine_level_mesh);
    fclose(cgfid);
    printf("Coarse grid loaded for lvl=%d...\n",lvl);
    /*-- Set boundary flags on the coarse mesh --*/
    mgl[0].set_bdry_flags(mgl[lvl+1].fine_level_mesh);
  }

  return status;
}

/***********************************************************************************************/
/**
 * \fn SHORT gmg_build_coarse_FE_spaces (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Build FE spaces on coarse grids (mostly for use with BC flags)
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Peter Ohm
 * \date   07/22/2019
 *
 */
SHORT gmg_build_coarse_FE_spaces( MG_blk_data *mgl,
                                  AMG_param *param)
{
  // local variables
  const SHORT prtlvl     = param->print_level;
  const SHORT cycle_type = param->cycle_type;
  const SHORT csolver    = param->coarse_solver;
  const SHORT max_levels = param->max_levels;

  INT lvl = 0;
  INT status = SUCCESS;
  INT dim = mgl[0].fine_level_mesh->dim;
  INT i,j,nf1d,nc1d,csize;

  INT nspaces = mgl[0].FE->nspaces;

  // Allocate and create FE spaces if needed.
  for( lvl = 1; lvl < max_levels; lvl++){
    mgl[lvl].FE = (block_fespace *) calloc( 1, sizeof(block_fespace));
    mgl[lvl].FE->ndof = 0;
    mgl[lvl].FE->nspaces = nspaces;
    mgl[lvl].FE->var_spaces = (fespace **) calloc( nspaces, sizeof(fespace *));
    for( i = 0; i < nspaces; i++){
      mgl[lvl].FE->var_spaces[i] = (fespace *) calloc(1, sizeof(fespace));
      create_fespace(mgl[lvl].FE->var_spaces[i], mgl[lvl].fine_level_mesh, mgl[0].FE->var_spaces[i]->FEtype);
      mgl[lvl].FE->ndof += mgl[lvl].FE->var_spaces[i]->ndof;
    }
  }

  return status;
}

/***********************************************************************************************/
/**
 * \fn SHORT gmg_apply_periodic_BC (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief 
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Peter Ohm
 * \date   07/22/2019
 *
 */
SHORT gmg_apply_periodic_BC( MG_blk_data *mgl,
                             AMG_param *param,
                             INT NoBBL)
{
  // local variables
  const SHORT prtlvl     = param->print_level;
  const SHORT cycle_type = param->cycle_type;
  const SHORT csolver    = param->coarse_solver;
  const SHORT max_levels = param->max_levels;

  INT lvl = 0;
  INT status = SUCCESS;
  INT dim = mgl[0].fine_level_mesh->dim;
  INT i,j,nf1d,nc1d,csize;
  INT ndof, cnt, dof_shift;
  INT brow = mgl[0].A.brow;

  dCSRmat tempRA;

  block_fespace FE_blk;
  FE_blk.nspaces = mgl[0].A.brow;
  FE_blk.var_spaces = (fespace **) calloc( FE_blk.nspaces, sizeof(fespace*));

  for( lvl = 0; lvl < max_levels; lvl++){
    set_periodic_bdry(mgl[lvl].FE->var_spaces[0], mgl[lvl].fine_level_mesh,0.0,1.0,0.0,1.0,0.0,1.0);
    set_periodic_bdry(mgl[lvl].FE->var_spaces[1], mgl[lvl].fine_level_mesh,0.0,1.0,0.0,1.0,0.0,1.0);
    set_periodic_bdry(mgl[lvl].FE->var_spaces[2], mgl[lvl].fine_level_mesh,0.0,1.0,0.0,1.0,0.0,1.0);
    set_periodic_bdry(mgl[lvl].FE->var_spaces[3], mgl[lvl].fine_level_mesh,0.0,1.0,0.0,1.0,0.0,1.0);
    set_periodic_bdry(mgl[lvl].FE->var_spaces[4], mgl[lvl].fine_level_mesh,0.0,1.0,0.0,1.0,0.0,1.0);

    // Create fake blockFE space that matches the 3x3 block matrix (put all displacements together)

    // Displacement Block
    FE_blk.var_spaces[0] = (fespace *) calloc(1, sizeof(fespace));
    ndof = 0;
    for(i=NoBBL; i<dim+1; i++){ ndof += mgl[lvl].FE->var_spaces[i]->ndof; }
    FE_blk.var_spaces[0]->periodic = (INT *) calloc(ndof, sizeof(INT));
    FE_blk.var_spaces[0]->ndof = ndof;
    cnt = 0;
    dof_shift = 0;
    for(i=NoBBL; i<dim+1; i++){
      for(j=0; j<mgl[lvl].FE->var_spaces[i]->ndof; j++){
        if( mgl[lvl].FE->var_spaces[i]->periodic[j] > -1 ){
          FE_blk.var_spaces[0]->periodic[cnt] = mgl[lvl].FE->var_spaces[i]->periodic[j]+dof_shift;
        } else {
          FE_blk.var_spaces[0]->periodic[cnt] = mgl[lvl].FE->var_spaces[i]->periodic[j];
        }
        cnt++;
      }
      dof_shift =cnt;
    }
    // Darcy Block
    FE_blk.var_spaces[1] = mgl[lvl].FE->var_spaces[3];
    // Pressure Block
    FE_blk.var_spaces[2] = mgl[lvl].FE->var_spaces[4];

    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].A_periodic));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].P_periodic));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R_periodic));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R_periodic_scaled));

    printf("\tGenerating periodic P\n");
    generate_periodic_P_blockFE( &FE_blk, &(mgl[lvl].P_periodic) );
    printf("\tGenerating periodic R scaled\n");
    generate_periodic_R_scaled_blockFE( &(mgl[lvl].P_periodic), &(mgl[lvl].R_periodic_scaled) );
    printf("\tGenerating periodic R\n");
    for(i=0; i<brow; i++){
      dcsr_trans( mgl[lvl].P_periodic.blocks[i+i*brow], mgl[lvl].R_periodic.blocks[i+i*brow] );
    }

    // Form triple matrix product on level for A_periodic
    printf("\tEliminating PeriodicBC...\n");
    eliminate_PeriodicBC_blockFE( &(mgl[lvl].P_periodic), &(mgl[lvl].A), NULL);
    printf("\tEliminated PeriodicBC...\n");
    // Form triple matrix product on level for P_periodic (for lvl-1)
    printf("\tForming Prolongation and Restriction matrices for periodic problem\n");
    if(lvl>0){
      for(i=0; i< brow; i++){
        // Take advantage of the fact that R and P are block diagonal
        // Create R with periodic BC
        // TODO
        //dcsr_mxm(mgl[lvl].R_periodic_scaled.blocks[i+i*brow],mgl[lvl-1].R.blocks[i+i*brow],&tempRA);
        //dcsr_mxm(&tempRA,mgl[lvl-1].P_periodic.blocks[i+i*brow],mgl[lvl-1].R.blocks[i+i*brow]);
        //dcsr_free(&tempRA);
        // Create P with periodic BC
        dcsr_mxm(mgl[lvl-1].R_periodic_scaled.blocks[i+i*brow],mgl[lvl-1].P.blocks[i+i*brow],&tempRA);
        dcsr_mxm(&tempRA,mgl[lvl].P_periodic.blocks[i+i*brow],mgl[lvl-1].P.blocks[i+i*brow]);
        dcsr_free(&tempRA);
        /////////////////////
        dcsr_trans( mgl[lvl-1].P.blocks[i+i*brow], mgl[lvl-1].R.blocks[i+i*brow] );
      }
    }

    dof_shift = 0;
    cnt =  0;
    cnt = NoBBL*mgl[lvl].FE->var_spaces[0]->ndof;
    for(i=0; i<brow; i++){
      for(j=0; j<mgl[lvl].P_periodic.blocks[i+i*brow]->nnz; j++){
        mgl[lvl].FE->dirichlet[ mgl[lvl].P_periodic.blocks[i+i*brow]->JA[j] + dof_shift ] = mgl[lvl].FE->dirichlet[cnt];
        cnt++;
      }
      dof_shift += mgl[lvl].P_periodic.blocks[i+i*brow]->col;
    }

  }

  return status;
}

/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/**
 * \fn static SHORT gmg_blk_setup (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of gmg assuming that each FE space relates one-to-one with a block in A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Peter Ohm
 * \date   07/22/2019
 *
 */
SHORT gmg_blk_setup_generic(MG_blk_data *mgl,
                            AMG_param *param)
{
  // local variables
  const SHORT prtlvl     = param->print_level;
  const SHORT cycle_type = param->cycle_type;
  const SHORT csolver    = param->coarse_solver;
  SHORT       max_levels = param->max_levels;

  INT   lvl = 0;
  INT   status = SUCCESS;
  INT   dim = mgl[0].fine_level_mesh->dim;
  INT   nspaces = mgl[0].FE->nspaces;
  INT   i,j,nf1d,nc1d,csize;
  REAL  setup_start, setup_end;

  dCSRmat tempRA;
  dCSRmat Rmerge;
  block_dCSRmat tempRblk;// To merge bbl and P1 blocks
  
  get_time(&setup_start);// Timing

  status = gmg_load_coarse_grids_from_file( mgl, param );
  status = gmg_build_coarse_FE_spaces( mgl, param );

  while ( lvl < max_levels-1 ) {
    /*-- Allocate for R A P --*/
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].P));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A_noBC));
    printf("Allocation of R A P for lvl=%d is finished...\n",lvl);

    /*-- Form restriction matrix R for each FE space--*/
    for( i = 0; i<nspaces; i++ ){
      // Check gmg type, this could be replaced with checking FE type
      switch( mgl[lvl].gmg_type[i] ){
        case 0:// P0
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nelm/2);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nelm/2);
          build_constant_R( mgl[lvl].R.blocks[i+i*nspaces], nf1d, nc1d);
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[i], mgl[lvl+1].fine_level_mesh, 1,1);
          break;
        case 1:// P1
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
          build_linear_R( mgl[lvl].R.blocks[i+i*nspaces], nf1d, nc1d);
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[i], mgl[lvl+1].fine_level_mesh, 1,1);
          break;
        case 30:// RT0
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
          nc1d = (nf1d)/2;
          build_face_R( mgl[lvl].R.blocks[i+i*nspaces], mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[i], mgl[lvl+1].fine_level_mesh, 1,1);
          break;
        case 999:// P1 + bubbles combined as a single block matrix
          bdcsr_alloc(dim+1,dim+1,&tempRblk);
          tempRblk.blocks[3] = NULL; tempRblk.blocks[5] = NULL;
          tempRblk.blocks[6] = NULL; tempRblk.blocks[7] = NULL;
          //P1 both x and Y components
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
          for(j=1; j<dim+1; j++) build_linear_R( tempRblk.blocks[j+j*tempRblk.brow], nf1d, nc1d);
          //Bubble
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
          nc1d = (nf1d)/2;
          build_bubble_R( tempRblk.blocks[0], tempRblk.blocks[1], tempRblk.blocks[2],
                          mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
          Rmerge = bdcsr_2_dcsr(&tempRblk);
          dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
          dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[0], mgl[lvl+1].fine_level_mesh, 1, 10); //bbl
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, 1, 10); // Ux
          set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, 1, 10); // Uy
          dcsr_free(&Rmerge);
          break;
      }// switch gmg_type
    }// for i<nspaces
    printf("Built R for lvl=%d...\n",lvl);

    /*-- Form Prolongation --*/
    for(i=0; i<nspaces; i++){
      dcsr_trans(mgl[lvl].R.blocks[i+i*nspaces], mgl[lvl].P.blocks[i+i*nspaces]);
    }
    printf("Built P for lvl=%d...\n",lvl);

    /*-- Form coarse level stiffness matrix --*/
    for(i=0; i<nspaces; i++){
      for(j=0; j<nspaces; j++){
        if(mgl[lvl].A.blocks[j+i*nspaces]){
          dcsr_mxm(mgl[lvl].R.blocks[i+i*nspaces], mgl[lvl].A_noBC.blocks[j+i*nspaces], &tempRA);
          dcsr_mxm(&tempRA, mgl[lvl].P.blocks[j+j*nspaces], mgl[lvl+1].A_noBC.blocks[j+i*nspaces]);
          dcsr_free(&tempRA);
        } else {
          mgl[lvl+1].A_noBC.blocks[j+i*nspaces] = NULL;
        }
      }//j
    }//i
    printf("Built RAP for lvl=%d...\n",lvl);

    /*-- Eliminate dirichlet boundaries from stiffness matrix --*/
    // Need to write this: Need to fill dirichlet BC array for each level.
    bdcsr_cp( &mgl[lvl+1].A_noBC, &mgl[lvl+1].A);
    set_dirichlet_bdry_block(mgl[lvl+1].FE, mgl[lvl+1].fine_level_mesh);
    eliminate_DirichletBC_blockFE_blockA(NULL, mgl[lvl+1].FE, mgl[lvl+1].fine_level_mesh,NULL,&mgl[lvl+1].A,0.0);

    lvl++;
  }

  // Setup coarse level systems for direct solvers
  switch (csolver) {

#if WITH_SUITESPARSE
      case SOLVER_UMFPACK: {
          printf("Setting up coarse solve: Using UMFPACK...\n");
          // Need to sort the matrix A for UMFPACK to work
          // merge blocks
          mgl[lvl].Ac = bdcsr_2_dcsr(&mgl[lvl].A);
          dCSRmat Ac_tran;
          dcsr_trans(&mgl[lvl].Ac, &Ac_tran);
          printf("Ac stats: row=%d col=%d nnz=%d\n",mgl[lvl].Ac.row,mgl[lvl].Ac.col,mgl[lvl].Ac.nnz);
          // It is equivalent to do transpose and then sort
          //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
          //     fasp_dcsr_sort(&Ac_tran);
          dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
          dcsr_free(&Ac_tran);
          mgl[lvl].Numeric = umfpack_factorize(&mgl[lvl].Ac, 0);
        break;
      }
#endif
      default:
          printf("We are not using SUITESPARSE!\n");
          // Do nothing!
          break;
  }

  // setup total level number and current level
  INT m=0;
  for(i=0; i<nspaces; i++) m += mgl[0].A.blocks[i+i*nspaces]->row;
  mgl[0].num_levels = max_levels = lvl+1;
  mgl[0].w          = dvec_create(m);

  for ( lvl = 1; lvl < max_levels; ++lvl) {
      INT mm = 0;
      for(i=0;i<nspaces;i++) mm += mgl[lvl].A.blocks[i+i*nspaces]->row;
      mgl[lvl].num_levels = max_levels;
      mgl[lvl].b          = dvec_create(mm);
      mgl[lvl].x          = dvec_create(mm);

      mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

      if ( cycle_type == NL_AMLI_CYCLE )
          mgl[lvl].w = dvec_create(3*mm);
      else
          mgl[lvl].w = dvec_create(2*mm);
  }

  if ( prtlvl > PRINT_NONE ) {
      get_time(&setup_end); // Timing
      print_cputime("geometric multigrid setup", setup_end - setup_start);
  }
  return status;
}


/**
 * \fn static SHORT gmg_blk_setup (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of gmg assuming that each FE space relates one-to-one with a block in A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Peter Ohm
 * \date   07/22/2019
 *
 */
SHORT gmg_blk_setup_biot_bubble(MG_blk_data *mgl,
                                AMG_param *param)
{
  // local variables
  const SHORT prtlvl     = param->print_level;
  const SHORT cycle_type = param->cycle_type;
  const SHORT csolver    = param->coarse_solver;
  SHORT       max_levels = param->max_levels;

  INT lvl     = 0;
  INT status  = SUCCESS;
  INT dim     = mgl[0].fine_level_mesh->dim;
  INT nspaces = mgl[0].FE->nspaces;
  INT brow    = mgl[0].A.brow;
  INT i,j,nf1d,nc1d,csize,m;
  REAL  setup_start, setup_end;

  mgl[0].dirichlet_blk = (INT**)calloc(brow,sizeof(INT*));
  INT bm = 0;
  INT NoBBL = 0;
  if( mgl[0].gmg_type[0] == 111 ) NoBBL = 1; // If bubble is eliminated skip that FE space

  for(i=NoBBL;i<mgl[0].FE->nspaces;i++){
    if(i==NoBBL) mgl[0].dirichlet_blk[0] = &(mgl[0].FE->dirichlet[bm]); // Displacement
    if(i==3)     mgl[0].dirichlet_blk[1] = &(mgl[0].FE->dirichlet[bm]); // Darcy (RT0)
    if(i==4)     mgl[0].dirichlet_blk[2] = &(mgl[0].FE->dirichlet[bm]); // Pressure
    bm += mgl[0].FE->var_spaces[i]->ndof;
  }

  block_dCSRmat tempRblk;// To merge bbl and P1 blocks
  dCSRmat Rmerge;// merged bbl and P1 blocks
  dCSRmat tempRA;
  
  get_time(&setup_start);// Timing

  status = gmg_load_coarse_grids_from_file( mgl, param );
  printf("Coarse grids loaded from file\n");
  status = gmg_build_coarse_FE_spaces( mgl, param );
  printf("Coarse finite element spaces built.\n");

  while ( lvl < max_levels-1 ) {
    /*-- Allocate for R A P --*/
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].P));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A));
    bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A_noBC));

    mgl[lvl+1].dirichlet_blk = (INT**)calloc(brow,sizeof(INT*));//remove eventually

    /*-- Form restriction matrix R for each FE space--*/
    for( i = 0; i<brow; i++ ){
      // Check gmg type, this could be replaced with checking FE type
      switch( mgl[lvl].gmg_type[i] ){
        case 0:// P0
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nelm/2);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nelm/2);
          build_constant_R( mgl[lvl].R.blocks[i+i*brow], nf1d, nc1d);

          if(mgl[0].periodic_BC){
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[4], mgl[lvl+1].fine_level_mesh, -1,-1); // P
//            mgl[lvl+1].FE->var_spaces[4]->dirichlet[0] = 1;
          } else {
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[4], mgl[lvl+1].fine_level_mesh, -1,-1); // P
          }

          break;
        case 111:// P1 combined as a vector space to form single block matrix
          bdcsr_alloc(dim,dim,&tempRblk);
          tempRblk.blocks[1] = NULL;
          tempRblk.blocks[2] = NULL;
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
          for(j=0; j<dim; j++) build_linear_R( tempRblk.blocks[j+j*tempRblk.brow], nf1d, nc1d);
          Rmerge = bdcsr_2_dcsr(&tempRblk);
          dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
          dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
          dcsr_free(&Rmerge);

          if( mgl[0].periodic_BC ){
            nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, -1,-1); // Ux
            mgl[lvl+1].FE->var_spaces[1]->dirichlet[nc1d+1] = 1;
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, -1,-1); // Uy
            mgl[lvl+1].FE->var_spaces[2]->dirichlet[nc1d+1] = 1;
          } else {
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, 5,5); // Ux
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, 1,4); // Uy
          }

          break;
        case 999:// P1 + bubbles combined as a single block matrix
          bdcsr_alloc(dim+1,dim+1,&tempRblk);
          tempRblk.blocks[3] = NULL;
          tempRblk.blocks[5] = NULL;
          tempRblk.blocks[6] = NULL;
          tempRblk.blocks[7] = NULL;
          //P1 both x and Y components
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
          nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
          for(j=1; j<dim+1; j++) build_linear_R( tempRblk.blocks[j+j*tempRblk.brow], nf1d, nc1d);
          //Bubble
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
          nc1d = (nf1d)/2;
          build_bubble_R( tempRblk.blocks[0], tempRblk.blocks[1], tempRblk.blocks[2],
                          mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
          Rmerge = bdcsr_2_dcsr(&tempRblk);
          dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
          dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
          dcsr_free(&Rmerge);

          printf("\n\n****************************************************************************************************\n");
          printf("Dirichlet boundary conditions for GMG have been butchered to make periodic BC work in GMG\n");
          printf("****************************************************************************************************\n\n");
          if( mgl[0].periodic_BC ){
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[0], mgl[lvl+1].fine_level_mesh, -1,-1); //bbl
//            mgl[lvl+1].FE->var_spaces[0]->dirichlet[3] = 1;

            nc1d = sqrt(mgl[lvl+1].fine_level_mesh->nv);
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, -1,-1); // Ux
//            mgl[lvl+1].FE->var_spaces[1]->dirichlet[nc1d+1] = 1;
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, -1,-1); // Uy
//            mgl[lvl+1].FE->var_spaces[2]->dirichlet[nc1d+1] = 1;
          } else {
            //set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[0], mgl[lvl+1].fine_level_mesh, 1, 5); //bbl
            //set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, 4, 5); // Ux
            //set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, 1, 4); // Uy
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[0], mgl[lvl+1].fine_level_mesh, 1, 10); //bbl
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[1], mgl[lvl+1].fine_level_mesh, 1, 10); // Ux
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[2], mgl[lvl+1].fine_level_mesh, 1, 10); // Uy
          }

          break;
        case 30:// RT0
          nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
          nc1d = (nf1d)/2;
          build_face_R( mgl[lvl].R.blocks[i+i*brow], mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);

          if( mgl[0].periodic_BC ){
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[3], mgl[lvl+1].fine_level_mesh, -1,-1); // w
//            mgl[lvl+1].FE->var_spaces[3]->dirichlet[3] = 1;
          } else {
            //set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[3], mgl[lvl+1].fine_level_mesh, 1, 5); // w
            set_dirichlet_bdry(mgl[lvl+1].FE->var_spaces[3], mgl[lvl+1].fine_level_mesh, 1, 10); // w
          }

          break;
        default:
          printf("### ERROR: Unknown geometric multigrid type: %d!\n",mgl[lvl].gmg_type[i]);
          break;
      }// switch gmg_type
    }// for i<brow
    printf("\tBuilt R for lvl=%d...\n",lvl);

    /*-- Form Prolongation --*/
    for(i=0; i<brow; i++){
      dcsr_trans(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].P.blocks[i+i*brow]);
    }
    printf("\tBuilt P for lvl=%d...\n",lvl);

    /*-- Form coarse level stiffness matrix --*/
    for(i=0; i<brow; i++){
      for(j=0; j<brow; j++){
        if(mgl[lvl].A.blocks[j+i*brow]){
          dcsr_mxm(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].A_noBC.blocks[j+i*brow], &tempRA);
          dcsr_mxm(&tempRA, mgl[lvl].P.blocks[j+j*brow], mgl[lvl+1].A_noBC.blocks[j+i*brow]);
          dcsr_free(&tempRA);
        } else {
          mgl[lvl+1].A_noBC.blocks[j+i*brow] = NULL;
        }
      }//j
    }//i
    printf("\tBuilt RAP for lvl=%d...\n",lvl);

    /*-- Eliminate dirichlet boundaries from stiffness matrix --*/
    bdcsr_cp( &mgl[lvl+1].A_noBC, &mgl[lvl+1].A);
    set_dirichlet_bdry_block(mgl[lvl+1].FE, mgl[lvl+1].fine_level_mesh);//TODO: this will not work in elim case. bbl in FE struct
    mgl[lvl+1].FE->dirichlet = mgl[lvl+1].FE->dirichlet + NoBBL*mgl[lvl+1].FE->var_spaces[0]->ndof;
    eliminate_DirichletBC_blockFE_blockA(NULL, mgl[lvl+1].FE, mgl[lvl+1].fine_level_mesh,NULL,&mgl[lvl+1].A,0.0);
    mgl[lvl+1].FE->dirichlet = mgl[lvl+1].FE->dirichlet - NoBBL*mgl[lvl+1].FE->var_spaces[0]->ndof;

    bm = 0;
    mgl[lvl+1].dirichlet_blk = (INT**)calloc(brow,sizeof(INT*));
    for(i=NoBBL; i<nspaces; i++){
      if(i==NoBBL) mgl[lvl+1].dirichlet_blk[0] = &mgl[lvl+1].FE->dirichlet[bm]; // Displacement
      if(i==3)     mgl[lvl+1].dirichlet_blk[1] = &mgl[lvl+1].FE->dirichlet[bm]; // Darcy (RT0)
      if(i==4)     mgl[lvl+1].dirichlet_blk[2] = &mgl[lvl+1].FE->dirichlet[bm]; // Pressure
      bm += mgl[lvl+1].FE->var_spaces[i]->ndof;
    }

    if( !mgl[0].periodic_BC ){
      mgl[lvl+1].FE->dirichlet = mgl[lvl+1].FE->dirichlet + NoBBL*mgl[lvl+1].FE->var_spaces[0]->ndof;
    }

    lvl++;
  }

  /*-- Apply periodic boundary condition operators --*/
  if( mgl[0].periodic_BC ){
    printf("Applying periodic BC\n");
    status = gmg_apply_periodic_BC( mgl, param, NoBBL);
    printf("Finished applying periodic BC\n");
  }


  /*-- Setup coarse level systems for direct solvers --*/
  switch (csolver) {

#if WITH_SUITESPARSE
      case SOLVER_UMFPACK: {
          printf("Setting up coarse solve: Using UMFPACK...\n");
          // Need to sort the matrix A for UMFPACK to work
          // merge blocks
          // // Add a small I to each diagonal block
          mgl[lvl].Ac = bdcsr_2_dcsr(&mgl[lvl].A);
          dCSRmat Ac_tran;
          if(0){
            dCSRmat I1 = dcsr_create_identity_matrix(mgl[lvl].Ac.row,0);
            dCSRmat ApE;
            dcsr_add( &mgl[lvl].Ac, 1.0, &I1, 1e-10, &ApE);
            dcsr_trans(&ApE, &Ac_tran);
          } else {
            dcsr_trans(&mgl[lvl].Ac, &Ac_tran);
          }
          // It is equivalent to do transpose and then sort
          //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
          //     fasp_dcsr_sort(&Ac_tran);
          dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
          dcsr_free(&Ac_tran);
          mgl[lvl].Numeric = umfpack_factorize(&mgl[lvl].Ac, 0);
        break;
      }
#endif
      default:
          printf("We are not using SUITESPARSE!\n");
          // Do nothing!
          break;
  }

  // setup total level number and current level
  m=0;
  for(i=0; i<brow; i++) m += mgl[0].A.blocks[i+i*brow]->row;
  mgl[0].num_levels = max_levels = lvl+1;
  mgl[0].w          = dvec_create(m);

  for ( lvl = 1; lvl < max_levels; ++lvl) {
      INT mm = 0;
      for(i=0;i<brow;i++) mm += mgl[lvl].A.blocks[i+i*brow]->row;
      mgl[lvl].num_levels = max_levels;
      mgl[lvl].b          = dvec_create(mm);
      mgl[lvl].x          = dvec_create(mm);

      mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!

      if ( cycle_type == NL_AMLI_CYCLE )
          mgl[lvl].w = dvec_create(3*mm);
      else
          mgl[lvl].w = dvec_create(2*mm);
  }

  if ( prtlvl > PRINT_NONE ) {
      get_time(&setup_end);
      print_cputime("geometric multigrid setup", setup_end - setup_start);
  }

  return status;
}

/**
 * \fn void smoother_block_setup( MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup of block smoothers
 *
 * \param mgl       Pointer to MG_blk_data
 * \param param     Pointer to AMG_param
 *
 */
void smoother_setup_biot_monolithic( MG_blk_data *bmgl, AMG_param *param)
{
  // Initialize Schwarz parameters
  Schwarz_param swzparam;
  bmgl->Schwarz_levels = param->Schwarz_levels;
  swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
  swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
  swzparam.Schwarz_type   = param->Schwarz_type;
  swzparam.Schwarz_blksolver = 32;
  bmgl[0].Schwarz.blk_solver = 32;

  dCSRmat Amerge = bdcsr_2_dcsr(&bmgl[0].A);
  bmgl[0].Schwarz.A = dcsr_sympat( &Amerge );
  Schwarz_setup_geometric( &bmgl[0].Schwarz, &swzparam, bmgl[0].fine_level_mesh);
}
