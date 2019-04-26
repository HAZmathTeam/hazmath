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
 *                                  trimesh *fmesh,
 *                                  trimesh *cmesh,
 *                                  INT     nf1d,
 *                                  INT     nc1d)
 * \brief Build tentative P for RT0 elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param fmesh              Fine Trimesh
 * \param cmesh              Coarse Trimesh
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
                             trimesh  *fmesh,
                             trimesh  *cmesh,
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
          rowa = cmesh->el_f->IA[celm]-1;
          rowb = cmesh->el_f->IA[celm+1]-1;
          locFaceId = 0;
          for(j=rowa;j<rowb;j++){
            cface = cmesh->el_f->JA[j]-1;
            // Get Fine DOF
            rowaf = fmesh->el_f->IA[felm]-1;
            rowbf = fmesh->el_f->IA[felm+1]-1;
            for(jj=rowaf;jj<rowbf;jj++){
              fface = fmesh->el_f->JA[jj]-1;
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

/***********************************************************************************************/
INT find_the_fine_vertex_the_hard_way(REAL* midpoint, trimesh* fmesh, INT fface)
{
  // Finds the fine vertex that corresponds with the midpoint of the coarse edge
  INT dim = fmesh->dim;
  INT i;
  INT vertex=-1;
  REAL diff;
  REAL diffmin=999;
  REAL x,y,z=0.0;
  // Allocate
  INT* v_on_f = (INT*)calloc(dim,sizeof(INT));
  get_incidence_row(fface,fmesh->f_v,v_on_f);

  for(i=0;i<dim;i++){
    x = fmesh->cv->x[v_on_f[i]-1] - midpoint[0];
    y = fmesh->cv->y[v_on_f[i]-1] - midpoint[1];
    if(dim==3) z = fmesh->cv->z[v_on_f[i]] - midpoint[2];
    diff = ABS(x+y+z);
    if(diff<diffmin){
      vertex = v_on_f[i]-1;
    }
  }
  if(vertex==-1) printf("\n\nERROR: NO VERTEX MATCHING MIDPOINT FOUND\n\n");
  return vertex;
}
/***********************************************************************************************/
/**
 * \fn void build_bubble_R (dCSRmat *R,
 *                          dCSRmat *Rblx,
 *                          dCSRmat *Rbly,
 *                          trimesh *fmesh,
 *                          trimesh *cmesh,
 *                          INT     nf1d,
 *                          INT     nc1d)
 * \brief Build tentative P for RT0 elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param fmesh              Fine Trimesh
 * \param cmesh              Coarse Trimesh
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
                     trimesh  *fmesh,
                     trimesh  *cmesh,
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
          rowa = cmesh->el_f->IA[celm]-1;
          rowb = cmesh->el_f->IA[celm+1]-1;
          locFaceId = 0;
          for(j=rowa;j<rowb;j++){
            cface = cmesh->el_f->JA[j]-1;
            // Get Fine DOF
            rowaf = fmesh->el_f->IA[felm]-1;
            rowbf = fmesh->el_f->IA[felm+1]-1;
            for(jj=rowaf;jj<rowbf;jj++){
              fface = fmesh->el_f->JA[jj]-1;
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
                x[0] = fmesh->cv->x[v_on_f[ed]-1];
                x[1] = fmesh->cv->y[v_on_f[ed]-1];
                // Evaluate Coarse Bubble at endpoints
                bubble_face_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
                // Save in R linear
                vertex = v_on_f[ed]-1;

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
                if( kk == 16 ){ printf("\n\n\n--------------------------------------------------\nNODE NOT STORED!!!!\n");}
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

              value = ( Fphi[locFaceId*dim+i] - lval ) / fmesh->f_norm[fface*dim+i];

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
INT gmg_setup_RT0(trimesh* fine_level_mesh)
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

  trimesh** meshHeirarchy = (trimesh**)calloc(max_levels, sizeof(trimesh*));
  //meshHeirarchy[0] = mgl->fine_level_mesh;
  meshHeirarchy[0] = fine_level_mesh;

  dCSRmat R;

  // Pre-loop stuff
  for(i=0;i<max_levels-1;i++){
    fSize = sqrt(meshHeirarchy[i]->nv)-1;// Is there a type problem here?
    cSize = (fSize/2);
    //sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",cSize+1);//Need to customize this to specific directory
    sprintf(cgridfile,"/home/pohm01/HAZMAT/hazmath/examples/grids/2D/unitSQ_n%d.haz",cSize+1);//Need to customize this to specific directory
    cgfid = HAZ_fopen(cgridfile,"r");
    meshHeirarchy[i+1] = (trimesh*)calloc(1, sizeof(trimesh));
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
/***********************************************************************************************/
/***********************************************************************************************/
/**
 * \fn static SHORT gmg_blk_setup (MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 *
 * \param mgl    Pointer to AMG_data
 * \param param  Pointer to AMG_param
 *
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 */
//static SHORT gmg_blk_setup(MG_blk_data *mgl,
SHORT gmg_blk_setup(MG_blk_data *mgl,
                            AMG_param *param)
{
    /*
     * TODO:
     * - Can't use brow to determine number of FE spaces since linear and bubble will be merged for A.
     *   Can't use A to determine nf1d for linears and bubbles. (Replace to use mesh?)
     *
     */
    // local variables
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;

    // local variables
    INT           nf1d, nc1d, csize;
    INT           brow,m;
    INT           bm;
    INT           i,j;
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           dim = mgl[lvl].fine_level_mesh->dim;
    REAL          setup_start, setup_end;
    fespace       FE_loc;
    block_fespace FE_blk;

    FE_blk.nspaces = mgl[0].A.brow;

    get_time(&setup_start);

    // dirichlet
    mgl[0].dirichlet = (INT*)calloc(mgl[0].FE->ndof,sizeof(INT));
    mgl[0].dirichlet_blk = (INT**)calloc(FE_blk.nspaces,sizeof(INT*));
    bm = 0;
    INT NoBBL = 0;
    if( mgl[0].gmg_type[0] == 111 ) NoBBL = 1;
    for(i=NoBBL;i<mgl[0].FE->nspaces;i++){
    //for(i=1;i<mgl[0].FE->nspaces;i++){
      //if(i==0) mgl[0].dirichlet_blk[0] = &mgl[0].dirichlet[bm]; // Bubbles + u1 + u2
      if(i==NoBBL) mgl[0].dirichlet_blk[0] = &mgl[0].dirichlet[bm]; // Displacement
      if(i==3) mgl[0].dirichlet_blk[1] = &mgl[0].dirichlet[bm]; // Darcy (RT0)
      if(i==4) mgl[0].dirichlet_blk[2] = &mgl[0].dirichlet[bm]; // Pressure
      for(j=0;j<mgl[0].FE->var_spaces[i]->ndof;j++){
        mgl[0].dirichlet[bm] = mgl[0].FE->var_spaces[i]->dirichlet[j];
        ++bm;
      }
    }

    // To Peter:  this is how to use this fuctions to delete the boundarys -- Xiaozhe 
    block_dCSRmat *B;
    B = (block_dCSRmat *)calloc(1, sizeof(block_dCSRmat));
    bdcsr_delete_rowcol(&mgl[0].A, mgl[0].dirichlet, mgl[0].dirichlet, B);

    // To read from file
    FILE* cgfid;
    char cgridfile[500];

    // Alloc temp block matrices here
    block_dCSRmat tempRblk;// needed to merge bubbles and linears, somehow...
    dCSRmat Rmerge;
    dCSRmat tempRA;

    printf("Beginning Main GMG setup loop...\n");
    /*-- Main GMG setup loop --*/
    while ( lvl < max_levels-1 ) {
      printf("A blocks: %d %d\n",mgl[lvl].A.brow,mgl[lvl].A.bcol);
      brow = mgl[lvl].A.brow;

      // Track dirichlet bdry
      m=0;
      for(i=0; i<brow; i++) m += mgl[lvl].A.blocks[i+i*brow]->row;
      mgl[lvl+1].dirichlet = (INT*)calloc(m,sizeof(INT));
      mgl[lvl+1].dirichlet_blk = (INT**)calloc(brow,sizeof(INT*));

      /*-- Build coarse level mesh --*/
      csize = (sqrt(mgl[lvl].fine_level_mesh->nv)-1)/2 + 1;
      //Need to customize this to specific directory
      sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",csize);
      //sprintf(cgridfile,"/home/xiaozhehu/Work/Projects/HAZMATH/hazmath/examples/grids/2D/unitSQ_n%d.haz",csize);
      //sprintf(cgridfile,"/home/pohm01/HAZMAT/hazmath/examples/grids/2D/unitSQ_n%d.haz",csize);
      cgfid = HAZ_fopen(cgridfile,"r");
      mgl[lvl+1].fine_level_mesh = (trimesh*)calloc(1, sizeof(trimesh));
      initialize_mesh(mgl[lvl+1].fine_level_mesh);
      creategrid_fread(cgfid,0,mgl[lvl+1].fine_level_mesh);
      fclose(cgfid);
      printf("Mesh Loaded for lvl=%d...\n",lvl);
      mgl[0].set_bdry_flags(mgl[lvl+1].fine_level_mesh);

      /*-- Allocate for R A P --*/
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R));
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].P));
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A));
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A_noBC));
      printf("Allocation of R A P for lvl=%d is finished...\n",lvl);

      bm=0; //indexing of dirichlet
      /*-- Form restriction --*/
      for(i=0; i<brow; i++){
        //Check gmg type
        switch( mgl[lvl].gmg_type[i] ) {
          case 0://P0
            //TODO: check nf1d and nc1d calculations
            nf1d = sqrt(mgl[lvl].A.blocks[i+i*brow]->row/2);
            nc1d = (nf1d)/2;
            printf("\tBuilding constant R...\n");
            build_constant_R( mgl[lvl].R.blocks[i+i*brow], nf1d, nc1d);
            printf("\tBuilt constant R...\n");
            // dirichlet
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 0);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, -1,-1); // p
            mgl[lvl+1].dirichlet_blk[i] = &mgl[lvl+1].dirichlet[bm];
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
          break;
          case 1://P1
            nf1d = sqrt(mgl[lvl].A.blocks[i+i*brow]->row);
            nc1d = (nf1d-1)/2 + 1;
            printf("+++++++++++++++++++++++++++++++++++++++++++++++++++ %d %d\n",nf1d,nc1d);
            build_linear_R( mgl[lvl].R.blocks[i+i*brow], nf1d, nc1d);
          break;
          case 30://RT0
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            // Build Here
            printf("\tBuilding RT0 R...\n");
            build_face_R( mgl[lvl].R.blocks[i+i*brow], mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
            //csr_print_matlab(stdout,mgl[lvl].R.blocks[i+i*brow]);
            printf("\tBuilt RT0 R...\n");
            // dirichlet
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 30);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 1,5); // w
            mgl[lvl+1].dirichlet_blk[i] = &mgl[lvl+1].dirichlet[bm];
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
          break;
          case 61://Bubble
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            // Build Here
          break;
          case 111://P1 for each dim then merge
            bdcsr_alloc(dim,dim,&tempRblk);
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
            nc1d = (nf1d-1)/2 + 1;
            printf("+++++++++++++++++++++++++++++++++++++++++++++++++++ %d %d\n",nf1d,nc1d);
            for(j=0; j<dim; j++) build_linear_R( tempRblk.blocks[j+j*tempRblk.brow], nf1d, nc1d);
            tempRblk.blocks[1] = NULL;
            tempRblk.blocks[2] = NULL;
            Rmerge = bdcsr_2_dcsr(&tempRblk);
            printf("\tMerged Displacement R... storage in %d\n",i+i*brow);
            dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
            dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
            dcsr_free(&Rmerge);

            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 1);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 5,5); // Ux
            mgl[lvl+1].dirichlet_blk[i] = &mgl[lvl+1].dirichlet[bm];
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 1);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 1,4); // Uy
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
          break;
          case 999:// P1+bubble
            bdcsr_alloc(dim+1,dim+1,&tempRblk);
            tempRblk.blocks[3] = NULL;
            tempRblk.blocks[5] = NULL;
            tempRblk.blocks[6] = NULL;
            tempRblk.blocks[7] = NULL;
            //P1
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
            nc1d = (nf1d-1)/2 + 1;
            for(j=1; j<dim+1; j++) build_linear_R( tempRblk.blocks[j+j*tempRblk.brow], nf1d, nc1d);
            printf("\tBuilt Linear R...\n");
            //Bubble
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            build_bubble_R( tempRblk.blocks[0], tempRblk.blocks[1], tempRblk.blocks[2],
                            mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
            printf("\tBuilt Bubble R...\n");
            Rmerge = bdcsr_2_dcsr(&tempRblk);
            printf("\tMerged Displacement R... storage in %d\n",i+i*brow);
            dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
            dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
            dcsr_free(&Rmerge);
            // BC flag dirichlet elim stuff
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 61);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 1,5); //bbl
            mgl[lvl+1].dirichlet_blk[i] = &mgl[lvl+1].dirichlet[bm];
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 1);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 5,5); // Ux
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
            create_fespace(&FE_loc, mgl[lvl+1].fine_level_mesh, 1);
            set_dirichlet_bdry(&FE_loc, mgl[lvl+1].fine_level_mesh, 1,4); // Uy
            for(j=0;j<FE_loc.ndof;j++){
              mgl[lvl+1].dirichlet[bm] = FE_loc.dirichlet[j];
              ++bm;
            }
            free_fespace(&FE_loc);
          break;
          default:
            printf("### ERROR: Unknown geometric multigrid type: %d!\n",mgl[lvl].gmg_type[i]);
          break;
        }// switch gmg_type
      }//for i<brow
      printf("Built R for lvl=%d...\n",lvl);

      /*-- Form Prolongation --*/
      for(i=0; i<brow; i++){
          dcsr_trans(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].P.blocks[i+i*brow]);
      }
      printf("Built P for lvl=%d...\n",lvl);
          //csr_print_matlab(stdout,mgl[lvl].A_noBC.blocks[0]);

      /*-- Form coarse level stiffness matrix --*/
      for(i=0; i<brow; i++){
        for(j=0; j<brow; j++){
          if(mgl[lvl].A.blocks[j+i*brow]){
            printf("RAP on block[%d,%d]: TODO: FIX NNZ FOR R(i think)\n",i,j);
            printf("Matrix:\n\tR.nnz=%d\tR.row=%d\tR.col=%d\n\tA.nnz=%d\tA.row=%d\tA.col=%d\n\tP.nnz=%d\tP.row=%d\tP.col=%d\n\n",
                    mgl[lvl].R.blocks[i+i*brow]->nnz,
                    mgl[lvl].R.blocks[i+i*brow]->row,
                    mgl[lvl].R.blocks[i+i*brow]->col,
                    mgl[lvl].A.blocks[j+i*brow]->nnz,
                    mgl[lvl].A.blocks[j+i*brow]->row,
                    mgl[lvl].A.blocks[j+i*brow]->col,
                    mgl[lvl].P.blocks[j+j*brow]->nnz,
                    mgl[lvl].P.blocks[j+j*brow]->row,
                    mgl[lvl].P.blocks[j+j*brow]->col);
            if(i==j){
              dcsr_rap(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].A_noBC.blocks[j+i*brow], mgl[lvl].P.blocks[j+j*brow], mgl[lvl+1].A_noBC.blocks[j+i*brow]);
              printf("RAP finished on block...\n");
            } else {
              dcsr_mxm(mgl[lvl].R.blocks[i+i*brow],mgl[lvl].A_noBC.blocks[j+i*brow],&tempRA);
              dcsr_mxm(&tempRA,mgl[lvl].P.blocks[j+j*brow],mgl[lvl+1].A_noBC.blocks[j+i*brow]);
              dcsr_free(&tempRA);
            }
          } else { mgl[lvl+1].A_noBC.blocks[j+i*brow] = NULL; }
        }//j
      }//i
      printf("Built RAP for lvl=%d...\n",lvl);
      /*-- Eliminate dirichlet boundaries from stiffness matrix --*/
      bdcsr_cp( &mgl[lvl+1].A_noBC, &mgl[lvl+1].A);
      FE_blk.dirichlet = mgl[lvl+1].dirichlet;
      printf("eliminating BC on coarse level\n");
      bdcsr_shift(&mgl[lvl+1].A,  1);
      //eliminate_DirichletBC_blockFE_blockA(NULL,&FE_blk,&mgl[lvl+1].fine_level_mesh,NULL,&mgl[lvl+1].A,0.0);
      eliminate_DirichletBC_blockFE_blockA(NULL,&FE_blk,mgl[lvl+1].fine_level_mesh,NULL,&mgl[lvl+1].A,0.0);
      bdcsr_shift(&mgl[lvl+1].A, -1);


      // PRINT MATRIX
      if(lvl==0){
          FILE* matid = HAZ_fopen("Acoarse.dat","w");
          csr_print_matlab(matid,mgl[lvl+1].A.blocks[0]);
          fclose(matid);
      }

      ++lvl;
    } // lvl

    // Setup coarse level systems for direct solvers
    switch (csolver) {

#if WITH_SUITESPARSE
        printf("Setting up coarse solve: Using UMFPACK...\n");
        case SOLVER_UMFPACK: {
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
//        print_amg_complexity(mgl,prtlvl);
        print_cputime("geometric multigrid setup", setup_end - setup_start);
    }

    return status;
}

/***********************************************************************************************/
//SHORT gmg_setup (AMG_data *mgl,
//                 GMG_param *param)
//{
//    SHORT status;
//    //SHORT status = gmg_setup_P1(mgl,param);
//
//    return status;
//}
