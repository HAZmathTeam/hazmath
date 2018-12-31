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

    //Test Print
    //INT ja,jb;
    //for(i=0; i<R->row; i++){
    //    ja = IA[i];
    //    jb = IA[i+1];
    //    for(j=ja; j<jb; j++){
    //        printf("%3d | %3d, %5f\n",i,R->JA[j],R->val[j]);
    //    }
    //}
    //printf("nnz: %d\t jstart: %d\n",R->nnz,jstart);
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

        //printf("C: %4d\tF:%4d, %4d, %4d, %4d\n",cdof,fdof,fdof+1,fdof+2,fdof+nf1d*2);

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

        //printf("C: %4d\tF:%4d, %4d, %4d, %4d\n",cdof,fdof+3,fdof+nf1d*2+1,fdof+nf1d*2+2,fdof+nf1d*2+3);
      }
    }
    //printf("Row: %d\tCol: %d\tnnz: %d\tIA.end = %d\t%d\n",R->row,R->col,R->nnz,cdof+1,jstart);
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
    INT i,j,ii,jj;
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
    R->nnz = R->row*9; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    I = (INT *)calloc(R->nnz,sizeof(INT));
    J = (INT *)calloc(R->nnz,sizeof(INT));
    V = (REAL *)calloc(R->nnz,sizeof(REAL));
    IJfilled = (INT *)calloc(R->col*R->row,sizeof(INT));//TODO:BAD!!!!!
    for(i=0; i<R->row*R->col; i++) IJfilled[i] = -1; // Set all to -1 //TODO:BAD!!!!
//    printf("We allocated all that. The loop is next.\n");

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
          //printf("coarse: %d\tfine: %d\n",celm,felm);

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

              //printf("%5d\t%5d\t%f\t\t%f, %f\t%d\t%f, %f \t%f, %f\n",cface,fface,value,
              //       fmesh->f_norm[fface*dim+0],fmesh->f_norm[fface*dim+1],locFaceId,
              //       cmesh->f_norm[cface*dim+0],cmesh->f_norm[cface*dim+1],
              //       phi[locFaceId*dim+0],phi[locFaceId*dim+1]);

              if( IJfilled[cface*R->row + fface] == -1 ) {
                IJfilled[cface*R->row + fface] = 1;
                if(ABS(value) > 1e-8){
                  //printf("Indexing|%3d\t%3d\n",cface,fface);
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

    //printf("index: %d\tnnz: %d\n",index,R->nnz);//TODO: Correct number of nnz in the matrix

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
    x = fmesh->cv->x[v_on_f[i]] - midpoint[0];
    y = fmesh->cv->y[v_on_f[i]] - midpoint[1];
    if(dim==3) z = fmesh->cv->z[v_on_f[i]] - midpoint[2];
    diff = ABS(x+y+z);
    if(diff<diffmin){
      vertex = v_on_f[i];
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
    INT nq1d = 4;
    INT quad;
    REAL qval;
    qcoordinates *cqface = get_quadrature_boundary(fmesh,nq1d,2);//TODO:free

    // Basis stuff
    INT dim = fmesh->dim;
    REAL* x = (REAL *) calloc(dim,sizeof(REAL));
    REAL* phi = (REAL *) calloc(dim*cmesh->f_per_elm,sizeof(REAL));
    REAL* dphi = (REAL *) calloc(dim*dim*cmesh->f_per_elm,sizeof(REAL));
    REAL value;
    REAL lval;

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
    R->nnz = R->row*9; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));
    // Rename for convenience
    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;
    // allocate helper arrays.
    I = (INT *)calloc(R->nnz,sizeof(INT));
    J = (INT *)calloc(R->nnz,sizeof(INT));
    V = (REAL *)calloc(R->nnz,sizeof(REAL));
    IJfilled = (INT *)calloc(R->col*R->row,sizeof(INT));//TODO:BAD!!!!!
    for(i=0; i<R->row*R->col; i++) IJfilled[i] = -1;

    // allocate memory for Rl1
    // allocate memory for Rl2
    INT vertex;
    INT* VertFilled = (INT *)calloc(fmesh->nv,sizeof(INT));//TODO:BAD!!!!!
    INT indexl=0;
    INT*  Il1 = (INT*)calloc(fmesh->nv,sizeof(INT));
    INT*  Jl1 = (INT*)calloc(fmesh->nv,sizeof(INT));
    REAL* Vl1 = (REAL*)calloc(fmesh->nv,sizeof(REAL));
    INT*  Il2 = (INT*)calloc(fmesh->nv,sizeof(INT));
    INT*  Jl2 = (INT*)calloc(fmesh->nv,sizeof(INT));
    REAL* Vl2 = (REAL*)calloc(fmesh->nv,sizeof(REAL));


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
//          printf("ELEMENT | coarse: %d\tfine: %d\n",celm,felm);

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
              for(quad=0;quad<cqface->nq_per_elm;quad++){
                qval = 0.0;
                x[0] = cqface->x[fface*cqface->nq_per_elm + quad];
                x[1] = cqface->y[fface*cqface->nq_per_elm + quad];
                bubble_face_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
                for(i=0;i<dim;i++){
                  qval += fmesh->f_norm[fface*dim+i]*phi[locFaceId*dim+i];
                }
                value += qval * cqface->w[fface*cqface->nq_per_elm+quad];
              }
              //value = value / cmesh->f_area[cface];
              //value = value / fmesh->f_area[fface];

              // Get Linear Part
              // note: 2D ONLY: only nonzero value of bubble on vertices is on the vertex in the middle of the coarse edge.
              x[0] = cmesh->f_mid[cface*dim];  //TODO: is cface the same as given by locFaceId?
              x[1] = cmesh->f_mid[cface*dim+1];
              if(dim==3) x[2] = cmesh->f_mid[cface*dim+2];
              bubble_face_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
              // Integrate linear part over fine edge (midpoint rule)
              lval = 0.0;
              for(i=0;i<dim;i++) lval += phi[locFaceId*dim+i]*fmesh->f_norm[fface*dim+i];
              lval = (1.0/dim)*fmesh->f_area[fface] * lval;

              // Subtract off linear part (Bubbles used to enrich P1)
              value = value - lval;


              value = value / fmesh->f_area[fface];

              // Save in R linear
              vertex = find_the_fine_vertex_the_hard_way(x, fmesh, fface);
              if( VertFilled[vertex] == 0 ){
//                printf("Vertex | (%d, %d) : %7.4f, %7.4f\n",cface,vertex,phi[locFaceId*dim],phi[locFaceId*dim+1]);
                VertFilled[vertex] = 1;
                Il1[indexl] = cface;
                Jl1[indexl] = vertex;
                Vl1[indexl] = phi[locFaceId*dim+0];
                Il2[indexl] = cface;
                Jl2[indexl] = vertex;
                Vl2[indexl] = phi[locFaceId*dim+1];
                indexl++;
              }

              if( IJfilled[cface*R->row + fface] == -1 ) {
                IJfilled[cface*R->row + fface] = 1;
                if(ABS(value) > 1e-8){
//                  printf("Indexing | (%3d, %3d) : %7.4f\n",cface,fface,value);
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

//    printf("index: %d\tnnz: %d\n",index,R->nnz);//TODO: Correct number of nnz in the matrix

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

    // allocate memory for R
    Rblx->row = cmesh->nface;
    Rblx->col = fmesh->nv;
    Rblx->nnz = R->row*9; //TODO: FIX THIS
    Rblx->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    Rblx->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    Rblx->val = (REAL *)calloc(R->nnz, sizeof(REAL));
    // allocate memory for R
    Rbly->row = cmesh->nface;
    Rbly->col = fmesh->nv;
    Rbly->nnz = R->row*9; //TODO: FIX THIS
    Rbly->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    Rbly->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    Rbly->val = (REAL *)calloc(R->nnz, sizeof(REAL));
    // Fill CSR matrix for Bubble Linear
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
    sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",cSize+1);//Need to customize this to specific directory
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
    const SHORT min_cdof   = MAX(param->coarse_dof,50);

    // local variables
    INT           nf1d, nc1d, csize;
    INT           brow,m;
    INT           i,j;
    SHORT         max_levels = param->max_levels, lvl = 0, status = SUCCESS;
    INT           dim = mgl[lvl].fine_level_mesh->dim;
    REAL          setup_start, setup_end;
    Schwarz_param swzparam;

    INT regionType = 0;
    // To read from file
    FILE* cgfid;
    char cgridfile[500];

    // Alloc temp block matrices here
    block_dCSRmat tempRblk;// needed to merge bubbles and linears, somehow...
    dCSRmat Rmerge;

    get_time(&setup_start);

    printf("Beginning Main GMG setup loop...\n");
    /*-- Main GMG setup loop --*/
    while ( lvl < max_levels-1 ) {
      brow = mgl[lvl].A.brow;
      /*-- Build coarse level mesh --*/
      csize = (sqrt(mgl[lvl].fine_level_mesh->nv)-1)/2 + 1;
      //Need to customize this to specific directory
      sprintf(cgridfile,"/Users/Yggdrasill/Research/HAZMAT/hazmat/examples/grids/2D/unitSQ_n%d.haz",csize);
      cgfid = HAZ_fopen(cgridfile,"r");
      mgl[lvl+1].fine_level_mesh = (trimesh*)calloc(1, sizeof(trimesh));
      initialize_mesh(mgl[lvl+1].fine_level_mesh);
      creategrid_fread(cgfid,0,mgl[lvl+1].fine_level_mesh);
      fclose(cgfid);
      printf("Mesh Loaded for lvl=%d...\n",lvl);
      /*-- Allocate for R A P --*/
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].R));
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl].P));
      bdcsr_alloc(mgl[lvl].A.brow, mgl[lvl].A.bcol, &(mgl[lvl+1].A));
      printf("Allocation of R A P for lvl=%d is finished...\n",lvl);
      /*-- Form restriction --*/
      for(i=0; i<brow; i++){
        //Check gmg type
        switch( mgl[lvl].gmg_type[i] ) {
          case 0://P0
              //TODO: check nf1d and nc1d calculations
            nf1d = sqrt(mgl[lvl].A.blocks[i+i*brow]->row/2);
            nc1d = (nf1d)/2;
            printf("Building constant R...\n");
            build_constant_R( mgl[lvl].R.blocks[i+i*brow], nf1d, nc1d);
            printf("Built constant R...\n");
            break;
          case 1://P1
            nf1d = sqrt(mgl[lvl].A.blocks[i+i*brow]->row);
            nc1d = (nf1d-1)/2 + 1;
            build_linear_R( mgl[lvl].R.blocks[i+i*brow], nf1d, nc1d);
            break;
          case 30://RT0
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            // Build Here
            printf("Building RT0 R...\n");
            build_face_R( mgl[lvl].R.blocks[i+i*brow], mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
            printf("Built RT0 R...\n");
            break;
          case 61://Bubble
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            // Build Here
            //build_bubble_R( mgl[lvl].R.blocks[i+i*brow], mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
            break;
          case 999:// P1+bubble
            bdcsr_alloc(dim+1,dim+1,&tempRblk);
            //P1
            // TODO:do for each dim
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv);
            nc1d = (nf1d-1)/2 + 1;
            for(j=1; j<dim+1; j++) build_linear_R( tempRblk.blocks[i+i*tempRblk.brow], nf1d, nc1d);
            printf("Built Linear R...\n");
            //Bubble
            nf1d = sqrt(mgl[lvl].fine_level_mesh->nv)-1;
            nc1d = (nf1d)/2;
            build_bubble_R( tempRblk.blocks[0], tempRblk.blocks[tempRblk.brow], tempRblk.blocks[2*tempRblk.brow],
                            mgl[lvl].fine_level_mesh, mgl[lvl+1].fine_level_mesh, nf1d, nc1d);
            printf("Built Bubble R...\n");
            Rmerge = bdcsr_2_dcsr(&tempRblk);
            printf("Merged Displacement R... storage in %d\n",i+i*brow);
            dcsr_alloc(Rmerge.row,Rmerge.col,Rmerge.nnz,mgl[lvl].R.blocks[i+i*brow]);
            dcsr_cp(&Rmerge,mgl[lvl].R.blocks[i+i*brow]);
            printf("Stored Displacement R...\n");
            break;
          default:
            printf("### ERROR: Unknown geometric multigrid type: %d!\n",mgl[lvl].gmg_type[i]);
            break;
        }// switch
      }//for brow

      printf("Built R for lvl=%d...\n",lvl);

      /*-- Form Prolongation --*/
      for(i=0; i<brow; i++){
          dcsr_trans(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].P.blocks[i+i*brow]);
      }
      printf("Built P for lvl=%d...\n",lvl);
      /*-- Form coarse level stiffness matrix --*/
      for(i=0; i<brow; i++){
        for(j=0; j<brow; j++){
          if(mgl[lvl].A.blocks[i+j*brow]){
            printf("RAP on block[%d,%d]: TODO: FIX NNZ FOR R(i think)\n",i,j);
            dcsr_rap(mgl[lvl].R.blocks[i+i*brow], mgl[lvl].A.blocks[i+j*brow], mgl[lvl].P.blocks[j+j*brow], mgl[lvl+1].A.blocks[i+j*brow]);
          }
        }//j
      }//i
      printf("Built RAP for lvl=%d...\n",lvl);
                             
      ++lvl;
    } // lvl

    // Setup coarse level systems for direct solvers
    switch (csolver) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            // merge blocks
            mgl[lvl].Ac = bdcsr_2_dcsr(&mgl[lvl].A)
            dCSRmat Ac_tran;
            dcsr_trans(&mgl[lvl].Ac, &Ac_tran);
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
SHORT gmg_setup (AMG_data *mgl,
                 GMG_param *param)
{
    SHORT status;
    //SHORT status = gmg_setup_P1(mgl,param);

    return status;
}
