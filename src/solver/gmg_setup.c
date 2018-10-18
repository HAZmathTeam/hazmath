/*! \file src/solver/gmg_setup.c
 *
 *  Geometric Multigrid: SETUP phase
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/24/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/***********************************************************************************************/
/**
 * \fn static void form_linear_R (dCSRmat *tentp,
 *                                ivector *vertices,
 *                                INT      levelNum)
 * \brief Build tentative R for piecewise linear elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param vertices           Pointer to the aggregation of vertices
 * \param levelNum           Level number
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
static void build_linear_R (dCSRmat *R,
                           GMG_data mgl,
                           GMG_param param,
                           INT      nf1d,
                           INT      nc1d)
{
    INT i,j;
    INT ci, cj;
    INT fdof, cdof;
    INT jstart=0;

    // allocate memory for R
    R->row = nc1d;
    R->col = nf1d;
    R->nnz = nc1d*nc1d; //TODO: FIX THIS
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
        }else if( cj == nc1d ){
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
}

/***********************************************************************************************/
/**
 * \fn static void build_constant_R (dCSRmat *tentp,
 *                                  INT     ndof
 *                                  INT     levelNum)
 * \brief Build tentative P for piecewise constant elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param ndof               number of degrees of freedom
 * \param levelNum           Level number
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
static void build_constant_R (dCSRmat *R,
                             GMG_data mgl,
                             GMG_param param,
                             INT      nf1d,
                             INT      nc1d)
{
    int i,j;
    int ci, cj;
    int fdof, cdof;
    int jstart=0;
    // nf1d is number of elements we would have in 1d (number of vertices in 1d minus 1)

    // allocate memory for R
    R->row = nc1d;
    R->col = nf1d;
    R->nnz = nc1d*nc1d; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;

    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d; ci=ci+2){
        cdof = ci + cj*nc1d*2; // TODO: check this

        // Lower Triangle
        fdof = ci*4 + cj*nf1d*4; // TODO: check this
        //Fill
        JA[0+jstart] = fdof;
        JA[1+jstart] = fdof+1;
        JA[2+jstart] = fdof+2;
        JA[3+jstart] = fdof+nf1d;
        for(i=0; i<4; i++){
          val[i+jstart] = 1.0;
        }
        jstart = jstart+4;
        IA[cdof+1] = jstart;

        // Upper Triangle
        cdof = cdof+1;
        fdof = fdof+nf1d+3;
        //Fill
        JA[0+jstart] = fdof-nf1d-2;
        JA[1+jstart] = fdof-2;
        JA[2+jstart] = fdof-1;
        JA[3+jstart] = fdof;
        for(i=0; i<4; i++){
          val[i+jstart] = 1.0;
        }
        jstart = jstart+4;
        IA[cdof+1] = jstart;
      }
    }
}

/***********************************************************************************************/
/**
 * \fn static void build_face_R (dCSRmat *tentp,
 *                                  INT     ndof
 *                                  INT     levelNum)
 * \brief Build tentative P for piecewise constant elements
 *
 * \param tentp              Pointer to the prolongation operators
 * \param ndof               number of degrees of freedom
 * \param levelNum           Level number
 *
 * \author Peter Ohm
 * \date   10/12/2018
 *
 * \note Modified by Peter Ohm on 10/12/2018
 *
 */
static void build_face_R (dCSRmat *R,
                             trimesh  *fmesh,
                             trimesh  *cmesh,
                             INT      nf1d,
                             INT      nc1d)
{
    INT i,j,ii,jj;
    INT ci, cj;
    INT fdof, cdof;
    INT felm, celm;
    INT fface, cface;
    INT felmList[8];
    INT jstart=0;
    INT index;
    // nf1d is number of elements we would have in 1d (number of vertices in 1d minus 1)

    INT rowa;
    INT rowb;
    INT rowaf;
    INT rowbf;

    // Basis stuff
    INT dim = fmesh->dim;
    REAL* x = (REAL *) calloc(dim,sizeof(REAL));
    printf("Made it here, now we want to know if cmesh actually exists\n");
    REAL* phi = (REAL *) calloc(dim*cmesh->f_per_elm,sizeof(REAL));
    printf("Made it here, so cmesh probably exists\n");
    REAL* dphi = (REAL *) calloc(cmesh->f_per_elm,sizeof(REAL));
    REAL value;

    INT* v_on_elm = (INT*)calloc(cmesh->v_per_elm,sizeof(INT));
    printf("Doesn't look like a problem with v_per_elm\n");
    INT* f_on_elm = (INT*)calloc(cmesh->f_per_elm,sizeof(INT));
    printf("Doesn't look like a problem with f_per_elm\n");


    //Garbage
    INT* I;
    INT* J;
    INT* V;

    // allocate memory for R
    R->row = nc1d;
    R->col = nf1d;
    R->nnz = nc1d*nc1d; //TODO: FIX THIS
    R->IA  = (INT *)calloc(R->row+1,sizeof(INT));
    R->JA  = (INT *)calloc(R->nnz, sizeof(INT));
    R->val = (REAL *)calloc(R->nnz, sizeof(REAL));

    printf("We allocated all that crap. The loop is next.\n");

    REAL *val = R->val;
    INT *JA  = R->JA;
    INT *IA  = R->IA;
    // This loops over the coarse and fine elements. Don't question it.
    for(cj = 0; cj<nc1d; cj++){
      for(ci = 0; ci<nc1d*2; ci=ci+2){
        // Lower Triangle and upper triangle box
        celm = ci + cj*nc1d*2; // TODO: check this
        felm = ci*2 + cj*nf1d*4; // TODO: check this
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
          printf("coarse: %d\tfine: %d\n",celm,felm);

          // Get Coarse DOF
          get_incidence_row(celm,cmesh->el_v,v_on_elm);
          get_incidence_row(celm,cmesh->el_f,f_on_elm);
          rowa = cmesh->el_f->IA[celm]-1;
          rowb = cmesh->el_f->IA[celm+1]-1;
          for(j=rowa;j<rowb;j++){
            cface = cmesh->el_f->JA[j]-1;
            // Get Fine DOF
            rowaf = fmesh->el_f->IA[felm]-1;
            rowbf = fmesh->el_f->IA[felm+1]-1;
            for(jj=rowaf;jj<rowbf;jj++){
              fface = fmesh->el_f->JA[jj]-1;
              // Fill P
              value = 0.0;
              x[0] = fmesh->f_mid[fface*dim];
              x[1] = fmesh->f_mid[fface*dim+1];
              rt_basis(phi,dphi,x,v_on_elm,f_on_elm,cmesh);
              for(i=0;i<dim;i++) value += fmesh->f_norm[fface*dim+i]*phi[i];
              //printf("This is dumb, the value is %f\n",value);
              printf("%5d\t%5d\t%f\n",cface,fface,value);
//              I[index] = cface;
//              J[index] = fface;
//              V[index] = value;
            }
          }
        }
      }
    }
}

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

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
  SHORT lvl = 0;

  // 
  INT i,j;
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
    build_face_R(&R, meshHeirarchy[i], meshHeirarchy[i+1], nf1d, nc1d);

    // Setup for next loop
  }

  return 0;

}


/***********************************************************************************************/
SHORT gmg_setup (GMG_data *mgl,
                 GMG_param *param)
{
    SHORT status;
    //SHORT status = gmg_setup_P1(mgl,param);

    return status;
}
