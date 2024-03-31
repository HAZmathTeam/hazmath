/*! \file src/amr/nested_iteration.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20240328.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note contains functions for setting up nested iteration runs.
 *
 */
#include "hazmath.h"

/*!
 * \fn void scfinalize_temp(scomplex *sc,const INT set_bndry_codes)
 *
 * \brief Removes all hierachy and make sc to represent only the final
 *        grid. computes connected components and connected components
 *        on the boundary.
 *
 * \param sc: simplicial complex
 * \param set_bndry_codes: if 0 then create the sparse matrix for all vertices;
 *
 * \note This is a modified version of the real scfinalize avoiding some memory issues with the boundary codes until that is fixed.
 *
 */
void scfinalize_temp(scomplex *sc,const INT set_bndry_codes)
{
  // INT n=sc->n;
  INT ns,j=-10,k=-10;
  INT n1=sc->n+1;
  /*
      store the finest mesh in sc structure.
      on input sc has all the hierarchy, on return sc only has the final mesh.
  */
  //  free(sc->parent_v->val);  sc->parent_v->val=NULL;
  ns=0;
  for (j=0;j<sc->ns;j++){
    /*
      On the last grid are all simplices that were not refined, so
      these are the ones for which child0 and childn are not set.
    */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++) {
	      sc->nodes[ns*n1+k]=sc->nodes[j*n1+k];
      }
      sc->child0[ns]=-1;
      sc->childn[ns]=-1;
      sc->gen[ns]=sc->gen[j];
      sc->flags[ns]=sc->flags[j];
      ns++;
    }
  }
  sc->ns=ns;
  sc->nodes=realloc(sc->nodes,n1*sc->ns*sizeof(INT));
  sc->nbr=realloc(sc->nbr,n1*sc->ns*sizeof(INT));
  sc->vols=realloc(sc->vols,sc->ns*sizeof(REAL));
  sc->child0=realloc(sc->child0,sc->ns*sizeof(INT));
  sc->childn=realloc(sc->childn,sc->ns*sizeof(INT));
  sc->gen=realloc(sc->gen,sc->ns*sizeof(INT));
  sc->flags=realloc(sc->flags,sc->ns*sizeof(INT));
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  // this also can be called separately
  // set_bndry_codes should always be set to 1.
  //  set_bndry_codes=1;
  find_cc_bndry_cc(sc,(INT )1); //set_bndry_codes);
  //
  /* if(set_bndry_codes){ */
  /*   for(j=0;j<sc->nv;++j){ */
  /*     if(sc->bndry[j]>128) sc->bndry[j]-=128; */
  /*   } */
  /* } */
  // clean up: TODO: DO NOT FREE ANYTHING UNTIL LUDMIL FIXES!
  // icsr_free(sc->bndry_v);
  // free(sc->bndry_v);
  // sc->bndry_v=NULL;
  // icsr_free(sc->parent_v);
  // free(sc->parent_v);
  // sc->parent_v=NULL;
  return;
}

/*!
 * \fn void get_initial_mesh_ni(nested_it* ni,INT dim,INT init_ref_levels)
 *
 * \brief Gets an initial uniform mesh to start the nested iteration process
 *
 * \param ni: nested iteration struct containing current mesh and simplicial complexes
 * \param dim:  Dimension of problem
 * \param init_ref_levels: Number of initial uniform refinements to get to first mesh
 *
 * \return ni: creates the initial simplicial complex and mesh for nested iteration
 *
 */
void get_initial_mesh_ni(nested_it* ni,INT dim,INT init_ref_levels)
{

  INT i;
  // Get the coarsest mesh on the cube in dimension dim and set the refinement type.
  ni->sc_all=mesh_cube_init(dim,1,11); // check why is this 11
  scomplex* sc=ni->sc_all[0];
  
  // uniform refine
  if(dim==2){
    for(i=0;i<init_ref_levels;++i){
      uniformrefine2d(sc);
      sc_vols(sc);
    }
  } else if(dim==3){
    for(i=0;i<init_ref_levels;++i){
      uniformrefine3d(sc);
      sc_vols(sc);
    }
  } else {
    check_error(ERROR_DIM, __FUNCTION__);
  }
  // Get boundary codes TODO: LTZ check on scfinalize version
  scfinalize_temp(sc,(INT )1);
  sc_vols(sc);

  // Convert to mesh_struct for FEM assembly
  ni->mesh=malloc(sizeof(mesh_struct));
  ni->mesh[0]=sc2mesh(sc);
  build_mesh_all(ni->mesh);

  return;
}