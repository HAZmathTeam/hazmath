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
 * \fn void initialize_ni(nested_it* ni)
 *
 * \brief Initialize the nested iteration struct
 *
 * \param ni: nested iteration struct containing current mesh, simplicial complexes, and FE info
 *
 * \return ni: allocate the components of the nested iteration struct
 *
 */
void initialize_ni(nested_it* ni)
{
  // Mesh
  ni->mesh=malloc(sizeof(mesh_struct));

  // FE space
  ni->FE=malloc(sizeof(block_fespace));

  // Quadrature
  ni->cq=malloc(sizeof(qcoordinates));

  // Assume at least 1 level of nested iteration
  ni->ni_levels=1;

  // Assume uniform refinement by default
  ni->mark_type=0;
  ni->mark_param=0.0;

  // Error estimator for AMR
  ni->err_est=malloc(sizeof(REAL));
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
  // Get the coarsest mesh on the cube in dimension dim and set the refinement type.
  // At the moment there are two types of refinement: > 10 uniform refinement; < 10 longest edge
  ni->sc_all=mesh_cube_init(dim,1,8); // we avoid > 10 because it does not complete sc structure
  scomplex* sc=ni->sc_all[0]; // Grab the finest level 
  
  // Perform init_ref_levels of refinement to get initial mesh (mark all elements to obtain cris-cross grid)
  ivector marked;
  marked=ivec_create(sc->ns);
  ivec_set(sc->ns,&marked,1);
  refine(init_ref_levels,sc,NULL);

 // OLD STUFF using uniformrefine which has boundary code issues
 // if(dim==2){
   //for(i=0;i<init_ref_levels;++i){
      //uniformrefine2d(sc);
      //sc_vols(sc); // compute volumes of simplices
   // }
  // } else if(dim==3){
  //   for(i=0;i<init_ref_levels;++i){
  //     uniformrefine3d(sc);
  //     sc_vols(sc);
  //   }
  // } else {
  //   check_error(ERROR_DIM, __FUNCTION__);
  // }

  // Get boundary codes and compute volumes of simplices TODO: LTZ check on scfinalize version
  scfinalize_nofree(sc,(INT )1);
  sc_vols(sc);

  // Convert to mesh_struct for FEM assembly
  //ni->mesh=malloc(sizeof(mesh_struct));
  ni->mesh[0]=sc2mesh(sc);
  build_mesh_all(ni->mesh);

  return;
}

/*!
 * \fn void ni_refine_mesh(nested_it* ni)
 *
 * \brief Refine mesh for next level in the nested iteration process
 *
 * \param ni: nested iteration struct containing current mesh, simplicial complexes, and FE info
 *
 * \return ni: updates the mesh for next NI level
 *
 */
void ni_refine_mesh(nested_it* ni)
{
  // Get to the finest grid in simplicial complex hierarchy (current mesh)
  scomplex *sc=NULL,*sctop=NULL; // Simplicial complex (all, current, finest)
  sc = ni->sc_all[0];
  sctop=scfinest(sc);

  // Mark element using error estimator if needed
  ivector marked;
  switch(ni->mark_type) {
  case 1: // maximal mark
    marked = amr_maximal_mark(sctop, ni->err_est, ni->mark_param); // elements with error > then mark_param*100% of largest error
    break;
  default:
    marked=ivec_create(sc->ns);
    ivec_set(sc->ns,&marked,1);
    break;
  }

  // Refine the grid once based on who's marked and convert to a mesh again
  refine(1,sc,&marked);
  scfinalize_nofree(sc,(INT )1);
  sc_vols(sc);

  // Convert to mesh_struct for FEM assembly
  ni->mesh[0]=sc2mesh(sc);
  build_mesh_all(ni->mesh);

  // Update the quadrature on the mesh using same order as previous level
  ni->cq = get_quadrature(ni->mesh,ni->cq->nq1d);

  return;
}

/*!
 * \fn void ni_update_sol(nested_it* ni,INT dim)
 *
 * \brief Updates solutions to the next level in the nested iteration process
 *
 * \param ni: nested iteration struct containing current mesh, simplicial complexes, and FE info
 * \param dim:  Dimension of problem
 *
 * \return ni: updates solution for next NI level
 * 
 * \note This assumes mesh is updated and FE spaces, as these are problem-dependent, are updated by user
 *
 */
void next_update_sol(nested_it* ni,INT dim)
{
 
  // Update solution
  return;
}

/*!
* \fn void free_ni(nested_it* ni)
*
* \brief Frees memory of ni struct
*
* \param ni     Pointer to ni struct to be freed
*
*/
void free_ni(nested_it* ni)
{
  if(ni==NULL) return;

  // Free Mesh
  if(ni->mesh) {
    free_mesh(ni->mesh);
    free(ni->mesh);
    ni->mesh=NULL;
  }

  // Figure out how to free simplicial complex (sc_all)
  
  // Free quadrature
  if(ni->cq) {
    free_qcoords(ni->cq);
    free(ni->cq);
    ni->cq=NULL;
  }

  // Free FE space
  if(ni->FE) {
    free_blockfespace(ni->FE);
    free(ni->FE);
    ni->FE=NULL;
  }

  // Free error estimate array
  if(ni->err_est) {
    free(ni->err_est);
    ni->err_est=NULL;
  }

  return;
}