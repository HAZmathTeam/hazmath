/*! \file examples/elliptic_p1/mesh_cube_init.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note meshes the unit cube in R(d), d=2,3,4,5,....
 *
 *  \note: modified by ltz on 20210601
 *
 */
/**********************************************************************/
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn scomplex *meshes_init(const INT dim)
 *
 * \brief creates initial meshes for the cube in 2,3,4,5 dimensions. 
 *
 * \param dim:   the spatial dimension
 *
 * \return the simplicial complex with the initial mesh.  
 *
 * \note
 *
 */
/**********************************************************************/
/*!
 * \fn scomplex *mesh_cube_init(const INT dim, const INT ref_type)
 *
 * \brief creates initial meshes for the cube in 2,3,4,5 dimensions. 
 *
 * \param dim: the spatial dimension \param ref_type: the refinement
 *             type: if >10, then uniform refinement; if <10, then
 *             newest vertex bisection. Check include/amr.h scomplex
 *             structure for more info.
 *
 * \return the simplicial complex with the initial mesh.  
 *
 * \note
 *
 */
/****************************************************************/
scomplex *mesh_cube_init(const INT dim, const INT ref_type)
{
  INT i,j,in;
  scomplex *sc=NULL;
  cube2simp *c2s=cube2simplex(dim);
  sc=haz_scomplex_init(dim,c2s->ns,c2s->nvcube,dim);
  memcpy(sc->nodes,c2s->nodes,(dim+1)*sc->ns*sizeof(INT));
  for(i=0;i<sc->nv;++i){
    in=i*sc->n;
    for(j=0;j<sc->n;j++){
      sc->x[in+j]=(REAL )c2s->bits[in+j];
    }
  }
  cube2simp_free(c2s);
  find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
  sc_vols(sc);
  if(sc->n > 3)
    sc->ref_type=0;
  else
    sc->ref_type=ref_type;
  return sc;
}
