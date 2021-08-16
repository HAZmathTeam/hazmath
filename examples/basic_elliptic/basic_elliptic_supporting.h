/*! \file examples/basic_elliptic/basic_elliptic_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program generates simplicial grids in 2,3,4... dimension.
 *
 * \note mesh making for basic elliptic
 */
/*********************************************************************/
/*!
 * \fn static mesh_struct make_uniform_mesh(const INT dim,
 *                                          const INT mesh_ref_levels,
 *                                          const INT mesh_ref_type,
 *       				    const INT set_bndry_codes)
 *
 * \brief makes a mesh_ref_levels of refined mesh on the unit cube in dimension dim. 
 *
 * \param dim               input array
 * \param mesh_ref_levels   number of refinement levels
 *
 * \param mesh_ref_type     if > 10: uniform refinement ;
 *                          if < 10: nearest vertex bisection ; 
 *
 * \param set_bndry_codes   if .eq. 0: the boundary codes come from sc->bndry[];
 *                          if .ne. 0  the boundary codes are set 
 *
 * \return mesh_struct with the mesh. 
 *
 * \note 20210815 (ltz)
 *
 */
static mesh_struct  make_uniform_mesh(const INT dim,			\
				      const INT mesh_ref_levels,	\
				      const INT mesh_ref_type,		\
				      const INT set_bndry_codes)
{
  INT i;
  INT jlevel,k;
  scomplex *sc=NULL,*sctop=NULL,*dsc=NULL;
  /* 
   * get the coarsest mesh on the cube in dimension dim and set the
   * refinement type.
   */
  sc=mesh_cube_init(dim,mesh_ref_type);
  if(sc->ref_type>10){
    /* Uniform refinement only for dim=2 or dim=3*/
    if(dim==3){
      for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
	uniformrefine3d(sc);
	sc_vols(sc);
      }
    } else if(dim==2){
      for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
	//	fprintf(stdout,".%d.",jlevel+1);fflush(stdout);
	uniformrefine2d(sc);
	sc_vols(sc);
      }
    } else {
      fprintf(stderr,"\n%%%% ERROR: WRONG spatial dimension for uniform refinement");
    }
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    sc_vols(sc);
  } else {
    ivector marked; marked.val=NULL;
    for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
      /* choose the finest grid */
      sctop=scfinest(sc);
      marked.row=sctop->ns;
      marked.val=realloc(marked.val,marked.row*sizeof(INT));
      /*mark everything*/
      for(k=0;k<marked.row;k++) marked.val[k]=TRUE;
      // now we refine.
      refine(1,sc,&marked);
      /* free */
      haz_scomplex_free(sctop);
    }
    ivec_free(&marked);
    scfinalize(sc);
    sc_vols(sc);
  }  
  find_cc_bndry_cc(sc,set_bndry_codes);
  for(i=0;i<sc->nv;++i){
    if(sc->bndry[i]>128) sc->bndry[i]-=128;
  }
  mesh_struct mesh0=sc2mesh(sc);
  haz_scomplex_free(sc);  
  return mesh0;
}
