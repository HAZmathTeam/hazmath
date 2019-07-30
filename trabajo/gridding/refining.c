/*! \file src/amr/refining.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 * \note: routines used to refine a simplicial grid grid ref_level times.

 */
#include "hazmath.h"
//#include "grid_defs.h"
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void n_refine(INT amr_ref_type, INT ref_levels, scomplex *sc,	\
	      void *anything, 
	      void (*solving)(INT , scomplex *, void *),		\
	      void (*estimating)(INT , scomplex *, void *),	\
	      void (*marking)(INT , scomplex *, void *))
{
/* 
 *
 * refine ref_levels. Follows the algorithm:
 * SOLVE-->ESTIMATE-->MARK-->REFINE. 
 *
 */
  if(ref_levels<=0) return;
  /* 
     "anything" is a void pointer to structure or array or bunch of
     arrays of mixed type which are used in solve, estimate and mark
     phases.
  */
  /**/
  INT j=-1,nsold,print_level=0;
  if(!sc->level){
    /*form neighboring list; */
    find_nbr(sc->ns,sc->nv,sc->n,sc->nodes,sc->nbr);
    //    haz_scomplex_print(sc,0,__FUNCTION__);  fflush(stdout);
    INT *wrk=calloc(5*(sc->n+2),sizeof(INT));
    /* construct bfs tree for the dual graph */
    abfstree(0,sc,wrk,print_level);
    //    haz_scomplex_print(sc,0,__FUNCTION__);fflush(stdout);
    //    exit(100);
    if(wrk) free(wrk);
  }
  //INT  n=sc->n;
  fprintf(stdout,"refine: ");
  while(sc->level < ref_levels && TRUE){
    //    if(dxstar->row && dxstar->val) markstar(sc,dxstar);
    /* SOLVE */
    (*solving)(amr_ref_type, sc, anything);
    /* ESTIMATE */
    (*estimating)(amr_ref_type,sc, anything);
    /* MARK */
    (*marking)(amr_ref_type, sc, anything);
      /* for(i = 0;i<sc->ns;i++){ */
      /* //    if(sc->gen[i] < level) continue;     */
      /* fprintf(stdout,"\n%s: ZZZZZZZZZZZzsimplex %d (mark=%d, gen=%d, c0=%d)", \ */
      /* 	      __FUNCTION__,i,sc->marked[i],sc->gen[i],sc->child0[i]); */
      /* } */
    /* REFINE FOLLOWS : */
    nsold=sc->ns;
    /* 
     * refine everything that is marked on the finest level and is
     * not yet refined: (marked>0 and child<0)
     */
    for(j = 0;j < nsold;j++)
      if(sc->marked[j] && (sc->child0[j]<0||sc->childn[j]<0))
	haz_refine_simplex(sc, j, -1);
    /* new mesh */
    sc->level++;
    fprintf(stdout,".%d.",sc->level);//,nsold,ns,nv);
  }
  fprintf(stdout,"\n");
  return;
}
/***********************************************************************/
void sc2mesh(scomplex *sc,mesh_struct *mesh)
{
  /* copy the final grid at position 1*/
  INT ns,n1=sc->n+1,jk=-10,k=-10,j=-10;
  /*  
      store the finest mesh in mesh_struct structure.  sc has all the
      hierarchy, mesh_struct will have only the last mesh.
  */
  ns=0;
  for (j=0;j<sc->ns;j++){
    /*  On the last grid are all simplices that were not refined, so
      these are the ones for which child0 and childn are not set.  */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      mesh->el_flag[ns]=sc->flags[j];
      mesh->el_vol[ns]=sc->vols[j];
      ns++;
    }
  }
  mesh->nv=sc->nv;
  for (j=0;j<mesh->nv;j++){
    mesh->v_flag[j]=sc->bndry[j];
  }
  mesh->el_v->IA = (INT *) calloc(ns+1,sizeof(INT)); 
  mesh->el_v->JA = (INT *) calloc(ns*n1,sizeof(INT));
  mesh->nelm=ns;
  mesh->el_v->IA[0]=0;   
  for (j=0;j<mesh->nelm;j++){
    mesh->el_v->IA[j+1]=mesh->el_v->IA[j]+n1;       
  }
  jk=0;
  for (j=0;j<sc->ns;j++){
    /*  copy el_v map using only the top grid;    */
    if(sc->child0[j]<0 || sc->childn[j]<0){
      for (k=0;k<n1;k++)
	memcpy(mesh->el_v->JA+jk*n1,sc->nbr+j*n1,n1*sizeof(INT));
      jk++;
  }
  fprintf(stdout,"\n%%After %d levels of refinement:\tsimplices=%d ; vertices=%d\n",sc->level,sc->nv,ns); fflush(stdout);
  }
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
