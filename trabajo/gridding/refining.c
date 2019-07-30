/*! \file src/amr/refining.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 * \note: routines used to refine a simplicial grid grid ref_level times.

 */
#include "hazmath.h"
/***********************************************************************/
void sc2mesh(scomplex *sc,mesh_struct *mesh)
{
  /* converts sc to a mesh struct and returns an ivector which has the
     correspondence between simplices*/
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
