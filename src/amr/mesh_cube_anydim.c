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
 * \fn scomplex **mesh_cube_init(const INT dim, const INT ref_type)
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
scomplex **mesh_cube_init(const INT dim, const INT ref_type)
{
  INT i,j,in;
  scomplex **sc_all=NULL;
  cube2simp *c2s=cube2simplex(dim);
  input_grid *g=malloc(sizeof(input_grid));
  // set up initial grid;
  g->title=strdup("Cube in n-d");
  g->fgrid=strdup("mesh.haz");
  g->fvtu=strdup("mesh.vtu");
  //
  g->ncsys=1;
  g->dim=c2s->n;
  INT nvface=c2s->nvface,nvcube=c2s->nvcube;
  g->nf=c2s->nf;
  g->nv=c2s->nvcube;
  g->ne=c2s->ne;
  g->nel=1;
  //
  g->ox=(REAL *)calloc(g->dim*abs(g->ncsys),sizeof(REAL));  
  g->systypes=(INT *)calloc(g->ncsys,sizeof(INT));
  g->syslabels=(INT *)calloc(g->ncsys,sizeof(INT));
  g->csysv=(INT *)calloc(g->nv,sizeof(INT));
  g->labelsv=(INT *)calloc(g->nv,sizeof(INT));
  g->bcodesv=(INT *)calloc(g->nv,sizeof(INT));
  // init these as they are not needed
  memset(g->systypes,0,g->ncsys*sizeof(INT));
  memset(g->syslabels,0,g->ncsys*sizeof(INT));
  memset(g->csysv,0,g->nv*sizeof(INT));
  memset(g->labelsv,0,g->nv*sizeof(INT));
  memset(g->bcodesv,0,g->nv*sizeof(INT));
  //nodes
  g->mnodes=(INT *)calloc(g->nel*(nvcube+1),sizeof(INT));
  for(i=0;i<nvcube;++i){
    g->mnodes[i]=i;
  }
  g->mnodes[nvcube]=1;
  // faces
  g->mfaces=(INT *)calloc(abs(g->nf)*(nvface+1),sizeof(INT));
  for(i=0;i<g->nf;++i){
    memcpy((g->mfaces+i*(nvface+1)),(c2s->faces+i*nvface),c2s->nvface*sizeof(INT));
    g->mfaces[nvface+i*(nvface+1)]=1;
  }
  // vertex coords:
  g->xv=(REAL *)calloc(g->dim*g->nv,sizeof(REAL));
  for(i=0;i<g->nv;++i){
    in=i*g->dim;
    for(j=0;j<g->dim;j++){
      g->xv[in+j]=(REAL )c2s->bits[in+j];
    }
  }
  // edges
  g->seg=(INT *)calloc(3*g->ne,sizeof(INT));
  g->xe=(REAL *)calloc(g->dim*g->ne,sizeof(REAL));
  INT j1,j2;
  for(i=0;i<g->ne;++i){
    j1=c2s->edges[2*i];
    j2=c2s->edges[2*i+1];
    g->seg[3*i]   =  j1;
    g->seg[3*i+1] =  j2;
    g->seg[3*i+2] =  1;
    for(j=0;j<dim;++j)
      g->xe[j+i*dim]=0.5*(g->xv[j+j1*dim]+g->xv[j+j2*dim]);
  }
  cube2simp_free(c2s);
  if(g->dim > 3)
    g->ref_type=0;
  else
    g->ref_type=-2;
  g->nref=0;   /* number of refinements (AMR)*/    
  g->mark_type=0;   /* AMR marking type (0)uniform; nonzero: user defined.. */
  g->err_stop=1e-20;   /* stop tolerance for AMR */
  g->num_refine_points=0;
  g->print_level=0;
  g->data_refine_points=NULL;
  //  input_grid_print(g);
  sc_all=generate_initial_grid(g);
  input_grid_free(g);
  //////////////////////////////////////////////////////
  if(sc_all[0]->n > 3)
    sc_all[0]->ref_type=0;
  else
    sc_all[0]->ref_type=ref_type;
  //////////////////////////////////////////////////////
  //  if(dim <4) vtkw(g->fvtu,sc,0,1.);
  //  haz_scomplex_print(sc_all[0],0,"XXX"); fflush(stdout);
  return sc_all;
}
