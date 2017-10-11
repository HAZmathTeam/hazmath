/*! \file src/grid/trans_grids.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 1/9/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  adds a line (a material) to every simplex in the grid. 
 *
 */

#include "hazmath.h"

void read_grid_haz0(FILE *gfid,trimesh *mesh) 
{
  /* old way without material (or flag) associated with every element */
  // Loop indices
  INT i,j,k,dummy;  
  INT nelm,nv,dim,nholes;
  i=fscanf(gfid,"%d %d %d %d",&nelm,&nv,&dim,&nholes);
  fprintf(stdout,"\ndim=%i; nelm=%i; nv=%i; nholes=%i\n",dim,nelm,nv,nholes);fflush(stdout);
  // Get number of vertices per element
  INT v_per_elm = dim+1;
  mesh->el_v=malloc(sizeof(iCSRmat));
  // Element-Vertex Map
  mesh->el_v->row=nelm+1;
  mesh->el_v->col=nv;
  mesh->el_v->IA = (INT *)calloc(nelm+1, sizeof(INT));
  mesh->el_v->JA = (INT *)calloc(nelm*v_per_elm, sizeof(INT));
  for(i=0;i<nelm+1;i++) {
    mesh->el_v->IA[i] = v_per_elm*i + 1;
      /* fprintf(stdout,"%d; ", *(mesh->el_v->IA+i));   fflush(stdout); */
  }  
  /* fprintf(stdout,"\n");   fflush(stdout); */
  for (i=0;i<v_per_elm;i++) {
    for (j=0;j<nelm;j++){
      k=v_per_elm*j+i;
      fscanf(gfid,"%d", (mesh->el_v->JA+k));
      /* fprintf(stdout,"%d; ", *(mesh->el_v->JA+k));   fflush(stdout); */
    }
  }
  mesh->el_flag = (INT *) calloc(nelm,sizeof(INT));
  for (j=0;j<nelm;j++){
    mesh->el_flag[j]=1;
  }    
  //new  rveci_(gfid,mesh->el_flag,&nelm);
  // Get next 2-3 lines for coordinate map
  mesh->cv = allocatecoords(nv,dim);
  INT nvdim=nv*dim;
  rvecd_(gfid,mesh->cv->x,&nvdim); // this is the only thing we need to
			     // read. It reads all coordinates by
			     // rows.
  // Get next 1-2 lines for boundary flags
  INT nbv = 0;
  mesh->v_flag = (INT *) calloc(nv,sizeof(INT));
  rveci_(gfid,mesh->v_flag,&nv);
  for(i=0;i<nv;i++) {
    if(mesh->v_flag[i]>0) {
      nbv++;
    }
  }
  // This format only allows for 1 connected region and
  // up to 2 connected boundaries (1 hole).
  INT nconn_reg = 0;
  INT nconn_bdry = 0;
  if(nholes==0) {
    nconn_reg = 1;
    nconn_bdry = 1;
  } else if(nholes==1) {
    nconn_reg = 1;
    nconn_bdry = 2;
  }

  // Update mesh with known quantities
  mesh->dim = dim;
  mesh->nelm = nelm;
  mesh->nv = nv;
  mesh->nbv = nbv;
  mesh->nconn_reg = nconn_reg;
  mesh->nconn_bdry = nconn_bdry;
  //  fprintf(stdout,"\nnv*dimension=%d\n",nvdim);fflush(stdout);
  mesh->v_per_elm = v_per_elm;
  return;
}
void hazw1(FILE *fmesh,				\
	   const trimesh mesh,			\
	   const INT nholes, const int shift)
{
  INT n=mesh.nv,ns=mesh.nelm, dim=mesh.dim,ndl=mesh.dim+1;
  INT *je = (mesh.el_v)->JA, *ib=mesh.v_flag;
  REAL *x = (mesh.cv)->x;
  INT k=-10,j=-10,kndl=-10;
  /* *******************************************
     HAZMATH way of writing mesh file. 
     *******************************************    */
  fprintf(fmesh,"%i %i %i %i\n",ns,n,dim,nholes);
  /* fprintf(stdout,"%i %i %li\n",n,ns,sizeof(ib)/sizeof(INT)); */
  for (j=0;j<ndl;j++) {
    for (k=0;k<ns;k++){
      kndl=ndl*k+j;
      /* shift if needed */
      fprintf(fmesh," %d ", je[kndl]+shift);
    }
    fprintf(fmesh,"\n");
  }  
  for (k=0;k<ns;k++){
    fprintf(fmesh," %d ", mesh.el_flag[k]);
  }
  fprintf(fmesh,"\n");
  INT jn;
  if(dim>1) {
    for(j=0;j<dim;j++){
      jn=j*n;
      for(k=0;k<n;k++){
	fprintf(fmesh," %23.16e",x[jn+k]);
      }
      fprintf(fmesh,"\n");
    }
    for(k=0;k<n;k++){
      fprintf(fmesh," %i ", ib[k]);
    }
    fprintf(fmesh,"\n");
  } else {
    for(k=0;k<n;k++){
      fprintf(fmesh,"%23.16e ",x[0*n+k]);
    }
    fprintf(fmesh,"\n");
    for(k=0;k<n;k++){
      fprintf(fmesh,"%i ", ib[k]);
    }
    fprintf(fmesh,"\n");
  }
  return;
}

int main(int argc,char *argv[])
{
  INT k;
  FILE *fmesh;
  /* for(k=0;k<argc;k++){ */
  /*   fprintf(stdout,"\narg %i: %s\n",k,argv[k]); */
  /* } */
  trimesh *mesh0=malloc(sizeof(trimesh));
  fmesh=HAZ_fopen(argv[1],"r");
  fprintf(stdout,"\n%%Reading %s\n",argv[1]);
  read_grid_haz0(fmesh,mesh0);
  fclose(fmesh);
  /* for(k=0;k<argc;k++){ */
  /*   fprintf(stdout,"\narg %i: %s\n",k,argv[k]); */
  /* } */
  fmesh=HAZ_fopen(argv[2],"w");
  hazw1(fmesh,*mesh0,mesh0->nconn_bdry-1,0);
  fprintf(stdout,"\n%%Output (hazmath) written on:%s\n",argv[2]);
  fclose(fmesh);
  //
  free(mesh0->el_v->IA);
  free(mesh0->el_v->JA);
  free(mesh0->el_flag);
  free(mesh0->el_v);
  free(mesh0->v_flag);
  free(mesh0->cv->x);
  free(mesh0->cv);
  free(mesh0);
  return 0;
}
