/*
 *  fespace.c
 *  
 *
 *  Created by James Adler on 2/17/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

/* Creates and destroys the structres for the finite-element spaces.  
 * Also contains routines for interpolating functions onto the FE space. */

// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"

/****************************************************************************************/
void create_fespace(fespace *FE,trimesh* mesh,INT FEtype)
{
  /* allocates memory and properties of fespace struct */

  fespace FEtmp;
  FEtmp.FEtype = FEtype;
  switch (FEtype)   
    {
    case 1: // Linears - P1
      FEtmp.cdof = mesh->cv;
      FEtmp.ndof = mesh->nv;
      FEtmp.nbdof = mesh->nbv;
      FEtmp.dof_per_elm = mesh->v_per_elm;
      FEtmp.el_dof = mesh->el_v;
      FEtmp.ed_dof = mesh->ed_v;
      FEtmp.f_dof = mesh->f_v;
      FEtmp.dof_bdry = mesh->v_bdry;
      break;
    case 2: // Quadratics - P2
      FEtmp.ndof = mesh->nv + mesh->nedge;
      FEtmp.nbdof = mesh->nbv + mesh->nbedge;
      FEtmp.dof_per_elm = mesh->v_per_elm + mesh->ed_per_elm;
      get_P2(&FEtmp,mesh);
      FEtmp.f_dof = NULL;
      break;
    case -1: // Nedelec Elements 
      FEtmp.cdof = NULL;
      FEtmp.ndof = mesh->nedge;
      FEtmp.nbdof = mesh->nbedge;
      FEtmp.dof_per_elm = mesh->ed_per_elm;
      FEtmp.el_dof = mesh->el_ed;
      FEtmp.ed_dof = NULL;
      FEtmp.f_dof = NULL;
      FEtmp.dof_bdry = mesh->ed_bdry;
      break;
    case -2: // Raviart-Thomas Elements
      FEtmp.cdof = NULL;
      FEtmp.ndof = mesh->nface;
      FEtmp.nbdof = mesh->nbface;
      FEtmp.dof_per_elm = mesh->f_per_elm;
      FEtmp.el_dof = mesh->el_f;
      FEtmp.ed_dof = NULL;
      FEtmp.f_dof = NULL;
      FEtmp.dof_bdry = mesh->f_bdry;
      break;
    default:
      printf("You haven't defined a finite-element space that I understand :-(");
      exit(0);
    }

  *FE = FEtmp;
  
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void free_fespace(fespace* FE)
{
  /* frees memory of arrays of fespace struct */

  free_coords(FE->cdof);
  icsr_free(FE->el_dof);
  icsr_free(FE->ed_dof);
  icsr_free(FE->f_dof);
  if(FE->dof_bdry) free(FE->dof_bdry);
  
  return;
}
/****************************************************************************************/

/***********************************************************************************************/
void get_P2(fespace* FE,trimesh* mesh) 
{
	
  /* Converts Node Coordinate and Element to Node maps to higher orders (just P2 for now) */
	
  INT i,j,k,s,jcntr;  /* Loop Indices */
  INT n1,n2,va,vb,ea,eb,v1,v2,mye,ed,na;

  fespace FEtmp;
  INT ndof = FE->ndof;
  INT dof_per_elm = FE->dof_per_elm;
  INT dim = mesh->dim;
  INT nv = mesh->nv;
  INT nedge = mesh->nedge;
  INT nelm = mesh->nelm;
  INT v_per_elm = mesh->v_per_elm;
  iCSRmat *ed_v = mesh->ed_v;
  iCSRmat *el_v = mesh->el_v;
  iCSRmat *el_ed = mesh->el_ed;
	
  INT* ipv = (INT *) calloc(mesh->v_per_elm,sizeof(INT));
	
  // Get Coordinates
  coordinates cdof;
  allocatecoords(&cdof,ndof,dim);	
  // First go through all vertices.
  for (i=0; i<nv; i++) {
    cdof.x[i] = mesh->cv->x[i];
    cdof.y[i] = mesh->cv->y[i];
    if (dim>2) { cdof.z[i] = mesh->cv->z[i]; }
  }
  // Now, go through and add extra nodes
  // These are simply the midpoints of the edges
  s = nv;
  for (i=0; i<nedge; i++) {
    ed = ed_v->IA[i]-1;
    n1 = ed_v->JA[ed]-1;
    n2 = ed_v->JA[ed+1]-1;
    cdof.x[s] = 0.5*(mesh->cv->x[n1]+mesh->cv->x[n2]);
    cdof.y[s] = 0.5*(mesh->cv->y[n1]+mesh->cv->y[n2]);
    if (dim>2) { cdof.z[s] = 0.5*(mesh->cv->z[n1]+mesh->cv->z[n2]); }
    s++;
  }
	
  // Get Element to Node map
  iCSRmat el_n = icsr_create(nelm,ndof,dof_per_elm*nelm);
  // Rows of Element to Node map
  for(i=0;i<nelm+1;i++) {
    el_n.IA[i] = dof_per_elm*i + 1;
  }	
  // Columns
  for(i=0;i<nelm;i++) {
    va = el_v->IA[i]-1;
    vb = el_v->IA[i+1]-1;
    ea = el_ed->IA[i]-1;
    eb = el_ed->IA[i+1]-1;
    na = el_n.IA[i]-1;
    jcntr = 0;
    for(j=va;j<vb;j++) {
      ipv[jcntr] = el_v->JA[j];
      el_n.JA[na+jcntr] = ipv[jcntr];
      jcntr++;
    }
    for (j=0;j<v_per_elm;j++) {
      for (k=j+1;k<v_per_elm;k++) {
	v1 = ipv[j];
	v2 = ipv[k];
	for(s=ea;s<eb;s++) {
	  mye = el_ed->JA[s];
	  n1 = ed_v->JA[ed_v->IA[mye-1]-1];
	  n2 = ed_v->JA[ed_v->IA[mye-1]];
	  if ((n1==v1 && n2==v2) || (n2==v1 && n1==v2)) {
	    el_n.JA[v_per_elm+na] = mye+nv;
	    na++;
	  }
	}
      }
			
    }
  }
	
  // Fix Boundaries
  INT* dof_bdry = (INT *) calloc(ndof,sizeof(INT));
  // First set of nodes are vertices
  for (i=0; i<nv; i++) {
    dof_bdry[i] = mesh->v_bdry[i];
  }
  // Rest are edges
  for(i=0;i<nedge;i++) {
    dof_bdry[nv+i] = mesh->ed_bdry[i];
  }

  // Get edge to DOF map
  iCSRmat ed_n = icsr_create(nedge,ndof,3*nedge);
  s = nv+1;
  for (i=0; i<nedge; i++) {
    ed = ed_v->IA[i]-1;
    ed_n.IA[i] = ed+1+i;
    n1 = ed_v->JA[ed]-1;
    n2 = ed_v->JA[ed+1]-1;
    ed_n.JA[ed+i] = n1+1;
    ed_n.JA[ed+1+i] = n2+1;
    ed_n.JA[ed+2+i] = s;
    s++;
  }
  ed_n.IA[nedge] = 3*nedge+1;
	
  FEtmp.dof_bdry = dof_bdry;
  FEtmp.el_dof = &el_n;
  FEtmp.ed_dof = &ed_n;
  FEtmp.cdof = &cdof;
  
  *FE = FEtmp;

  if(ipv) free(ipv);
	
  return;
}
/******************************************************************************************/

/****************************************************************************************/
void dump_el_dof(FILE* fid,iCSRmat *el_dof) 
{
  /* Dump the element to DOF map to file for plotting purposes
   *
   * Input:		
   *          el_dof:      Element to DOF map
   *
   * Output:		
   *          coord.dat:    coord(nelm,dof_per_elm)  coordinates of nodes
   *
   */
	
  INT i,j,acol,bcol; /* Loop Indices */
		
  for (i=0; i<el_dof->row; i++) {
    acol = el_dof->IA[i]-1;
    bcol = el_dof->IA[i+1]-1;
    for (j=acol; j<bcol; j++) {
      fprintf(fid,"%d\t",el_dof->JA[j]);
    }
    fprintf(fid,"\n");
  }
   		
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void dump_fespace(fespace *FE) 
{
  /* Dump the mesh data to file for plotting purposes
   *
   * Input:		
   *	      FE:   Finite-element space
   *
   * Output:		
   *   	      el_dof.dat    el_nd(nelm,dof_per_elm) Elements for each node
   *          bdry.dat      bdry(ndof) Indicates whether dof is on boundary
   *
   */
	
  INT i;
  INT totdof = FE->ndof;
  FILE* fid1 = fopen("output/el_dof.dat","w");
  FILE* fid2 = fopen("output/bdry.dat","w");

  // Dump Element to DOF map
  dump_el_dof(fid1,FE->el_dof);

  // Dump boundary data
  for(i=0;i<totdof;i++) {
    fprintf(fid2,"%d\n",FE->dof_bdry[i]);
  }

  fclose(fid1);
  fclose(fid2);
	
  return;
}
/****************************************************************************************/
