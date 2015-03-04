/*
 *  grid.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

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

/**************************************************************************************************************************************************/
/********* Read in arbritray grid and create FEM structure for the grid ***************************************************************************/
/**************************************************************************************************************************************************/
void creategrid(FILE *gfid,INT dim,INT nholes,trimesh* mesh) 
{
  /* Reads in gridfile of the following form:
   *
   * First line:			nelms nnodes nboundaryedges nboundaryedges
   * Next four lines:		       	Each element will have 3 vertices in 2D and 4 in 3D
   *			line 1:	       	First node of all elements
   *			line 2:		Second node of all elements
   *			line 3:	        Third node of all elements
   *			line 4:		Dummy line (ignore) (or 4th node for 3D)
   *					Columns are node-element map not sparse
   * Next two (or 3 in 3D) lines:	x,y,z coordinates
   *			line 1:		x coordinates of all nodes
   *			line 2:		y coordinates of all nodes
   *			(line 3:	z coordinates of all nodes (only in 3D))
   * Next two lines in 2D:		List of nodes that land on boundary
   *			line 1:		Corresponds to first part of edge on bondary
   *			line 2:		Corresponds to second part of edge on boundary
   *					Columns of this are edge_node map of boundary
   *
   * OR Next line in 3D:                List of nodes and whether they land on boundary (binary)
   *
   * Rest of lines are dummy and not used
   *
   * INPUT: gfid    Grid File ID
   *        dim     Dimension of Problem 
   *        nholes  Number of Holes in the domain (usually 0)
   *        TODO: change mesh files so that dim and nholes is automatically read in
   */
	
  INT i,j,k; /* Loop indices */

  // Define dummy mesh to load data
  trimesh mesh_temp;

  // Get Number of elements, nodes and boundary edges first
  INT* line1 = calloc(4,sizeof(INT));
  INT lenhead = 4;
  rveci_(gfid,line1,&lenhead);
  INT nelm = line1[0];
  INT nv = line1[1];
  INT nbedge = line1[2];
  free(line1);

  // Get number of vertices, edges, and faces per element
  INT v_per_elm = dim+1;
  INT ed_per_elm = 3*(dim-1);
  INT f_per_elm = v_per_elm;
			
  // Allocate arrays to read in other data such as coordinate information
  INT* element_vertex = (INT *) calloc(nelm*v_per_elm,sizeof(INT));
  coordinates cv;
  allocatecoords(&cv,nv,dim);
  INT* bdry_v = (INT *) calloc(nbedge*2,sizeof(INT));
		
  // Get next 3-4 lines Element-Vertex Map
  INT* line2 = (INT *) calloc(nelm,sizeof(INT));
  for (i=0; i<v_per_elm; i++) {
    rveci_(gfid,line2,&nelm);
    for (j=0; j<nelm; j++) {
      element_vertex[j*v_per_elm+i] = line2[j];
    }
  }
  // Next line is dummy in 2D
  if(dim==2) {
    rveci_(gfid,line2,&nelm);
  }
  free(line2);
	
  // Get next 2-3 lines for coordinate map
  rvecd_(gfid,cv.x,&nv);
  rvecd_(gfid,cv.y,&nv);
  if(dim==3)
    rvecd_(gfid,cv.z,&nv);
	
  // Get next 1-2 lines for boundary nodes
  INT nbv = 0;
  INT* v_bdry = (INT *) calloc(nv,sizeof(INT));
  if(dim==2) {
    INT* line3 = (INT *) calloc(nbedge,sizeof(INT));
    INT* line4 = (INT *) calloc(nbedge,sizeof(INT));
    rveci_(gfid,line3,&nbedge);
    rveci_(gfid,line4,&nbedge);
    for (i=0; i<nbedge; i++) {
      bdry_v[i*2+0] = line3[i];
      bdry_v[i*2+1] = line4[i];
    }
    free(line3);
    free(line4);
  } else if(dim==3) {
    rveci_(gfid,v_bdry,&nv);
    for(i=0;i<nv;i++)
      if(v_bdry[i])
	nbv++;
  }
    
    
  // if(bdry_v) free(bdry_v);  // I comment this out since otherwise it give seg fault later -- Xiaozhe
    

  printf("\nConverting Grid Maps to CSR and Computing Data Structures:\n ");
  /* Element Vertex Map */
  iCSRmat el_v = convert_elmnode(element_vertex,nelm,nv,v_per_elm);
  if(element_vertex) free(element_vertex);
    
  /* Edge to Vertex Map */
  INT nedge = 0;
  if(dim==2) {  
    nedge = nelm+nv-(dim-1);
  } else if(dim==3) {
    get_nedge(&nedge,el_v); 
  }


  iCSRmat ed_v = get_edge_v(nedge,el_v);

    
  /* Get Boundary Edges and Vertices */
  INT* ed_bdry = (INT *) calloc(nedge,sizeof(INT));
  if (dim==2) {
      
    isboundary_ed(ed_v,nedge,nbedge,bdry_v,ed_bdry);
    isboundary_v(nv,bdry_v,v_bdry,nbedge,&nbv);


  } else if(dim==3) {
    isboundary_ed3D(ed_v,nedge,cv,&nbedge,v_bdry,ed_bdry);
  }

  if(bdry_v) free(bdry_v);
	
  /* Element to Edge Map */
  iCSRmat el_ed = get_el_ed(el_v,ed_v);
    


  /* Get Edge Stats such as edge lengths, midpoints, and tangent vectors */
  REAL* ed_len = (REAL *) calloc(nedge,sizeof(REAL));
  REAL* ed_tau = (REAL *) calloc(nedge*dim,sizeof(REAL));
  REAL* ed_mid = (REAL *) calloc(nedge*dim,sizeof(REAL));
    

  edge_stats_all(ed_len,ed_tau,ed_mid,cv,ed_v,dim);
    

  /* Figure out Number of Faces */
  INT nface = 0;
  INT euler = -10; // Euler Number should be 1 for simplices
  if(dim==2) {
    nface = nedge;
    euler = nv - nedge + nelm;
  } else if (dim==3) {
    nface = 1 + nedge-nv+nelm;
    nholes=0;
    nface = nface + nholes; // add number of holes!
    euler = nv - nedge + nface - nelm;
  } else {
    baddimension();
  }
  if(euler!=1+nholes) {
    printf("Your simplices are all messed up.  Euler Characteristic doesn't equal 1+nholes!\teuler=%d\tnholes=%d",euler,nholes);
    exit(0);
  }


    
  // In order to get some of the face maps, we need to know how the faces are ordered on each element
  // This is done by the same ordering as the nodes, connecting the opposite face with each node.
  INT* fel_order= (INT *)calloc(f_per_elm*dim,sizeof(INT));
  get_face_ordering(v_per_elm,dim,f_per_elm,fel_order);

  // Next get the element to face map, the face to vertex map, and the face boundary information.
  iCSRmat el_f = icsr_create(nelm,nface,f_per_elm*nelm);
  iCSRmat f_v = icsr_create(nface,nv,nface*dim);
  INT* f_bdry = (INT *) calloc(nface,sizeof(INT));
  INT nbface;
  get_face_maps(el_v,v_per_elm,nface,dim,f_per_elm,&el_f,f_bdry,&nbface,&f_v,fel_order);
  // In case v_bdry has different types of boundaries, match the face boundary to the same:
  for(i=0;i<nface;i++) {
    if(f_bdry[i]==1) {
      j = f_v.IA[i]-1;
      k = f_v.JA[j]-1;
      f_bdry[i] = v_bdry[k];
    }
  }
    
  // Next get the face data, such as areas, midpoints, and normal vectors 
  REAL* f_area = (REAL *) calloc(nface,sizeof(REAL));
  REAL* f_mid = (REAL *) calloc(nface*dim,sizeof(REAL));
  REAL* f_norm = (REAL *) calloc(nface*dim,sizeof(REAL));
  // Need to update mesh_temp for a few things first
  mesh_temp.dim = dim;
  mesh_temp.cv = &cv;
  mesh_temp.f_per_elm = f_per_elm;
  mesh_temp.el_v = &el_v;
  mesh_temp.el_f = &el_f;
  mesh_temp.v_per_elm = v_per_elm;
  mesh_temp.f_v = &f_v;
    
    printf("hello\n");
    
  // TODO check if passing the mesh was done correctly
  face_stats(f_area,f_mid,f_norm,fel_order,mesh_temp);
    
    printf("hello1\n");


  // Finally get volumes/areas of elements and the midpoint (barycentric)
  REAL* el_mid = (REAL *) calloc(nelm*dim,sizeof(REAL));
  REAL* el_vol = (REAL *) calloc(nelm,sizeof(REAL));
  get_el_vol(el_vol,el_v,cv,dim,v_per_elm);
  get_el_mid(el_mid,el_v,cv,dim);
  
  // Assign components to the mesh
  mesh_temp.nelm = nelm;
  mesh_temp.nv = nv;
  mesh_temp.nedge = nedge;
  mesh_temp.ed_per_elm = ed_per_elm;
  mesh_temp.nface = nface;
  mesh_temp.nholes = nholes;
  mesh_temp.nbv = nbv;
  mesh_temp.nbedge = nbedge;
  mesh_temp.nbface = nbface;
  mesh_temp.el_ed = &el_ed;
  mesh_temp.ed_v = &ed_v;
  mesh_temp.el_vol = el_vol;
  mesh_temp.el_mid = el_mid;
  mesh_temp.ed_len = ed_len;
  mesh_temp.ed_tau = ed_tau;
  mesh_temp.ed_mid = ed_mid;
  mesh_temp.f_area = f_area;
  mesh_temp.f_mid = f_mid;
  mesh_temp.f_norm = f_norm;
  
  mesh_temp.v_bdry = v_bdry;
  mesh_temp.ed_bdry = ed_bdry;
  mesh_temp.f_bdry = f_bdry;
  
  *mesh = mesh_temp;
  return;
}
/**************************************************************************************************************************************************/

/**************************************************************************************************************************************************/
void initialize_mesh(trimesh* mesh) 
{
  /* Initializes all components of the structure trimesh. */
	
  mesh->dim = -666;
  mesh->nelm = -666;
  mesh->cv = NULL;
  mesh->f_per_elm = -666;
  mesh->el_v = NULL;
  mesh->el_f = NULL;
  mesh->v_per_elm = -666;
  mesh->f_v = NULL;
  mesh->nv = -666;
  mesh->nedge = -666;
  mesh->ed_per_elm = -666;
  mesh->nface = -666;
  mesh->nholes = -666;
  mesh->nbv = -666;
  mesh->nbedge = -666;
  mesh->nbface = -666;
  mesh->el_ed = NULL;
  mesh->ed_v = NULL;
  mesh->el_vol = NULL;
  mesh->el_mid = NULL;
  mesh->ed_len = NULL;
  mesh->ed_tau = NULL;
  mesh->ed_mid = NULL;
  mesh->f_area = NULL;
  mesh->f_mid = NULL;
  mesh->f_norm = NULL;
  mesh->v_bdry = NULL;
  mesh->ed_bdry = NULL;
  mesh->f_bdry = NULL;
  return;
}
/**************************************************************************************************************************************************/

/***********************************************************************************************/
iCSRmat convert_elmnode(INT *element_node,INT nelm,INT nv,INT nve) 
{
	
  /* Convert the input element to vertex map into sparse matrix form 
   *
   * Input: nelm:	            Number of elements
   *	    nv:		    Number of vertices
   *        nve:	            Number of vertices per element
   *	    element_node(nelm,nve): Each row is an element, each column is the corresponding vertices for that element 
   *										
   *
   * Output:	el_v	    element to vertex map in CSR format.  Rows are all elements, columns are all vertices.  Entries are 1 if there's a connection.
   */
	
  INT i,j,k; /* Loop Indices */

  iCSRmat el_v;
	
  /* Exactly nve vertices per row */
  if ( nelm > 0 ) {
    el_v.IA = (INT *)calloc(nelm+1, sizeof(INT));
  }
  else {
    el_v.IA = NULL;
  }
  for(i=0;i<nelm+1;i++) {
    el_v.IA[i] = nve*i + 1;
  }
	
  /* Columns are extracted directly from input array */
  if ( nv > 0 ) {
    el_v.JA = (INT *)calloc(nelm*nve, sizeof(INT));
  }
  else {
    el_v.JA = NULL;
  }
  k=0;
  for(i=0;i<nelm;i++) {
    for(j=0;j<nve;j++) {
      el_v.JA[k] = element_node[i*nve+j];
      k=k+1;
    }
  }
  
  el_v.val = NULL;
  el_v.row = nelm; el_v.col = nv; el_v.nnz = nelm*nve;
	
  return el_v;
}
/***********************************************************************************************/

/***********************************************************************************************/
void get_nedge(INT* nedge, iCSRmat el_v) 
{
	
  /* Gets the Number of Edges (NEEDED for 3D)
   *
   *	Input:  iel_n,j_eln:		Element to Node Map
   *			n					Total number of Nodes
   *			nelm				Total number of Elements
   *			element_order:		Number of nodes per element
   *
   *	Output: nedge				Total number of Edges
   *
   */
	
  //INT i,j,col_b,col_e,icntr; /* Loop indices and counters */
	
  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans_1(&el_v,&v_el);
	
  /* Create Vertex to Vertex Map by v_el*el_v */
  iCSRmat v_v;
  icsr_mxm_1(&v_el,&el_v,&v_v);
	
  // Now go through upper triangular (above diagonal) portion of vertex-vertex and count non-zeros
  // Each of these vertex-vertex connections are edges.  Upper so we don't count i->j and j->i as two
  // separate edges, and no diagonal since i->i is not an edge.
  /* INT* iv_v = v_v.IA; */
  /* INT* jv_v = v_v.JA; */
  /* icntr = 0; */
  /* for(i=0;i<nv-1;i++) { */
  /*   col_b = in_n[i]; */
  /*   col_e = in_n[i+1]-1; */
  /*   for(j=col_b;j<=col_e;j++) { */
  /*     if((jn_n[j-1]-1)>i) { */
  /* 	icntr++; */
  /*     } */
  /*   } */
  /* }	 */
  /* Free Node Node */
  /* free(in_n); */
  /* free(jn_n); */
  /* *nedge = icntr; */

  *nedge = (INT) (v_v.nnz - v_v.row)/2;

  /* Free matrices */
  icsr_free(&v_v);
  icsr_free(&v_el);
  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
iCSRmat get_edge_v(INT nedge,iCSRmat el_v) 
{
	
  /* Gets the Edge to Vertex mapping in CSR Fromat (Should be dim independent)
   * This is used to get Element to Edge map, but also for determining
   * things such as edge_length and the tangent vectors of edges 
   *
   *	Input: el_v:	  Element to Vertex Map
   *	       nv:	  Total number of Vertices
   *	       v_per_elm: Number of vertices per element
   *
   */

  INT nv = el_v.col;
  iCSRmat ed_v;
  if ( nedge > 0 ) {
    ed_v.IA = (INT *)calloc(nedge+1, sizeof(INT));
  }
  else {
    ed_v.IA = NULL;
  }
  if ( nv > 0 ) {
    ed_v.JA = (INT *)calloc(2*nedge, sizeof(INT));
  }
  else {
    ed_v.JA = NULL;
  }
	
  INT i,j,col_b,col_e,icntr,jcntr; /* Loop indices and counters */
	
  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans_1(&el_v,&v_el);
	
  /* Create Vertex to Vertex Map by v_el*el_v */
  iCSRmat v_v;
  icsr_mxm_symb_1(&v_el,&el_v,&v_v);
  INT* iv_v = v_v.IA;
  INT* jv_v = v_v.JA;
	  
  // Now go through upper triangular (above diagonal) portion of node-node and count non-zeros
  // Each of these node-node connections are edges.  Upper so we don't count i->j and j->i as two
  // separate edges, and no diagonal since i->i is not an edge.
  jcntr = 0;
  icntr = 0;
  for(i=0;i<nv-1;i++) {
    col_b = iv_v[i];
    col_e = iv_v[i+1]-1;
    for(j=col_b;j<=col_e;j++) {
      if((jv_v[j-1]-1)>i) {
	ed_v.IA[icntr] = jcntr+1;
	ed_v.JA[jcntr] = jv_v[j-1];
	ed_v.JA[jcntr+1] = i+1;
	jcntr=jcntr+2;
	icntr++;
      }
    }
  }
  ed_v.IA[icntr] = jcntr+1;		
	
  /* Free Node Node */
    //free(iv_v);
    //free(jv_v);
  icsr_free(&v_v);
  icsr_free(&v_el);

  ed_v.val=NULL;
  ed_v.row = nedge;
  ed_v.col = nv;
  ed_v.nnz = 2*nedge;
    
  return ed_v;
}
/***********************************************************************************************/

/***********************************************************************************************/
void isboundary_v(INT nv,INT *bdry_v,INT *v_bdry,INT nbedge,INT *nbv) 
{
	
  /* Indicates whether a vertex is on boundary or not
   *
   *	Input:  nv:	Number of vertices
   *		bdry_v	List of boundaries and corresponding vertices
   *		nbedge	Total Number of boundary edges
   *	Output: 
   *		v_bdry	List of vertices and whether or not they're a boundary
   *		nbv	Total number of boundary vertices
   */
	
  INT i,m,cnt;
    
	
  for (i=0; i<nv; i++) {
    v_bdry[i] = 0;
  }
    

	
  for (i=0; i<2*nbedge; i++) {
    m = bdry_v[i];
    v_bdry[m-1] = 1;
  }
    

	
  cnt = 0;
  for (i=0; i<nv; i++) {
    if(v_bdry[i]==1) {
      cnt++;
    }
  }
  *nbv = cnt;
    

  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
void isboundary_ed(iCSRmat ed_v,INT nedge,INT nbedge,INT *bdry_v,INT *ed_bdry) 
{
	
  /* Indicates whether an edge is on boundary or not
   *
   *	Input:  ed_v		Edge to Vertex Map
   *		nedge		Total number of Edges
   *		nbedge:		Number of boundary edges
   *		bdry_v		List of boundaries and corresponding vertices
   *	Output: 
   *		ed_bdry		List of Edges and whether or not they're a boundary
   */
	
  INT i,j,col_b,col_e,k,jcntr; /* Loop indices and counters */
  INT* n = calloc(2,sizeof(INT));
  INT* m = calloc(2,sizeof(INT));
	
  // Zero out ed_bdry
  for (i=0; i<nedge; i++) {
    ed_bdry[i] = 0;
  }
	
  for (i=0; i<nedge; i++) {
    col_b = ed_v.IA[i];
    col_e = ed_v.IA[i+1]-1;
    jcntr=0;
    for (j=col_b; j<=col_e; j++) {
      n[jcntr] = ed_v.JA[j-1];
      jcntr++;
    }
    for (k=0; k<nbedge; k++) {
      m[0] = bdry_v[k*2+0];
      m[1] = bdry_v[k*2+1];
      if (((n[0]==m[0]) && (n[1]==m[1])) || ((n[1]==m[0]) && (n[0]==m[1]))) {
	ed_bdry[i] = 1;
      }
    }
  }
  if(n) free(n);
  if(m) free(m);
	
  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
void isboundary_ed3D(iCSRmat ed_v,INT nedge,coordinates cv,INT *nbedge,INT *v_bdry,INT *ed_bdry) 
{
	
  /* Counts the number of boundary edges and indicates whether an edge is a boundary
   *
   *	Input:  ed_v:		Edge to Vertex Map
   *		nedge	    	Total number of edges
   *		cv		Vertex Coordinates
   *		v_bdry		List of vertices and whether they are on boundary
   *	Output: 
   *		ed_bdry		List of edges and whether they are on boundary
   *		nbedge	        Number of boundary edges
   */
	
  INT i,col_b,col_e,jcntr; /* Loop indices and counters */
  INT n;
  INT m;
	
  jcntr=0;
  // For every edge get nodes
  for (i=0; i<nedge; i++) {
    col_b = ed_v.IA[i];
    col_e = ed_v.IA[i+1]-1;
    n = ed_v.JA[col_b-1]-1;
    m = ed_v.JA[col_e-1]-1;
    //Check if two nodes are on boundary
    if (v_bdry[n]!=0 && v_bdry[m]!=0) {
      if(cv.x[n]==cv.x[m] || cv.y[n]==cv.y[m] || cv.z[n]==cv.z[m]) {
	ed_bdry[i] = v_bdry[n];
	jcntr++;
      } else {
	ed_bdry[i] = 0;
      }
    }
  }
  *nbedge = jcntr;
  return;
}
/***********************************************************************************************/

/***********************************************************************************************/
iCSRmat get_el_ed(iCSRmat el_v,iCSRmat ed_v) 
{
  /* Gets the Element to Edge mapping in CSR Fromat (Should be dim independent)
   *
   *	Input: el_v:	  Element to Vertex Map
   *	       ed_v:	  Edge to Vertex Map
   *	       v_per_elm: Number of vertices per element
   *
   */

  iCSRmat el_ed;
  
  // Get transpose of edge to vertex
  iCSRmat v_ed;
  icsr_trans_1(&ed_v,&v_ed);
  icsr_mxm_symb_max_1(&el_v,&v_ed,&el_ed,2);
  printf("nedge=%d\tnv=%d\n",ed_v.row,ed_v.col);
  pveci_(ed_v.IA,(ed_v.row)+1);
  pveci_(ed_v.JA,ed_v.nnz);
  pveci_(v_ed.IA,v_ed.row+1);
  pveci_(v_ed.JA,v_ed.nnz);
  pveci_(el_ed.IA,el_ed.row+1);
  pveci_(el_ed.JA,el_ed.nnz);
  icsr_free(&v_ed);

    
  return el_ed;
}
/***********************************************************************************************/

/****************************************************************************************/
void allocatecoords(coordinates *A,INT ndof,INT mydim)
{
  /* allocates memory and properties of coordinates struct */

  coordinates Atmp;
  Atmp.x = calloc(ndof,sizeof(REAL));
  Atmp.y = calloc(ndof,sizeof(REAL));
  if (mydim==3) {
    Atmp.z = calloc(ndof,sizeof(REAL));
  } else {
    Atmp.z = NULL;
  }
  Atmp.n = ndof;

  *A = Atmp;
  
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void free_coords(coordinates* A)
{
  /* fres memory of arrays of Incident matrix struct */

  if(A->x) free(A->x);
  if(A->y) free(A->y);
  if(A->z) free(A->z);
  
  return;
}
/****************************************************************************************/

/*********************************************************************************************************/
void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates cv,iCSRmat ed_v,INT dim) 
{

  /***************************************************************************
   * Get length, tangent vector (tau), and midpoint of every edge **********
   *
   *
   *    Input:   cv        Coordinates of Vertices
   *             ed_v      Edge to Vertex Map
   *             dim       Dimension of Problem
   *
   *    Output:  ed_len    Lengths of Edges
   *             ed_tau    Tangent Vector of Edges ordered (tx1,ty1,tz1,tx2,ty2,tz2,...)
   *             ed_mid    Midpoint of Edges ordered (mx1,my1,mz1,mx2,my2,mz2,...) 
   */
	
  INT i,jcnt,j,j_a,j_b; /* loop index */
  INT ip[2];
  REAL x[2],y[2],z[2];
  ip[0]=0; ip[1]=0; x[0]=0.0; x[1]=0.0; y[0]=0.0; y[1]=0.0; z[0]=0.0; z[1]=0.0;
	
  INT ihi1,ilo1; /* Orders nodes from largest to smallest for correct orientation */
  INT nedge = ed_v.row;

  for(i=0;i<nedge;i++) {

    /* Find Nodes in given Edge */
    j_a = ed_v.IA[i]-1;
    j_b = ed_v.IA[i+1]-1;
    jcnt = 0;
    for (j=j_a; j<j_b;j++) {
      ip[jcnt] = ed_v.JA[j];
      x[jcnt] = cv.x[ip[jcnt]-1];
      y[jcnt] = cv.y[ip[jcnt]-1];
      if (dim==3) {
	z[jcnt] = cv.z[ip[jcnt]-1];
      } else {
	z[jcnt] = 0;
      }
      jcnt++;
    }
	
    /* Compute Length */
    ed_len[i] = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + (z[1]-z[0])*(z[1]-z[0]));
	
    /* Get correct orientation for tangent vectors */
    if(MAX(ip[0],ip[1])==ip[0]) {
      ihi1 = 0;
      ilo1 = 1;
    } else {
      ihi1 = 1;
      ilo1 = 0;
    }
	
    /* Compute Tangent Vectors */
    ed_tau[i*dim] = (x[ihi1]-x[ilo1])/ed_len[i];
    ed_tau[i*dim+1] = (y[ihi1]-y[ilo1])/ed_len[i];
    if (dim==3) {
      ed_tau[i*dim+2] = (z[ihi1]-z[ilo1])/ed_len[i];		
    }
	
    /* Get midpoint of edges for later use */
    ed_mid[i*dim] = 0.5*(x[0]+x[1]);
    ed_mid[i*dim+1] = 0.5*(y[0]+y[1]);
    if (dim==3) {
      ed_mid[i*dim+2] = 0.5*(z[0]+z[1]);		
    }
  }
  return;
}
/*********************************************************************************************************/

/***********************************************************************************************/
void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order)
{
	
  /* Gets the face ordering on element (Should be dim independent and consistent for all elements)
   *
   *               3
   *             /   \
   *            2     1
   *           /       \
   *          1 ---3----2
   *
   *	Input:  
   *		el_order:	Number of nodes per element
   *            dim             Dimension
   *            f_order         Number of Faces per element 3 in 2D, 4 in 3D
   *
   *	Output:
   *            fel_order Indicates order of faces on an element
   *
   */
	
  INT i,j,jcntr,col_b; /* Loop Indices */

  // First order the faces per element.  Start with node 1 and indicate opposite face as face 1.  Then loop over nodes
  for(j=0;j<f_order;j++) {
    for(i=0;i<dim;i++) {
      jcntr = j+i+2;
      col_b = j*dim+i;
      if(jcntr>el_order) {
	fel_order[col_b] = jcntr-el_order;
      } else {
	fel_order[col_b] = jcntr;
      }
    }
  }
  return;
}
/***********************************************************************************************/


/***********************************************************************************************/
void get_face_maps(iCSRmat el_v,INT el_order,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,INT *fel_order)
{
	
  /* Gets the Element to Face, Face to Node, Face to Boundary mapping in CSR Fromat as well as face ordering on element (Should be dim independent)
   *
   *               3
   *             /   \
   *            2     1
   *           /       \
   *          1 ---3----2
   *
   *	Input:  el_v:		Element to Vertex Map
   *		el_order:	Number of nodes per element
   *            nface           Total Number of Faces
   *            dim             Dimension
   *            f_order         Number of Faces per element 3 in 2D, 4 in 3D
   *            fel_order Indicates order of faces on an element
   *
   *	Output:
   *		el_f	Element to Face Map where value represents face number in element
   *            f_bdry  Indicates whether given face is on boundary or not
   *            nbf     Number of boundary faces
   *            f_v     Face to Node Map
   *
   */
	
  INT i,j,k,m,p,jk,col_b,icntr,jcntr,kcntr; /* Loop Indices */
  INT ncol1,ncol2,ncol3,ncol4,ncol5,ncol6;  /* Node indices */
  INT el1=-1,el2=-1,el3=-1,el=-1; /* Element indices */
  INT iflag = 0;  /* Marks if we found a face */
  INT nbf = 0; /* hold number of boundary faces */
  INT nd[dim];
  INT f_num=-1;  /* Current face number of element */
  INT nelm = el_v.row;
  INT nv = el_v.col;

  // We will build face to element map first and then transpose it
  iCSRmat f_el = icsr_create (nface,nelm,nelm*f_order);
  
  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans_1(&el_v,&v_el);

  // Get interior faces first
  // Loop over All elements and each face on the element
  // Then populate face to element and face to vertex map
  // Also determines if a face is on the boundary.
  icntr=0;
  jcntr=0;
  kcntr=0;
  iflag=0;
  for(i=0;i<nelm;i++) {
    // Get first node of element
    col_b = el_v.IA[i]-1;
    // Now loop over all faces of element
    for(j=0;j<f_order;j++) {
      // Get appropriate nodes on face
      for(k=0;k<dim;k++) {
	jk = col_b+fel_order[j*dim+k]-1;
	nd[k] = el_v.JA[jk];
      }
      // Next Loop over all elements that touch the first node
      ncol1 = v_el.IA[nd[0]-1]-1;
      ncol2 = v_el.IA[nd[0]]-1;
      for(k=ncol1;k<ncol2;k++) {
	el1 = v_el.JA[k]-1;
	if(el1!=i) { // If not the same element we started in then check other nodes elements
	  ncol3 = v_el.IA[nd[1]-1]-1;
	  ncol4 = v_el.IA[nd[1]]-1;
	  for(m=ncol3;m<ncol4;m++) {
	    el2 = v_el.JA[m]-1;
	    if(el2==el1) { // Our two nodes share the element!  In 2D this is the face! In 3D, we check third guy
	      if(dim==2) {
		// Mark the element that shares the face as well as the nodes that are shared
		iflag = 1;
		el = el1;
	      } else if(dim==3) {
		ncol5 = v_el.IA[nd[2]-1]-1;
		ncol6 = v_el.IA[nd[2]]-1;
		for(p=ncol5;p<ncol6;p++) {
		  el3 = v_el.JA[p]-1;
		  if(el3==el2) {
		    // Mark the element that shares the face as well as the nodes that are shared
		    iflag = 1;
		    el = el1;
		  }
		}
	      } else {
		baddimension();
	      }
	    }
	  }
	}
      }
      if(iflag && el>i) {
	f_el.IA[icntr] = jcntr+1;
	f_el.JA[jcntr] = i+1;
	f_el.val[jcntr] = j+1;
	f_el.JA[jcntr+1] = el+1;
	// Find face number for other element
	find_facenumber(el_v,fel_order,el+1,nd,dim,&f_num);
	f_el.val[jcntr+1] = f_num;
	f_bdry[icntr] = 0;
	f_v->IA[icntr] = kcntr+1;
	f_v->JA[kcntr] = nd[0];
	f_v->JA[kcntr+1] = nd[1];
	kcntr+=2;
	if(dim==3) {
	  f_v->JA[kcntr] = nd[2];
	  kcntr++;
	}
	icntr++;
	jcntr+=2;
      } else if(!iflag){  // this must be a boundary face!
	nbf++;
	f_el.IA[icntr] = jcntr+1;
	f_el.JA[jcntr] = i+1;
	f_el.val[jcntr] = j+1;
	f_bdry[icntr] = 1;
	f_v->IA[icntr] = kcntr+1;
	f_v->JA[kcntr] = nd[0];
	f_v->JA[kcntr+1] = nd[1];
	kcntr+=2;
	if(dim==3) {
	  f_v->JA[kcntr] = nd[2];
	  kcntr++;
	}
	icntr++;
	jcntr++;
      }
      iflag=0;
    }
  }

  f_el.IA[icntr] = jcntr+1;
  f_v->IA[icntr] = kcntr+1;
  
  /* Get Transpose of f_el -> el_f */
  icsr_trans_1(&f_el,el_f);

  /* Transpose face_node and back again to order nodes in increasing order (from global structure) */
  iCSRmat v_f = icsr_create(nv,nface,nface*dim);
  icsr_trans_1(f_v,&v_f);
  icsr_trans_1(&v_f,f_v);
  
  *nbface = nbf;

  icsr_free(&v_el);
  icsr_free(&f_el);
  icsr_free(&v_f);

  return;
}
/***********************************************************************************************/

/****************************************************************************************/
void find_facenumber(iCSRmat el_v,INT* fel_order,INT elm,INT* nd,INT dim,INT *f_num)          
{

  /* Find the face number of the given element, using the element it shares the face with
   *
   * Input:
   *        el_v              Element to Vertex Map (all elements and vertices)
   *        fel_order         Ordering of faces on element with corresponding nodes
   *        elm               Current Elment we consider
   *        nd                Nodes of current face
   *        dim               Dimension
   *
   * Output:
   *       f_num             Face Number for that element
   *
   */

  INT cola,colb,mark,icnt,i,j,mynd;
  INT fntmp=-1;

  cola = el_v.IA[elm-1]-1;
  colb = el_v.IA[elm]-1;
  mark = 0;
  icnt=0;
  for(i=cola;i<colb;i++) {
    icnt++;
    mynd = el_v.JA[i];
    for(j=0;j<dim;j++) {
      if(mynd!=nd[j]) {
	mark++;
      }
    }
    if(mark==dim) {
      fntmp = icnt;
    }
    mark=0;
  }

  *f_num = fntmp;

  return;
}
/****************************************************************************************/

/*********************************************************************************************************/
void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,INT* fel_order,trimesh mesh) 
{
  /***************************************************************************
   * Get area, normal vector for all faces **********
   *
   *
   *    Input:   
   *             fel_order            Ordering of Faces on a given element
   *             mesh                 Information needed for mesh
   *
   *    Output:  f_area               Area of each face (length in 2D)
   *             f_norm(nface,dim)    Normal vectors of each face
   *             f_mid(nface,dim)     Midpoints of face
   *            
   */
    
    printf("hello-in\n");
	
  INT i,jcnt,j,j_a,j_b; /* loop index */
  INT nface = mesh.el_f->col;
  INT nelm = mesh.el_f->row;
  INT dim = mesh.dim;
  INT el_order = mesh.v_per_elm;
  INT f_order = mesh.f_per_elm;

  coordinates *cv = mesh.cv;

  iCSRmat *el_f = mesh.el_f;
  iCSRmat *f_v = mesh.f_v; 
  iCSRmat *el_v = mesh.el_v;

  // Face Node Stuff
  INT* ipf = calloc(dim,sizeof(INT));
  REAL* xf = calloc(dim,sizeof(REAL));
  REAL* yf = calloc(dim,sizeof(REAL));
  REAL* zf;
  if(dim==3) {
    zf = calloc(dim,sizeof(REAL));
  }

  // Face Element Stuff  
  INT notbdry=-666;
  INT myel,myopn;
  INT* ie = (INT *) calloc(2,sizeof(INT));
  INT* op_n = (INT *) calloc(2,sizeof(INT));

  // Element Node Stuff
  INT* myel_n = (INT *) calloc(el_order,sizeof(INT));
  REAL* p = (REAL *) calloc(el_order,sizeof(REAL));
  REAL* dpx = (REAL *) calloc(el_order,sizeof(REAL));
  REAL* dpy = (REAL *) calloc(el_order,sizeof(REAL));
  REAL* dpz=NULL;
  if(dim==3) { dpz = (REAL *) calloc(el_order,sizeof(REAL)); }
  REAL grad_mag,x,y,z,e1x,e1y,e1z,e2x,e2y,e2z;
  
  /* Get Face to Element Map */
  /* Get Transpose of f_el -> el_f */
    
    printf("hello-in1\n");
    
  iCSRmat f_el = icsr_create(nface,nelm,f_order*nelm);
  icsr_trans_1(el_f,&f_el);
    
    printf("hello-in2\n");


  // Loop over all Faces
  for(i=0;i<nface;i++) {
    /* Find Vertices in given Face */
    j_a = f_v->IA[i]-1;
    j_b = f_v->JA[i+1]-1;
    jcnt = 0;
    for (j=j_a; j<j_b;j++) {
      ipf[jcnt] = f_v->JA[j];
      xf[jcnt] = cv->x[ipf[jcnt]-1];
      yf[jcnt] = cv->y[ipf[jcnt]-1];
      if (dim==3) {
	zf[jcnt] = cv->z[ipf[jcnt]-1];
      }
      jcnt++;
    }

      printf("hello-in3\n");

      
    // Find Corresponding Elements and order in element 
    // Also picks correct opposite node to form vector 
    // normal vectors point from lower number element to higher one
    // or outward from external boundary
    j_a = f_el.IA[i]-1;
    j_b = f_el.IA[i+1]-1;
    jcnt=0;
    for (j=j_a; j<j_b; j++) {
      notbdry = j_b-j_a-1;
      ie[jcnt] = f_el.JA[j];
      op_n[jcnt] = f_el.val[j]-1;
      jcnt++;
    }
    if(notbdry && (ie[1]<ie[0])) {
      myopn = op_n[1];
      myel = ie[1];
    } else {
      myopn = op_n[0];
      myel = ie[0];
    }
    // Get Nodes of this chosen element
    j_a = el_v->IA[myel-1]-1;
    j_b = el_v->IA[myel]-1;
    jcnt=0;
    for(j=j_a;j<j_b;j++) {
      myel_n[jcnt] = el_v->JA[j];
      jcnt++;
    }
    x = cv->x[myel_n[myopn]];
    y = cv->y[myel_n[myopn]];
    z = -666.66;
    if(dim==3) {
      z = cv->z[myel_n[myopn]];
    }
      
      printf("hello-in3\n");


    /* Compute Area (length if 2D) and get midpt of face */
    if(dim==2) {
      f_area[i] = sqrt(pow(fabs(xf[1]-xf[0]),2)+pow(fabs(yf[1]-yf[0]),2));
      f_mid[i*dim] = (xf[0]+xf[1])/2.0;
      f_mid[i*dim+1] = (yf[0]+yf[1])/2.0;
    } else if(dim==3) {
      e1x = xf[1]-xf[0];
      e1y = yf[1]-yf[0];
      e1z = zf[1]-zf[0];
      e2x = xf[2]-xf[0];
      e2y = yf[2]-yf[0];
      e2z = zf[2]-zf[0];
      f_area[i] = 0.5*sqrt(pow(e1y*e2z-e2y*e1z,2)+pow(e1z*e2x-e2z*e1x,2)+pow(e1x*e2y-e2x*e1y,2));
      f_mid[i*dim] = (xf[0]+xf[1]+xf[2])/3.0;
      f_mid[i*dim+1] = (yf[0]+yf[1]+yf[2])/3.0;
      f_mid[i*dim+2] = (zf[0]+zf[1]+zf[2])/3.0;
    } else {
      baddimension();
    }
	
      printf("hello-in4\n");

      
    //Compute Normal Vectors based on opposite node
    //Get Linear Basis Functions for particular element
    PX_H1_basis(p,dpx,dpy,dpz,x,y,z,myel_n,1,mesh);
      printf("hello-in5\n");
    grad_mag = dpx[myopn]*dpx[myopn]+dpy[myopn]*dpy[myopn];
      printf("hello-in6\n");
    if(dim==3) {
      grad_mag = grad_mag + dpz[myopn]*dpz[myopn];
    }
      printf("hello-in7\n");
    grad_mag = -sqrt(grad_mag);
    f_norm[i*dim] = dpx[myopn]/grad_mag;
    f_norm[i*dim+1] = dpy[myopn]/grad_mag;
      printf("hello-in8\n");
    if(dim==3) {
      f_norm[i*dim+2] = dpz[myopn]/grad_mag;
    }
      printf("hello-in9\n");

  }


    
  icsr_free(&f_el);
  if(ipf) free(ipf);
  if(xf) free(xf);
  if(yf) free(yf);
  if(p) free(p);
  if(dpx) free(dpx);
  if(dpy) free(dpy);
  if(op_n) free(op_n);
  if(ie) free(ie);
  if(myel_n) free(myel_n);
  if(dim==3) { 
    if(dpz) free(dpz);
    if(zf) free(zf);
  }

  return;
}
/*********************************************************************************************************/

/*********************************************************************************************************/
void get_el_mid(REAL *el_mid,iCSRmat el_v,coordinates cv,INT dim) 
{
	
	
  /*** Compute the midpts of triangluar element using the vertices
   *
   *    Input:   
   *             el_v		Element to Vertex Map
   *             cv		Coordinates of Vertices
   *
   *    Output: el_mid		Midpoint of element
   */
	
  INT i,j,cnt,nd,acol,bcol; /* Loop Index */
  INT nelm = el_v.row;
	
  if (dim==2) {
		
    for (i=0; i<nelm; i++) {
      acol = el_v.IA[i]-1;
      bcol = el_v.IA[i+1]-1;
      cnt=0;
      el_mid[i*dim]=0;
      el_mid[i*dim+1]=0;
      for (j=acol; j<bcol; j++) {
	nd = el_v.JA[j]-1;
	el_mid[i*dim] = el_mid[i*dim]+cv.x[nd];
	el_mid[i*dim+1] = el_mid[i*dim+1]+cv.y[nd];
	cnt++;
      }
      el_mid[i*dim]=el_mid[i*dim]/3.0;
      el_mid[i*dim+1]=el_mid[i*dim+1]/3.0;
    }
  } else if (dim==3) {
    for (i=0; i<nelm; i++) {
      acol = el_v.IA[i]-1;
      bcol = el_v.IA[i+1]-1;
      cnt=0;
      el_mid[i*dim]=0;
      el_mid[i*dim+1]=0;
      el_mid[i*dim+2]=0;
      for (j=acol; j<bcol; j++) {
	nd = el_v.JA[j]-1;
	el_mid[i*dim] = el_mid[i*dim]+cv.x[nd];
	el_mid[i*dim+1] = el_mid[i*dim+1]+cv.y[nd];
	el_mid[i*dim+2] = el_mid[i*dim+2]+cv.z[nd];;
	cnt++;
      }
      el_mid[i*dim]=0.25*el_mid[i*dim];
      el_mid[i*dim+1]=0.25*el_mid[i*dim+1];
      el_mid[i*dim+2]=0.25*el_mid[i*dim+2];
    }
  }
	
  return;
}
/*********************************************************************************************************/

/*********************************************************************************************************/
void get_el_vol(REAL *el_vol,iCSRmat el_v,coordinates cv,INT dim,INT v_per_elm) 
{
  /*** Compute the area/volume of ALL triangluar/tetrahedral elements using the vertices
   *
   *    Input:   
   *             el_v		Element to Vertices Map
   *             coordinates	Coordinates of Vertices
   *             dim            Dimension of Problem
   *             v_per_elm      Number of vertices per element
   *
   *    Output:  el_vol         Area/Volume of each element
   */
	
  INT i,j_a,j_b,j,jcnt; /* Loop Index */
  REAL x2,x3,x4,y2,y3,y4,z2,z3,z4;
  REAL* x = (REAL *) calloc(v_per_elm,sizeof(REAL));		/* Coordinates of nodes on elements */
  REAL* y = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* z=NULL;
  INT nelm = el_v.row;

  if(dim==2) {
    for (i=0;i<nelm;i++) {
      j_a = el_v.IA[i];
      j_b = el_v.IA[i+1]-1;
      jcnt=0;
      for (j=j_a; j<=j_b; j++) {
	x[jcnt] = cv.x[el_v.JA[j-1]-1];
	y[jcnt] = cv.y[el_v.JA[j-1]-1];
	jcnt++;
      }
      el_vol[i] = 0.5*fabs(y[0]*(x[1]-x[2]) + y[1]*(x[2]-x[0]) + y[2]*(x[0]-x[1]));
    }
  } else if(dim==3) {
    z = (REAL *) calloc(v_per_elm,sizeof(REAL));
	
    for (i=0;i<nelm;i++) {
      j_a = el_v.IA[i];
      j_b = el_v.IA[i+1]-1;
      jcnt=0;
      for (j=j_a; j<=j_b; j++) {
	x[j] = cv.x[el_v.JA[j-1]-1];
	y[j] = cv.y[el_v.JA[j-1]-1];
	z[j] = cv.z[el_v.JA[j-1]-1];
	jcnt++;
      }
      x2 = x[1]-x[0];
      x3 = x[2]-x[0];
      x4 = x[3]-x[0];
      y2 = y[1]-y[0];
      y3 = y[2]-y[0];
      y4 = y[3]-y[0];
      z2 = z[1]-z[0];
      z3 = z[2]-z[0];
      z4 = z[3]-z[0];
	
      el_vol[i] = (REAL) (fabs(x2*(y3*z4-y4*z3) - y2*(x3*z4-x4*z3) + z2*(x3*y4-x4*y3)))/6.0;
    }	
  }
	
  if(x) free(x);
  if(y) free(y);
  if(z) free(z);
  return;
}
/*********************************************************************************************************/

/****************************************************************************************/
void free_mesh(trimesh* mesh)
{
  /* frees memory of arrays in mesh struct */

  free_coords(mesh->cv);
  icsr_free(mesh->el_v);
  icsr_free(mesh->el_ed);
  icsr_free(mesh->el_f);
  icsr_free(mesh->ed_v);
  icsr_free(mesh->f_v);
  if(mesh->el_vol) free(mesh->el_vol);
  if(mesh->el_mid) free(mesh->el_mid);
  if(mesh->ed_len) free(mesh->ed_len);
  if(mesh->ed_tau) free(mesh->ed_tau);
  if(mesh->ed_mid) free(mesh->ed_mid);
  if(mesh->f_area) free(mesh->f_area);
  if(mesh->f_norm) free(mesh->f_norm);
  if(mesh->f_mid) free(mesh->f_mid);
  if(mesh->v_bdry) free(mesh->v_bdry);
  if(mesh->ed_bdry) free(mesh->ed_bdry);
  if(mesh->f_bdry) free(mesh->f_bdry);
  
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void get_incidence_row(INT row,iCSRmat *fem_map,INT* thisrow)
{
  /* Gets single row of an incidence map (i.e., Gets vertices of given element from el_v) */

  INT j;
  INT rowa = fem_map->IA[row];
  INT rowb = fem_map->IA[row+1]-1;
  INT jcntr = 0;
  for (j=rowa; j<=rowb; j++) {
    thisrow[jcntr] = fem_map->JA[j-1];
    jcntr++;
  }
  
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void dump_coords(FILE* fid,coordinates *c) 
{
  /* Dump the coordinate data to file for plotting purposes
   *
   * Input:		
   *          c:      Coordinate data.
   *
   * Output:		
   *          coord.dat:    coord(n,dim)  coordinates of nodes
   *
   */
	
  INT i; /* Loop Indices */
		
  if (c->z) {
    for (i=0; i<c->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\t%25.16e\n",c->x[i],c->y[i],c->z[i]);
    }
  } else {
    for (i=0; i<c->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\n",c->x[i],c->y[i]);
    }
  }
   		
  return;
}
/****************************************************************************************/


