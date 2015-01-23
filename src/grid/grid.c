/*
 *  grid.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
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
   */
	
  INT i,j; /* Loop indices */

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
  if(bdry_v) free(bdry_v);

  printf("\nConverting Grid Maps to CSR and Computing Data Structures:\n ");	
  /* Element Vertex Map */
  iCSRmat el_v = convert_elmnode(element_vertex,nelm,nv,v_per_elm);
  if(element_vertex) free(element_vertex);
	
  /* Edge to Node Map */
  INT nedge = 0;
  if(dim==2) {  
    nedge = nelm+nv-(dim-1);
  } else if(dim==3) {
    get_nedge(&nedge,el_v); 
  }
  iCSRmat ed_v = get_edge_v(nedge,el_v);
	
  /* Get Boundary Edges and Nodes */
  INT* ed_bdry = (INT *) calloc(nedge,sizeof(INT));
  if (dim==2) {
    isboundary_ed(ed_v,nedge,nbedge,bdry_v,ed_bdry);
    isboundary_v(nv,bdry_v,v_bdry,nbedge,&nbv);
    if(bdry_v) free(bdry_v);
  } else if(dim==3) {
    isboundary_ed3D(ed_v,nedge,cv,&nbedge,v_bdry,ed_bdry);
  }
	
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
  

  // Assign components to the mesh
  mesh_temp.dim = dim;
  mesh_temp.nv = nv;
  mesh_temp.v_per_elm = v_per_elm;
  mesh_temp.nedge = nedge;
  mesh_temp.ed_per_elm = ed_per_elm;
  mesh_temp.nface = nface;
  mesh_temp.f_per_elm = f_per_elm;
  mesh_temp.x = cv.x;
  mesh_temp.y = cv.y;
  mesh_temp.z = cv.z;
  mesh_temp.nbv = nbv;
  mesh_temp.nbedge = nbedge;
  //mesh_temp.nbface = nbface;
  mesh_temp.el_v = &el_v;
  mesh_temp.el_ed = &el_ed;
  //mesh_temp.el_f = &el_f;
  mesh_temp.ed_v = &ed_v;
  //mesh_temp.f_v = &f_v;
  //mesh_temp.el_vol = &el_vol;
  //mesh_temp.el_mid = &el_mid;
  mesh_temp.ed_len = ed_len;
  mesh_temp.ed_tau = ed_tau;
  mesh_temp.ed_mid = ed_mid;
  //mesh_temp.f_area = &f_area;
  //mesh_temp.f_norm = &f_norm;
  
  mesh_temp.v_bdry = v_bdry;
  mesh_temp.ed_bdry = ed_bdry;
  //mesh_temp.f_bdry = f_bdry;
  
  *mesh = mesh_temp;
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
  icsr_mxm_1(&v_el,&el_v,&v_v);
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
  free(iv_v);
  free(jv_v);
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
#define max(a,b) ((a) > (b) ? (a) : (b))
    if(max(ip[0],ip[1])==ip[0]) {
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

/* /\***********************************************************************************************\/ */
/* void get_face_maps(iCSRmat el_v,INT el_order,INT nface,INT dim,INT f_order,iCSRmat el_f,INT *f_bdry,INT *nbf,iCSRmat f_n,INT *fel_order)  */
/* { */
	
/*   /\* Gets the Element to Face, Face to Node, Face to Boundary mapping in CSR Fromat as well as face ordering on element (Should be dim independent) */
/*    *     */
/*    *               3 */
/*    *             /   \ */
/*    *            2     1 */
/*    *           /       \ */
/*    *          1 ---3----2 */
/*    * */
/*    *	Input:  el_v:		Element to Vertex Map */
/*    *		el_order:	Number of nodes per element */
/*    *            nface           Total Number of Faces  */
/*    *            dim             Dimension */
/*    *            f_order         Number of Faces per element 3 in 2D, 4 in 3D */
/*    * */
/*    *	Output:  */
/*    *		el_f	Element to Face Map where value represents face number in element */
/*    *            f_bdry  Indicates whether given face is on boundary or not */
/*    *            nbf     Number of boundary faces */
/*    *            f_n     Face to Node Map */
/*    *            fel_order Indicates order of faces on an element */
/*    * */
/*    *\/ */
	
/*   INT i,j,k,m,p,jk,col_b,icntr,jcntr,kcntr; /\* Loop Indices *\/ */
/*   INT ncol1,ncol2,ncol3,ncol4,ncol5,ncol6;  /\* Node indices *\/ */
/*   INT el1,el2,el3,el; /\* Element indices *\/ */
/*   INT iflag = 0;  /\* Marks if we found a face *\/ */
/*   INT nbf = 0; /\* hold number of boundary faces *\/ */
/*   INT nd[dim]; */
/*   INT f_num;  /\* Current face number of element *\/ */
/*   INT* in_f; */
/*   INT* jn_f; */

/*   // We will build face to element map first and then transpose it */
/*   INT* iface_el = (INT *) calloc(nface+1,sizeof(INT)); */
/*   INT* jface_el = (INT *) calloc(nelm*f_order,sizeof(INT)); */
/*   INT* face_el = (INT *) calloc(nelm*f_order,sizeof(INT)); */

/*   // First order the faces per element.  Start with node 1 and indicate opposite face as face 1.  Then loop over nodes */
/*   for(j=0;j<f_order;j++) { */
/*     for(i=0;i<dim;i++) { */
/*       jcntr = j+i+2; */
/*       col_b = j*dim+i; */
/*       if(jcntr>el_order) { */
/* 	fel_order[col_b] = jcntr-el_order; */
/*       } else { */
/* 	fel_order[col_b] = jcntr; */
/*       } */
/*     } */
/*   } */
/*   /\* Get Transpose of el_n -> n_el *\/ */
/*   INT* in_el = calloc(n+1,sizeof(INT)); */
/*   INT* jn_el = calloc(element_order*nelm,sizeof(INT)); */
/*   atransps(iel_n,jel_n,nelm,n,in_el,jn_el); */

/*   // Get interior faces first */
/*   // Loop over All elements and each face on the element */
/*   icntr=0; */
/*   jcntr=0; */
/*   kcntr=0; */
/*   iflag=0; */
/*   for(i=0;i<nelm;i++) { */
/*     // Get first node of element */
/*     col_b = iel_n[i]-1; */
/*     // Now loop over all faces of element */
/*     for(j=0;j<face_order;j++) { */
/*       // Get appropriate nodes on face */
/*       for(k=0;k<mydim;k++) { */
/* 	jk = col_b+fel_order[j*mydim+k]-1; */
/* 	nd[k] = jel_n[jk]; */
/*       } */
/*       // Next Loop over all elements that touch the first node */
/*       ncol1 = in_el[nd[0]-1]-1; */
/*       ncol2 = in_el[nd[0]]-1; */
/*       for(k=ncol1;k<ncol2;k++) { */
/* 	el1 = jn_el[k]-1; */
/* 	if(el1!=i) { // If not the same element we started in then check other nodes elements */
/* 	  ncol3 = in_el[nd[1]-1]-1; */
/* 	  ncol4 = in_el[nd[1]]-1; */
/* 	  for(m=ncol3;m<ncol4;m++) { */
/* 	    el2 = jn_el[m]-1; */
/* 	    if(el2==el1) { // Our two nodes share the element!  In 2D this is the face! In 3D, we check third guy */
/* 	      if(mydim==2) { */
/* 		// Mark the element that shares the face as well as the nodes that are shared */
/* 		iflag = 1; */
/* 		el = el1; */
/* 	      } else if(mydim==3) { */
/* 		ncol5 = in_el[nd[2]-1]-1; */
/* 		ncol6 = in_el[nd[2]]-1; */
/* 		for(p=ncol5;p<ncol6;p++) { */
/* 		  el3 = jn_el[p]-1; */
/* 		  if(el3==el2) { */
/* 		    // Mark the element that shares the face as well as the nodes that are shared */
/* 		    iflag = 1; */
/* 		    el = el1; */
/* 		  } */
/* 		} */
/* 	      } else { */
/* 		baddimension(); */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*       if(iflag && el>i) { */
/* 	iface_el[icntr] = jcntr+1; */
/* 	jface_el[jcntr] = i+1; */
/* 	face_el[jcntr] = j+1; */
/* 	jface_el[jcntr+1] = el+1; */
/* 	// Find face number for other element */
/* 	find_facenumber(iel_n,jel_n,fel_order,el+1,nd,mydim,&f_num); */
/* 	face_el[jcntr+1] = f_num; */
/* 	face_bdry[icntr] = 0; */
/* 	if_n[icntr] = kcntr+1; */
/* 	jf_n[kcntr] = nd[0]; */
/* 	jf_n[kcntr+1] = nd[1]; */
/* 	//	fprintf(stdout,"\n1:kcntr=%d\n",kcntr+1); */
/* 	kcntr+=2; */
/* 	if(mydim==3) {  */
/* 	  jf_n[kcntr] = nd[2]; */
/* 	  //	  fprintf(stdout,"\n2:kcntr=%d\n",kcntr); */
/* 	  kcntr++; */
/* 	} */
/* 	icntr++; */
/* 	jcntr+=2; */
/*       } else if(!iflag){  // this must be a boundary face! */
/* 	nbf++; */
/* 	iface_el[icntr] = jcntr+1; */
/* 	jface_el[jcntr] = i+1; */
/* 	face_el[jcntr] = j+1; */
/* 	face_bdry[icntr] = 1; */
/* 	if_n[icntr] = kcntr+1; */
/* 	jf_n[kcntr] = nd[0]; */
/* 	jf_n[kcntr+1] = nd[1]; */
/* 	//	fprintf(stdout,"\n3:kcntr=%d\n",kcntr+1); */
/* 	kcntr+=2; */
/* 	if(mydim==3) {  */
/* 	  //	  fprintf(stdout,"\n4:kcntr=%d\n",kcntr); */
/* 	  jf_n[kcntr] = nd[2]; */
/* 	  kcntr++; */
/* 	} */
/* 	icntr++; */
/* 	jcntr++; */
/*       }	 */
/*       iflag=0; */
/*     } */
/*   } */

/*   iface_el[icntr] = jcntr+1; */
/*   if_n[icntr] = kcntr+1; */
/*   //  fprintf(stdout,"\nicntr=%d,kcntr=%d,jcntr=%d\n",icntr,kcntr,jcntr); */
/*   //  fflush(stdout); */
/*   /\* Get Transpose of face_el -> el_face *\/ */
/*   atransps(iface_el,jface_el,nface,nelm,iel_face,jel_face); */
/*   atransp_int(iface_el,jface_el,face_el,nface,nelm,iel_face,jel_face,el_face); */

/*   /\* Transpose face_node and back again to order nodes in increasing order (from global structure) *\/ */
/*   in_f = calloc(n+1,sizeof(INT)); */
/*   jn_f = calloc(nface*mydim,sizeof(INT)); */
/*   atransps(if_n,jf_n,nface,n,in_f,jn_f); */
/*   atransps(in_f,jn_f,n,nface,if_n,jf_n); */

/*   *nbface = nbf; */

/*   if(in_el) free(in_el); */
/*   if(jn_el) free(jn_el); */
/*   if(iface_el) free(iface_el); */
/*   if(jface_el) free(jface_el); */
/*   if(face_el) free(face_el); */
/*   if(in_f) free(in_f); */
/*   if(jn_f) free(jn_f); */

/*   return; */
/* } */
/* /\***********************************************************************************************\/ */
