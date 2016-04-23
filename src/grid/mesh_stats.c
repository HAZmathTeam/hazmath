/*
 *  mesh_stats.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 *  Obtains routines for generating properties of the mesh,
 *  including all incidence matrices, element volumes, etc.
 *
 */

#include "hazmat.h"

/******************************************************************************/
iCSRmat convert_elmnode(INT *element_vertex,INT nelm,INT nv,INT nve) 
{	
  /* Convert the input element to vertex map into sparse matrix form 
   *
   * Input: 
   *  nelm:	                Number of elements
   *  nv:		        Number of vertices
   *  nve:	                Number of vertices per element
   *  element_vertex(nelm,nve): Each row is an element, each column is
   *                            the corresponding vertices for that element 
   *										
   * Output:	
   *  el_v:	                Element to vertex map in CSR format.  Rows  
   *                            are all elements, columns are all vertices. 
   *                            Entries are 1 if there's a connection.
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
      el_v.JA[k] = element_vertex[i*nve+j];
      k=k+1;
    }
  }
  
  el_v.val = NULL;
  el_v.row = nelm; el_v.col = nv; el_v.nnz = nelm*nve;
	
  return el_v;
}
/*******************************************************************************/

/*******************************************************************************/
void get_nedge(INT* nedge, iCSRmat el_v) 
{	
  /* Gets the Number of Edges 
   *
   * Input:  
   *  el_v:	Element to Node Map
   *
   * Output: 
   *  nedge:	Total number of Edges
   *
   */
		
  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans_1(&el_v,&v_el);
	
  /* Create Vertex to Vertex Map by v_el*el_v */
  iCSRmat v_v;
  icsr_mxm_symb_1(&v_el,&el_v,&v_v);

	
  // Now go through upper triangular (above diagonal) portion of vertex-vertex 
  // and count non-zeros
  // Each of these vertex-vertex connections are edges.  Upper so we don't 
  // count i->j and j->i as two separate edges, and no diagonal since i->i 
  // is not an edge.

  *nedge = (INT) (v_v.nnz - v_v.row)/2;

  /* Free matrices */
  icsr_free(&v_v);
  icsr_free(&v_el);
  return;
}
/*******************************************************************************/

/*******************************************************************************/
iCSRmat get_edge_v(INT nedge,iCSRmat el_v) 
{
	
  /* Gets the Edge to Vertex mapping in CSR Fromat (Should be dim independent)
   * This is used to get Element to Edge map, but also for determining
   * things such as edge_length and the tangent vectors of edges 
   *
   * Input: 
   *  el_v:    Element to Vertex Map
   *  nedge:   Total number of Edges
   *	       
   * Output:
   *  ed_v:    Edge to vertex map
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
	  
  // Now go through upper triangular (above diagonal) portion of node-node and 
  // count non-zeros.
  // Each of these node-node connections are edges.  Upper so we don't count
  // i->j and j->i as two separate edges, and no diagonal since i->i is not 
  // an edge.
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
  icsr_free(&v_v);
  icsr_free(&v_el);

  ed_v.val=NULL;
  ed_v.row = nedge;
  ed_v.col = nv;
  ed_v.nnz = 2*nedge;
    
  return ed_v;
}
/*******************************************************************************/

/*******************************************************************************/
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
/*******************************************************************************/

/*******************************************************************************/
void isboundary_ed(iCSRmat ed_v,INT nedge,INT nbedge,INT *bdry_v,INT *ed_bdry) 
{	
  /* Indicates whether an edge is on boundary or not
   *
   *	Input:  ed_v		Edge to Vertex Map
   *		nedge		Total number of Edges
   *		nbedge:		Number of boundary edges
   *		bdry_v		List of boundaries and corresponding vertices
   *	Output: 
   *		ed_bdry		List of Edges whether boundary or not 
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
/*******************************************************************************/

/*******************************************************************************/
void isboundary_ed3D(iCSRmat ed_v,INT nedge,coordinates *cv,INT *nbedge,INT *v_bdry,INT *ed_bdry) 
{
	
  /* Counts the number of boundary edges and indicates whether an edge is a
   * boundary
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
    col_b = ed_v.IA[i]-1;
    col_e = ed_v.IA[i+1]-1;
    n = ed_v.JA[col_b]-1;
    m = ed_v.JA[col_b+1]-1;
    //Check if two nodes are on boundary
    if (v_bdry[n]!=0 && v_bdry[m]!=0) {
      if(cv->x[n]==cv->x[m] || cv->y[n]==cv->y[m] || cv->z[n]==cv->z[m]) {
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
/*******************************************************************************/

/*******************************************************************************/
void isboundary_ed_f(iCSRmat* f_ed,iCSRmat* ed_v,INT nedge,INT nface,INT *f_bdry,INT *v_bdry,INT *nbedge,INT *ed_bdry) 
{	
  /* Counts the number of boundary edges and indicates whether an edge is 
   * a boundary, using the f_bdry map.
   *
   *	Input:  f_ed:		Edge to Vertex Map
   *		nedge	    	Total number of edges
   *		nface		Total number of faces
   *		f_bdry		List of faces and whether they are on boundary
   *	Output: 
   *		ed_bdry		List of edges and whether they are on boundary
   *		nbedge	        Number of boundary edges
   */
	
  INT i,j,ed,v1,v2,col_b,col_e,jcntr; /* Loop indices and counters */
	
  // For every face get edges
  // All edges on a face should have the same property as the face.
  for (i=0; i<nface; i++) {
    col_b = f_ed->IA[i]-1;
    col_e = f_ed->IA[i+1]-1;
    for(j=col_b;j<col_e;j++) {
      ed = f_ed->JA[j]-1;
      v1 = ed_v->JA[ed_v->IA[ed]-1]-1;
      v2 = ed_v->JA[ed_v->IA[ed]]-1;
      if(v_bdry[v1]==v_bdry[v2]) {
	if(f_bdry[i]==v_bdry[v1]) {
	  ed_bdry[ed] = f_bdry[i];
	}
      }
    }
  }

  // Now go back through edges and count how many are on boundary.
  jcntr=0;
  for (i=0; i<nedge; i++) {
    if(ed_bdry[i]!=0) {
      jcntr++;
    }
  }

  *nbedge = jcntr;

  return;
}
/*******************************************************************************/

/*******************************************************************************/
iCSRmat get_el_ed(iCSRmat el_v,iCSRmat ed_v) 
{
  /* Gets the Element to Edge mapping in CSR Fromat (Should be dim independent)
   *
   * Input: 
   *  el_v:	  Element to Vertex Map
   *  ed_v:	  Edge to Vertex Map
   *
   * Output:
   *  el_ed:      Element to Edge map
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
/*******************************************************************************/

/*******************************************************************************/
void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates *cv,iCSRmat ed_v,INT dim) 
{
  /* Get length, tangent vector (tau), and midpoint of every edge 
   *
   * Input:   
   *  cv        Coordinates of Vertices
   *  ed_v      Edge to Vertex Map
   *  dim       Dimension of Problem
   *
   * Output:  
   *  ed_len    Lengths of Edges
   *  ed_tau    Tangent Vector of Edges ordered (tx1,ty1,tz1,tx2,ty2,tz2,...)
   *  ed_mid    Midpoint of Edges ordered (mx1,my1,mz1,mx2,my2,mz2,...) 
   */
	
  INT i,jcnt,j,j_a,j_b; /* loop index */
  INT ip[2];
  REAL x[2],y[2],z[2];
  ip[0]=0; ip[1]=0; x[0]=0.0; x[1]=0.0; y[0]=0.0; y[1]=0.0; z[0]=0.0; z[1]=0.0;
	
  /* Orders nodes from largest to smallest for correct orientation */
  INT ihi1,ilo1; 
  INT nedge = ed_v.row;

  for(i=0;i<nedge;i++) {

    /* Find Nodes in given Edge */
    j_a = ed_v.IA[i]-1;
    j_b = ed_v.IA[i+1]-1;
    jcnt = 0;
    for (j=j_a; j<j_b;j++) {
      ip[jcnt] = ed_v.JA[j];
      x[jcnt] = cv->x[ip[jcnt]-1];
      y[jcnt] = cv->y[ip[jcnt]-1];
      if (dim==3) {
	z[jcnt] = cv->z[ip[jcnt]-1];
      } else {
	z[jcnt] = 0;
      }
      jcnt++;
    }
	
    /* Compute Length */
    ed_len[i] = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + \
		     (z[1]-z[0])*(z[1]-z[0]));
	
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
/*******************************************************************************/

/*******************************************************************************/
void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order)
{
  /* Gets the face ordering on element 
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

  // First order the faces per element.  Start with node 1 and indicate 
  // opposite face as face 1.  Then loop over nodes
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
/*******************************************************************************/

/*******************************************************************************/
void get_face_maps(iCSRmat el_v,INT el_order,iCSRmat ed_v,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,iCSRmat *f_ed,INT *fel_order)
{
  /* Gets the Element to Face, Face to Node, Face to Boundary mapping in CSR 
   * format as well as face ordering on element (Should be dim independent)
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
   *            fel_order       Indicates order of faces on an element
   *
   *	Output:
   *		el_f	Element to Face Map where value represents face
   *                     number in element
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
  INT edpf = 2*dim - 3;

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
	// If not the same element we started in then check other nodes elements
	if(el1!=i) { 
	  ncol3 = v_el.IA[nd[1]-1]-1;
	  ncol4 = v_el.IA[nd[1]]-1;
	  for(m=ncol3;m<ncol4;m++) {
	    el2 = v_el.JA[m]-1;
	    // Our two nodes share the element!  In 2D this is the face! 
	    // In 3D, we check third guy
	    if(el2==el1) { 
	      if(dim==2) {
		// Mark the element that shares the face as well as the
		// nodes that are shared
		iflag = 1;
		el = el1;
	      } else if(dim==3) {
		ncol5 = v_el.IA[nd[2]-1]-1;
		ncol6 = v_el.IA[nd[2]]-1;
		for(p=ncol5;p<ncol6;p++) {
		  el3 = v_el.JA[p]-1;
		  if(el3==el2) {
		    // Mark the element that shares the face as well as 
		    // the nodes that are shared
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
	find_facenumber(el_v,el+1,nd,dim,&f_num);
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

  /* Transpose face_node and back again to order nodes in increasing order 
   * (from global structure) */
  iCSRmat v_f;
  icsr_trans_1(f_v,&v_f);
  icsr_free(f_v);
  icsr_trans_1(&v_f,f_v);

  /* Get Face to Edge Map */
  iCSRmat face_ed;
  iCSRmat v_ed;
  icsr_trans_1(&ed_v,&v_ed);
  icsr_mxm_symb_max_1(f_v,&v_ed,&face_ed,2);
  for(i=0;i<nface+1;i++)
    f_ed->IA[i] = face_ed.IA[i];
  for(i=0;i<nface*edpf;i++)
    f_ed->JA[i] = face_ed.JA[i];
  icsr_free(&face_ed);
  icsr_free(&v_ed);
  
  *nbface = nbf;

  icsr_free(&v_el);
  icsr_free(&f_el);
  icsr_free(&v_f);

  return;
}
/*******************************************************************************/

/*******************************************************************************/
void find_facenumber(iCSRmat el_v,INT elm,INT* nd,INT dim,INT *f_num)          
{

  /* Find the face number of the given element, using the element it shares the face with
   *
   * Input:
   *        el_v              Element to Vertex Map (all elements and vertices)
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
/*******************************************************************************/

/*******************************************************************************/
void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,iCSRmat *f_v,trimesh *mesh) 
{
  /* Get area, normal vector for all faces
   *
   *    Input:   
   *             mesh                 Information needed for mesh
   *
   *    Output:  f_area               Area of each face (length in 2D)
   *             f_norm(nface,dim)    Normal vectors of each face
   *             f_mid(nface,dim)     Midpoints of face
   *            
   */
    	
  INT i,jcnt,j,j_a,j_b; /* loop index */
  INT nface = mesh->el_f->col;
  INT dim = mesh->dim;
  INT el_order = mesh->v_per_elm;

  coordinates *cv = mesh->cv;

  iCSRmat *el_f = mesh->el_f;
  iCSRmat *el_v = mesh->el_v;

  // Face Node Stuff
  INT* ipf = (INT *) calloc(dim,sizeof(INT));
  REAL* xf = (REAL *) calloc(dim,sizeof(REAL));
  REAL* yf = (REAL *) calloc(dim,sizeof(REAL));
  REAL* zf;
  if(dim==3) {
    zf = calloc(dim,sizeof(REAL));
  }
  REAL* myx = (REAL *) calloc(dim,sizeof(REAL));

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
  REAL grad_mag,e1x,e1y,e1z,e2x,e2y,e2z;
  
  /* Get Face to Element Map */
  /* Get Transpose of f_el -> el_f */    
  iCSRmat f_el;
  icsr_trans_1(el_f,&f_el);

  // Loop over all Faces
  for(i=0;i<nface;i++) {
    /* Find Vertices in given Face */
    j_a = f_v->IA[i]-1;
    j_b = f_v->IA[i+1]-1;
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
    myx[0] = cv->x[myel_n[myopn]-1];
    myx[1] = cv->y[myel_n[myopn]-1];
    if(dim==3) {
      myx[2] = cv->z[myel_n[myopn]-1];
    }

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
      f_area[i] = 0.5*sqrt(pow(e1y*e2z-e2y*e1z,2)+pow(e1z*e2x-e2z*e1x,2)+ \
			   pow(e1x*e2y-e2x*e1y,2));
      f_mid[i*dim] = (xf[0]+xf[1]+xf[2])/3.0;
      f_mid[i*dim+1] = (yf[0]+yf[1]+yf[2])/3.0;
      f_mid[i*dim+2] = (zf[0]+zf[1]+zf[2])/3.0;
    } else {
      baddimension();
    }
      
    // Compute Normal Vectors based on opposite node
    // Get Linear Basis Functions for particular element
    PX_H1_basis(p,dpx,dpy,dpz,myx,myel_n,1,mesh);
    grad_mag = dpx[myopn]*dpx[myopn]+dpy[myopn]*dpy[myopn];
    if(dim==3) {
      grad_mag = grad_mag + dpz[myopn]*dpz[myopn];
    }
    grad_mag = -sqrt(grad_mag);
    f_norm[i*dim] = dpx[myopn]/grad_mag;
    f_norm[i*dim+1] = dpy[myopn]/grad_mag;
    if(dim==3) {
      f_norm[i*dim+2] = dpz[myopn]/grad_mag;
    }
  }

  icsr_free(&f_el);
  if(ipf) free(ipf);
  if(xf) free(xf);
  if(yf) free(yf);
  if(myx) free(myx);
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
/********************************************************************************/

/********************************************************************************/
void sync_facenode(iCSRmat *f_v,REAL* f_norm,trimesh *mesh)           
{
  /* Reorder the Face-Node mapping so it has positive orientation with 
   * respect to the face's normal vector
   *
   * Input:
   *        if_n,jf_n         Face to Node Map
   *        ndpf              Number of Nodes per Face
   *        mydim             Dimension
   *        nface             Total Number of Faces
   *
   * Output:
   *       if_n,jf_n         Reordered Face to Node Map
   *
   */
  
  INT i,j; /* loop index */
  INT nface = mesh->el_f->col;
  INT dim = mesh->dim;
  INT ndpf = dim;
  coordinates *cv = mesh->cv;
  
  REAL nx,ny,nz,tx,ty,tz,mysign;
  INT nd,rowa,rowb,jcnt,nf1,nf2,nf3;
  REAL* xf = calloc(ndpf,sizeof(REAL));
  REAL* yf = calloc(ndpf,sizeof(REAL));
  REAL* zf = calloc(ndpf,sizeof(REAL));

  if(dim==2) {
    for(i=0;i<nface;i++) {
      // Get normal vector of face
      nx = f_norm[(i)*dim];
      ny = f_norm[(i)*dim+1];

      // Get Coordinates of Nodes
      rowa = f_v->IA[i]-1;
      rowb = f_v->IA[i+1]-1;
      jcnt=0;
      for(j=rowa;j<rowb;j++) {
	nd = f_v->JA[j]-1;
	xf[jcnt] = cv->x[nd];
	yf[jcnt] = cv->y[nd];
	jcnt++;
      }
      // Determine proper orientation of basis vectors  Compute n^(\perp)*t. 
      // If + use face_node ordering, if - switch sign
      tx = xf[1]-xf[0];
      ty = yf[1]-yf[0];
      mysign = -ny*tx + nx*ty;
      if(mysign<0) {
	nf2 = f_v->JA[rowa];
	nf1 = f_v->JA[rowa+1];
	f_v->JA[rowa+1] = nf2;
	f_v->JA[rowa] = nf1;
      } 
    }
  } else if (dim==3) {
    for(i=0;i<nface;i++) {
      // Get normal vector of face
      nx = f_norm[(i)*dim];
      ny = f_norm[(i)*dim+1];
      nz = f_norm[(i)*dim+2];

      // Get Coordinates of Nodes
      rowa = f_v->IA[i]-1;
      rowb = f_v->IA[i+1]-1;
      jcnt=0;
      for(j=rowa;j<rowb;j++) {
	nd = f_v->JA[j]-1;
	xf[jcnt] = cv->x[nd];
	yf[jcnt] = cv->y[nd];
	zf[jcnt] = cv->z[nd];
	jcnt++;
      }
      // Determine proper orientation of basis vectors  Compute n^(\perp)*t. 
      // If + use face_node ordering, if - switch sign
      tx = (yf[1]-yf[0])*(zf[2]-zf[0]) - (zf[1]-zf[0])*(yf[2]-yf[0]);
      ty = (zf[1]-zf[0])*(xf[2]-xf[0]) - (xf[1]-xf[0])*(zf[2]-zf[0]);
      tz = (xf[1]-xf[0])*(yf[2]-yf[0]) - (yf[1]-yf[0])*(xf[2]-xf[0]);
      mysign = nx*tx + ny*ty + nz*tz;
      if(mysign<0) {
	nf3=f_v->JA[rowa+1];
	nf2=f_v->JA[rowa+2];
	f_v->JA[rowa+1]=nf2;
	f_v->JA[rowa+2]=nf3;
      } 
    }
  }

  if(xf) free(xf);
  if(yf) free(yf);
  if(zf) free(zf);
  return;
}
/********************************************************************************/

/********************************************************************************/
void get_el_mid(REAL *el_mid,iCSRmat el_v,coordinates *cv,INT dim) 
{	
  /* Compute the midpts of triangluar element using the vertices
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
	el_mid[i*dim] = el_mid[i*dim]+cv->x[nd];
	el_mid[i*dim+1] = el_mid[i*dim+1]+cv->y[nd];
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
	el_mid[i*dim] = el_mid[i*dim]+cv->x[nd];
	el_mid[i*dim+1] = el_mid[i*dim+1]+cv->y[nd];
	el_mid[i*dim+2] = el_mid[i*dim+2]+cv->z[nd];;
	cnt++;
      }
      el_mid[i*dim]=0.25*el_mid[i*dim];
      el_mid[i*dim+1]=0.25*el_mid[i*dim+1];
      el_mid[i*dim+2]=0.25*el_mid[i*dim+2];
    }
  }
	
  return;
}
/********************************************************************************/

/********************************************************************************/
void get_el_vol(REAL *el_vol,iCSRmat el_v,coordinates *cv,INT dim,INT v_per_elm) 
{
  /* Compute the area/volume of ALL triangluar/tetrahedral elements using the vertices
   *
   *    Input:   
   *             el_v		Element to Vertices Map
   *             coordinates	Coordinates of Vertices
   *             dim            Dimension of Problem
   *             v_per_elm      Number of vertices per element
   *
   *    Output:  el_vol         Area/Volume of each element
   */
	
  /* Loop Index */
  INT i,j_a,j_b,j,jcnt; 
  REAL x2,x3,x4,y2,y3,y4,z2,z3,z4;

  /* Coordinates of nodes on elements */
  REAL* x = (REAL *) calloc(v_per_elm,sizeof(REAL));		
  REAL* y = (REAL *) calloc(v_per_elm,sizeof(REAL));
  REAL* z=NULL;
  INT nelm = el_v.row;

  if(dim==2) {
    for (i=0;i<nelm;i++) {
      j_a = el_v.IA[i];
      j_b = el_v.IA[i+1]-1;
      jcnt=0;
      for (j=j_a; j<=j_b; j++) {
	x[jcnt] = cv->x[el_v.JA[j-1]-1];
	y[jcnt] = cv->y[el_v.JA[j-1]-1];
	jcnt++;
      }
      el_vol[i] = 0.5*fabs(y[0]*(x[1]-x[2]) + y[1]*(x[2]-x[0]) + \
			   y[2]*(x[0]-x[1]));
    }
  } else if(dim==3) {
    z = (REAL *) calloc(v_per_elm,sizeof(REAL));
	
    for (i=0;i<nelm;i++) {
      j_a = el_v.IA[i];
      j_b = el_v.IA[i+1]-1;
      jcnt=0;
      for (j=j_a; j<=j_b; j++) {
	x[jcnt] = cv->x[el_v.JA[j-1]-1];
	y[jcnt] = cv->y[el_v.JA[j-1]-1];
	z[jcnt] = cv->z[el_v.JA[j-1]-1];
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
	
      el_vol[i] = (REAL) (fabs(x2*(y3*z4-y4*z3) - y2*(x3*z4-x4*z3) + \
			       z2*(x3*y4-x4*y3)))/6.0;
    }	
  }
	
  if(x) free(x);
  if(y) free(y);
  if(z) free(z);
  return;
}
/********************************************************************************/


