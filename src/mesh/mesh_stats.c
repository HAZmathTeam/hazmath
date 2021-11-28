/*! \file src/mesh/mesh_stats.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 1/9/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \brief Obtains routines for generating properties of the mesh,
 *  including all incidence matrices, element volumes, etc.
 * \note updated by James Adler 07/25/2018
 *
 */

#include "hazmath.h"

/******************************************************************************/
/*!
 * \fn iCSRmat convert_elmnode(INT *element_vertex,INT nelm,INT nv,INT nve)
 *
 * \brief Convert the input element to vertex map into sparse matrix form
 *
 * \param nelm	                   Number of elements
 * \param nv		                   Number of vertices
 * \param nve	                       Number of vertices per element
 * \param element_vertex(nelm,nve)   Each row is an element, each column is
 *                                   the corresponding vertices for that element
 *
 * \return el_v:	                   Element to vertex map in CSR format.  Rows
 *                                   are all elements, columns are all vertices.
 *
 * \note Not currently used
 *
 */
iCSRmat convert_elmnode(INT *element_vertex,INT nelm,INT nv,INT nve)
{
  // Loop indices
  INT i,j,k;

  iCSRmat el_v;

  /* Exactly nve vertices per row */
  if ( nelm > 0 ) {
    el_v.IA = (INT *)calloc(nelm+1, sizeof(INT));
  }
  else {
    el_v.IA = NULL;
  }
  for(i=0;i<nelm+1;i++) {
    el_v.IA[i] = nve*i;
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
/*!
 * \fn iCSRmat get_edge_v(INT* nedge,iCSRmat* el_v)
 *
 * \brief Gets the Edge to Vertex mapping in CSR Fromat (Should be dim independent)
 *        This is used to get Element to Edge map, but also for determining
 *        things such as edge_length and the tangent vectors of edges
 *        The routine also counts the number of edges.
 *
 * \param nedge	  Number of edges
 * \param el_v	  Element to vertex map in CSR format.
 *
 * \return ed_v:	  Edge to vertex map in CSR format.
 *
 */
iCSRmat get_edge_v(INT* nedge,iCSRmat* el_v)
{
  INT nv = el_v->col;
  INT i,j,col_b,col_e,icntr,jcntr; /* Loop indices and counters */

  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans(el_v,&v_el);

  /* Create Vertex to Vertex Map by v_el*el_v */
  iCSRmat v_v;
  icsr_mxm_symb(&v_el,el_v,&v_v);
  INT* iv_v = v_v.IA;
  INT* jv_v = v_v.JA;

  // Get the number of edges.
  // This is the number of nonzeroes in the upper triangular
  // part of vertex-vertex map (not including diagonal).
  INT ned = (INT) (v_v.nnz - v_v.row)/2;
  *nedge = ned;

  iCSRmat ed_v;
  if ( ned > 0 ) {
    ed_v.IA = (INT *)calloc(ned+1, sizeof(INT));
  }
  else {
    ed_v.IA = NULL;
  }
  if ( nv > 0 ) {
    ed_v.JA = (INT *)calloc(2*ned, sizeof(INT));
  }
  else {
    ed_v.JA = NULL;
  }

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
      if(jv_v[j]>i) {
        ed_v.IA[icntr] = jcntr;
        ed_v.JA[jcntr] = jv_v[j];
        ed_v.JA[jcntr+1] = i;
        jcntr=jcntr+2;
        icntr++;
      }
    }
  }
  ed_v.IA[icntr] = jcntr;

  /* Free Node Node */
  icsr_free(&v_v);
  icsr_free(&v_el);
  ed_v.val=NULL;
  ed_v.row = ned;
  ed_v.col = nv;
  ed_v.nnz = 2*ned;

  return ed_v;
}
/*******************************************************************************/

/*******************************************************************************/
/*!
 * \fn void boundary_f_ed(iCSRmat* f_ed,iCSRmat* ed_v,INT nedge,INT nface,INT *f_flag,INT *v_flag,INT *nbedge,INT *ed_flag)
 *
 * \brief Counts the number of boundary edges and indicates whether an edge and a face is
 *        a boundary, using the f_flag and v_flag maps.
 *
 * \param f_ed                     Face to Edge Map
 * \param ed_v	                   Edge to vertex map in CSR format.
 * \param nedge	                   Number of edges
 * \param nface	                   Number of faces
 * \param f_flag		           Boundary flag array for faces
 * \param v_bdry                   Boundary flag array for vertices
 *
 * \return nbedge                  Number of boundary edges
 * \return ed_flag                 Boundary flag array for edges
 *
 */
void boundary_f_ed(iCSRmat* f_ed,iCSRmat* ed_v,INT nedge, INT nface,INT *f_flag,INT *v_flag,INT *nbedge,INT *ed_flag,INT dim)
{
  INT i,j,ed,v1,v2,jcntr; /* Loop indices and counters */

  // Edge flag is max of its two vertex flags (assuming both are on boundary)
  // Face flag is max of its edge flags
  // This way if one vertex is Neumann the whole edge/face will be Neumann
  INT ed_per_f = (dim*(dim-1))/2;
  INT* ed_on_f = (INT *) calloc(ed_per_f,sizeof(INT));
  INT maxe=0;
  for(i=0; i<nface; i++) {
    if( f_flag[i]==1 ) {
      get_incidence_row(i,f_ed,ed_on_f);
      maxe = -666;
      for(j=0; j<ed_per_f; j++) {
        ed = ed_on_f[j];
        v1 = ed_v->JA[ed_v->IA[ed]];
        v2 = ed_v->JA[ed_v->IA[ed]+1];
        ed_flag[ed] = MAX(v_flag[v1],v_flag[v2]);
        if(ed_flag[ed]>maxe)
          maxe = ed_flag[ed];
      }
      f_flag[i] = maxe;
    }
  }
  if(ed_on_f) free(ed_on_f);

  // Now go back through edges and count how many are on boundary.
  jcntr=0;
  for (i=0; i<nedge; i++) {
    if(ed_flag[i]!=0) {
      jcntr++;
    }
  }

  *nbedge = jcntr;

  return;
}
/*******************************************************************************/

/*******************************************************************************/
/*!
 * \fn iCSRmat get_el_ed(iCSRmat* el_v,iCSRmat* ed_v)
 *
 * \brief Gets the Element to Edge mapping in CSR Fromat (Should be dim independent)
 *
 * \param el_v                       Element to vertex map
 * \param ed_v	                   Edge to vertex map
 *
 * \return el_ed                     Element to edge map
 *
 */
iCSRmat get_el_ed(iCSRmat* el_v,iCSRmat* ed_v)
{
  iCSRmat el_ed;

  // Get transpose of edge to vertex
  iCSRmat v_ed;
  icsr_trans(ed_v,&v_ed);
  icsr_mxm_symb_max(el_v,&v_ed,&el_ed,2);
  //zicsr_mxm_symb_max(el_v,&v_ed,&el_ed,2);
  //fprintf(stdout,"\nXXXX: %d, %d %d\n",el_ed.nnz,el_ed.row,el_ed.col);
  //print_full_mat_int(el_ed.nnz,1,el_ed.JA,"eled");
  icsr_free(&v_ed);

  return el_ed;
}
/*******************************************************************************/

/*******************************************************************************/
/*!
 * \fn void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates *cv,iCSRmat* ed_v,INT dim)
 *
 * \brief Get length, tangent vector (tau), and midpoint of every edge
 *
 * \param cv                         Coordinates of vertices
 * \param ed_v	                   Edge to vertex map
 * \param dim                        Dimension of problem
 *
 * \return ed_len                    Length of each edge
 * \return ed_tau                    Tangent vector or each edge, ordered (tx1,ty1,tz1,tx2,ty2,tz2,...)
 * \return ed_mid                    Midpoint of each edge, ordered (mx1,my1,mz1,mx2,my2,mz2,...)
 *
 *
 */
void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates *cv,iCSRmat* ed_v,INT dim)
{
  // Loop indices
  INT i,jcnt,j,j_a,j_b;

  INT ip[2];
  REAL x[2],y[2],z[2];
  ip[0]=0; ip[1]=0; x[0]=0.0; x[1]=0.0; y[0]=0.0; y[1]=0.0; z[0]=0.0; z[1]=0.0;

  /* Orders nodes from largest to smallest for correct orientation */
  INT ihi1,ilo1;
  INT nedge = ed_v->row;

  for(i=0;i<nedge;i++) {

    /* Find Nodes in given Edge */
    j_a = ed_v->IA[i];
    j_b = ed_v->IA[i+1];
    jcnt = 0;
    for (j=j_a; j<j_b;j++) {
      ip[jcnt] = ed_v->JA[j];
      x[jcnt] = cv->x[ip[jcnt]];
      y[jcnt] = cv->y[ip[jcnt]];
      if (dim==3) {
        z[jcnt] = cv->z[ip[jcnt]];
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
/*!
 * \fn void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order)
 *
 * \brief Gets the face ordering on element
 *
 *               2
 *             /   \
 *            1     0
 *           /       \
 *          0 ---2----1
 *
 * \param el_order                   Number of nodes per element
 * \param dim                        Dimension of problem
 * \param f_order                    Number of faces per element: 3 in 2D, 4 in 3D
 *
 * \return fel_order                 Indicates order of faces on an element
 *
 *
 */
void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order)
{
  // Loop indices
  INT i,j,jcntr,col_b;

  // First order the faces per element.  Start with node 0 and indicate
  // opposite face as face 0.  Then loop over nodes
  for(j=0;j<f_order;j++) {
    for(i=0;i<dim;i++) {
      jcntr = j+i+2;
      col_b = j*dim+i;
      if(jcntr>el_order) {
        fel_order[col_b] = jcntr-el_order-1;
      } else {
        fel_order[col_b] = jcntr-1;
      }
    }
  }
  return;
}
/*******************************************************************************/

/*******************************************************************************/
/*!
 * \fn void get_face_maps(iCSRmat* el_v,INT el_order,iCSRmat* ed_v,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,iCSRmat *f_ed,INT *fel_order)
 *
 * \brief Gets the Element to Face, Face to Node, Face to Boundary mapping in CSR
 *        format as well as face ordering on element (Should be dim independent)
 *
 *               2
 *             /   \
 *            1     0
 *           /       \
 *          0 ---2----1
 *
 * \param el_v                       Element to vertex map
 * \param el_order                   Number of nodes per element
 * \param ed_v                       Edge to vertex map
 * \param nface                      Number of faces
 * \param dim                        Dimension of problem
 * \param f_order                    Number of faces per element: 3 in 2D, 4 in 3D
 * \param fel_order                  Indicates order of faces on an element
 *
 * \return el_f                      Element to face map
 * \return f_bdry                    Binary boundary array for faces
 * \return nbf                       Number of boundary faces
 * \return f_v                       Face to vertex map
 *
 */
void get_face_maps(iCSRmat* el_v,INT el_order,iCSRmat* ed_v,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,iCSRmat *f_ed,INT *fel_order)
{
  // Flag for errors
  SHORT status;

  INT i,j,k,m,p,jk,col_b,icntr,jcntr,kcntr; /* Loop Indices */
  INT ncol1,ncol2,ncol3,ncol4,ncol5,ncol6;  /* Node indices */
  INT el1=-1,el2=-1,el3=-1,el=-1; /* Element indices */
  INT iflag = 0;  /* Marks if we found a face */
  INT nbf = 0; /* hold number of boundary faces */
  INT nd[dim];
  INT f_num=-1;  /* Current face number of element */
  INT nelm = el_v->row;
  INT edpf = 2*dim - 3;

  // We will build face to element map first and then transpose it
  iCSRmat f_el = icsr_create (nface,nelm,nelm*f_order);

  /* Get Transpose of el_v -> v_el */
  iCSRmat v_el;
  icsr_trans(el_v,&v_el);

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
    col_b = el_v->IA[i];
    // Now loop over all faces of element
    for(j=0;j<f_order;j++) {
      // Get appropriate nodes on face
      for(k=0;k<dim;k++) {
        jk = col_b+fel_order[j*dim+k];
        nd[k] = el_v->JA[jk];
      }
      // Next Loop over all elements that touch the first node
      ncol1 = v_el.IA[nd[0]];
      ncol2 = v_el.IA[nd[0]+1];
      for(k=ncol1;k<ncol2;k++) {
        el1 = v_el.JA[k];
        // If not the same element we started in then check other nodes elements
        if(el1!=i) {
          ncol3 = v_el.IA[nd[1]];
          ncol4 = v_el.IA[nd[1]+1];
          for(m=ncol3;m<ncol4;m++) {
            el2 = v_el.JA[m];
            // Our two nodes share the element!  In 2D this is the face!
            // In 3D, we check third guy
            if(el2==el1) {
              if(dim==2) {
                // Mark the element that shares the face as well as the
                // nodes that are shared
                iflag = 1;
                el = el1;
              } else if(dim==3) {
                ncol5 = v_el.IA[nd[2]];
                ncol6 = v_el.IA[nd[2]+1];
                for(p=ncol5;p<ncol6;p++) {
                  el3 = v_el.JA[p];
                  if(el3==el2) {
                    // Mark the element that shares the face as well as
                    // the nodes that are shared
                    iflag = 1;
                    el = el1;
                  }
                }
              } else {
                status = ERROR_DIM;
                check_error(status, __FUNCTION__);
              }
            }
          }
        }
      }
      if(iflag && el>i) {
        f_el.IA[icntr] = jcntr;
        f_el.JA[jcntr] = i;
        f_el.val[jcntr] = j;
        f_el.JA[jcntr+1] = el;
        // Find face number for other element
        find_facenumber(el_v,el,nd,dim,&f_num); 
        f_el.val[jcntr+1] = f_num;
        f_bdry[icntr] = 0;
        f_v->IA[icntr] = kcntr;
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
        f_el.IA[icntr] = jcntr;
        f_el.JA[jcntr] = i;
        f_el.val[jcntr] = j;
        f_bdry[icntr] = 1;
        f_v->IA[icntr] = kcntr;
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

  f_el.IA[icntr] = jcntr;
  f_v->IA[icntr] = kcntr;

  /* Get Transpose of f_el -> el_f */
  icsr_trans(&f_el,el_f);

  /* Transpose face_node and back again to order nodes in increasing order
   * (from global structure) */
  iCSRmat v_f;
  icsr_trans(f_v,&v_f);
  icsr_free(f_v);
  icsr_trans(&v_f,f_v);

  /* Get Face to Edge Map */
  iCSRmat face_ed;
  iCSRmat v_ed;
  icsr_trans(ed_v,&v_ed);
  icsr_mxm_symb_max(f_v,&v_ed,&face_ed,2);
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
/*!
 * \fn void find_facenumber(iCSRmat* el_v,INT elm,INT* nd,INT dim,INT *f_num)
 *
 * \brief Find the face number of the given element, using the element it shares the face with
 *
 * \param el_v                       Element to vertex map
 * \param elm                        Current element we consider
 * \param nd                         Nodes of current face
 * \param dim                        Dimension of problem
 *
 * \return f_num                     Face number for that element
 *
 *
 */
void find_facenumber(iCSRmat* el_v,INT elm,INT* nd,INT dim,INT *f_num)
{
  INT cola,colb,mark,icnt,i,j,mynd;
  INT fntmp=-1;

  cola = el_v->IA[elm];
  colb = el_v->IA[elm+1];
  mark = 0;
  icnt=-1;
  for(i=cola;i<colb;i++) {
    icnt++;
    mynd = el_v->JA[i];
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
/*!
 * \fn void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,iCSRmat *f_v,mesh_struct *mesh)
 *
 * \brief Get area, normal vector, and midpoints for all faces
 *
 * \param mesh                       Mesh struct
 *
 * \return f_area                    Area of each face (length in 2D)
 * \return f_norm                    Normal vector or each face
 * \return f_mid                     Midpoint of each face
 *
 *
 */
void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,iCSRmat *f_v,mesh_struct *mesh)
{
  // Flag for errors
  SHORT status;

  // Loop indices
  INT i,jcnt,j,j_a,j_b;

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
  REAL* zf=NULL;
  REAL dimover=(REAL )1./((REAL )dim);
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
  REAL* dp = (REAL *) calloc(el_order*dim,sizeof(REAL));
  REAL grad_mag,e1x,e1y,e1z,e2x,e2y,e2z;

  /* Get Face to Element Map */
  /* Get Transpose of f_el -> el_f */
  iCSRmat f_el;
  icsr_trans(el_f,&f_el);
  
  // Loop over all Faces
  for(i=0;i<nface;i++) {
    /* Find Vertices in given Face */
    j_a = f_v->IA[i];
    j_b = f_v->IA[i+1];
    jcnt = 0;
    for (j=j_a; j<j_b;j++) {
      ipf[jcnt] = f_v->JA[j];
      xf[jcnt] = cv->x[ipf[jcnt]];
      yf[jcnt] = cv->y[ipf[jcnt]];
      if (dim==3) {
        zf[jcnt] = cv->z[ipf[jcnt]];
      }
      jcnt++;
    }

    // Find Corresponding Elements and order in element
    // Also picks correct opposite node to form vector
    // normal vectors point from lower number element to higher one
    // or outward from external boundary
    j_a = f_el.IA[i];
    j_b = f_el.IA[i+1];
    jcnt=0;
    for (j=j_a; j<j_b; j++) {
      notbdry = j_b-j_a-1;
      ie[jcnt] = f_el.JA[j];
      op_n[jcnt] = f_el.val[j];
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
    j_a = el_v->IA[myel];
    j_b = el_v->IA[myel+1];
    jcnt=0;
    for(j=j_a;j<j_b;j++) {
      myel_n[jcnt] = el_v->JA[j];
      jcnt++;
    }
    myx[0] = cv->x[myel_n[myopn]];
    myx[1] = cv->y[myel_n[myopn]];
    if(dim==3) {
      myx[2] = cv->z[myel_n[myopn]];
    }

    /* Compute Area (length if 2D) and get midpt of face */
    if(dim==2) {
      f_area[i] = sqrt(pow(fabs(xf[1]-xf[0]),2)+pow(fabs(yf[1]-yf[0]),2));
      f_mid[i*dim] = (xf[0]+xf[1])*dimover;
      f_mid[i*dim+1] = (yf[0]+yf[1])*dimover;
    } else if(dim==3) {
      e1x = xf[1]-xf[0];
      e1y = yf[1]-yf[0];
      e1z = zf[1]-zf[0];
      e2x = xf[2]-xf[0];
      e2y = yf[2]-yf[0];
      e2z = zf[2]-zf[0];
      f_area[i] = 0.5*sqrt(pow(e1y*e2z-e2y*e1z,2)+pow(e1z*e2x-e2z*e1x,2)+ \
                           pow(e1x*e2y-e2x*e1y,2));
      f_mid[i*dim] = (xf[0]+xf[1]+xf[2])*dimover;
      f_mid[i*dim+1] = (yf[0]+yf[1]+yf[2])*dimover;
      f_mid[i*dim+2] = (zf[0]+zf[1]+zf[2])*dimover;
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }

    // Compute Normal Vectors based on opposite node
    // Get Linear Basis Functions for particular element
    PX_H1_basis(p,dp,myx,myel_n,1,mesh);
    grad_mag = dp[myopn*dim]*dp[myopn*dim]+dp[myopn*dim+1]*dp[myopn*dim+1];
    if(dim==3) {
      grad_mag += dp[myopn*dim+2]*dp[myopn*dim+2];
    }
    grad_mag = -sqrt(grad_mag);
    f_norm[i*dim] = dp[myopn*dim]/grad_mag;
    f_norm[i*dim+1] = dp[myopn*dim+1]/grad_mag;
    if(dim==3) {
      f_norm[i*dim+2] = dp[myopn*dim+2]/grad_mag;
    }
  }

  icsr_free(&f_el);
  if(ipf) free(ipf);
  if(xf) free(xf);
  if(yf) free(yf);
  if(myx) free(myx);
  if(p) free(p);
  if(dp) free(dp);
  if(op_n) free(op_n);
  if(ie) free(ie);
  if(myel_n) free(myel_n);
  if(dim==3) {
    if(zf) free(zf);
  }

  return;
}
/********************************************************************************/

/********************************************************************************/
/*!
 * \fn void sync_facenode(iCSRmat *f_v,REAL* f_norm,mesh_struct *mesh)
 *
 * \brief Reorder the Face-Node mapping so it has positive orientation with
 *        respect to the face's normal vector
 *
 * \param mesh                       Mesh struct
 * \param f_norm                     Normal vector or each face
 *
 * \return f_v                       Reordered face to vertex map
 *
 */
void sync_facenode(iCSRmat *f_v,REAL* f_norm,mesh_struct *mesh)
{
  // Loop indices
  INT i,j;

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
      rowa = f_v->IA[i];
      rowb = f_v->IA[i+1];
      jcnt=0;
      for(j=rowa;j<rowb;j++) {
        nd = f_v->JA[j];
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
      rowa = f_v->IA[i];
      rowb = f_v->IA[i+1];
      jcnt=0;
      for(j=rowa;j<rowb;j++) {
        nd = f_v->JA[j];
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
/*!
 * \fn void get_el_mid(REAL *el_mid,iCSRmat* el_v,coordinates *cv,INT dim)
 *
 * \brief Compute the midpoints of a triangluar element using the vertices
 *
 * \param el_v                   Element to vertex map
 * \param cv                     Coordinates of vertices
 *
 * \return el_mid                Midpoint of elements
 *
 *
 */
void get_el_mid(REAL *el_mid,iCSRmat* el_v,coordinates *cv,INT dim)
{
  // Flag for errors
  SHORT status;

  INT i,j,cnt,nd,acol,bcol; /* Loop Index */
  INT nelm = el_v->row;

  if (dim==1) {
    for (i=0; i<nelm; i++) {
      acol = el_v->IA[i];
      bcol = el_v->IA[i+1];
      cnt=0;
      el_mid[i]=0;
      for (j=acol; j<bcol; j++) {
        nd = el_v->JA[j];
        el_mid[i] += cv->x[nd];
        cnt++;
      }
      el_mid[i]=el_mid[i]/2.0;
    }
  } else if (dim==2) {
    for (i=0; i<nelm; i++) {
      acol = el_v->IA[i];
      bcol = el_v->IA[i+1];
      cnt=0;
      el_mid[i*dim]=0;
      el_mid[i*dim+1]=0;
      for (j=acol; j<bcol; j++) {
        nd = el_v->JA[j];
        el_mid[i*dim] += cv->x[nd];
        el_mid[i*dim+1] += cv->y[nd];
        cnt++;
      }
      el_mid[i*dim]=el_mid[i*dim]/3.0;
      el_mid[i*dim+1]=el_mid[i*dim+1]/3.0;
    }
  } else if (dim==3) {
    for (i=0; i<nelm; i++) {
      acol = el_v->IA[i];
      bcol = el_v->IA[i+1];
      cnt=0;
      el_mid[i*dim]=0;
      el_mid[i*dim+1]=0;
      el_mid[i*dim+2]=0;
      for (j=acol; j<bcol; j++) {
        nd = el_v->JA[j];
        el_mid[i*dim] += cv->x[nd];
        el_mid[i*dim+1] += cv->y[nd];
        el_mid[i*dim+2] += cv->z[nd];;
        cnt++;
      }
      el_mid[i*dim]=0.25*el_mid[i*dim];
      el_mid[i*dim+1]=0.25*el_mid[i*dim+1];
      el_mid[i*dim+2]=0.25*el_mid[i*dim+2];
    }
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  return;
}
/********************************************************************************/

/********************************************************************************/
/*!
 * \fn void get_el_vol(REAL *el_vol,iCSRmat *el_v,coordinates *cv,INT dim,INT v_per_elm)
 *
 * \brief Compute the area/volume of ALL triangluar/tetrahedral elements using the vertices
 *
 * \param el_v                   Element to vertex map
 * \param cv                     Coordinates of vertices
 * \param dim                    Dimension of problem
 * \param v_per_elm              Number of vertices per element
 *
 * \return el_vol                Area/Volume of each element
 *
 *
 */
void get_el_vol(REAL *el_vol,iCSRmat *el_v,coordinates *cv,INT dim,INT v_per_elm)
{
  //  old_get_el_vol(el_vol,el_v,cv,dim,v_per_elm);
  //  return;
  // Loop Indices
  INT i,j_a,j_b,j,k,l,jcnt,lv,lnv;

  REAL x2,x3,x4,y2,y3,y4,z2,z3,z4;

  /* Coordinates of nodes on elements */
  REAL* x = (REAL *) calloc(dim*v_per_elm,sizeof(REAL));
  REAL *y=NULL, *z=NULL;
  INT nelm = el_v->row;
  INT nv = el_v->col;
  for (i=0;i<nelm;i++) {
    j_a = el_v->IA[i];
    j_b = el_v->IA[i+1];
    for (l=0;l<dim;l++){
      jcnt=0;
      lv=l*v_per_elm;
      lnv=l*nv;
      for (j=j_a; j<j_b; j++) {
        k=el_v->JA[j];
        x[lv+jcnt] = cv->x[lnv+k];
        jcnt++;
      }
    }
    if(dim==1) {
      el_vol[i]=fabs(x[1]-x[0]);
      //      fprintf(stdout,"dim=%d; diff=%e\n",dim,el_voli-el_vol[i]);
    } else if(dim==2) {
      y=x+v_per_elm;
      el_vol[i] = 0.5*fabs(y[0]*(x[1]-x[2]) + y[1]*(x[2]-x[0]) +	\
			   y[2]*(x[0]-x[1]));
      //      fprintf(stdout,"dim=%d; diff=%e\n",dim,el_voli-el_vol[i]);
    } else if(dim==3) {
      y=x+v_per_elm;z=y+v_per_elm;
      x2 = x[1]-x[0];
      x3 = x[2]-x[0];
      x4 = x[3]-x[0];
      y2 = y[1]-y[0];
      y3 = y[2]-y[0];
      y4 = y[3]-y[0];
      z2 = z[1]-z[0];
      z3 = z[2]-z[0];
      z4 = z[3]-z[0];
      el_vol[i] = fabs(x2*(y3*z4-y4*z3) - y2*(x3*z4-x4*z3) +	\
		       z2*(x3*y4-x4*y3))/6.0;
    }
  }

  if(x) free(x);

  return;
}
/**************************************************************************/
