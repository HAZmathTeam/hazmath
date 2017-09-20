/*! \file src/grid/grid.c
 *
 * \brief Obtains generic routines for initializing, building, and creating mesh data.
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 1/9/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 11/14/2016
 */

#include "hazmath.h"

/******************************************************************************/
/*!
 * \fn void initialize_mesh(trimesh* mesh)
 *
 * \brief Initializes all components of the structure trimesh.
 *
 * \return mesh     Struct for Mesh
 *
 */
void initialize_mesh(trimesh* mesh) 
{
  mesh->dim = -666;
  mesh->nelm = -666;
  mesh->cv = NULL;
  mesh->f_per_elm = -666;
  mesh->el_v = NULL;
  mesh->el_f = NULL;
  mesh->v_per_elm = -666;
  mesh->f_v = NULL;
  mesh->f_ed = NULL;
  mesh->nv = -666;
  mesh->nedge = -666;
  mesh->ed_per_elm = -666;
  mesh->nface = -666;
  mesh->nconn_reg = -666;
  mesh->nconn_bdry = -666;
  mesh->v_component = NULL;
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
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void creategrid_fread(FILE *gfid,INT file_type,trimesh* mesh)
 *
 * \brief Creates grid by reading in from file.
 *
 * \param gfid      Grid FILE ID
 * \param file_type Type of File Input: 0 - haz format; 1 - vtk format
 *
 * \return mesh     Struct for Mesh
 *
 */
void creategrid_fread(FILE *gfid,INT file_type,trimesh* mesh) 
{
  // Initialize mesh for read in.
  initialize_mesh(mesh);
  mesh->el_v = malloc(sizeof(struct iCSRmat));

  // READ FILE
  if(file_type==0) {
    read_grid_haz(gfid,mesh);
  } else if(file_type==1) {
    read_grid_vtk(gfid,mesh);
  } else {
    fprintf(stderr, \
            "Unknown mesh file type, %d. Try using vtk format. -Exiting\n", \
            file_type);
    exit(255);
  }

  // Build rest of mesh
  build_mesh(mesh);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void build_mesh(trimesh* mesh)
 *
 * \brief Builds the mesh.
 * \note Assumes basic info has already been read in or created.
 *       including dim, el_v, cv, nv, and nbv.
 *
 * \return mesh     Struct for Mesh
 *
 */
void build_mesh(trimesh* mesh) 
{
  // Flag for errors
  SHORT status;

  // Loop Indices
  INT i,j,k;

  // Get data already built
  INT dim = mesh->dim;
  INT nelm = mesh->nelm;
  INT nv = mesh->nv;
  INT nconn_bdry = mesh->nconn_bdry;
  INT v_per_elm = dim+1;

  // Stuff for Edges and Faces.  None for 1D
  if(dim==2 || dim==3) {
    INT ed_per_elm = 3*(dim-1);
    INT f_per_elm = v_per_elm;

    // Edge to Vertex Map
    INT nedge = 0;
    iCSRmat ed_v = get_edge_v(&nedge,mesh->el_v);

    /* Element to Edge Map */
    iCSRmat el_ed = get_el_ed(mesh->el_v,&ed_v);

    /* Get Edge Stats such as edge lengths, midpoints, and tangent vectors */
    REAL* ed_len = (REAL *) calloc(nedge,sizeof(REAL));
    REAL* ed_tau = (REAL *) calloc(nedge*dim,sizeof(REAL));
    REAL* ed_mid = (REAL *) calloc(nedge*dim,sizeof(REAL));
    edge_stats_all(ed_len,ed_tau,ed_mid,mesh->cv,&ed_v,dim);

    /* Figure out Number of Faces */
    INT nholes = nconn_bdry - 1;
    INT nface = 0;
    INT euler = -10; // Euler Number should be 1 for simplices
    if(dim==2) {
      nface = nedge;
      euler = nv - nedge + nelm + nholes;
    } else if (dim==3) {
      nface = 1 + nedge-nv+nelm;
      nface = nface + nholes; // add number of holes!
      euler = nv - nedge + nface - nelm - nholes;
    } else {
      status = ERROR_DIM;
      check_error(status, __FUNCTION__);
    }
    if(euler!=1) {
      printf("ERROR HAZMATH DANGER: in function %s, your simplices are all messed up.  "
             "Euler Characteristic doesn't equal 1+nholes!\teuler=%d\tnholes=%d.\n\n",__FUNCTION__,euler,nholes);
      exit(ERROR_DIM);
    }

    // In order to get some of the face maps, we need to know how the faces are
    // ordered on each element
    // This is done by the same ordering as the nodes, connecting the opposite
    // face with each node.
    INT* fel_order= (INT *) calloc(f_per_elm*dim,sizeof(INT));
    get_face_ordering(v_per_elm,dim,f_per_elm,fel_order);

    // Next get the element to face map, the face to vertex map, and the face
    //  boundary information.
    iCSRmat f_v = icsr_create(nface,nv,nface*dim);
    iCSRmat f_ed = icsr_create(nface,nedge,nface*(2*dim-3));
    INT* f_bdry = (INT *) calloc(nface,sizeof(INT));
    INT nbface;
    mesh->el_f = malloc(sizeof(struct iCSRmat));
    get_face_maps(mesh->el_v,v_per_elm,&ed_v,nface,dim,f_per_elm,mesh->el_f, \
                  f_bdry,&nbface,&f_v,&f_ed,fel_order);

    // In case v_bdry has different types of boundaries, match the
    // face boundary to the same:
    /*(ltz) correction: we need to get the max boundary code for a
      face, because we can say that a face is on the Neumann boundary
      (open set) iff the face is on the boundary and at least one of
      its vertices is on the Neumann boundary */
    INT maxflag,jk;
    for(i=0;i<nface;i++) {
      if(f_bdry[i]==1) {	
        j = f_v.IA[i]-1;
        k = f_v.JA[j]-1;
	maxflag = mesh->v_bdry[k];
	for(jk=f_v.IA[i];jk<f_v.IA[i+1]-1;jk++){
	  if(mesh->v_bdry[k] > maxflag) maxflag = mesh->v_bdry[k];
	}
        f_bdry[i] = maxflag; /*(old)mesh->v_bdry[k];*/
	/* 	(debugging) fprintf(stdout,"\nyyy%i %i",i,f_bdry[i]);fflush(stdout); */
	/* (debugging) } else { */
	/*(debugging) 	fprintf(stdout,"\nxxx%i %i",i,f_bdry[i]);fflush(stdout); */
      }
    }
    
    // Next get the face data, such as areas, midpoints, and normal vectors
    REAL* f_area = (REAL *) calloc(nface,sizeof(REAL));
    REAL* f_mid = (REAL *) calloc(nface*dim,sizeof(REAL));
    REAL* f_norm = (REAL *) calloc(nface*dim,sizeof(REAL));
    // Need to update mesh for a few things first
    mesh->f_per_elm = f_per_elm;

    // Get Statistics of the faces (midpoint, area, normal vector, ordering, etc.)
    face_stats(f_area,f_mid,f_norm,&f_v,mesh);

    // Reorder appropriately
    sync_facenode(&f_v,f_norm,mesh);
    mesh->f_v = malloc(sizeof(struct iCSRmat));
    mesh->f_ed = malloc(sizeof(struct iCSRmat));
    *(mesh->f_v) = f_v;
    *(mesh->f_ed) = f_ed;

    /* Get Boundary Edges */
    INT nbedge=0;
    INT* ed_bdry = (INT *) calloc(nedge,sizeof(INT));
    isboundary_ed(&f_ed,&ed_v,nedge,nface,f_bdry,mesh->v_bdry,&nbedge,ed_bdry);

    // Assign components to the mesh
    mesh->nedge = nedge;
    mesh->ed_per_elm = ed_per_elm;
    mesh->nface = nface;
    mesh->nbedge = nbedge;
    mesh->nbface = nbface;
    mesh->el_ed = malloc(sizeof(struct iCSRmat));
    *(mesh->el_ed) = el_ed;
    mesh->ed_v = malloc(sizeof(struct iCSRmat));
    *(mesh->ed_v) = ed_v;
    mesh->ed_len = ed_len;
    mesh->ed_tau = ed_tau;
    mesh->ed_mid = ed_mid;
    mesh->f_area = f_area;
    mesh->f_mid = f_mid;
    mesh->f_norm = f_norm;
    mesh->ed_bdry = ed_bdry;
    mesh->f_bdry = f_bdry;

    if(fel_order) free(fel_order);
  } else if (dim==1) {
    mesh->nedge = 0;
    mesh->ed_per_elm = 0;
    mesh->nface = 0;
    mesh->nbedge = 0;
    mesh->nbface = 0;
  } else {
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  // Finally get volumes/areas of elements and the midpoint (barycentric)
  REAL* el_mid = (REAL *) calloc(nelm*dim,sizeof(REAL));
  REAL* el_vol = (REAL *) calloc(nelm,sizeof(REAL));
  get_el_vol(el_vol,mesh->el_v,mesh->cv,dim,v_per_elm);
  get_el_mid(el_mid,mesh->el_v,mesh->cv,dim);
  mesh->el_vol = el_vol;
  mesh->el_mid = el_mid;

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn struct coordinates *allocatecoords(INT ndof,INT mydim)
 *
 * \brief Allocates memory and properties of coordinates struct
 *
 * \param ndof      Total number of coordinates
 * \param mydim     Dimension of coordinates
 *
 * \return A        Coordinate struct
 *
 */
struct coordinates *allocatecoords(INT ndof,INT mydim)
{  
  // Flag for errors
  SHORT status;

  struct coordinates *A = malloc(sizeof(struct coordinates));
  assert(A != NULL);

  switch (mydim)
  {
  case 1:
    A->x = (REAL *) calloc(ndof,sizeof(REAL));
    A->y = NULL;
    A->z = NULL;
    break;
  case 2:
    A->x = (REAL *) calloc(ndof,sizeof(REAL));
    A->y = (REAL *) calloc(ndof,sizeof(REAL));
    A->z = NULL;
    break;
  case 3:
    A->x = (REAL *) calloc(ndof,sizeof(REAL));
    A->y = (REAL *) calloc(ndof,sizeof(REAL));
    A->z = (REAL *) calloc(ndof,sizeof(REAL));
    break;
  default:
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }

  A->n = ndof;
  
  return A;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void free_coords(coordinates* A)
 *
 * \brief Frees memory of arrays of coordinate struct
 *
 * \param A         Pointer to coordinates struct to be freed
 *
 */
void free_coords(coordinates* A)
{
  if (A==NULL) return;

  if(A->x) {
    free(A->x);
    A->x = NULL;
  }

  if(A->y) {
    free(A->y);
    A->y = NULL;
  }

  if(A->z) {
    free(A->z);
    A->z = NULL;
  }
  
  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_coords(FILE* fid,coordinates *c)
 *
 * \brief Dump the coordinate data to file for plotting purposes
 *
 * \param fid         Output FILE ID
 * \param c           Coordinate struct
 *
 * \return coord.dat  Coordinates of node in format: coord(n,dim)
 *
 */
void dump_coords(FILE* fid,coordinates *c) 
{
  INT i;

  if (c->z) {
    for (i=0; i<c->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\t%25.16e\n",c->x[i],c->y[i],c->z[i]);
    }
  } else if(c->y) {
    for (i=0; i<c->n; i++) {
      fprintf(fid,"%25.16e\t%25.16e\n",c->x[i],c->y[i]);
    }
  } else {
    for (i=0; i<c->n; i++) {
      fprintf(fid,"%25.16e\n",c->x[i]);
    }
  }
  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void free_mesh(trimesh* mesh)
 *
 * \brief Frees memory of arrays of mesh struct
 *
 * \param mesh     Pointer to mesh struct to be freed
 *
 */
void free_mesh(trimesh* mesh)
{  
  if(mesh==NULL) return;

  if (mesh->cv){
    free_coords(mesh->cv);
    free(mesh->cv);
    mesh->cv = NULL;
  }
  if (mesh->el_v) {
    icsr_free(mesh->el_v);
    free(mesh->el_v);
    mesh->el_v = NULL;
  }
  if(mesh->el_ed) {
    icsr_free(mesh->el_ed);
    free(mesh->el_ed);
    mesh->el_ed = NULL;
  }

  if(mesh->el_f) {
    icsr_free(mesh->el_f);
    free(mesh->el_f);
    mesh->el_f = NULL;
  }
  
  if(mesh->ed_v) {
    icsr_free(mesh->ed_v);
    free(mesh->ed_v);
    mesh->ed_v = NULL;
  }

  if(mesh->f_v) {
    icsr_free(mesh->f_v);
    free(mesh->f_v);
    mesh->f_v = NULL;
  }

  if(mesh->f_ed) {
    icsr_free(mesh->f_ed);
    free(mesh->f_ed);
    mesh->f_ed = NULL;
  }

  if(mesh->el_vol) {
    free(mesh->el_vol);
    mesh->el_vol = NULL;
  }
  if(mesh->el_mid) {
    free(mesh->el_mid);
    mesh->el_mid = NULL;
  }

  if(mesh->ed_len) {
    free(mesh->ed_len);
    mesh->ed_len = NULL;
  }

  if(mesh->ed_tau) {
    free(mesh->ed_tau);
    mesh->ed_tau = NULL;
  }

  if(mesh->ed_mid) {
    free(mesh->ed_mid);
    mesh->ed_mid = NULL;
  }

  if(mesh->f_area) {
    free(mesh->f_area);
    mesh->f_area = NULL;
  }

  if(mesh->f_norm) {
    free(mesh->f_norm);
    mesh->f_norm = NULL;
  }

  if(mesh->f_mid) {
    free(mesh->f_mid);
    mesh->f_mid = NULL;
  }

  if(mesh->v_bdry) {
    free(mesh->v_bdry);
    mesh->v_bdry = NULL;
  }
  
  if(mesh->ed_bdry) {
    free(mesh->ed_bdry);
    mesh->ed_bdry = NULL;
  }

  if(mesh->f_bdry) {
    free(mesh->f_bdry);
    mesh->f_bdry = NULL;
  }

  if(mesh->v_component) {
    free(mesh->v_component);
    mesh->v_component = NULL;
  }
  
  return;
}
/******************************************************************************/
