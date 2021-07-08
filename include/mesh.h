//
//  mesh.h
//
//
//  Created by James Adler, Xiaozhe Hi, and Ludmil Zikatanov 2015-01-09.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _mesh_h
#define _mesh_h

#include "sparse.h"
#include "vec.h"

/**
 * \struct coords
 * \brief Returns coordinates of nodes
 */
typedef struct coords{

  //! coordinate values (ordered x(0) y(0) z(0) x(1) y(1) z(1) etc...))
  REAL* x;

  //! Size of arrays (number of nodes)
  INT n;

  //! dimension
  INT dim;

} coords;

/**
 * \struct coordinates
 * \brief Returns coordinates of nodes
 */
typedef struct coordinates{

  //! x values
  REAL* x;

  //! y values
  REAL* y;

  //! z values (if in 3D)
  REAL* z;

  //! Size of arrays (number of nodes)
  INT n;

} coordinates;


/**
 * \struct mesh_struct
 * \brief Builds a triangular/tetrahedral mesh, including all its
 * properties and mappings between vertices, edges, and faces
 */
typedef struct mesh_struct{

  //! dimension
  INT dim;

  //! Number of elements
  INT nelm;

  //! number of vertices
  INT nv;

  //! number of vertices per element
  INT v_per_elm;

  //! number of edges
  INT nedge;

  //! number of edges per element
  INT ed_per_elm;

  //! number of faces (in 2D faces=edges)
  INT nface;

  //! number of faces per element
  INT f_per_elm;

  //! Number of connected regions in domain - usually 1
  INT nconn_reg;

  //! Number of connected boundaries - usually 1
  INT nconn_bdry;

  //! Array indicating which vertices are on which component
  INT* v_component;

  //! coordinates of vertices
  coordinates* cv;

  //! number of vertices on boundary
  INT nbv;

  //! number of edges on boundary
  INT nbedge;

  //! number of faces on boundary
  INT nbface;

  //! element to vertex map (CSR Format)
  iCSRmat* el_v;

  //! element to edge map (CSR Format)
  iCSRmat* el_ed;

  //! element to face map (CSR Format)
  iCSRmat* el_f;

  //! edge to vertex map (CSR Format)
  iCSRmat* ed_v;

  //! face to vertex map (CSR Format)
  iCSRmat* f_v;

  //! face to edge map (CSR Format)
  iCSRmat* f_ed;

  //! element volumes (areas in 2D)
  REAL* el_vol;

  //! barycenter of element
  REAL* el_mid;

  //! edge lengths
  REAL* ed_len;

  //! tangent vectors on edges
  REAL* ed_tau;

  //! midpoint of edges
  REAL* ed_mid;

  //! area of faces
  REAL* f_area;

  //! normal vector on face
  REAL* f_norm;

  //! midpoint of face
  REAL* f_mid;

  //! indicates a flag for vertex such as whether a vertex is on boundary
  INT* v_flag;

  //! indicates a flag for edge such as whether an edge is on boundary
  INT* ed_flag;

  //! indicates a flag for face such as whether a face is on boundary
  INT* f_flag;

  //! indicates a flag for element such as in what domain the element is.
  INT* el_flag;

  //! extra double working array for random things
  REAL* dwork;

  //! extra int working array for random things
  INT* iwork;

} mesh_struct;


#endif
