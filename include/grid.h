//
//  grid.h
//  
//
//  Created by Adler, James on 1/9/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _grid_h
#define _grid_h

#include "sparse.h"
#include "vec.h"

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
 * \struct trimesh
 * \brief Builds a triangular/tetrahedral mesh, including all its properties and mappings between vertices, edges, and faces
 */
typedef struct trimesh{
    
  //! dimension
  INT dim;

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

  //! Number of holes in domain - usually 0
  INT nholes;

  //! coordinates of vertices
  REAL* x;
  REAL* y;
  REAL* z;
    
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

  //! indicates whether a vertex is on boundary
  INT* v_bdry;

  //! indicates whether an edge is on boundary
  INT* ed_bdry;

  //! indicates whether a face is on boundary
  INT* f_bdry;
    
} trimesh; 


#endif
