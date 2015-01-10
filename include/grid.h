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

/**
 * \struct coordinates
 * \brief Returns coordinates of nodes
 */
typedef struct coordinates{

  //! x values
  dvector x;

  //! y values
  dvector y;

  //! z values (if in 3D)
  dvector z;

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
    
  //! number of edges
  INT nedge;
    
  //! number of faces (in 2D faces=edges)
  INT nface;
    
  //! number of vertices on boundary
  INT nbv;

  //! number of edges on boundary
  INT nbedge;

  //! number of faces on boundary
  INT nbface;

  //! element to vertex map (CSR Format)
  iCSRmat el_v;                

  //! element to edge map (CSR Format)
  iCSRmat el_ed;

  //! element to face map (CSR Format)
  iCSRmat el_f;

  //! edge to vertex map (CSR Format)
  iCSRmat ed_n;

  //! face to vertex map (CSR Format)
  iCSRmat f_n;

  //! element volumes (areas in 2D)
  dvector el_vol;

  //! barycenter of element
  dvector el_mid;

  //! edge lengths 
  dvector ed_len;

  //! tangent vectors on edges
  dvector ed_tau;

  //! midpoint of edges
  dvector ed_mid;

  //! area of faces
  dvector f_area;

  //! normal vector on face
  dvector f_norm;

  //! coordinates of vertices
  coordinates cv;

  

  
  coordinates cv;               /* Coordinates of vertices */
  
  CSRinc ed_n;                  /* Edge to Vertex Map */
  INT* ed_bdry=NULL;	       	/* Indicates whether an edge is a boundary */
  INT* v_bdry=NULL;	       	/* Indicates whether a vertex is a boundary */
  INT* n_bdry=NULL; 	      	/* Indicates whether a node is a boundary */
  INT* bdry_v=NULL; 		/* Indicates nodes of an edge on boundary */
    
} dCSRmat; /**< Sparse matrix of REAL type in CSR format */


#endif
