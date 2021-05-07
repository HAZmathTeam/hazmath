//
//  fem.h
//
//
//  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2015-02-01
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _fem_h
#define _fem_h

#include "sparse.h"
#include "vec.h"
#include "mesh.h"

/* markers for boundary conditions:
   in the mesh structure these will be the values of the array
   BOUNDARY FACES ARE MARKED WITH
   mesh_struct.f_bdry[i]=0 then (i) is an interior face.
   1 <= mesh_struct.f_bdry[i] <= 16 (i) is on the DIRICHLET boundary;
   17 <= mesh_struct.f_bndry[i] <=32 (i) is  on the NEUMANN boundary;
   33 <= mesh_struct.f_bndry[i] <=64 (i) is  on the ROBIN boundary;
*/
#define MARKER_DIRICHLET 1
#define MARKER_NEUMANN  17
#define MARKER_ROBIN  33
#define MARKER_BOUNDARY_NO  65

/**
 * \struct qcoordinates
 * \brief Returns coordinates of quadrature nodes
 * \note TODO - will remove this and replace with the structure qcoords below.
 */
typedef struct qcoordinates{

  //! x values
  REAL* x;

  //! y values
  REAL* y;

  //! z values (if in 3D)
  REAL* z;

  //! weights
  REAL* w;

  //! Size of arrays (number of quadrature nodes)
  INT n;

  //! Number of quadrature nodes on 1 element
  INT nq_per_elm;

  //! Number of quadrature nodes in one direction
  INT nq1d;

} qcoordinates;

// Newer version
/**
 * \struct quadrature
 * \brief Returns coordinates of quadrature nodes
 */
typedef struct quadrature{

  //! Dimension of problem
  INT dim;

  //! quadrature points (ordered x(q1) y(q1) z(q1) x(q2) y(q2) z(q2) etc...)
  REAL* x;

  //! weights
  REAL* w;

  //! Size of arrays (total number of quadrature nodes)
  INT nq;

  //! Number of quadrature nodes on 1 entity (face or element or edge)
  INT nq_simplex;

  //! Number of quadrature nodes in one direction
  INT nq1d;

} quadrature;

/**
 * \struct fespace
 * \brief Returns properties of the finite-element space on entire mesh
 * \note TODO - will remove some unnecessary components
 */
typedef struct fespace{

  //! Type of finite element: 0-9 PX | 10-19 QX (not yet) | 20 Ned | 30 RT | -9 - -1 DGX (not yet) | 61 - face bubbles | 99 - single DOF for constraints
  // TODO: Rename to fe_type
  INT FEtype;

  //! Indicates if this is a space of scalara functions, 0, or a space of vector functions, 1.
  INT scal_or_vec;

  //! Where is the DoF defined: 0 - vertex; 1 - edge; 2 - face; 3 - tetrahedra; etc...
  INT dof_form;

  //! Number of Elements
  // TODO: Remove cause it's redundant
  INT nelm;

  //! Coordinates of DOF
  // TODO: Remove -> not needed
  coordinates* cdof;

  //! number of DOF
  INT ndof;

  //! number of DOF on boundary
  // TODO: Remove -> not needed
  INT nbdof;

  //! number of DOF per element
  INT dof_per_elm;

  //! Element to DOF map
  iCSRmat* el_dof;

  //! Dof per edge
  INT dof_per_edge;

  //! Edge to DOF map
  iCSRmat* ed_dof;

  //! number of DOF per face
  INT dof_per_face;

  //! Face to DOF map
  iCSRmat* f_dof;

  //! Dirichlet Boundaries (1 if Dirichlet; 0 if not)
  INT* dirichlet;

  //! DOF flags - indicates if the DOF is a special DOF (i.e. on certain boundary)
  INT* dof_flag;

  //! Perioidc Boundaries (For each DOF indicate if it is periodic with another DOF.  Mark -1 for non-periodic)
  INT* periodic;

  //! Basis Functions and Derivatives at a single quad point
  // Note: fe_local_data will point to these directly so it is allocated once
  REAL* phi;
  REAL* dphi;
  REAL* ddphi;

} fespace;

/**
 * \struct simplex_local_data
 * \brief Contains all relevant local data on a given element/face/edge
 *
 */
typedef struct simplex_local_data {

  //! Dimension of mesh
  INT dim;

  // Element Info
  //! Index of current entity (simplex)
  INT sindex;
  //! Volume (area/length) of ximplex
  REAL vol;
  //! Flag of element
  INT flag;

  // Vertex Info
  //! Number of vertices locally
  INT n_v;
  //! vertices on element/face/edge
  INT* local_v;
  //! coordinates of local vertices
  REAL* xv;

  // Edge Info
  //! Number of edges locally
  INT n_ed;
  //! edges on element/face/edge
  INT* local_ed;
  //! vertices on edges
  INT* v_on_ed;
  //! Lengths of edges
  REAL* ed_len;
  //! Tangent vector on edges
  REAL* ed_tau;
  //! Midpoint of edges
  REAL* ed_mid;

  // Face Info
  //! Number of faces locally
  INT n_f;
  //! Faces on element/face/edge
  INT* local_f;
  //! Edges on faces
  INT* ed_on_f;
  //! Vertices on Faces
  INT* v_on_f;
  //! Area of face
  REAL* f_area;
  //! Normal vector on faces
  REAL* f_norm;
  //! Barycenter midpoint of face
  REAL* f_mid;

  // Assembly data
  //! Quadrature on entity
  quadrature* quad_local;
  //! Reference Element maps
  REAL* ref_map; // The map x = B*xr + x0
  REAL* lams; // P1 basis at quad points (actually on ref element)
  REAL* gradlams; // Inverse map plus row sum which gives gradients of P1 basis

  // Space for extra stuff if needed
  REAL* dwork;
  INT* iwork;

} simplex_local_data;


/**
 * \struct local_fe_data
 * \brief Contains the local data of fe spaces and DoF on local entity
 */
 typedef struct fe_local_data {

   //! Number of spaces
   INT nspaces;

   //! Number of scalar unknowns
   INT nun;

   //! Type of FE spaces
   INT* fe_types;

   //! Indicates whether the spaces are scalar or vector-valued
   INT* scal_or_vec;

   //! Total number of DoF locally for all spaces
   INT n_dof;

   //! Array continaing number of DoF locally per space
   INT* n_dof_per_space;

   //! DoF on element/face/edge for all spaces
   INT* local_dof;

   //! Local Flags of DoF
   INT* local_dof_flags;

   //! Solution at local DoF for all spaces
   REAL* u_local;

   //! Basis functions for all spaces at a single quadrature point
   // Points to array allocated inside each fespace struct
   REAL** phi;
   REAL** dphi;
   REAL** ddphi;

 } fe_local_data;

/**
 * \struct block_fespace
 * \brief Block of fespaces for coupled problems -> holds all data
 *
 * \note TODO: Will change this to fe_system and remove nbdof
 */
typedef struct block_fespace {

  //! number of FEM spaces in system
  INT nspaces;

  //! number of unknowns in system (includes # of components for vectors)
  INT nun;

  //! total number of elements
  INT nelm;

  //! total number of dof
  INT ndof;

  //! total number of boundary dof
  // TODO: Remove -> unnecessary
  INT nbdof;

  //! blocks of fespaces
  fespace **var_spaces;

  //! Dirichlet Boundaries (1 if Dirichlet; 0 if not) for ALL unknowns
  INT* dirichlet;

  //! All DOF flags - indicates if the DOF is a special DOF (i.e. on certain boundary)
  INT* dof_flag;

  //! Local Mesh Data - stuff needed on a given element (or face or edge)
  simplex_local_data *simplex_data;

  //! Local FE Data - stuff needed on a given element for DoF Info
  fe_local_data *fe_data;

} block_fespace;


//**************** NEW STUFF **********************************//




#endif
