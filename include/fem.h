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

/**
 * \struct fespace
 * \brief Returns properties of the finite-element space
 * \note TODO - will replace with fe_space setup defined below soon
 */
typedef struct fespace{

  //! Type of finite element: 0-9 PX | 10-19 QX (not yet) | 20 Ned | 30 RT | -9 - -1 DGX (not yet) | 61 - face bubbles | 99 - single DOF for constraints
  INT FEtype;

  //! Indicates if this is a space of scalara functions, 0, or a space of vector functions, 1.
  INT scal_or_vec;

  //! Number of Elements
  INT nelm;

  //! Coordinates of DOF
  coordinates* cdof;

  //! number of DOF
  INT ndof;

  //! number of DOF on boundary
  INT nbdof;

  //! number of DOF per element
  INT dof_per_elm;

  //! Element to DOF map
  iCSRmat* el_dof;

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

  //! Basis Functions and Derivatives
  REAL* phi;
  REAL* dphi;

} fespace;

/**
 * \struct local_data
 * \brief Contains all relevant data on a given element/face/edge
 *
 */
typedef struct fe_local_data {

  //! Number of vertices locally
  INT nlocal_vert;

  //! Number of DoF locally for all spaces
  INT nlocal_dof;

  //! Numer of DoF locally per space
  INT* nlocal_dof_space;

  //! Dimension of mesh
  INT dim;

  //! DoF on element/face/edge
  INT* local_dof;

  //! Local Flags of DoF
  INT* local_dof_flags;

  //! vertices on element/face/edge
  INT* local_vert;

  //! coordinates of local vertices
  REAL* xv;

  //! Number of vertices per dof for each space
  INT* nvert_per_dof_space;

  //! vertices on dof per space (nspaces * nlocal_dof_space[i] * nver_per_dof_space[i])
  INT* v_on_dof_space;

  //! Solution at local DoF
  REAL* u_local;

  //! Quadrature on entity
  qcoordinates* quad_local;

  //! Reference Element maps
  REAL* ref_map; // The map x = B*xr + x0
  REAL* lams; // P1 basis at quad points (actually on ref element)
  REAL* gradlams; // Inverse map plus row sum which gives gradients of P1 basis

  // Space for extra stuff if needed
  REAL* dwork;
  INT* iwork;

} fe_local_data;

/**
 * \struct block_fespace
 * \brief Block of fespaces for coupled problems
 *
 * \note TODO: Will change this to fe_system and remove nbdof
 */
typedef struct block_fespace {

  //! number of FEM spaces in system
  INT nspaces;

  //! number of unknowns in system (includes # of components for vectors)
  INT nun;

  //! total number of dof
  INT ndof;

  //! total number of boundary dof
  INT nbdof;

  //! blocks of fespaces
  fespace **var_spaces;

  //! Dirichlet Boundaries (1 if Dirichlet; 0 if not) for ALL unknowns
  INT* dirichlet;

  //! All DOF flags - indicates if the DOF is a special DOF (i.e. on certain boundary)
  INT* dof_flag;

  //! Local Data - stuff needed on a given element (or face or edge)
  fe_local_data *loc_data;

} block_fespace;


//**************** NEW STUFF **********************************//

/**
 * \struct qcoords
 * \brief Returns coordinates of quadrature nodes
 */
typedef struct qcoords{

  //! Dimension of problem
  INT dim;

  //! quadrature points (ordered x(q1) y(q1) z(q1) x(q2) y(q2) z(q2) etc...)
  REAL* x;

  //! weights
  REAL* w;

  //! Size of arrays (number of quadrature nodes)
  INT n;

  //! Number of quadrature nodes on 1 entity (face or element or edge)
  INT nq_per_region;

  //! Number of quadrature nodes in one direction
  INT nq1d;

} qcoords;

/**
 * \struct fe_space
 * \brief Returns properties of the finite-element space
 */
typedef struct fe_space{

  //! Type of finite element:
  // Implemented:  0-2: PX | 20: Ned-0 | 30: RT-0 | 62: P1 + facebubble | 99: single DoF for constraints
  // TBI: -9--1: DGX | 10-19: QX | 21-29: Ned-X | 31:39: RT-X | 40-59: Ned/RT-X on quads | 60: vector of scalars | 61: MINI element
  INT fe_type;

  //! Indicates if this is a space of scalar functions, 0, or a space of vector functions, 1.
  INT scal_or_vec;

  //! Where is the DoF defined: 0 - vertex; 1 - edge; 2 - face; 3 - tetrahedra; etc...
  INT dof_form;

  //! number of DOF
  INT ndof;

  //! number of DOF per element
  INT dof_per_elm;

  //! number of DOF per face
  INT dof_per_face;

  //! number of DOF per edge
  INT dof_per_edge;

  //! Element to DOF map
  iCSRmat* el_dof;

  //! Face to DOF map
  iCSRmat* f_dof;

  //! Edge to DOF map
  iCSRmat* ed_dof;

  //! Dirichlet Boundaries (1 if Dirichlet; 0 if not)
  INT* dirichlet;

  //! DOF flags - indicates if the DOF is a special DOF (i.e. on certain boundary)
  INT* dof_flag;

  //! Perioidc Boundaries (For each DOF indicate if it is periodic with another DOF.  Mark -1 for non-periodic)
  INT* periodic;

  //! Basis Functions and Derivatives on quadrature of reference element
  REAL* phi;
  REAL* dphi;
  REAL* ddphi;

} fe_space;


#endif
