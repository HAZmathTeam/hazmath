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
 */
typedef struct qcoordinates{

  //! quadrature points (ordered x(q1) y(q1) z(q1) x(q2) y(q2) z(q2) etc...)
  REAL* x;

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
 * \struct block_fespace
 * \brief Block of fespaces for coupled problems
 *
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

} block_fespace; /**< Matrix of REAL type in Block CSR format */


#endif
