//
//  fem.h
//  
//
//  Created by Adler, James on 2/1/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _fem_h
#define _fem_h

#include "sparse.h"
#include "vec.h"
#include "grid.h"

/**
 * \struct qcoordinates
 * \brief Returns coordinates of quadrature nodes
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
	
} qcoordinates;

/**
 * \struct fespace
 * \brief Returns properties of the finite-element space
 */
typedef struct fespace{

  //! Type of finite element: 0-9 PX | 10-19 QX (not yet) | 20 Ned | 30 RT | -9 - -1 DGX (not yet)
  INT FEtype;

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

  //! Face to DOF map
  iCSRmat* f_dof;

  //! Dirichlet Boundary map
  INT* dof_bdry;
	
} fespace;

/**
 * \struct block_fespace
 * \brief Block of fespaces for coupled problems
 *
 */
typedef struct block_fespace {
    
  //! number of FEM spaces in system
  INT nspaces;

  //! total number of dof
  INT ndof;

  //! total number of boundary dof
  INT nbdof;
  
  //! blocks of dCSRmat, point to blocks[brow][bcol]
  fespace **var_spaces;

  //! Boundary DOF array for ALL unknowns
  INT* dof_bdry;
    
} block_fespace; /**< Matrix of REAL type in Block CSR format */


#endif
