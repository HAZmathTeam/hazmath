//
//  nonlinear.h
//
//
//  Created by Adler, James on 10/18/16.
//
//  Contains Structs for nonlinear iterations
//  For now just assumes Newton, and FAS in the works
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _nonlinear_h
#define _nonlinear_h

/**
 * \struct newton
 * \brief Returns data for Newton Stepping
 */
typedef struct newton{

  //! Indicate if matrices are in block form 1=yes 0=no
  INT isblock;

  //! Max number of Newton Steps
  INT max_steps;

  //! Current Newton step
  INT current_step;

  //! Tolerance Type: 0 - ||nonlinear residual||<tol | 1 - ||update|| < tol
  INT tol_type;

  //! Stopping Tolerance
  REAL tol;

  //! Step Length: sol = sol_prev + step_length*update
  REAL step_length;

  // Assume the form A(sol) = f gets linearized to
  // Jacobian(sol_prev)[update] = f - A(sol_prev)
  //! Jacobian-matrix
  dCSRmat* Jac;

  //! Jacobian-matrix Block CSR
  block_dCSRmat* Jac_block;

  //! Solution at previous Newton step
  dvector* sol_prev;

  //! Current solution
  dvector* sol;

  //! Newton update
  dvector* update;

  //! RHS of Newton Iteration (nonlinear residual f- (A(sol_prev))
  dvector* rhs;

  //! Norm of nonlinear residual (combined total if in block form)
  REAL res_norm;

  //! Norm of update (combined total if in block form)
  REAL update_norm;

} newton;

#endif
