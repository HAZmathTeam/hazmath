//
//  timestep.h
//  
//
//  Created by Adler, James on 10/15/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _timestep_h
#define _timesetp_h

/**
 * \struct timestepper
 * \brief Returns timestepping data
 */
// Assume the form: a*M du/dt + A u = f
// After timestepping we get: A_time*u = rhs_time
// For now we assume the following time-steppers:
// CN: (aM + 0.5*dt*A)u = 0.5*dt*(fprev+f) + (aM - 0.5*dt*A)uprev
// BDF1: (aM + dt*A)u = dt*f + aM*uprev
typedef struct timestepper{

  //! Number of Time Steps
  INT tsteps;

  //! Size of Time Step
  REAL dt;

  //! Timestepping Scheme: 0 = Crank-Nicolson; X = BDF-X
  INT time_scheme;

  //! Timestepping Scheme Name (string)
  char time_scheme_str[30];

  //! Indicates if RHS is time-dependent
  INT rhs_timedep;

  //! Indicates if Boundary Conditions are time-dependent
  INT bc_timedep;

  //! Current time step
  INT current_step;

  //! Current time
  REAL time;

  // Assume the PDE is always of the form
  // du/dt + Lu = f
  // -> d/dt Mx + Ax = f
  // -> At*x = ft
  //! M-matrix
  dCSRmat* M;

  //! A-matrix
  dCSRmat* A;

  //! Time-propagator matrix (A_time)
  dCSRmat* At;

  //! Time-propagator matrix with no BC
  dCSRmat* At_noBC;

  //! Solution at previous time step (eventually this should be array of solutions)
  dvector* sol_prev;

  //! Current solution
  dvector* sol;

  //! Current RHS (f)
  dvector* rhs;

  //! RHS at previous step (fprev)
  dvector* rhs_prev;

  //! RHS of Time Propagator
  dvector* rhs_time;
	
} timestepper;


#endif
