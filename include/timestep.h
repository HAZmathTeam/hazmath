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
 *
 * \note Assume the form: d(Mu)/dt + L(u) = f
 *        After timestepping we get: A_time*u = rhs_time
 * For now we assume the following time-steppers:
 *   CN: (M + 0.5*dt*A)u = 0.5*dt*(fprev+f) + M*uprev - 0.5*dt*L(uprev)
 *   BDF-1: (M + dt*A)u = dt*f + M*uprev
 *   BDF-k: (ak*M + bk*dt*A)u = bk*dt*f + \sum_{s=0:k-1} as*M*u_s
 * Here, A is the discretization of L.
 * Note that L can be nonlinear, so A represents the discrete Gateaux derivative of L, L'
 * In the linear case, L(uprev) = A*uprev (=> L = mat-vec and Ldata = A)
 */
typedef struct timestepper{

  //! Number of Time Steps
  INT tsteps;

  //! Size of Time Step
  REAL dt;

  //! Timestepping Scheme: 0 = Crank-Nicolson; X = BDF-X
  INT time_scheme;

  //! Timestepping Scheme Name (string)
  char time_scheme_str[30];

  //! Number of previous time steps to store
  INT old_steps;

  //! Indicates if RHS is time-dependent
  INT rhs_timedep;

  //! Current time step
  INT current_step;

  //! Current time
  REAL time;

  //! M-matrix (dcsr)
  dCSRmat* M;

  //! A-matrix (in nonlinear case this would be L'(u)[du]
  dCSRmat* A;

  //! L(u)-vector (in linear case this is the same as Au)
  void (*L)(void *,REAL *,REAL *);

  //! Data for L operator function (Linear Case: data=A)
  void *Ldata;

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

/**
 * \struct block_timestepper
 * \brief Returns timestepping data
 *
 * \note Assume the form: d(Mu)/dt + L(u) = f
 *        After timestepping we get: A_time*u = rhs_time
 * For now we assume the following time-steppers:
 *   CN: (M + 0.5*dt*A)u = 0.5*dt*(fprev+f) + M*uprev - 0.5*dt*L(uprev)
 *   BDF-1: (M + dt*A)u = dt*f + M*uprev
 *   BDF-k: (ak*M + bk*dt*A)u = bk*dt*f + \sum_{s=0:k-1} as*M*u_s
 * Here, A is the discretization of L.
 * Note that L can be nonlinear, so A represents the discrete Gateaux derivative of L, L'
 * In the linear case, L(uprev) = A*uprev (=> L = mat-vec and Ldata = A)
 * Assumes the matrices are all block_dcsr
 */
typedef struct block_timestepper{

  //! Number of Time Steps
  INT tsteps;

  //! Size of Time Step
  REAL dt;

  //! Timestepping Scheme: 0 = Crank-Nicolson; X = BDF-X
  INT time_scheme;

  //! Timestepping Scheme Name (string)
  char time_scheme_str[30];

  //! Number of previous time steps to store
  INT old_steps;

  //! Indicates if RHS is time-dependent
  INT rhs_timedep;

  //! Current time step
  INT current_step;

  //! Current time
  REAL time;

  //! M-matrix (block_dcsr)
  block_dCSRmat* M;

  //! A-matrix (in nonlinear case this would be L'(u)[du]
  block_dCSRmat* A;

  //! L(u)-vector (in linear case this is the same as Au)
  void (*L)(void *,REAL *,REAL *);

  //! Data for L operator function (Linear Case: data=A)
  void *Ldata;

  //! Time-propagator matrix (A_time)
  block_dCSRmat* At;

  //! Time-propagator matrix with no BC
  block_dCSRmat* At_noBC;

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

} block_timestepper;


#endif
