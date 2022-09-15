/*! \file src/timestepping/timestep.c
 *
 * \brief This code will contain all the tools needed to perform timestepping
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 2/18/16.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 * \note modified by James Adler 02/22/2019 for 0-1 fix
 */

#include "hazmath.h"

//**** NON-BLOCK STUFF ******/
/******************************************************************************************************/
/*!
 * \fn void initialize_timestepper(timestepper *tstepper,input_param *inparam,INT rhs_timedep,INT ndof)
 *
 * \brief Initialize the timestepping struct.
 *
 * \param inparam       Input from input parameter list
 * \param ndof          Number of DOF in the system
 * \param rhs_timedep   Indicates if the RHS (f) is time-dependent (1 or 0)
 *
 * \return tstepper     Struct for Timestepping
 *
 * \note Assume the form: d(Mu)/dt + L(u) = f
 *       After timestepping we get: A_time*u = rhs_time
 *       For now we assume the following time-steppers:
 *        CN: (M + 0.5*dt*A)u = 0.5*dt*(fprev+f) + M*uprev - 0.5*dt*L(uprev)
 *        BDF-1: (M + dt*A)u = dt*f + M*uprev
 *        BDF-k: (ak*M + bk*dt*A)u = bk*dt*f + \sum_{s=0:k-1} as*M*u_s
 *       Here, A is the discretization of L.
 *       Note that L can be nonlinear, so A represents the discrete Gateaux derivative of L, L'
 *       In the linear case, L(uprev) = A*uprev (=> L = mat-vec and Ldata = A)
 *
 */
void initialize_timestepper(timestepper *tstepper,input_param *inparam,INT rhs_timedep,INT ndof)
{

  // Number of time steps
  tstepper->tsteps = inparam->time_steps;

  // Time step size
  tstepper->dt = inparam->time_step_size;

  // Time step Scheme
  tstepper->time_scheme = inparam->time_step_type;
  if(tstepper->time_scheme==0) {
    sprintf(tstepper->time_scheme_str,"Crank-Nicolson");
    tstepper->old_steps = 1;
  } else {
    sprintf(tstepper->time_scheme_str,"BDF-%lld",(long long )tstepper->time_scheme);
    tstepper->old_steps = tstepper->time_scheme;
  }

  // Indicator if rhs or boundaries are time-dependent
  tstepper->rhs_timedep = rhs_timedep;

  // Current Time and Time Step
  tstepper->current_step = 0;
  if( inparam->time_start){
    tstepper->time = inparam->time_start;
  } else {
    tstepper->time = 0.0;
  }

  // Matrices and Vectors
  tstepper->M=NULL;
  tstepper->A=NULL;
  tstepper->At=malloc(sizeof(struct dCSRmat)); /* A_time */
  tstepper->At_noBC=malloc(sizeof(struct dCSRmat)); /* A_time with no boundary elimination */
  tstepper->sol_prev=malloc(sizeof(struct dvector)); /* uprev */
  tstepper->sol=NULL;      /* u */
  tstepper->rhs=NULL;     /* f */
  tstepper->rhs_prev=malloc(sizeof(struct dvector)); /* fprev */
  tstepper->rhs_time=malloc(sizeof(struct dvector));

  dvec_alloc(tstepper->old_steps*ndof,tstepper->sol_prev);
  dvec_alloc(ndof,tstepper->rhs_prev);
  dvec_alloc(ndof,tstepper->rhs_time);

  // Set L operator for RHS (L=A if linear, otherwise call an assembly)
  if(inparam->nonlinear_itsolver_type==0) { // Linear
    tstepper->L=dcsr_mxv_forts;
  } else {
    tstepper->L=NULL;
  }
  tstepper->Ldata=NULL;

  return;
}
/******************************************************************************************************/

/****************************************************************************************/
/*!
 * \fn void free_timestepper(timestepper* ts)
 *
 * \brief Frees memory of arrays of timestepping struct
 *
 * \return n_it         Freed struct for Timestepping
 *
 */
void free_timestepper(timestepper* ts)
{
  if(ts->A) {
    dcsr_free(ts->A);
    ts->A=NULL;
  }

  if(ts->M) {
    dcsr_free(ts->M);
    ts->M=NULL;
  }

  if(ts->At) {
    dcsr_free(ts->At);
    free(ts->At);
    ts->At=NULL;
  }

  if(ts->At_noBC) {
    dcsr_free(ts->At_noBC);
    free(ts->At_noBC);
    ts->At_noBC=NULL;
  }

  if(ts->sol) {
    dvec_free(ts->sol);
    ts->sol=NULL;
  }

  if(ts->sol_prev) {
    dvec_free(ts->sol_prev);
      free(ts->sol_prev);
      ts->sol_prev=NULL;
  }

  if(ts->rhs) {
      dvec_free(ts->rhs);
      ts->rhs=NULL;
  }

  if(ts->rhs_prev) {
      dvec_free(ts->rhs_prev);
      free(ts->rhs_prev);
      ts->rhs_prev=NULL;
  }

  if(ts->rhs_time) {
      dvec_free(ts->rhs_time);
      free(ts->rhs_time);
      ts->rhs_time=NULL;
  }

  return;
}
/****************************************************************************************/

/****************************************************************************************/
/*!
 * \fn void update_timestep(timestepper *tstepper)
 *
 * \brief Updates the Timestepping data at each step.
 *
 * \return tstepper     Updated timestepping struct
 *
 */
void update_timestep(timestepper *tstepper)
{
  INT i,j;
  INT ndof = tstepper->sol->row;

  // Counters and Physical Time
  tstepper->current_step++;
  tstepper->time = tstepper->current_step*tstepper->dt;

  // If using BDF-k store last k timesteps
  INT k = tstepper->old_steps;

  // Solution
  for(i=1;i<k;i++) {
    for(j=0;j<ndof;j++) {
      tstepper->sol_prev->val[(k-i)*ndof+j] = tstepper->sol_prev->val[(k-i-1)*ndof+j];
    }
  }
  for(j=0;j<ndof;j++) {
    tstepper->sol_prev->val[j] = tstepper->sol->val[j];
  }

  // RHS
  for(j=0;j<ndof;j++) {
    tstepper->rhs_prev->val[j] = tstepper->rhs->val[j];
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_timeoperator(timestepper* ts,INT first_visit,INT cpyNoBC)
 *
 * \brief Gets the matrix to solve for timestepping scheme
 *        Assumes we have: M du/dt + L(u) = b
 *
 * \param ts            Timestepping struct
 * \param first_visit   Indicates if this is the first visit in order to do allocation.
 * \param cpyNoBC       Indicates if you'd like to store the time matrix without BC eliminated
 *
 * \return ts.Atime     Matrix to solve with
 *
 */
void get_timeoperator(timestepper* ts,INT first_visit,INT cpyNoBC)
{
  // Flag for errors
  SHORT status;

  REAL dt = ts->dt;
  INT time_scheme = ts->time_scheme;

  switch (time_scheme) {
  case 0: // Crank-Nicolson: (M + 0.5*dt*A)u = (M*uprev - 0.5*dt*L(uprev) + 0.5*dt*(b_old + b)
    dcsr_add(ts->M,1.0,ts->A,0.5*dt,ts->At);
    break;
  case 1: // Backward Euler: (M + dt*A)u = M*uprev + dt*b
    dcsr_add(ts->M,1.0,ts->A,dt,ts->At);
    break;
  case 2: // BDF-2: (M + (2/3)*dt*A)u = (4/3)*M*uprev - (1/3)*M*uprevprev+ (2/3)*dt*b
    dcsr_add(ts->M,1.0,ts->A,2.0*dt/3.0,ts->At);
    break;
  default:
    status = ERROR_TS_TYPE;
    check_error(status, __FUNCTION__);
  }

  if(first_visit)
    dcsr_alloc(ts->At->row,ts->At->col,ts->At->nnz,ts->At_noBC);
  if(cpyNoBC)
    dcsr_cp(ts->At,ts->At_noBC);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void update_time_rhs(timestepper *ts)
 *
 * \brief Updates the right-hand side for a timestepping scheme
 *
 * \param ts            Timestepping struct
 *
 * \return ts.rhstime   RHS to solve with
 *
 */
void update_time_rhs(timestepper *ts)
{
  // Flag for errors
  SHORT status;

  if(ts->time_scheme==0) { // Crank-Nicolson: (M + 0.5*dt*A)u = (M - 0.5*dt*L(uprev) + 0.5*dt*(b_old + b)

    dvector btmp = dvec_create(ts->rhs->row);
    dvector Mu = dvec_create(ts->rhs->row);
    dvector Lu = dvec_create(ts->rhs->row);

    // Add new and old RHS
    dvec_axpyz(1.0,ts->rhs_prev,ts->rhs,&btmp);
    // btmp = 0.5*dt*(b_old+b)
    dvec_ax(0.5*ts->dt,&btmp);

    // Obtain M*uprev
    dcsr_mxv(ts->M,ts->sol_prev->val,Mu.val);

    // M*uprev + 0.5*dt(b_old+b)
    dvec_axpy(1.0,&Mu,&btmp);

    // Obtain L(uprev) (could be nonlinear)
    ts->L(ts->Ldata,ts->sol_prev->val,Lu.val);

    // Compute updated RHS
    dvec_axpyz(-0.5*ts->dt,&Lu,&btmp,ts->rhs_time);

    dvec_free(&btmp);
    dvec_free(&Mu);
    dvec_free(&Lu);

   } else if(ts->time_scheme==1) { // Backward Euler: (M + dt*A)u = M*uprev + dt*b

    dvector btmp = dvec_create(ts->sol->row);

    // Get M*uprev
    dcsr_mxv(ts->M,ts->sol_prev->val,btmp.val);

    // Compute updated RHS
    dvec_axpyz(ts->dt,ts->rhs,&btmp,ts->rhs_time);

    dvec_free(&btmp);
  } else if(ts->time_scheme==2) { // BDF-2: (M + (2/3)*dt*A)u = (4/3)*M*uprev - (1/3)*M*uprevprev + (2/3)*dt*b

    dvector btmp1 = dvec_create(ts->sol->row);
    dvector btmp2 = dvec_create(ts->sol->row);
    REAL* solprevptr;

    // Get (4/3)*M*uprev
    dcsr_mxv(ts->M,ts->sol_prev->val,btmp1.val);
    dvec_ax(4.0/3.0,&btmp1);

    // Get -(1/3)*M*uprevprev
    solprevptr = ts->sol_prev->val + ts->sol->row;
    dcsr_mxv(ts->M,solprevptr,btmp2.val);
    dvec_ax(-1.0/3.0,&btmp2);

    // Add first two components
    dvec_axpy(1.0,&btmp1,&btmp2);

    // Compute updated RHS
    dvec_axpyz((2.0/3.0)*ts->dt,ts->rhs,&btmp2,ts->rhs_time);

    dvec_free(&btmp1);
    dvec_free(&btmp2);
  } else {
    status = ERROR_TS_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
/******************************************************************************************************/

//**** BLOCK Versions ******/
/******************************************************************************************************/
/*!
 * \fn void initialize_blktimestepper(block_timestepper *tstepper,input_param *inparam,INT rhs_timedep,INT ndof,INT blksize)
 *
 * \brief Initialize the BLOCK timestepping struct.
 *
 * \param inparam       Input from input parameter list
 * \param rhs_timedep   Indicates if the RHS (f) is time-dependent (1 or 0)
 * \param ndof          Number of DOF in the system
 * \param blksize       Number of block rows in matrix (assume rows=cols)
 *
 * \return tstepper     Struct for Block Timestepping
 *
 * \note Assume the form: d(Mu)/dt + L(u) = f
 *       After timestepping we get: A_time*u = rhs_time
 *       For now we assume the following time-steppers:
 *        CN: (M + 0.5*dt*A)u = 0.5*dt*(fprev+f) + M*uprev - 0.5*dt*L(uprev)
 *        BDF-1: (M + dt*A)u = dt*f + M*uprev
 *        BDF-k: (ak*M + bk*dt*A)u = bk*dt*f + \sum_{s=0:k-1} as*M*u_s
 *       Here, A is the discretization of L.
 *       Note that L can be nonlinear, so A represents the discrete Gateaux derivative of L, L'
 *       In the linear case, L(uprev) = A*uprev (=> L = mat-vec and Ldata = A)
 */
void initialize_blktimestepper(block_timestepper *tstepper,input_param *inparam,INT rhs_timedep,INT ndof,INT blksize)
{
  // Number of time steps
  tstepper->tsteps = inparam->time_steps;

  // Time step size
  tstepper->dt = inparam->time_step_size;

  // Time step Scheme
  tstepper->time_scheme = inparam->time_step_type;
  if(tstepper->time_scheme==0) {
    sprintf(tstepper->time_scheme_str,"Crank-Nicolson");
    tstepper->old_steps = 1;
  } else {
    sprintf(tstepper->time_scheme_str,"BDF-%lld",(long long )tstepper->time_scheme);
    tstepper->old_steps = tstepper->time_scheme;
  }

  // Indicator if rhs or boundaries are time-dependent
  tstepper->rhs_timedep = rhs_timedep;

  // Current Time and Time Step
  tstepper->current_step = 0;
  if( inparam->time_start){
    tstepper->time = inparam->time_start;
  } else {
    tstepper->time = 0.0;
  }

  // Matrices and Vectors
  tstepper->M=NULL;
  tstepper->A=NULL;
  tstepper->At=malloc(sizeof(struct block_dCSRmat)); /* A_time */
  tstepper->At_noBC=malloc(sizeof(struct block_dCSRmat)); /* A_time with no boundary elimination */
  tstepper->sol_prev=malloc(sizeof(struct dvector)); /* uprev */
  tstepper->sol=NULL;      /* u */
  tstepper->rhs=NULL;     /* f */
  tstepper->rhs_prev=malloc(sizeof(struct dvector)); /* fprev */
  tstepper->rhs_time=malloc(sizeof(struct dvector));

  bdcsr_alloc(blksize,blksize,tstepper->At);
  bdcsr_alloc(blksize,blksize,tstepper->At_noBC);
  dvec_alloc(tstepper->old_steps*ndof,tstepper->sol_prev);
  dvec_alloc(ndof,tstepper->rhs_prev);
  dvec_alloc(ndof,tstepper->rhs_time);

  // Set L operator for RHS (L=A if linear, otherwise call an assembly)
  if(inparam->nonlinear_itsolver_type==0) { // Linear
    tstepper->L=bdcsr_mxv_forts;
  } else {
    tstepper->L=NULL;
  }
  tstepper->Ldata=NULL;

  return;
}
/******************************************************************************************************/

/****************************************************************************************/
/*!
 * \fn free_blktimestepper(block_timestepper* ts, INT flag)
 *
 * \brief Frees memory of arrays of BLOCK timestepping struct
 *
 * \param ts      point to the block_timestepper
 * \param flag    flag of how the block dcsr matrices are allocated (0: minimal | 1: standard)
 *
 * \return n_it   Freed struct for block Timestepping
 *
 */
void free_blktimestepper(block_timestepper* ts, INT flag)
{

  if(ts->A) {
    if (flag == 0){
        bdcsr_free_minimal(ts->A);
    }
    else {
        bdcsr_free(ts->A);
    }
    ts->A=NULL;
  }

  if(ts->M) {
      if (flag == 0){
          bdcsr_free_minimal(ts->M);
      }
      else {
        bdcsr_free(ts->M);
      }
    ts->M=NULL;
  }

  if(ts->At) {
    bdcsr_free(ts->At);
    free(ts->At);
    ts->At=NULL;
  }

  if(ts->At_noBC) {
    bdcsr_free(ts->At_noBC);
    free(ts->At_noBC);
    ts->At_noBC=NULL;
  }

  if(ts->sol) {
    dvec_free(ts->sol);
    ts->sol=NULL;
  }

  if(ts->sol_prev) {
    dvec_free(ts->sol_prev);
      free(ts->sol_prev);
      ts->sol_prev=NULL;
  }

  if(ts->rhs) {
      dvec_free(ts->rhs);
      ts->rhs=NULL;
  }

  if(ts->rhs_prev) {
      dvec_free(ts->rhs_prev);
      free(ts->rhs_prev);
      ts->rhs_prev=NULL;
  }

  if(ts->rhs_time) {
      dvec_free(ts->rhs_time);
      free(ts->rhs_time);
      ts->rhs_time=NULL;
  }

  return;
}
/****************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void update_blktimestep(block_timestepper *tstepper)
 *
 * \brief Updates the BLOCK Timestepping data at each step.
 *
 * \return tstepper Updated block timestepping struct
 *
 */
void update_blktimestep(block_timestepper *tstepper)
{
  INT i,j;
  INT ndof = tstepper->sol->row;

  // Counters and Physical Time
  tstepper->current_step++;
  //tstepper->time = tstepper->current_step*tstepper->dt;
  tstepper->time = tstepper->time + tstepper->dt;

  // If using BDF-k store last k timesteps
  INT k = tstepper->old_steps;

  // Solution
  for(i=1;i<k;i++) {
    for(j=0;j<ndof;j++) {
      tstepper->sol_prev->val[(k-i)*ndof+j] = tstepper->sol_prev->val[(k-i-1)*ndof+j];
    }
  }
  for(j=0;j<ndof;j++) {
    tstepper->sol_prev->val[j] = tstepper->sol->val[j];
  }

  // RHS
  for(j=0;j<ndof;j++) {
    tstepper->rhs_prev->val[j] = tstepper->rhs->val[j];
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
* \fn void get_blktimeoperator(block_timestepper* ts,INT cpyNoBC)
*
* \brief Gets the matrix to solve for BLOCK timestepping scheme
*        Assumes we have: M du/dt + L(u) = b
*
* \param ts            block Timestepping struct
* \param first_visit   Indicates if this is the first visit in order to do allocation.
* \param cpyNoBC       Indicates if you'd like to store the time matrix without BC eliminated
*
* \return ts.Atime     Matrix to solve with
*
*/
void get_blktimeoperator(block_timestepper* ts,INT first_visit,INT cpyNoBC)
{
    // Flag for errors
    SHORT status;

    REAL dt = ts->dt;
    INT time_scheme = ts->time_scheme;

    switch (time_scheme) {
      case 0: // Crank-Nicolson: (M + 0.5*dt*A)u = (M*uprev - 0.5*dt*L(uprev) + 0.5*dt*(b_old + b)
        bdcsr_add(ts->M,1.0,ts->A,0.5*dt,ts->At);
        break;
      case 1: // Backward Euler: (M + dt*A)u = M*uprev + dt*b
        bdcsr_add(ts->M,1.0,ts->A,dt,ts->At);
        break;
      case 2: // BDF-2: (M + (2/3)*dt*A)u = (4/3)*M*uprev - (1/3)*M*uprevprev+ (2/3)*dt*b
        bdcsr_add(ts->M,1.0,ts->A,2.0*dt/3.0,ts->At);
        break;
      default:
        status = ERROR_TS_TYPE;
        check_error(status, __FUNCTION__);
    }
    if(cpyNoBC)
      bdcsr_cp(ts->At,ts->At_noBC);

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void update_blktime_rhs(block_timestepper *ts)
 *
 * \brief Updates the right-hand side for a BLOCK timestepping scheme
 *
 * \param ts            block Timestepping struct
 *
 * \return ts.rhstime   RHS to solve with
 *
 */
void update_blktime_rhs(block_timestepper *ts)
{
  // Flag for errors
  SHORT status;
  if(ts->time_scheme==0) { // Crank-Nicolson: (M + 0.5*dt*A)u = (M - 0.5*dt*L(uprev) + 0.5*dt*(b_old + b)

    dvector btmp = dvec_create(ts->rhs->row);
    dvector Mu = dvec_create(ts->rhs->row);
    dvector Lu = dvec_create(ts->rhs->row);

    // Add new and old RHS
    dvec_axpyz(1.0,ts->rhs_prev,ts->rhs,&btmp);
    // btmp = 0.5*dt*(b_old+b)
    dvec_ax(0.5*ts->dt,&btmp);

    // Obtain M*uprev
    bdcsr_mxv(ts->M,ts->sol_prev->val,Mu.val);

    // M*uprev + 0.5*dt(b_old+b)
    dvec_axpy(1.0,&Mu,&btmp);

    // Obtain L(uprev) (could be nonlinear)
    ts->L(ts->Ldata,ts->sol_prev->val,Lu.val);

    // Compute updated RHS
    dvec_axpyz(-0.5*ts->dt,&Lu,&btmp,ts->rhs_time);

    dvec_free(&btmp);
    dvec_free(&Mu);
    dvec_free(&Lu);

   } else if(ts->time_scheme==1) { // Backward Euler: (M + dt*A)u = M*uprev + dt*b

    dvector btmp = dvec_create(ts->sol->row);

    // Get M*uprev
    bdcsr_mxv(ts->M,ts->sol_prev->val,btmp.val);

    // Compute updated RHS
    dvec_axpyz(ts->dt,ts->rhs,&btmp,ts->rhs_time);

    dvec_free(&btmp);
  } else if(ts->time_scheme==2) { // BDF-2: (M + (2/3)*dt*A)u = (4/3)*M*uprev - (1/3)*M*uprevprev + (2/3)*dt*b

    dvector btmp1 = dvec_create(ts->sol->row);
    dvector btmp2 = dvec_create(ts->sol->row);
    REAL* solprevptr;

    // Get (4/3)*M*uprev
    bdcsr_mxv(ts->M,ts->sol_prev->val,btmp1.val);
    dvec_ax(4.0/3.0,&btmp1);

    // Get -(1/3)*M*uprevprev
    solprevptr = ts->sol_prev->val + ts->sol->row;
    bdcsr_mxv(ts->M,solprevptr,btmp2.val);
    dvec_ax(-1.0/3.0,&btmp2);

    // Add first two components
    dvec_axpy(1.0,&btmp1,&btmp2);

    // Compute updated RHS
    dvec_axpyz((2.0/3.0)*ts->dt,ts->rhs,&btmp2,ts->rhs_time);

    dvec_free(&btmp1);
    dvec_free(&btmp2);
  } else {
    status = ERROR_TS_TYPE;
    check_error(status, __FUNCTION__);
  }

  return;
}
/******************************************************************************************************/

// OLD STUFF NEEDED FOR MAXWELL RUNS

/******************************************************************************************************/
/*!
 * \fn void fixrhs_time(dvector* b,dvector* b_old,dCSRmat* M,dCSRmat* A,dvector* uprev,INT time_scheme,REAL dt,dvector* b_update)
 *
 * \brief Updates the right-hand side for a timestepping scheme
 *        Assume the form: a*M du/dt + A u = f
 *        After timestepping we get: A_time*u = rhs_time
 *        For now we assume the following time-steppers:
 *         CN:   (aM + 0.5*dt*A)u = 0.5*dt*(fprev+f) + (aM - 0.5*dt*A)uprev
 *         BDF1: (aM + dt*A)u = dt*f + aM*uprev
 *
 * \param b            Original right-hand side from current time-step
 * \param b_old        Original RHS from previous time-step (only needed if RHS is time-dependent.  If not just use b twice)
 * \param A            Spatial Matrix
 * \param M            Mass Matrix
 * \param uprev        Previous solution
 * \param timescheme   What type of timestepping to use (0->CN 1->Backward Euler (BDF1) etc...)
 * \param dt           Time step size
 *
 * \return b_update     Updated rhs
 *
 */
void fixrhs_time(dvector* b,dvector* b_old,dCSRmat* M,dCSRmat* A,dvector* uprev,INT time_scheme,REAL dt,dvector* b_update)
{
  INT i;

    for(i=0;i<b->row;i++) {

      b_update->val[i] = 0.0;

    }

  // Crank-Nicolson (alpha = 2/dt): (alpha M + A)u = (alpha M - A)uprev + (b_old + b)
  if(time_scheme==0) {

      dCSRmat Atemp;

    // Add new and old RHS
    dvec_axpyz(1.0,b_old,b,b_update);

    // Obtain alpha M - A

    dcsr_add(M,2.0/dt,A,-1.0,&Atemp);

    // Compute updated RHS
    dcsr_aAxpy(1.0,&Atemp,uprev->val,b_update->val);

    // Free Atemp
    dcsr_free(&Atemp);

  // Backward Euler (alpha = 1/dt): (alpha M + A)u = alpha M uprev + b
  } else if(time_scheme==1) {

    dcsr_aAxpy(1.0,M,uprev->val,b_update->val);

  } else { // Unknown
    printf("Not sure how to do timestepping, so I have decided to give up.\n\n");
    exit(2);
  }

  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
/*!
 * \fn void get_timeoperator_old(dCSRmat* M,dCSRmat* A,INT time_scheme,REAL dt,dCSRmat* Atime)
 *
 * \brief Gets the matrix to solve for timestepping scheme
 *        Assumes we have: M du/dt + Au = b
 *        After timestepping we get: A_time*u = rhs_time
 *        For now we assume the following time-steppers:
 *         CN:   (aM + 0.5*dt*A)u = 0.5*dt*(fprev+f) + (aM - 0.5*dt*A)uprev
 *         BDF1: (aM + dt*A)u = dt*f + aM*uprev
 *
 * \param A            Spatial Matrix
 * \param M            Mass Matrix
 * \param timescheme   What type of timestepping to use (0->CN 1->Backward Euler (BDF1) etc...)
 * \param dt           Time step size
 *
 * \return Atime       Updated propagation matrix
 *
 */
void get_timeoperator_old(dCSRmat* M,dCSRmat* A,INT time_scheme,REAL dt,dCSRmat* Atime)
{
  // Crank-Nicolson (alpha = 2/dt): (alpha M + A)u = (alpha M - A)uprev + (b_old + b)
  if(time_scheme==0) {

    dcsr_add(M,2.0/dt,A,1.0,Atime);

  // Backward Euler (alpha = 1/dt): (alpha M + A)u = alpha M uprev + b
  } else if(time_scheme==1) {

    dcsr_add(M,1.0/dt,A,1.0,Atime);

  } else { // Unknown
    printf("Not sure how to do timestepping, so I have decided to give up.\n\n");
    exit(2);
  }

  return;
}
/******************************************************************************************************/
