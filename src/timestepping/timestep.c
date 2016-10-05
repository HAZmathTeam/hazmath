/*! \file timestep.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 2/18/16.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 */

/* This code will contain all the tools needed to perform timestepping */

#include "hazmat.h"

/******************************************************************************************************/
void fixrhs_time(dvector* b,dvector* b_old,dCSRmat* M,dCSRmat* A,dvector* uprev,INT time_scheme,REAL dt,dvector* b_update)
{
  /********* Updates the right-hand side for a timestepping scheme *********************
   *
   *     Assumes we have: M du/dt + Au = b
   * 
   *	Input:		
   *            b            Original right-hand side from current time-step
   *            b_old        Original RHS from previous time-step (only needed if RHS is time-dependent.  If not just use b twice)
   *            A            Spatial Matrix
   *            M            Mass Matrix
   *            uprev        Previous solution
   *            dof_bdry     Indicates which DOF are on boundary
   *            timescheme   What type of timestepping to use (0->CN 1->Backward Euler (BDF1) etc...)
   *            dt           Time step size
   *
   *	Output:		
   *            b_update     Updated rhs
   *
   */
	
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
      
    dcsr_add_1(M,2.0/dt,A,-1.0,&Atemp);
      
    // Compute updated RHS
    dcsr_aAxpy_1(1.0,&Atemp,uprev->val,b_update->val);
      
    // Free Atemp
    dcsr_free(&Atemp);

  // Backward Euler (alpha = 1/dt): (alpha M + A)u = alpha M uprev + b
  } else if(time_scheme==1) { 

    dcsr_aAxpy_1(1.0,M,uprev->val,b_update->val);
    
  } else { // Unknown
    printf("Not sure how to do timestepping, so I have decided to give up.\n\n");
    exit(2);
  }
	
  return;
}
/******************************************************************************************************/

/******************************************************************************************************/
void get_timeoperator(dCSRmat* M,dCSRmat* A,INT time_scheme,REAL dt,dCSRmat* Atime)
{
  /********* Gets the matrix to solve for timestepping scheme *********************
   *
   *     Assumes we have: M du/dt + Au = b
   * 
   *	Input:		
   *            A            Spatial Matrix
   *            M            Mass Matrix
   *            timescheme   What type of timestepping to use (0->CN 1->Backward Euler (BDF1) etc...)
   *            dt           Time step size
   *
   *	Output:		
   *            Atime        Matrix to solve with
   *
   */
	
  // Crank-Nicolson (alpha = 2/dt): (alpha M + A)u = (alpha M - A)uprev + (b_old + b) 
  if(time_scheme==0) { 
    
    dcsr_add_1(M,2.0/dt,A,1.0,Atime);

  // Backward Euler (alpha = 1/dt): (alpha M + A)u = alpha M uprev + b
  } else if(time_scheme==1) { 

    dcsr_add_1(M,1.0/dt,A,1.0,Atime);

  } else { // Unknown
    printf("Not sure how to do timestepping, so I have decided to give up.\n\n");
    exit(2);
  }
	
  return;
}
/******************************************************************************************************/

