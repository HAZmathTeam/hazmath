//
//  solver.h
//  
//
//  Created by Hu, Xiaozhe on 5/13/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _solver_h
#define _solver_h

/**
 * \struct precond
 * \brief Preconditioner data and action
 *
 * \note This is the preconditioner structure for preconditioned iterative methods.
 */
typedef struct {
    
    //! data for preconditioner, void pointer
    void *data;
    
    //! action for preconditioner, void function pointer
    void (*fct)(REAL *, REAL *, void *);
    
} precond; /**< Data for general preconditioner passed to iterative solvers */



#endif
