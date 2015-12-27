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

/**
 * \struct Link
 * \brief Struct for Links
 */
typedef struct
{
    
    //! previous node in the linklist
    INT prev;
    
    //! next node in the linklist
    INT next;
    
} Link; /**< General data structure for Links */

/**
 * \struct linked_list
 * \brief A linked list node
 *
 * \note This definition is adapted from hypre 2.0.
 */
typedef struct linked_list
{
    
    //! data
    INT data;
    
    //! starting of the list
    INT head;
    
    //! ending of the list
    INT tail;
    
    //! next node
    struct linked_list *next_node;
    
    //! previous node
    struct linked_list *prev_node;
    
} ListElement; /**< Linked element in list */

/**
 * List of links
 */
typedef ListElement *LinkList; /**< linked list */

#endif
