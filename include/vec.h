//
//  vec.h
//  
//
//  Created by Hu, Xiaozhe on 1/9/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _vec_h
#define _vec_h

/**
 * \struct dvector
 * \brief Vector with n entries of REAL type
 */
typedef struct dvector{
    
    //! number of rows
    INT row;
    
    //! actual vector entries
    REAL *val;
    
} dvector; /**< Vector of REAL type */

/**
 * \struct ivector
 * \brief Vector with n entries of INT type
 */
typedef struct ivector{
    
    //! number of rows
    INT row;
    
    //! actual vector entries
    INT *val;
    
} ivector; /**< Vector of INT type */



#endif
