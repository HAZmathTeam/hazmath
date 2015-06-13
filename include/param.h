//
//  param.h
//  
//
//  Created by Hu, Xiaozhe on 06/13/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _param_h
#define _param_h

/**
 * \struct input_param
 * \brief Input parameters
 *
 * Input parameters, reading from disk file
 */
typedef struct {
    
    //----------------
    // output flags
    //----------------
    SHORT print_level;   /**< print level */
    SHORT output_type;   /**< type of output stream */
    
    //----------------
    // files
    //----------------
    char workdir[256];   /**< working directory for data files */
    char inifile[256];   /**< ini file name */
    char gridfile[256];  /**< grid file name */
    
    //--------------------------
    // finite element parameters
    //--------------------------
    // genearal parameters
    INT dim;            /**< dimension */
    INT nquad;          /**< quadature nodes in each direction */
    
    // parameters for H(D) equations
    INT FE_type;        /**< finite element type */
    
    // paramters for Stokes/NS equations
    INT FE_type_velocity;   /**< finite element type for velocity */
    INT FE_type_pressure;    /**< finite element type for pressure */
    
    //----------------------------
    // time steppng paramters
    //----------------------------
    INT time_step_type;  /**< time step type */
    INT time_steps;      /**< time steps */
    REAL time_step_size; /**< time step type */
    
    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    INT nonlinear_itsolver_maxit;   /**< maximal iterations of nonlinear solver*/
    REAL nonlinear_itsolver_tol;    /**< tolerance for nonlinear solver */
    
    //-------------------------
    // linear solver parameters
    //-------------------------
    // Iterative solver
    INT linear_itsolver_type;         /**< type of iteative linear solver */
    INT linear_itsolver_maxit;        /**< maximal iterations of linear iterative solver*/
    REAL linear_itsolver_tol;         /**< tolerance for linear iterative solver */
    
} input_param; /**< Input parameters */



#endif
