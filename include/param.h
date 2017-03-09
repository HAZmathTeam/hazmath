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
    
    //----------------
    // files
    //----------------
    char inifile[256];   /**< input parameter file name */
    char gridfile[256];  /**< grid file name */
    char output_dir[256];   /**< output directory */
    
    //--------------------------
    // finite element parameters
    //--------------------------
    // genearal parameters
    INT nquad;          /**< quadature nodes in each direction */
    
    // parameters for H(D) equations (examples only)
    INT FE_type;        /**< finite element type */

    //----------------------------
    // time steppng paramters
    //----------------------------
    INT  time_step_type;  /**< time step type */
    INT  time_steps;      /**< time steps */
    REAL time_step_size; /**< time step type */
    INT  rhs_time_dep;    /**< indicates if rhs is time-dependent */
    
    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    INT nonlinear_itsolver_type;      /**< type of nonlinear solver */
    INT  nonlinear_itsolver_maxit;   /**< maximal iterations of nonlinear solver*/
    REAL nonlinear_itsolver_tol;    /**< tolerance for nonlinear solver */
    INT  nonlinear_itsolver_toltype; /**< type of stopping tolerance for nonlinear solver */

    //-------------------------
    // linear solver parameters
    //-------------------------
    // Iterative solver
    INT   linear_itsolver_type;         /**< type of iteative linear solver */
    INT   linear_itsolver_maxit;        /**< maximal iterations of linear iterative solver*/
    REAL  linear_itsolver_tol;         /**< tolerance for linear iterative solver */
    SHORT linear_stop_type;           /**< stop type of linear iterative solver */
    INT   linear_restart;                      /**< restart number used in GMRES */
    
    // Preconditioner
    INT linear_precond_type;                 /**< type of preconditioner for iterative solvers */
    
    // Algebraic Multigrid
    SHORT AMG_type;                /**< Type of AMG */
    SHORT AMG_levels;              /**< maximal number of levels */
    SHORT AMG_cycle_type;          /**< type of cycle */
    SHORT AMG_smoother;            /**< type of smoother */
    SHORT AMG_smooth_order;        /**< order for smoothers */
    REAL  AMG_relaxation;           /**< over-relaxation parameter for SOR */
    SHORT AMG_presmooth_iter;      /**< number of presmoothing */
    SHORT AMG_postsmooth_iter;     /**< number of postsmoothing */
    SHORT AMG_polynomial_degree;   /**< degree of polynomial smoothers  */
    INT   AMG_coarse_dof;            /**< max number of coarsest level DOF */
    REAL  AMG_tol;                  /**< tolerance for AMG if used as preconditioner */
    INT   AMG_maxit;                 /**< number of iterations for AMG used as preconditioner */
    SHORT AMG_coarse_solver;       /**< coarse solver type */
    SHORT AMG_coarse_scaling;      /**< switch of scaling of the coarse grid correction */
    SHORT AMG_amli_degree;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type; /**< type of Krylov method used by nonlinear AMLI cycle */
    INT AMG_Schwarz_levels;        /**< number of levels use Schwarz smoother */

    // Classsical AMG
    SHORT AMG_coarsening_type;     /**< coarsening type */
    SHORT AMG_interpolation_type;  /**< interpolation type */
    REAL AMG_strong_threshold;     /**< strong threshold for coarsening */
    REAL AMG_truncation_threshold; /**< truncation factor for interpolation */
    REAL AMG_max_row_sum;          /**< maximal row sum */
    INT AMG_aggressive_level;      /**< number of levels use aggressive coarsening */
    INT AMG_aggressive_path;       /**< number of paths used to determine strongly coupled C-set */

    // Aggregation AMG
    SHORT AMG_aggregation_type;    /**< aggregation type */
    REAL AMG_strong_coupled;       /**< strong coupled threshold for aggregate */
    INT AMG_max_aggregation;       /**< max size of each aggregate */
    INT AMG_pair_number;           /**< number of pairs in matching algorithm */
    REAL AMG_quality_bound;        /**< threshold for pair wise aggregation */
    
    // HX preconditioner
    SHORT HX_smooth_iter;            /**< number of smoothing */
    
} input_param; /**< Input parameters */

/**
 * \struct linear_itsolver_param
 * \brief Parameters passed to linear iterative solvers
 */
typedef struct {
    
    SHORT linear_itsolver_type; /**< solver type: see message.h */
    SHORT linear_precond_type;  /**< preconditioner type: see message.h */
    SHORT linear_stop_type;     /**< stopping criteria type */
    INT   linear_maxit;         /**< max number of iterations */
    REAL  linear_tol;           /**< convergence tolerance */
    INT   linear_restart;       /**< number of steps for restarting: for GMRES etc */
    SHORT linear_print_level;   /**< print level: 0--10 */
    
    // HX preconditioner
    SHORT HX_smooth_iter;            /**< number of smoothing */
    
} linear_itsolver_param; /**< Parameters for iterative solvers */


/**
 * \struct AMG_param
 * \brief Parameters for AMG solver
 *
 * \note This is needed for the AMG solver/preconditioner.
 */
typedef struct {
    
    //! type of AMG method
    SHORT AMG_type;
    
    //! print level for AMG
    SHORT print_level;
    
    //! max number of iterations of AMG
    INT maxit;
    
    //! stopping tolerance for AMG solver
    REAL tol;
    
    //! max number of levels of AMG
    SHORT max_levels;
    
    //! max number of coarsest level DOF
    INT coarse_dof;
    
    //! type of AMG cycle
    SHORT cycle_type;
    
    //! quality threshold for pairwise aggregation
    REAL quality_bound;
    
    //! smoother type
    SHORT smoother;
    
    //! smoother order
    SHORT smooth_order; // 1: nature order 2: C/F order (both are symmetric)
    
    //! number of presmoothers
    SHORT presmooth_iter;
    
    //! number of postsmoothers
    SHORT postsmooth_iter;
    
    //! relaxation parameter for SOR smoother
    REAL relaxation;
    
    //! degree of the polynomial smoother
    SHORT polynomial_degree;
    
    //! coarse solver type
    SHORT coarse_solver;
    
    //! switch of scaling of the coarse grid correction
    SHORT coarse_scaling;
    
    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;
    
    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;
    
    //! type of Krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
    
    //! coarsening type
    SHORT coarsening_type;
    
    //! aggregation type
    SHORT aggregation_type;
    
    //! interpolation type
    SHORT interpolation_type;
    
    //! strong connection threshold for coarsening
    REAL strong_threshold;
    
    //! maximal row sum parameter
    REAL max_row_sum;
    
    //! truncation threshold
    REAL truncation_threshold;
    
    //! number of levels use aggressive coarsening
    INT aggressive_level;
    
    //! number of paths use to determine strongly coupled C points
    INT aggressive_path;
    
    //! number of pairwise matchings
    INT pair_number;
    
    //! strong coupled threshold for aggregate
    REAL strong_coupled;
    
    //! max size of each aggregate
    INT max_aggregation;

} AMG_param; /**< Parameters for AMG */


#endif
