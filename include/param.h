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
    SHORT linear_stop_type;           /**< stop type of linear iterative solver */
    INT linear_restart;                      /**< restart number used in GMRES */
    
    // Preconditioner
    INT linear_precond_type;                 /**< type of preconditioner for iterative solvers */
    
    // parameters for ILU
    SHORT ILU_type;      /**< ILU type for decomposition*/
    INT ILU_lfil;        /**< level of fill-in */
    REAL ILU_droptol;    /**< drop tolerance */
    REAL ILU_relax;      /**< scaling factor: add the sum of dropped entries to diagonal */
    REAL ILU_permtol;    /**< permutation tolerance */
    
    // Algebraic Multigrid
    SHORT AMG_type;                /**< Type of AMG */
    SHORT AMG_levels;              /**< maximal number of levels */
    SHORT AMG_cycle_type;          /**< type of cycle */
    SHORT AMG_smoother;            /**< type of smoother */
    SHORT AMG_smooth_order;        /**< order for smoothers */
    REAL AMG_relaxation;           /**< over-relaxation parameter for SOR */
    //SHORT AMG_polynomial_degree;   /**< degree of the polynomial smoother */
    SHORT AMG_presmooth_iter;      /**< number of presmoothing */
    SHORT AMG_postsmooth_iter;     /**< number of postsmoothing */
    INT AMG_coarse_dof;            /**< max number of coarsest level DOF */
    REAL AMG_tol;                  /**< tolerance for AMG if used as preconditioner */
    INT AMG_maxit;                 /**< number of iterations for AMG used as preconditioner */
    SHORT AMG_ILU_levels;          /**< how many levels use ILU smoother */
    SHORT AMG_coarse_solver;       /**< coarse solver type */
    SHORT AMG_coarse_scaling;      /**< switch of scaling of the coarse grid correction */
    SHORT AMG_amli_degree;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type; /**< type of Krylov method used by nonlinear AMLI cycle */
    //INT AMG_Schwarz_levels;        /**< number of levels use Schwarz smoother */
    
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
    //REAL AMG_tentative_smooth;     /**< relaxation factor for smoothing the tentative prolongation */
    //SHORT AMG_smooth_filter;       /**< use filter for smoothing the tentative prolongation or not */
    
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
 * \struct ILU_param
 * \brief Parameters for ILU
 */
typedef struct {
    
    //! print level
    SHORT print_level;
    
    //! ILU type for decomposition
    SHORT ILU_type;
    
    //! level of fill-in for ILUk
    INT ILU_lfil;
    
    //! drop tolerance for ILUt
    REAL ILU_droptol;
    
    //! add the sum of dropped elements to diagonal element in proportion relax
    REAL ILU_relax;
    
    //! permuted if permtol*|a(i,j)| > |a(i,i)|
    REAL ILU_permtol;
    
} ILU_param; /**< Parameters for ILU */

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
    
    //! relaxation parameter for smoothing the tentative prolongation
    REAL tentative_smooth;
    
    //! switch for filtered matrix used for smoothing the tentative prolongation
    SHORT smooth_filter;
    
    //! number of levels use ILU smoother
    SHORT ILU_levels;
    
    //! ILU type for smoothing
    SHORT ILU_type;
    
    //! level of fill-in for ILUs and ILUk
    INT ILU_lfil;
    
    //! drop tolerance for ILUt
    REAL ILU_droptol;
    
    //! relaxation for ILUs
    REAL ILU_relax;
    
    //! permuted if permtol*|a(i,j)| > |a(i,i)|
    REAL ILU_permtol;
    
    //! number of levels use Schwarz smoother
    INT Schwarz_levels;
    
    //! maximal block size
    INT Schwarz_mmsize;
    
    //! maximal levels
    INT Schwarz_maxlvl;
    
    //! type of Schwarz method
    INT Schwarz_type;
    
    //! type of Schwarz block solver
    INT Schwarz_blksolver;
    
} AMG_param; /**< Parameters for AMG */


#endif
