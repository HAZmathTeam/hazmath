//
//  param.h
//
//
//  Created by Hu, Xiaozhe on 06/13/15.
//
//
// \note added frac. exponent in input_param and AMG_param (Ana Budisa, 2020-05-13)

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
    // mesh parameters
    //--------------------------
    INT read_mesh_from_file; /**< read from file or use built-in generator */
    INT spatial_dim; /**< dimension of computational domain */
    INT refinement_type; /**< Type of refinement >10 -> uniform */
    INT refinement_levels; /** < Since of grid = 2^{refinement_levels+1} vertices in each direction */
    INT boundary_codes; //** 1 to be default

    //--------------------------
    // finite element parameters
    //--------------------------
    // genearal parameters
    INT nquad;          /**< quadature nodes in each direction */

    // parameters for H(D) equations (examples only)
    INT FE_type;        /**< finite element type */

    // Mass lumping
    INT Mass_lump;     /** < whether to do mass lumping or not */

    //----------------------------
    // time steppng paramters
    //----------------------------
    REAL time_start;      /**< start time */
    INT  time_step_type;  /**< time step type */
    INT  time_steps;      /**< time steps */
    REAL time_step_size; /**< time step type */
    INT  rhs_time_dep;    /**< indicates if rhs is time-dependent */

    //----------------------------
    // nonlinear solver parameters
    //----------------------------
    INT  nonlinear_itsolver_type;      /**< type of nonlinear solver */
    INT  nonlinear_itsolver_maxit;   /**< maximal iterations of nonlinear solver*/
    REAL nonlinear_itsolver_tol;    /**< tolerance for nonlinear solver */
    INT  nonlinear_itsolver_toltype; /**< type of stopping tolerance for nonlinear solver */
    INT fas_presmoothers; /* Number of presmoothing steps for FAS */
    INT fas_postsmoothers; /* Number of postsmoothing steps for FAS */
    INT fas_smooth_tol; /* Stopping tolerance for nonlinear smoother */

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

    // parameter for Schwarz method
    INT Schwarz_mmsize;  /**< maximal block size */
    INT Schwarz_maxlvl;  /**< maximal levels */
    INT Schwarz_type;    /**< type of Schwarz method */
    INT Schwarz_blksolver; /**< type of Schwarz block solver */

    // Algebraic Multigrid
    SHORT AMG_type;                /**< Type of AMG */
    SHORT AMG_levels;              /**< maximal number of levels */
    SHORT AMG_cycle_type;          /**< type of cycle */
    SHORT AMG_smoother;            /**< type of smoother */
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
    REAL  AMG_fpwr;                 /**< fractional exponent for fractional smoothers */

    // Unsmoothed Aggregation AMG (UA AMG)
    SHORT AMG_aggregation_type;    /**< aggregation type */
    REAL  AMG_strong_coupled;       /**< strong coupled threshold for aggregate */
    INT   AMG_max_aggregation;       /**< max size of each aggregate */

    // Smoothed Aggregation AMG (SA AMG)
    SHORT AMG_smooth_filter;       /**< use filter for smoothing the tentative */
    REAL AMG_tentative_smooth;     /**< relaxation factor for smoothing the tentative prolongation */

    // HX preconditioner
    SHORT HX_smooth_iter;            /**< number of smoothing */

    // BSR preconditioner
    REAL BSR_alpha;                 /**< weight on diagonal matrix alpha*D approx of A */
    REAL BSR_omega;                 /**< weight on update x = x + omega*Binv*(Ax-b) */

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

    // scaling parameter used in Argumented Lagrange type block preconditioners
    REAL AL_scaling_param;

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

    //! smoother type
    SHORT smoother;

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

    //! aggregation type
    SHORT aggregation_type;

    //! strong coupled threshold for aggregate
    REAL strong_coupled;

    //! max size of each aggregate
    INT max_aggregation;

    //! switch for filtered matrix used for smoothing the tentative prolongation
    SHORT smooth_filter;

    //! relaxation parameter for smoothing the tentative prolongation
    REAL tentative_smooth;

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

    /* Hacking in parameters for gmg smoothers */
    //! HAZMATH install dir
    char* HAZDIR;
    //! Track if schwarz should be used as a relaxation method on a block
    INT* Schwarz_on_blk;
    INT* Schwarz_patch_type;

    //! damping parameter for relaxation
    REAL damping_param;

    // BSR preconditioner
    REAL BSR_alpha;                 /**< weight on diagonal matrix alpha*D approx of A */
    REAL BSR_omega;                 /**< weight on update x = x + omega*Binv*(Ax-b) */

    // Fractional exponent for fractional smoothers
    REAL fpwr;

    // User defined smoother
    void *smoother_function;

} AMG_param; /**< Parameters for AMG */

/**
 * \struct GMG_param
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

    //! smoother type
    SHORT smoother;

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

    //! aggregation type
    SHORT aggregation_type;

    //! strong coupled threshold for aggregate
    REAL strong_coupled;

    //! max size of each aggregate
    INT max_aggregation;

    //! switch for filtered matrix used for smoothing the tentative prolongation
    SHORT smooth_filter;

    //! relaxation parameter for smoothing the tentative prolongation
    REAL tentative_smooth;

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

} GMG_param; /**< Parameters for AMG */

/*!
 * \struct Schwarz_param
 * \brief Parameters for Schwarz method
 *
 */
typedef struct {

    //! print leve
    SHORT print_level;

    //! type for Schwarz method
    SHORT Schwarz_type;

    //! maximal level for constructing the blocks
    INT Schwarz_maxlvl;

    //! maximal size of blocks
    INT Schwarz_mmsize;

    //! type of Schwarz block solver
    INT Schwarz_blksolver;

    //! gmg stuff
    INT *patch_type_gmg;

} Schwarz_param; /**< Parameters for Schwarz method */

#endif
