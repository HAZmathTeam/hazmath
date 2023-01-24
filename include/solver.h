//
//  solver.h
//
//
//  Created by Hu, Xiaozhe on 5/13/15.
//
// \note added frac. exponent in precond_data (Ana Budisa, 2020-05-13)
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "param.h"

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

    //! setup time of the preconditioner
    REAL setup_time;

} precond; /**< Data for general preconditioner passed to iterative solvers */

/***********************************************************************************************/

/**
 * \struct Schwarz_data
 * \brief Data for Schwarz methods
 *
 * This is needed for the Schwarz preconditioner/smoother.
 */
typedef struct {

    /* matrix information */

    //! pointer to the matrix
    dCSRmat A;  //

    /* blocks information */
    //! number of blocks
    INT nblk;

    //! row index of blocks
    INT *iblock;

    //! column index of blocks
    INT *jblock;

    //! local right hand side
    dvector rhsloc1;

    //! local solution
    dvector xloc1;

    /* THIS SHOULD BE IN SCHWARZ_PARAM: */
    //! Schwarz method type
    //   1 forward, 2 backward, 3 symmetric, 4 additive (global LU of the block matrix);
    //  If local LU of the block matrices is done every iteration, then
    //  11 forward, 12 backward, 13 symmetric, 14 additive
    SHORT Schwarz_type;
    //! Schwarz block solver
    INT blk_solver;

    //! working space size
    INT memt;

    //! mask
    INT *mask;

    //! maximal block size
    INT maxbs;

    //! maximal nnz in a block
    INT maxbnnz;

    //! maxa
    INT *maxa;

    //! matrix for each block during iterations (if Schwarz_type>10)
    //! or for all blocks (if Schwarz_type < 10);
    dCSRmat *blk_data;

    //! UMFPACK or HAZmath factorization place
    void **numeric;

    //! other parameters for the  Schwarz method
    Schwarz_param *swzparam;

} Schwarz_data;

/***********************************************************************************************/

/**
 * \struct AMG_data
 * \brief Data for AMG solvers
 *
 * \note This is needed for the AMG solver/preconditioner.
 */
typedef struct {

    /* Level information */

    //! max number of levels
    SHORT max_levels;

    //! number of levels in use <= max_levels
    SHORT num_levels;

    /* Problem information */
    //! pointer to the matrix at level level_num
    dCSRmat A;

    //! restriction operator at level level_num
    dCSRmat R;

    //! prolongation operator at level level_num
    dCSRmat P;

    //! pointer to the right-hand side at level level_num
    dvector b;

    //! pointer to the iterative solution at level level_num
    dvector x;

    // pointer to the mass matrix at level level_num
    dCSRmat M;

    /* Extra information */
    //! pointer to the numerical factorization from UMFPACK
    void *Numeric;

    //! dimension of the near kernel for SAMG
    INT near_kernel_dim;

    //! basis of near kernel space for SAMG
    REAL **near_kernel_basis;

    //! number of levels use Schwarz smoother
    INT Schwarz_levels;

    //! data of Schwarz smoother
    Schwarz_data Schwarz;

    //! Temporary work space
    dvector w;

    //! cycle type
    INT cycle_type;

    // User defined smoother
    void *wdata;

} AMG_data; /**< Data for AMG */

/***********************************************************************************************/
/**
 * \struct AMG_data_bsr
 * \brief Data for multigrid levels in dBSRmat format
 *
 * \note This structure is needed for the AMG solver/preconditioner in BSR format
 */
typedef struct {

    //! max number of levels
    INT max_levels;

    //! number of levels in use <= max_levels
    INT num_levels;

    //! pointer to the matrix at level level_num
    dBSRmat A;

    //! restriction operator at level level_num
    dBSRmat R;

    //! prolongation operator at level level_num
    dBSRmat P;

    //! pointer to the right-hand side at level level_num
    dvector b;

    //! pointer to the iterative solution at level level_num
    dvector x;

    //! pointer to the diagonal inverse at level level_num
    dvector diaginv;

    //! pointer to the matrix at level level_num (csr format)
    dCSRmat Ac;

    //! pointer to the numerical dactorization from UMFPACK
    void *Numeric;

    //! pointer to the pressure block (only for reservoir simulation)
    dCSRmat PP;

    //! pointer to the auxiliary vectors for pressure block
    REAL *pw;

    //! pointer to the saturation block (only for reservoir simulation)
    dBSRmat SS;

    //! pointer to the auxiliary vectors for saturation block
    REAL *sw;

    //! pointer to the diagonal inverse of the saturation block at level level_num
    dvector diaginv_SS;

    //! ILU data for pressure block
    //ILU_data PP_LU;

    //! pointer to the CF marker at level level_num
    ivector cfmark;

    //! number of levels use ILU smoother
    //INT ILU_levels;

    //! ILU matrix for ILU smoother
    //ILU_data LU;

    //! dimension of the near kernel for SAMG
    INT near_kernel_dim;

    //! basis of near kernel space for SAMG
    REAL **near_kernel_basis;

    //-----------------------------------------
    // extra near kernal space for extra solve
    //! Matrix data for near kernal
    dCSRmat *A_nk;

    //! Prolongation for near kernal
    dCSRmat *P_nk;

    //! Resriction for near kernal
    dCSRmat *R_nk;
    //-----------------------------------------

    //! temporary work space
    dvector w;

} AMG_data_bsr; /**< AMG data for BSR matrices */

/**
 * \struct AMG_data_bdcsr
 * \brief Data for MG solvers
 *
 * \note This is needed for the AMG solver/preconditioner in block_dCSR format.
 */
typedef struct {

    //! max number of levels
    INT max_levels;

    //! number of levels in use <= max_levels
    INT num_levels;

    //! pointer to the matrix at level level_num
    block_dCSRmat A;

    //! restriction operator at level level_num
    block_dCSRmat R;

    //! prolongation operator at level level_num
    block_dCSRmat P;

    //! pointer to the right-hand side at level level_num
    dvector b;

    //! pointer to the iterative solution at level level_num
    dvector x;

    /* Solver information */
    //! Block diagional of A (use those diagonal blocks for coarsening)
    dCSRmat *A_diag;

    //! pointer to the composite matrix (for coarsest level only)
    dCSRmat Ac;

    /* Extra information */
    //! pointer to the numerical factorization from UMFPACK
    void *Numeric;

    //! dimension of the near kernel for SAMG
    INT near_kernel_dim;

    //! basis of near kernel space for SAMG
    REAL **near_kernel_basis;

    /* Smoother information */
    //! AMG data for A_diag blocks
    AMG_data **mgl_diag;

    //! number of levels use Schwarz smoother
    INT Schwarz_levels;

    //! data of Schwarz smoother
    Schwarz_data Schwarz;

    //! Temporary work space
    dvector w;

    //! cycle type
    INT cycle_type;

    //------------------------------------------
    // extra data for specific methods/problems
    // -----------------------------------------
    // data for metric AMG only (A = AD + M)
    // -----------------------------------------
    //! pointer to the AD matrix
    block_dCSRmat *AD;

    //! pointer to the M matrix
    block_dCSRmat *M;

    //! pointer to the interface dof matrix
    dCSRmat *interface_dof;

    //! pointer to the submatrix corresponding to the interface part (block_dCSRmat format)
    block_dCSRmat *A_gamma;

    //! pointer to the submatrix corresponding to the interface part (dBSRmat format)
    dBSRmat A_gamma_bsr;

    //! pointer to the diagonal inverse of the interface part (dBSRmat format only)
    dvector A_gamma_diaginv;

} AMG_data_bdcsr; /**< Data for block MG */

/**
 * \struct MG_blk_data
 * \brief Data for MG solvers
 *
 * \note This is needed for the MG solver/preconditioner.
 */
typedef struct {

    /* Geometric information */
    mesh_struct *fine_level_mesh;

    //! Geometric Type
    INT *gmg_type;

    /* Level information */

    //! max number of levels
    SHORT max_levels;

    //! number of levels in use <= max_levels
    SHORT num_levels;

    /* Problem information */
    //! number of FE spaces
    INT num_spaces;

    //! FE spaces
    block_fespace *FE;

    //! bdry flag stuff
    void (*set_bdry_flags)(mesh_struct*);
    INT *dirichlet;
    INT **dirichlet_blk;

    //! pointer to the matrix at level level_num
    block_dCSRmat A;

    //! pointer to the matrix without dirichlet boundary elimination at level_num
    block_dCSRmat A_noBC;

    //! restriction operator at level level_num
    block_dCSRmat R;

    //! prolongation operator at level level_num
    block_dCSRmat P;

    //! pointer to the right-hand side at level level_num
    dvector b;

    //! pointer to the iterative solution at level level_num
    dvector x;

    /* Solver information */
    //! pointer to the composite matrix (for coarsest level only)
    dCSRmat Ac;

    /* Extra information */
    //! pointer to the numerical factorization from UMFPACK
    void *Numeric;

    //! dimension of the near kernel for SAMG
    INT near_kernel_dim;

    //! basis of near kernel space for SAMG
    REAL **near_kernel_basis;

    /* Smoother information */
    //! Block diagional of A
    dCSRmat *A_diag;

    //! AMG data for A_diag blocks
    AMG_data **mgl;

    //! number of levels use Schwarz smoother
    INT Schwarz_levels;

    //! data of Schwarz smoother
    Schwarz_data Schwarz;

    //! Temporary work space
    dvector w;

    //! cycle type
    INT cycle_type;

    /* Periodic Boundary */
    //! Periodic
    bool periodic_BC;

    //! Ap = Pp^T (A_neuman) Pp
    block_dCSRmat A_periodic;

    //! prolongation to periodic BC
    block_dCSRmat P_periodic;

    //! prolongation to periodic BC
    block_dCSRmat R_periodic;

    //! restriction for going from r_coarse to r_coarse_periodic
    block_dCSRmat R_periodic_scaled;

    /* Other Info */
    //! Stokes system for Biot
    block_dCSRmat *As;
    block_dCSRmat *nAs;
    block_fespace *FES;

} MG_blk_data; /**< Data for block MG */

typedef struct {

    /*!
     * \struct precond_data
     *
     * \brief Data that need to be passed to the preconditioners
     */

    //! type of AMG method
    SHORT AMG_type;

    //! print level in AMG preconditioner
    SHORT print_level;

    //! max number of iterations of AMG preconditioner
    INT   maxit;

    //! max number of AMG levels
    SHORT max_levels;

    //! tolerance for AMG preconditioner
    REAL  tol;

    //! AMG cycle type
    SHORT cycle_type;

    //! AMG smoother type
    SHORT smoother;

    //! number of presmoothing
    SHORT presmooth_iter;

    //! number of postsmoothing
    SHORT postsmooth_iter;

    //! relaxation parameter for SOR smoother
    REAL relaxation;

    //! degree of the polynomial smoother
    SHORT polynomial_degree;

    //! coarse solver type for AMG
    SHORT coarse_solver;

    //! switch of scaling of the coarse grid correction
    SHORT coarse_scaling;

    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;

    //! type of Krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;

    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;

    //! AMG preconditioner data
    AMG_data *mgl_data;

    //! Fractional exponent for fractional smoothers in AMG
    REAL fpwr;

    //! Matrix data
    dCSRmat *A;

    /****************************/
    /*  extra near kernel space */
    /****************************/
    //! Matrix data for near kernel
    dCSRmat *A_nk;

    //! Prolongation for near kernel
    dCSRmat *P_nk;

    //! Restriction for near kernel
    dCSRmat *R_nk;

    /**************************/
    /*  temporary work space  */
    /****************************/
    //! temporary dvector used to store and restore the residual
    dvector *r;

    //! temporary work space for other usage
    REAL *w;


} precond_data; /*! Data for general preconditioner */


/***********************************************************************************************/

typedef struct {

    /*!
     * \struct HX_curl_data
     * \brief Data for HX preconditioner for H(curl) problems
     */

    //! Curl Matrix
    dCSRmat *A;

    /* ---------------------*/
    /* smoother information */
    /* ---------------------*/
    //! Smoother type
    SHORT smooth_type;

    //! number of smoothing
    SHORT smooth_iter;

    /* ---------------------*/
    /* vector Laplacian information */
    /* ---------------------*/
    //! P_curl operator
    dCSRmat *P_curl;

    //! transpose of P_curl operator
    dCSRmat *Pt_curl;

    //! vector Laplacian
    dCSRmat *A_vgrad;

    //! AMG parameters for vector Laplacian
    AMG_param *amgparam_vgrad;

    //! AMG data for vector Laplacian
    AMG_data *mgl_vgrad;

    /* ---------------------*/
    /* scalar Laplacian information */
    /* ---------------------*/
    //! Grad operator
    dCSRmat *Grad;

    //! transpose of Grad operator
    dCSRmat *Gradt;

    //! vector Laplacian
    dCSRmat *A_grad;

    //! AMG parameters for vector Laplacian
    AMG_param *amgparam_grad;

    //! AMG data for vector Laplacian
    AMG_data *mgl_grad;

    /* ---------------------*/
    /* HX preconditioner information */
    /* ---------------------*/
    //! backup residual space
    REAL *backup_r;

    //! temporary work space for other usage
    REAL *w;

} HX_curl_data;

/***********************************************************************************************/

typedef struct {

    /*!
     * \struct HX_div_data
     * \brief Data for HX preconditioner for H(div) problems
     */

    //! Div Matrix
    dCSRmat *A;

    /* ---------------------*/
    /* smoother information */
    /* ---------------------*/
    //! Smoother type
    SHORT smooth_type;

    //! number of smoothing
    SHORT smooth_iter;

    /* ---------------------*/
    /* vector Laplacian information */
    /* ---------------------*/
    //! P_curl operator
    dCSRmat *P_curl;

    //! transpose of P_curl operator
    dCSRmat *Pt_curl;

    //! P_div operator
    dCSRmat *P_div;

    //! transpose of P_div operator
    dCSRmat *Pt_div;

    //! vector Laplacian
    dCSRmat *A_curlgrad;

    //! vector Laplacian
    dCSRmat *A_divgrad;

    //! AMG parameters for vector Laplacian
    AMG_param *amgparam_curlgrad;

    //! AMG data for vector Laplacian
    AMG_data *mgl_curlgrad;

    //! AMG parameters for vector Laplacian
    AMG_param *amgparam_divgrad;

    //! AMG data for vector Laplacian
    AMG_data *mgl_divgrad;

    /* ---------------------*/
    /* scalar Laplacian information */
    /* ---------------------*/
    //! Grad operator
    //dCSRmat *Grad;

    //! transpose of Grad operator
    //dCSRmat *Gradt;

    //! Curl operator
    dCSRmat *Curl;

    //! transpose of Curl operator
    dCSRmat *Curlt;

    //! scalar Laplacian
    dCSRmat *A_grad;

    //! vecror Curl (CtAC)
    dCSRmat *A_curl;

    //! AMG parameters for vector Laplacian
    AMG_param *amgparam_grad;

    //! AMG data for vector Laplacian
    AMG_data *mgl_grad;

    /* ---------------------*/
    /* HX preconditioner information */
    /* ---------------------*/
    //! backup residual space
    REAL *backup_r;

    //! temporary work space for other usage
    REAL *w;

} HX_div_data;

/***********************************************************************************************/
/**
 * \struct precond_diag_bsr
 * \brief Data for diagnal preconditioners in dBSRmat format
 *
 * \note This is needed for the diagnal preconditioner.
 */
typedef struct {

    //! dimension of each sub-block
    INT nb;

    //! diagnal elements
    dvector diag;

} precond_diag_bsr; /**< Data for diagnal preconditioners in dBSRmat format */

/***********************************************************************************************/
/**
 * \struct precond_data_bsr
 * \brief Data for preconditioners in dBSRmat format
 *
 * \note This structure is needed for the AMG solver/preconditioner in BSR format
 */
typedef struct {

    //! type of AMG method
    SHORT AMG_type;

    //! print level in AMG preconditioner
    SHORT print_level;

    //! max number of iterations of AMG preconditioner
    INT maxit;

    //! max number of AMG levels
    INT max_levels;

    //! tolerance for AMG preconditioner
    REAL tol;

    //! AMG cycle type
    SHORT cycle_type;

    //! AMG smoother type
    SHORT smoother;

    //! AMG smoother ordering
    SHORT smooth_order;

    //! number of presmoothing
    SHORT presmooth_iter;

    //! number of postsmoothing
    SHORT postsmooth_iter;

    //! coarsening type
    SHORT coarsening_type;

    //! relaxation parameter for SOR smoother
    REAL relaxation;

    //! coarse solver type for AMG
    SHORT coarse_solver;

    //! switch of scaling of the coarse grid correction
    SHORT coarse_scaling;

    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;

    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;

    //! smooth factor for smoothing the tentative prolongation
    REAL tentative_smooth;

    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;

    //! AMG preconditioner data
    AMG_data_bsr *mgl_data;

    //! AMG preconditioner data for pressure block
    AMG_data *pres_mgl_data;

    //! ILU preconditioner data (needed for CPR type preconditioner)
    //ILU_data *LU;

    //! Matrix data
    dBSRmat *A;

    // extra near kernal space

    //! Matrix data for near kernal
    dCSRmat *A_nk;

    //! Prolongation for near kernal
    dCSRmat *P_nk;

    //! Resriction for near kernal
    dCSRmat *R_nk;

    //! temporary dvector used to store and restore the residual
    dvector r;

    //! temporary work space for other usage
    REAL *w;

} precond_data_bsr; /**< Data for preconditioners in dBSRmat format */

/***********************************************************************************************/
/**
 * \struct precond_data_bdcsr
 * \brief Data for preconditioners in block_dCSRmat format
 *
 * \note This structure is needed for the AMG solver/preconditioner in block_dCSRmat format
 */
typedef struct {

    //! type of AMG method
    SHORT AMG_type;

    //! print level in AMG preconditioner
    SHORT print_level;

    //! max number of iterations of AMG preconditioner
    INT maxit;

    //! max number of AMG levels
    INT max_levels;

    //! tolerance for AMG preconditioner
    REAL tol;

    //! AMG cycle type
    SHORT cycle_type;

    //! AMG smoother type
    SHORT smoother;

    //! AMG smoother ordering
    SHORT smooth_order;

    //! number of presmoothing
    SHORT presmooth_iter;

    //! number of postsmoothing
    SHORT postsmooth_iter;

    //! coarsening type
    SHORT coarsening_type;

    //! relaxation parameter for SOR smoother
    REAL relaxation;

    //! coarse solver type for AMG
    SHORT coarse_solver;

    //! switch of scaling of the coarse grid correction
    SHORT coarse_scaling;

    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;

    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;

    //! smooth factor for smoothing the tentative prolongation
    REAL tentative_smooth;

    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;

    //! AMG preconditioner data
    AMG_data_bdcsr *mgl_data;

    //! Schwarz method data
    Schwarz_data *schwarz_data;

    //! Schwarz method parameters
    Schwarz_param *schwarz_param;

    //! LU data
    void **LU_data;

    //! Matrix data
    block_dCSRmat *A;

    //! Total size of of the matrix
    INT total_row;
    INT total_col;

    // extra near kernal space
    //! Matrix data for near kernal
    //dCSRmat *A_nk;

    //! Prolongation for near kernal
    //dCSRmat *P_nk;

    //! Resriction for near kernal
    //dCSRmat *R_nk;

    //! temporary dvector used to store and restore the residual
    dvector r;

    //! temporary work space for other usage
    REAL *w;

    //! temporary vector space for permutation that assembles the block format
    ivector perm;

} precond_data_bdcsr; /**< Data for preconditioners in block_dCSRmat format */


/***********************************************************************************************/

/**
 * \brief Data passed to the preconditioner for block preconditioning for block_dCSRmat format
 *
 * This is needed for the block preconditioner.
 */
typedef struct {

    /*-------------------------------------*/
    /* Basic data for block preconditioner */
    /*-------------------------------------*/
    block_dCSRmat *Abcsr; /**< problem data, the blocks */

    dCSRmat *A_diag;      /**< data for each diagonal block*/

    dvector r;            /**< temp work space */

    /*------------------------------*/
    /* Data for the diagonal blocks */
    /*------------------------------*/
    INT *block_solve_type;  /**<  how to solve each block  */

    /*--- solve by direct solver ---*/
    void **LU_diag;       /**< LU decomposition for the diagonal blocks (for UMFpack) */

    /*--- solve by diagonal preconditioner ---*/
    dvector **diag;

    /*---  solve by AMG ---*/
    AMG_data **mgl;       /**< AMG data for the diagonal blocks */
    AMG_param *amgparam;  /**< parameters for AMG */

    /*---  solve by GMG ---*/
    MG_blk_data *bmgl;    /**< Block MG data for monolithic */

    /*--- solver by HX preconditioner */
    HX_curl_data **hxcurldata; /**< HX data for the diagonal CURL blocks */
    HX_div_data **hxdivdata; /**< HX data for the diagonal DIV blocks */

    /*------------------------------*/
    /* Data for mesh infomation */
    /*------------------------------*/
    dvector *el_vol;   /**< volume of each element */

    /*------------------------------*/
    /* Data for Maxwell problem only!! */
    /*------------------------------*/
    dCSRmat *G;         /**< scaled gradiend operator */
    dCSRmat *K;         /**< scaled curl operator */
    dCSRmat *Gt;        /**< scaled transpose gradiend operator */
    dCSRmat *Kt;        /**< scaled transpose gradiend operator */

    /*------------------------------*/
    /* Data for fractional problem only! */
    /*------------------------------*/
    //dCSRmat A;
    dCSRmat *scaled_M;   /**< scaled Mass matrix */
    dvector *diag_scaled_M; /**< diagonal of scaled mass matrix */
    REAL scaled_alpha;   /**< scaled alpha */
    REAL scaled_beta;    /**< scaled beta */
    dvector *poles;      /**< poles for rational approximation */
    dvector *residues;   /**< residues for rational approximation */

} precond_block_data; /**< Precond data for block matrices */


/***********************************************************************************************/

/*!
 * \struct precond_ra_data
 *
 * \brief Data on rational approximation
 *        that need to be passed to the preconditioners
 */
typedef struct {

    /*--- solve by direct solver ---*/
    void **LU_diag;    /**< LU decomposition for shifted Laplacians (UMFpack) */

    /*---  solve by AMG ---*/
    AMG_data **mgl;       /**< AMG data for shifted Laplacians */
    AMG_param *amgparam;  /**< parameters for AMG */

    /*-----------------------------*/
    /* Data for fractional problem */
    /*-----------------------------*/
    dCSRmat *scaled_A;      /**< scaled stiffness matrix */
    dCSRmat *scaled_M;      /**< scaled Mass matrix */
    dvector *diag_scaled_M; /**< diagonal of scaled mass matrix */
    REAL scaled_alpha;      /**< scaled alpha coeff (goes with s_power) */
    REAL scaled_beta;       /**< scaled beta coeff (goes with t_power) */
    REAL s_power;           /**< first fractionality (goes with alpha) */
    REAL t_power;           /**< second fractionality (goes with beta) */
    dvector *poles;         /**< poles for rational approximation */
    dvector *residues;      /**< residues for rational approximation */

    /*------------------------*/
    /*  temporary work space  */
    /*------------------------*/
    dvector *r;     /**< temporary dvector to store and restore the residual */
    REAL *w;        /**< temporary work space for other usage */

} precond_ra_data; /**< Data for RA preconditioner */


/***********************************************************************************************/


/**
 * \struct matvec
 * \brief Matrix-vector multiplication
 */
typedef struct {

    //! data for Matrix-vector multiplication
    void *data;

    //! action for Matrix-vector, should be a pointer to a function
    //void (*fct)(void *, REAL *, REAL *);
    void (*fct)(REAL *, REAL *, void *);

} matvec; /**< Data for general Matrix-vector multiplication */


/**
 * \struct solve_stats
 * \brief statistics about the solve
 */
typedef struct{

    INT iteration_count;

    REAL time_setup;
    REAL time_precondition_setup;
    REAL time_solve;

    REAL rho1;
    REAL rho2;

    // Lazy way to pass in variables that I need.
    block_dCSRmat *As;
    block_dCSRmat *nAs;
    block_fespace *FES;

} solve_stats; /**< statistics about solve */

/**
 * \struct smoother_matvec
 * \brief  Smoother matrix-vector multiplication
 */
typedef struct {
    //! smoother type
    SHORT type;

    //! data for the smoother (e.g. smoother_data type)
    void *data;

    //! action of smoother as some matrix-vector application
    //!! void pointer, but usually function pointer of type
    //!! void (*fct)(REAL *, REAL *, void *)
    //!! unless its a e.g. python function
    void *fct;


} smoother_matvec; /**< data for smoother matvec multiplication */


/**
 * \struct smoother_data
 * \brief  data for smoother application
 */
typedef struct {

    //! starting index
    INT istart;

    //! ending index
    INT iend;

    //! step size
    INT istep;

    //! number of smoother iterations
    INT nsweeps;

    //! optional relaxation parameter
    REAL relax;

    //! smoother matrix
    dCSRmat *A;

} smoother_data; /**< data for smoother application */


#endif
