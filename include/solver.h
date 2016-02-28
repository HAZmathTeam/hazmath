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
    
} precond; /**< Data for general preconditioner passed to iterative solvers */

/**
 * \struct ILU_data
 * \brief Data for ILU setup
 */
typedef struct {
    
    //! row number of matrix LU, m
    INT row;
    
    //! column of matrix LU, n
    INT col;
    
    //! number of nonzero entries
    INT nzlu;
    
    //! integer array of row pointers and column indexes, the size is nzlu
    INT *ijlu;
    
    //! nonzero entries of LU
    REAL *luval;
    
    //! block size for BSR type only
    INT nb;
    
    //! work space size
    INT nwork;
    
    //! work space
    REAL *work;
    
} ILU_data; /**< Data for ILU */

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
    
    /* Extra information */
    
    //! pointer to the numerical factorization from UMFPACK
    void *Numeric;
    
    //! pointer to the CF marker at level level_num
    ivector cfmark;
    
    //! number of levels use ILU smoother
    INT ILU_levels;
    
    //! ILU matrix for ILU smoother
    ILU_data LU;
    
    //! dimension of the near kernel for SAMG
    INT near_kernel_dim;
    
    //! basis of near kernel space for SAMG
    REAL **near_kernel_basis;
    
    // Smoother order information
    
    //! number of levels use Schwarz smoother
    INT Schwarz_levels;
    
    //! data of Schwarz smoother
    // Schwarz_data Schwarz;
    
    //! Temporary work space
    dvector w;
    
    //! data for MUMPS
    // Mumps_data mumps;
    
    //! cycle type
    INT cycle_type;
    
} AMG_data; /**< Data for AMG */

/**
 * \struct precond_data
 * \brief Data passed to the preconditioners
 */
typedef struct {
    
    //! type of AMG method
    SHORT AMG_type;
    
    //! print level in AMG preconditioner
    SHORT print_level;
    
    //! max number of iterations of AMG preconditioner
    INT maxit;
    
    //! max number of AMG levels
    SHORT max_levels;
    
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
    
    //! relaxation parameter for SOR smoother
    REAL relaxation;
    
    //! degree of the polynomial smoother
    SHORT polynomial_degree;
    
    //! switch of scaling of the coarse grid correction
    SHORT coarsening_type;
    
    //! coarse solver type for AMG
    SHORT coarse_solver;
    
    //! switch of scaling of the coarse grid correction
    SHORT coarse_scaling;
    
    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;
    
    //! type of Krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
    
    //! smooth factor for smoothing the tentative prolongation
    REAL tentative_smooth;
    
    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;
    
    //! AMG preconditioner data
    AMG_data *mgl_data;
    
    //! ILU preconditioner data (needed for CPR type preconditioner)
    ILU_data *LU;
    
    //! Matrix data
    dCSRmat *A;
    
    // extra near kernel space
    
    //! Matrix data for near kernel
    dCSRmat *A_nk;
    
    //! Prolongation for near kernel
    dCSRmat *P_nk;
    
    //! Restriction for near kernel
    dCSRmat *R_nk;
    
    // temporary work space
    
    //! temporary dvector used to store and restore the residual
    dvector r;
    
    //! temporary work space for other usage
    REAL *w;
    
} precond_data; /**< Data for general preconditioner */


/**
 * \struct HX_curl_data
 * \brief Data for HX preconditioner for H(curl) problems
 */
typedef struct {
    
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
    /*--- solve by direct solver ---*/
    void **LU_diag;       /**< LU decomposition for the diagonal blocks (for UMFpack) */
    
    /*---  solve by AMG ---*/
    AMG_data **mgl;       /**< AMG data for the diagonal blocks */
    AMG_param *amgparam;  /**< parameters for AMG */
    
    /*--- solver by HX preconditioner */
    HX_curl_data **hxcurldata; /**< HX data for the diagonal blocks */
    
} precond_block_data; /**< Precond data for block matrices */

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
