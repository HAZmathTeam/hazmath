/*! \file src/utilities/wrapper.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 09/02/16.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmath.h"

/*************************************************************************************/
/*!
 * \fn void python_wrapper_krylov_amg(INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                                     REAL *b, REAL *u, REAL *tol, INT *maxit,
 *                                     INT *ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by AMG (this is an interface with PYTHON)
 *
 * \param n             Number of cols of A
 * \param nnz           Number of nonzeros of A
 * \param ia            IA of A in CSR format
 * \param ja            JA of A in CSR format
 * \param a             VAL of A in CSR format
 * \param b             RHS vector
 * \param u             Solution vector
 * \param tol           Tolerance for iterative solvers
 * \param maxit         Max number of iterations
 * \param print_lvl     Print level for iterative solvers
 *
 * \author Xiaozhe Hu
 * \date   09/02/2016
 *
 */
void python_wrapper_krylov_amg(INT *n,
                               INT *nnz,
                               INT *ia,
                               INT *ja,
                               REAL *a,
                               REAL *b,
                               REAL *u,
                               REAL *tol,
                               INT *maxit,
                               INT *print_lvl)
{
    dCSRmat         mat;      // coefficient matrix
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    linear_itsolver_param  itparam;  // parameters for linear itsolver
    
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);
    
    // Set parameters for linear iterative methods
    param_linear_solver_init(&itparam);
    param_linear_solver_set(&itparam, &inparam);
    if (*print_lvl > PRINT_MIN) param_linear_solver_print(&itparam);
    
    // Set parameters for algebriac multigrid methods
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (*print_lvl > PRINT_MIN) param_amg_print(&amgparam);
        
    amgparam.print_level          = *print_lvl;
    itparam.linear_tol            = *tol;
    itparam.linear_print_level    = *print_lvl;
    itparam.linear_maxit          = *maxit;
    
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA = ia;  mat.JA  = ja; mat.val = a;
    
    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;
    
    linear_solver_dcsr_krylov_amg(&mat, &rhs, &sol, &itparam, &amgparam);

}

/*!
 * \fn void python_wrapper_direct(INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                                     REAL *b, REAL *u, INT *ptrlvl)
 *
 * \brief Solve Ax=b by direct method (this is an interface with PYTHON)
 *
 * \param n             Number of cols of A
 * \param nnz           Number of nonzeros of A
 * \param ia            IA of A in CSR format
 * \param ja            JA of A in CSR format
 * \param a             VAL of A in CSR format
 * \param b             RHS vector
 * \param u             Solution vector
 * \param print_lvl     Print level for iterative solvers
 *
 * \author Xiaozhe Hu
 * \date   08/03/2018
 *
 */
void python_wrapper_direct(INT *n,
                           INT *nnz,
                           INT *ia,
                           INT *ja,
                           REAL *a,
                           REAL *b,
                           REAL *u,
                           INT *print_lvl)
{
    dCSRmat         mat;      // coefficient matrix
    dvector         rhs, sol; // right-hand-side, solution

    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA = ia;  mat.JA  = ja; mat.val = a;

    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;

    directsolve_UMF(&mat, &rhs, &sol, *print_lvl);
}


/*!
 * \fn void python_wrapper_krylov_block_2(INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                                     REAL *b, REAL *u, REAL *tol, INT *maxit,
 *                                     INT *ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned in 2 by 2 block form (this is an interface with PYTHON)
 *
 * \param n             Number of cols of A
 * \param nnz           Number of nonzeros of A
 * \param ia            IA of A in CSR format
 * \param ja            JA of A in CSR format
 * \param a             VAL of A in CSR format
 * \param b             RHS vector
 * \param u             Solution vector
 * \param tol           Tolerance for iterative solvers
 * \param maxit         Max number of iterations
 * \param print_lvl     Print level for iterative solvers
 *
 * \author Xiaozhe Hu
 * \date   08/08/2018
 *
 */
void python_wrapper_krylov_block_2(INT *n,
                               INT *nnz,
                               INT *ia,
                               INT *ja,
                               REAL *a,
                               REAL *b,
                               REAL *u,
                               REAL *tol,
                               INT *maxit,
                               INT *print_lvl)
{
    dCSRmat         mat_csr;      // coefficient matrix in csr format
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    linear_itsolver_param  itparam;  // parameters for linear itsolver

    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);

    // Set parameters for linear iterative methods
    param_linear_solver_init(&itparam);
    param_linear_solver_set(&itparam, &inparam);
    if (*print_lvl > PRINT_MIN) param_linear_solver_print(&itparam);

    // Set parameters for algebriac multigrid methods
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (*print_lvl > PRINT_MIN) param_amg_print(&amgparam);

    amgparam.print_level          = *print_lvl;
    itparam.linear_tol            = *tol;
    itparam.linear_print_level    = *print_lvl;
    itparam.linear_maxit          = *maxit;

    // form CSR matrix
    mat_csr.row = *n; mat_csr.col = *n; mat_csr.nnz = *nnz;
    mat_csr.IA = ia;  mat_csr.JA  = ja; mat_csr.val = a;

    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;

    // convert into 2 by 2 block CSR matrix
    int *bsize;
    bsize = (int *)calloc(2, sizeof(int));
    bsize[0] = mat_csr.row/2; bsize[1] = mat_csr.row/2;

    block_dCSRmat mat_bdcsr = dcsr_2_bdcsr(&mat_csr, 2, bsize);

    // get diagonal blocks
    dCSRmat *mat_diag;
    mat_diag = (dCSRmat *)calloc(2, sizeof(dCSRmat));

    dcsr_alloc(mat_bdcsr.blocks[0]->row, mat_bdcsr.blocks[0]->col, mat_bdcsr.blocks[0]->nnz, &mat_diag[0]);
    dcsr_cp(mat_bdcsr.blocks[0], &mat_diag[0]);

    dcsr_alloc(mat_bdcsr.blocks[3]->row, mat_bdcsr.blocks[3]->col, mat_bdcsr.blocks[3]->nnz, &mat_diag[1]);
    dcsr_cp(mat_bdcsr.blocks[3], &mat_diag[1]);


    // solve in 2 by 2 block form
    linear_solver_bdcsr_krylov_block_2(&mat_bdcsr, &rhs, &sol, &itparam, &amgparam, mat_diag);

    // clean memory
    bdcsr_free(&mat_bdcsr);
    dcsr_free( &mat_diag[0]);
    dcsr_free( &mat_diag[1]);
    if(mat_diag) free(mat_diag);

}


/***************************** END ***************************************************/

