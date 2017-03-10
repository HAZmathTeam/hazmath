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
/***************************** END ***************************************************/

