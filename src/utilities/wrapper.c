/*! \file src/utilities/wrapper.c
 *
 *  Created by James Adler and Xiaozhe Hu on 09/02/16.
 *  Copyright 2016__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

void python_wrapper_krylov_amg(INT *n,
                                INT *nnz,
                                INT *ia,
                                INT *ja,
                                REAL *a,
                                REAL *b,
                                REAL *u,
                                REAL *tol,
                                INT *maxit,
                                INT *ptrlvl)
{
    /**
     * \fn void python_wrapper_krylov_amg(INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
     *                                     REAL *b, REAL *u, REAL *tol, INT *maxit,
     *                                     INT *ptrlvl)
     *
     * \brief Solve Ax=b by Krylov method preconditioned by AMG
     *
     * \param n       Number of cols of A
     * \param nnz     Number of nonzeros of A
     * \param ia      IA of A in CSR format
     * \param ja      JA of A in CSR format
     * \param a       VAL of A in CSR format
     * \param b       RHS vector
     * \param u       Solution vector
     * \param tol     Tolerance for iterative solvers
     * \param maxit   Max number of iterations
     * \param ptrlvl  Print level for iterative solvers
     *
     * \author Xiaozhe Hu
     * \date   09/02/2016
     */
    
    dCSRmat         mat;      // coefficient matrix
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    linear_itsolver_param  itparam;  // parameters for linear itsolver
    
    input_param inparam;
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);
    
    // Set parameters for linear iterative methods
    param_linear_solver_init(&itparam);
    param_solver_set(&itparam, &inparam);
    //param_linear_solver_print(&linear_itparam);
    
    // Set parameters for algebriac multigrid methods
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    //param_amg_print(&amgparam);
    
    /*
    param_amg_init(&amgparam);
    amgparam.AMG_type             = UA_AMG;
    amgparam.aggregation_type     = VMB;
    
    amgparam.coarse_dof           = 100;
    amgparam.presmooth_iter       = 1;
    amgparam.postsmooth_iter      = 1;
    
    amgparam.strong_coupled       = 0.00;
    amgparam.max_aggregation      = 100;
    
    amgparam.ILU_type             = ILUt;
    amgparam.ILU_levels           = 0;
    amgparam.ILU_lfil             = 0;
    amgparam.ILU_droptol          = 0.01;
    amgparam.ILU_relax            = 0;
    
    param_linear_solver_init(&itparam);
    itparam.linear_itsolver_type = SOLVER_VFGMRES;
    */
    
    amgparam.print_level          = *ptrlvl;
    itparam.linear_tol            = *tol;
    itparam.linear_print_level    = *ptrlvl;
    itparam.linear_maxit          = *maxit;
    
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA = ia;  mat.JA  = ja; mat.val = a;
    
    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;
    
    linear_solver_dcsr_krylov_amg(&mat, &rhs, &sol, &itparam, &amgparam);
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
