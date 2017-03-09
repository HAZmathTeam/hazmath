/*! \file src/utilities/message.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note modified by Xiaozhe Hu 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
void print_itsolver_info(const INT  print_lvl,
                         const INT  stop_type,
                         const INT  iter,
                         const REAL rel_res,
                         const REAL abs_res,
                         const REAL factor)
{
  /*!
     * \fn void print_itsolver_info (const INT print_lvl, const INT stop_type, const INT iter,
     *                        const REAL rel_res, const REAL abs_res, const REAL factor)
     *
     * \brief Print out iteration information for linear iterative solvers at each iteration
     *
     * \param print_lvl     how much information to print (higher number means print more information)
     * \param stop_type     Type of stopping criteria
     * \param iter          Number of iterations
     * \param rel_res       Relative residual (different for different stop_type)
     * \param abs_res       Absolute residual (different for different stop_type)
     * \param factor        Contraction factor at each iteration
     *
     */

  if ( print_lvl >= PRINT_SOME ) {

    if ( iter > 0 ) { // iter > 0: not the first iteration
      printf("%6d | %13.6e   | %13.6e  | %10.4f\n", iter, rel_res, abs_res, factor);
    }
    else { // iter = 0: first iteration
      printf("-----------------------------------------------------------\n");
      switch (stop_type) {
      case STOP_REL_RES:
        printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
        break;
      case STOP_REL_PRECRES:
        printf("It Num | ||r||_B/||b||_B |    ||r||_B     |  Conv. Factor\n");
        break;
      case STOP_MOD_REL_RES:
        printf("It Num |   ||r||/||x||   |     ||r||      |  Conv. Factor\n");
        break;
      }
      printf("-----------------------------------------------------------\n");
      printf("%6d | %13.6e   | %13.6e  |     -.-- \n", iter, rel_res, abs_res);
    } // end if iter

  } // end if print_lvl
}

/***********************************************************************************************/
void print_cputime (const char *message,
                    const REAL cputime)
{

  /*!
     * \fn void print_cputime (const char *message, const REAL cputime)
     *
     * \brief Print CPU walltime
     *
     * \param message   Pointer to the string to print out
     * \param cputime   Walltime since start to end
     *
     */

  printf("%s costs %.4f seconds.\n", message, cputime);
}

/***********************************************************************************************/
void print_message (const INT print_lvl,
                    const char *message)
{
  /*!
     * \fn void print_message (const INT ptrlvl, const char *message)
     *
     * \brief Print out the message if necessary
     *
     * \param print_lvl   Level for output
     * \param message     Pointer to the error message
     *
     */

  if ( print_lvl > PRINT_NONE ) printf("%s", message);
}

/***********************************************************************************************/
void print_amg_complexity (AMG_data *mgl,
                           const SHORT print_lvl)
{
  /*!
     * \fn void print_amg_complexity (AMG_data *mgl, const SHORT prtlvl)
     *
     * \brief Print grid and operator complexity of AMG method
     *
     * \param mgl      Multilevel structure for AMG
     * \param print_lvl   How much information to print
     *
     */

  const SHORT   max_levels = mgl->num_levels;
  SHORT         level;
  REAL          grid_complexity=0.0, operator_complexity=0.0;

  if ( print_lvl >= PRINT_SOME ) {

    printf("-----------------------------------------------------------\n");
    printf("  Level   Num of rows   Num of nonzeros   Avg. NNZ / row   \n");
    printf("-----------------------------------------------------------\n");

    for ( level = 0; level < max_levels; ++level) {
      REAL AvgNNZ = (REAL) mgl[level].A.nnz/mgl[level].A.row;
      printf("%5d %13d %17d %14.2f\n",
             level, mgl[level].A.row, mgl[level].A.nnz, AvgNNZ);
      grid_complexity    += mgl[level].A.row;
      operator_complexity += mgl[level].A.nnz;
    }
    printf("-----------------------------------------------------------\n");

    grid_complexity     /= mgl[0].A.row;
    operator_complexity /= mgl[0].A.nnz;
    printf("  Grid complexity = %.3f  |", grid_complexity);
    printf("  Operator complexity = %.3f\n", operator_complexity);

    printf("-----------------------------------------------------------\n");

  }
}


/***********************************************************************************************/
void check_error (const SHORT status,
                  const char *func_name)
{

  /*!
     * \fn void check_error (const SHORT status, const char *fctname)
     *
     * \brief Check error status and print out error messages before quit
     *
     * \param status   Error status
     * \param fctname  piinter to the function name where this routine is called
     *
     */

  if ( status >= SUCCESS ) return; // No error at all

  switch ( status ) {
  case ERROR_DIM:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- You have entered the Twilight Zone (dimension is not 1-3) !!! \n\n", func_name);
    break;
  case ERROR_FE_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- The FE Space you want is not implemented !!! \n\n", func_name);
    break;
  case ERROR_OPEN_FILE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Cannot open file !!!\n\n", func_name);
    break;
  case ERROR_WRONG_FILE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Wrong file format !!!\n\n", func_name);
    break;
  case ERROR_INPUT_PAR:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Wrong input arguments !!!\n\n", func_name);
    break;
  case ERROR_REGRESS:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Regression test failed !!!\n\n", func_name);
    break;
  case ERROR_ALLOC_MEM:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Cannot allocate memory !!!\n\n", func_name);
    break;
  case ERROR_NUM_BLOCKS:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Wrong number of blocks !!!\n\n", func_name);
    break;
  case ERROR_DATA_STRUCTURE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Data structure mismatch !!!\n\n", func_name);
    break;
  case ERROR_DATA_ZERODIAG:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Matrix has zero diagonal entries !!!\n\n", func_name);
    break;
  case ERROR_DUMMY_VAR:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unexpected input argument !!!\n\n", func_name);
    break;
  case ERROR_AMG_INTERP_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown AMG interpolation type !!!\n\n", func_name);
    break;
  case ERROR_AMG_COARSE_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown AMG coarsening type !!!\n\n", func_name);
    break;
  case ERROR_AMG_SMOOTH_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown AMG smoother type !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown solver type !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_PRECTYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown preconditioner type !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_STAG:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Solver stagnation error !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_SOLSTAG:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Solution is close to zero !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_TOLSMALL:
    printf("\n!!! ERROR HAZMAT DANGER: in function '%s' -- Tol is too small for the solver !!!\n\n", func_name);
      break;
  case ERROR_SOLVER_MAXIT:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Max iteration number reached !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_EXIT:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Solver exited unexpected !!!\n\n", func_name);
    break;
  case ERROR_SOLVER_MISC:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown solver runtime error !!!\n\n", func_name);
    break;
  case ERROR_MISC:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown error occurred !!!\n\n", func_name);
    break;
  case ERROR_QUAD_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown quadrature rules !!!\n\n", func_name);
    break;
  case ERROR_QUAD_DIM:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Num of quad points is not supported !!!\n\n", func_name);
    break;
  case ERROR_UNKNOWN:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Function does not exit successfully !!!\n\n", func_name);
    break;
  case ERROR_MAT_DOF:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Your matrix size doesn't match your number of DOF !!!\n\n", func_name);
    break;
  case ERROR_TS_TYPE:
    printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- The time-stepping scheme you want is not implemented !!!\n\n", func_name);
    break;
  default:
    break;
  }

  exit(status);
}
