/*! \file src/utilities/message.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 3/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note modified by Xiaozhe Hu 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016 & 08/28/2021
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
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
void print_itsolver_info(const INT  print_lvl,
                         const INT  stop_type,
                         const INT  iter,
                         const REAL rel_res,
                         const REAL abs_res,
                         const REAL factor)
{
  if ( print_lvl >= PRINT_SOME ) {

    // case iter > 0:  not the first iteration
    if ( iter > 0 ) {
      printf("%6lld | %13.6e   | %13.6e  | %10.4f\n",   (long long )iter, rel_res, abs_res, factor);
    }
    // case iter == 0: first iteration
    else {
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
      printf("%6lld | %13.6e   | %13.6e  |     -.-- \n",   (long long )iter, rel_res, abs_res);
    } // end if iter

  }

}

/***********************************************************************************************/
/*!
   * \fn void print_cputime (const char *message, const REAL cputime)
   *
   * \brief Print CPU walltime
   *
   * \param message   Pointer to the string to print out
   * \param cputime   Walltime since start to end
   *
   */
void print_cputime(const char *message,
                   const REAL cputime)
{
  printf("%s costs %.4f seconds.\n", message, cputime);
}

/***********************************************************************************************/
/*!
   * \fn void print_message (const INT print_lvl, const char *message)
   *
   * \brief Print out the message if necessary
   *
   * \param print_lvl   Level for output
   * \param message     Pointer to the error message
   *
   */
void print_message(const INT print_lvl,
                   const char *message)
{
  if ( print_lvl > PRINT_NONE ) printf("%s", message);
}

/***********************************************************************************************/
/*!
 * \fn void print_amg_complexity (AMG_data *mgl, const SHORT prtlvl)
 *
 * \brief Print grid and operator complexity of AMG method
 *
 * \param mgl      Multilevel structure for AMG
 * \param print_lvl   How much information to print
 *
 */
void print_amg_complexity(AMG_data *mgl,
                          const SHORT print_lvl)
{
  // local variables
  const SHORT   max_levels = mgl->num_levels;
  SHORT         level;
  REAL          grid_complexity=0.0, operator_complexity=0.0;
  REAL          AvgNNZ;

  if ( print_lvl >= PRINT_SOME ) {

    printf("-----------------------------------------------------------\n");
    printf("  Level   Num of rows   Num of nonzeros   Avg. NNZ / row   \n");
    printf("-----------------------------------------------------------\n");

    for ( level = 0; level < max_levels; ++level) {
        AvgNNZ = (REAL) mgl[level].A.nnz/mgl[level].A.row;
        printf("%5lld %13lld %17lld %14.2f\n", (long long )level, (long long )mgl[level].A.row, (long long )mgl[level].A.nnz, AvgNNZ);
        grid_complexity     += mgl[level].A.row;
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
/**
 * \fn void void print_amgcomplexity_bsr (const AMG_data_bsr *mgl,
 *                                       const SHORT prtlvl)
 *
 * \brief Print complexities of AMG method for BSR matrices
 *
 * \param mgl      Multilevel hierachy for AMG
 * \param prtlvl   How much information to print
 *
 */
void print_amgcomplexity_bsr(const AMG_data_bsr  *mgl,
                             const SHORT          prtlvl)
{
    const SHORT  max_levels = mgl->num_levels;
    SHORT        level;
    REAL         gridcom = 0.0, opcom = 0.0;

    if ( prtlvl >= PRINT_SOME ) {

        printf("-----------------------------------------------------------\n");
        printf("  Level   Num of rows   Num of nonzeros   Avg. NNZ / row   \n");
        printf("-----------------------------------------------------------\n");
        for ( level = 0; level < max_levels; ++level ) {
            const REAL AvgNNZ = (REAL) mgl[level].A.NNZ/mgl[level].A.ROW;
            printf("%5lld  %13lld  %17lld  %14.2f\n",
                   (long long )level,(long long )mgl[level].A.ROW, (long long )mgl[level].A.NNZ, AvgNNZ);
            gridcom += mgl[level].A.ROW;
            opcom   += mgl[level].A.NNZ;
        }
        printf("-----------------------------------------------------------\n");

        gridcom /= mgl[0].A.ROW;
        opcom   /= mgl[0].A.NNZ;
        printf("  Grid complexity = %.3f  |", gridcom);
        printf("  Operator complexity = %.3f\n", opcom);

        printf("-----------------------------------------------------------\n");

    }
}


/***********************************************************************************************/
/*!
 * \fn void check_error(const SHORT status, const char *func_name)
 *
 * \brief Check error status and print out error messages before quit
 *
 * \param status   Error status
 * \param fctname  piinter to the function name where this routine is called
 *
 */
void check_error(const SHORT status,
                 const char *func_name)
{
  if ( status >= SUCCESS ) return; // No error at all

  switch ( status ) {

  case ERROR_INPUT_PAR:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Wrong input arguments !!!\n\n", func_name);
      break;
  case ERROR_OPEN_FILE:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Cannot open file !!!\n\n", func_name);
      break;
  case ERROR_WRONG_FILE:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Wrong file format !!!\n\n", func_name);
      break;
  case ERROR_DIM:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- You have entered the Twilight Zone (dimension is not 1-3) !!! \n\n", func_name);
      break;
  case ERROR_FE_TYPE:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- The FE Space you want is not implemented !!! \n\n", func_name);
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
  case ERROR_AMG_SMOOTH_TYPE:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown AMG smoother type !!!\n\n", func_name);
      break;
  case ERROR_AMG_AGG_TYPE:
      printf("\n!!! ERROR HAZMATH DANGER: in function '%s' -- Unknown AMG aggregation type !!!\n\n", func_name);
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
/***********************************************************************************************/
/*!
 * \fn void error_lib(const SHORT status, const char *func_name, const char *libname)
 *
 * \brief Exit if a function is called from external library which has
 * not been compiled in the __HAZMATH__ library.
 *
 * \param status   Error status
 * \param f_name  piinter to the function name where this routine is called
 * \param library name External library name
 *
 */
void error_extlib(const SHORT status, const char *func_name, \
		  const char *lib_name)
{
  fprintf(stderr,"\n\n******** HAZMATH FATAL ERROR *******\n");
  fprintf(stderr,"   In function \"%s\":\n",func_name);
  fprintf(stderr,"   A call to a function from %s library occured,\n",lib_name);
  fprintf(stderr,"   but %s support was not compiled in the HAZmath library.\n",lib_name);
  fprintf(stderr,"   IF you have %s installed, THEN RECOMPILE \"libhazmath\"\n",lib_name);
  fprintf(stderr,"   with %s support.\n",lib_name);
  fprintf(stderr,"********\n\n");
  //  fprintf(stderr,"***by issuing the commands: make config %s=yes; make install\n\n",lib_name);
  exit(status);
}
