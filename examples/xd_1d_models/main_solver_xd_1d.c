/*! \file examples/solver/solver_bdcsr.c
 *
 *  Created by Xiaozhe Hu on 04/16/22.
 *  Copyright 2022_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a matrix (in block_dCSRmat format) and a right hand side and solve it by certain linear solvers
 *
 * \note
 *
 */

/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
#include "hazmath.h"
#include "definitions_xd_1d.h"
#include "supporting_xd_1d.h"
#include "solver_xd_1d.h"
/***********************************************************************/
INT main(int argc, char* argv[])
{
  /*************** ACTION *************************************/
  //  char *dir_matrices=strdup("./input/1d_matrices_3d/");
  char *dir_matrices=strdup("./input/1d_matrices_2d/");
  char *finput_solver=strdup("./input/solver.input");
  solver_xd_1d(finput_solver,dir_matrices);
  free(finput_solver);
  free(dir_matrices);
  /* Set Solver Parameters */
  /* ivec_free(idofs); */
}	/* End of Program */
/*******************************************************************/
