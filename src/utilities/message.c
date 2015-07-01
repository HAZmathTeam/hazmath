/*
 *  message.c
 *
 *  Created by James Adler and Xiaozhe Hu on 3/6/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/***********************************************************************************************/
void print_itsolver_info (const INT ptrlvl,
                   const INT stop_type,
                   const INT iter,
                   const REAL relres,
                   const REAL absres,
                   const REAL factor)
{
    /**
     * \fn void print_itsolver_info (const INT ptrlvl, const INT stop_type, const INT iter,
     *                        const REAL relres, const REAL absres, const REAL factor)
     *
     * \brief Print out iteration information for iterative solvers
     *
     * \param ptrlvl     Level for output
     * \param stop_type  Type of stopping criteria
     * \param iter       Number of iterations
     * \param relres     Relative residual of different kinds
     * \param absres     Absolute residual of different kinds
     * \param factor     Contraction factor
     *
     */
    
    if ( ptrlvl >= PRINT_SOME ) {
        
        if ( iter > 0 ) {
            printf("%6d | %13.6e   | %13.6e  | %10.4f\n", iter, relres, absres, factor);
        }
        else { // iter = 0: initial guess
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
            printf("%6d | %13.6e   | %13.6e  |     -.-- \n", iter, relres, absres);
        } // end if iter
        
    } // end if ptrlvl
}

/***********************************************************************************************/
void print_cputime (const char *message,
                    const REAL cputime)
{
    
    /**
     * \fn void void print_cputime (const char *message, const REAL cputime)
     *
     * \brief Print CPU walltime
     *
     * \param message   Some string to print out
     * \param cputime   Walltime since start to end
     *
     */
    
    printf("%s costs %.4f seconds.\n", message, cputime);
}

/***********************************************************************************************/
void print_message (const INT ptrlvl,
                    const char *message)
{
    /**
     * \fn void print_message (const INT ptrlvl, const char *message)
     *
     * \brief Print output information if necessary
     *
     * \param ptrlvl   Level for output
     * \param message  Error message to print
     *
     */
    
    if ( ptrlvl > PRINT_NONE ) printf("%s", message);
}

