//
//  macro.h
//  
//
//  Created by Hu, Xiaozhe on 1/10/15.
//
//

#ifndef _macro_h
#define _macro_h

/**
 * \brief Some global constants
 */
#define BIGREAL          1e+20 /**< A large real number */
#define SMALLREAL        1e-20 /**< A small real number */
#define MAX_REFINE_LVL   20    /**< Maximal refinement level */
#define MAX_AMG_LVL      20    /**< Maximal AMG coarsening level */
#define MIN_CDOF         20    /**< Minimal number of coarsest variables */
#define MIN_CRATE        0.9   /**< Minimal coarsening ratio */
#define MAX_CRATE        20.0  /**< Maximal coarsening ratio */
#define STAG_RATIO       1e-4  /**< Stagnation tolerance = tol*STAGRATIO */
#define MAX_STAG         20    /**< Maximal number of stagnation times */
#define MAX_RESTART      20    /**< Maximal number of restarting for BiCGStab */
#define OPENMP_HOLDS     2000  /**< Switch to sequence version when size is small */


/**
 * \brief Definition of return status and error messages
 */
#define SUCCESS            0  /**< return from function successfully */
//---------------------------------------------------------------------------------
#define ERROR_OPEN_FILE       -10  /**< fail to open a file */
#define ERROR_WRONG_FILE      -11  /**< input contains wrong format */
#define ERROR_INPUT_PAR       -13  /**< wrong input argument */
#define ERROR_REGRESS         -14  /**< regression test fail */
#define ERROR_MAT_SIZE        -15  /**< wrong problem size */
#define ERROR_NUM_BLOCKS      -18  /**< wrong number of blocks */
#define ERROR_MISC            -19  /**< other error */
//---------------------------------------------------------------------------------
#define ERROR_ALLOC_MEM       -20  /**< fail to allocate memory */
#define ERROR_DATA_STRUCTURE  -21  /**< problem with data structures */
#define ERROR_DATA_ZERODIAG   -22  /**< matrix has zero diagonal entries */
#define ERROR_DUMMY_VAR       -23  /**< unexpected input data */
//---------------------------------------------------------------------------------
#define ERROR_AMG_INTERP_TYPE -30  /**< unknown interpolation type */
#define ERROR_AMG_SMOOTH_TYPE -31  /**< unknown smoother type */
#define ERROR_AMG_COARSE_TYPE -32  /**< unknown coarsening type */
#define ERROR_AMG_COARSEING   -33  /**< coarsening step failed to complete */
//---------------------------------------------------------------------------------
#define ERROR_SOLVER_TYPE     -40  /**< unknown solver type */
#define ERROR_SOLVER_PRECTYPE -41  /**< unknown precond type */
#define ERROR_SOLVER_STAG     -42  /**< solver stagnates */
#define ERROR_SOLVER_SOLSTAG  -43  /**< solver's solution is too small */
#define ERROR_SOLVER_TOLSMALL -44  /**< solver's tolerance is too small */
#define ERROR_SOLVER_ILUSETUP -45  /**< ILU setup error */
#define ERROR_SOLVER_MISC     -46  /**< misc solver error during run time */
#define ERROR_SOLVER_MAXIT    -48  /**< maximal iteration number exceeded */
#define ERROR_SOLVER_EXIT     -49  /**< solver does not quit successfully */
//---------------------------------------------------------------------------------
#define ERROR_QUAD_TYPE       -60  /**< unknown quadrature type */
#define ERROR_QUAD_DIM        -61  /**< unsupported quadrature dim */
//---------------------------------------------------------------------------------
#define ERROR_LIC_TYPE        -80  /**< wrong license type */
//---------------------------------------------------------------------------------
#define ERROR_UNKNOWN         -99  /**< an unknown error type */

/**
 * \brief Definition of logic type
 */
#define TRUE                    1  /**< logic TRUE */
#define FALSE                   0  /**< logic FALSE */

/**
 * \brief Definition of switch
 */
#define ON                      1  /**< turn on certain parameter */
#define OFF                     0  /**< turn off certain parameter */

/**
 * \brief Print level for all subroutines -- not including DEBUG output
 */
#define PRINT_NONE              0  /**< silent: no printout at all */
#define PRINT_MIN               1  /**< quiet: min info, error, important warnings */
#define PRINT_SOME              2  /**< some: more info, less important warnings */
#define PRINT_MORE              4  /**< more: print some useful debug information */
#define PRINT_MOST              8  /**< most: maximal printouts, no files */
#define PRINT_ALL              10  /**< everything: all printouts, including files */

/**
 * \brief integer and floating point numbers
 */
#define SHORT            short      /**< short integer type */
#define INT              int        /**< regular integer type: int or long */
#define LONG             long       /**< long integer type */
#define LONGLONG         long long  /**< long integer type */
#define REAL             double     /**< float type */

/**
 * \brief Definition of iterative solver stopping criteria types
 */
#define STOP_REL_RES            1  /**< relative residual ||r||/||b|| */
#define STOP_REL_PRECRES        2  /**< relative B-residual ||r||_B/||b||_B */
#define STOP_MOD_REL_RES        3  /**< modified relative residual ||r||/||x|| */

/**
 * \brief Definition of max, min, abs
 */
#define MAX(a,b) (((a)>(b))?(a):(b)) /**< bigger one in a and b */
#define MIN(a,b) (((a)<(b))?(a):(b)) /**< smaller one in a and b */
#define ABS(a) (((a)>=0.0)?(a):-(a)) /**< absolute value of a */

/**
 * \brief Definition of >, >=, <, <=, and isnan
 */
#define GT(a,b) (((a)>(b))?(TRUE):(FALSE))   /**< is a > b? */
#define GE(a,b) (((a)>=(b))?(TRUE):(FALSE))  /**< is a >= b? */
#define LS(a,b) (((a)<(b))?(TRUE):(FALSE))   /**< is a < b? */
#define LE(a,b) (((a)<=(b))?(TRUE):(FALSE))  /**< is a <= b? */
#define ISNAN(a) (((a)!=(a))?(TRUE):(FALSE)) /**< is a == NAN? */

#endif
