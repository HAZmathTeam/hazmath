//
//  macro.h
//
//
//  Created by Hu, Xiaozhe on 1/10/15.
//  Modified 20190327 (ltz)
//  Modified 20190728 (ltz)
//

#ifndef _macro_h
#define _macro_h

/* pi=3.141592653589793238512809e+00; e= 2.718281828459045235428168e+00
 in long double, below is rounded off a bit
*/
#define PI 3.14159265358979323851281e+00
#define EXP_OF_1  2.71828182845904523542817e+00
#define MAX_CHARS_INPUT_GRID_FILE   ((1<<15) - 1)

/**
 * \brief Some global constants
 */
#define BIGREAL          1e+20 /**< A large real number */
#define SMALLREAL        1e-20 /**< A small real number */
#define SMALLREAL2       1e-40  /**< An extremely small real number */
#define MAX_REFINE_LVL   20    /**< Maximal refinement level */
#define MAX_AMG_LVL      20    /**< Maximal AMG coarsening level */
#define MIN_CDOF         1    /**< Minimal number of coarsest variables */
#define MIN_CRATE        0.9   /**< Minimal coarsening ratio */
#define MAX_CRATE        100.0  /**< Maximal coarsening ratio */
#define STAG_RATIO       1e-4  /**< Stagnation tolerance = tol*STAGRATIO */
#define MAX_STAG         20    /**< Maximal number of stagnation times */
#define MAX_RESTART      20    /**< Maximal number of restarting for BiCGStab */

/**
 * \brief Definition of return status and error messages
 */
#define SUCCESS                 0  /**< return from function successfully */
//---------------------------------------------------------------------------------
#define ERROR_DIM              -1 /**< wrong dimension */
#define ERROR_FE_TYPE          -2 /**< wrong type of FEM */
#define ERROR_QUAD_TYPE        -3 /**< unknown quadrature type */
#define ERROR_QUAD_DIM         -4  /**< unsupported quadrature dim */
#define ERROR_MAT_DOF          -5 /**< stiffness matrix size and dof not consistent */
#define ERROR_TS_TYPE          -6 /**< unknown time-stepping scheme */
//---------------------------------------------------------------------------------
#define ERROR_OPEN_FILE       -10  /**< fail to open a file */
#define ERROR_WRONG_FILE      -11  /**< input contains wrong format */
#define ERROR_INPUT_PAR       -13  /**< wrong input argument */
#define ERROR_MAT_SIZE        -15  /**< wrong problem size */
#define ERROR_BLKMAT_ZERO     -16  /**< block matrix is singular (at least one row or one column is zero) */
#define ERROR_NUM_BLOCKS      -18  /**< wrong number of blocks */
#define ERROR_MISC            -19  /**< other error */
//---------------------------------------------------------------------------------
#define ERROR_ALLOC_MEM       -20  /**< fail to allocate memory */
#define ERROR_DATA_STRUCTURE  -21  /**< problem with data structures */
#define ERROR_DATA_ZERODIAG   -22  /**< matrix has zero diagonal entries */
#define ERROR_DUMMY_VAR       -23  /**< unexpected input data */
//---------------------------------------------------------------------------------
#define ERROR_AMG_SMOOTH_TYPE -31  /**< unknown smoother type */
#define ERROR_AMG_AGG_TYPE    -32  /**< unknown aggregation type */
#define ERROR_AMG_COARSEING   -33  /**< coarsening step failed to complete */
//---------------------------------------------------------------------------------
#define ERROR_SOLVER_TYPE     -40  /**< unknown solver type */
#define ERROR_SOLVER_PRECTYPE -41  /**< unknown precond type */
#define ERROR_SOLVER_STAG     -42  /**< solver stagnates */
#define ERROR_SOLVER_SOLSTAG  -43  /**< solver's solution is too small */
#define ERROR_SOLVER_TOLSMALL -44  /**< solver's tolerance is too small */
#define ERROR_SOLVER_MISC     -47  /**< misc solver error during run time */
#define ERROR_SOLVER_MAXIT    -48  /**< maximal iteration number exceeded */
#define ERROR_SOLVER_EXIT     -49  /**< solver does not quit successfully */
//---------------------------------------------------------------------------------
#define ERROR_LIC_TYPE        -80  /**< wrong license type */
//---------------------------------------------------------------------------------
#define ERROR_UNKNOWN         -99  /**< an unknown error type */

/**
 * \brief Definition of type "logical"
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
#ifndef REAL16
#define REAL16 long double
#endif


/**
 * \brief Definition of solver types for nonlinear methods
 */
#define NONLINEAR_NEWTON          0  /**< Newton's Method */
#define NONLINEAR_PICARD          1  /**< Picard Iterations */

/**
 * \brief Definition of solver types for linear iterative methods
 */
#define SOLVER_DEFAULT          0  /**< Use default solver */
//---------------------------------------------------------------------------------
#define SOLVER_CG               1  /**< Conjugate Gradient */
#define SOLVER_MinRes           2  /**< Minimal Residual */
#define SOLVER_VGMRES           3  /**< Variable Restarting GMRES */
#define SOLVER_VFGMRES          4  /**< Variable Restarting Flexible GMRES */
#define SOLVER_GCG              5  /**< Generalized Conjugate Gradient */
#define SOLVER_GCR              6  /**< Generalized Conjugate Residual */
//---------------------------------------------------------------------------------
#define SOLVER_AMG             21  /**< AMG as an iterative solver */
//---------------------------------------------------------------------------------
#define SOLVER_UMFPACK         32  /**< UMFPack Direct Solver */

/**
 * \brief Definition of iterative solver stopping criteria types
 */
#define STOP_REL_RES            1  /**< relative residual ||r||/||b|| */
#define STOP_REL_PRECRES        2  /**< relative B-residual ||r||_B/||b||_B */
#define STOP_MOD_REL_RES        3  /**< modified relative residual ||r||/||x|| */

/**
 * \brief Definition of preconditioner type for iterative methods
 */
#define PREC_NULL               0  /**< with no precond */
#define PREC_DIAG               1  /**< with diagonal precond */
#define PREC_AMG                2  /**< with AMG precond */
#define PREC_HX_CURL_A          6  /**< with additive HX preconditioner for H(curl) problem */
#define PREC_HX_CURL_M          7  /**< with multiplicative HX preconditioner for H(curl) problem */
#define PREC_HX_DIV_A           8  /**< with additive HX preconditioner for H(div) problem */
#define PREC_HX_DIV_M           9  /**< with multiplicative HX preconditioner for H(div) problem */

/**
 * \brief Definition of AMG types
 */
#define UA_AMG                  1  /**< unsmoothed aggregation AMG */
#define SA_AMG                  2  /**< smoothed aggregation AMG */
#define C_AMG                   3  /**< classical AMG */
#define MUA_AMG                 10 /**< Metric Unsmoothed aggregation AMG for interface problem */
#define MSA_AMG                 11 /**< Metric Smoothed aggregation AMG for interface problem */


/**
 * \brief Definition of aggregation types
 */
#define VMB                     1  /**< VMB aggregation */
#define MIS                     2  /**< MIS aggregation */
#define MWM                     3  /**< Maximal Weighted Matching aggregation */
#define HEC                     4  /**< Hearvy Edge Coarsening aggregation */

/**
 * \brief Definition of cycle types
 */
#define V_CYCLE                 1  /**< V-cycle */
#define W_CYCLE                 2  /**< W-cycle */
#define AMLI_CYCLE              3  /**< AMLI-cycle */
#define NL_AMLI_CYCLE           4  /**< Nonlinear AMLI-cycle */
#define ADD_CYCLE               5  /**< additive cycle */

/**
 * \brief Type of Schwarz smoother
 */
#define SCHWARZ_FORWARD         1  /**< Forward ordering */
#define SCHWARZ_BACKWARD        2  /**< Backward ordering */
#define SCHWARZ_SYMMETRIC       3  /**< Symmetric smoother */

/**
 * \brief Definition of standard smoother types
 */
#define SMOOTHER_JACOBI         1  /**< Jacobi smoother */
#define SMOOTHER_GS             2  /**< Gauss-Seidel smoother */
#define SMOOTHER_SGS            3  /**< Symmetric Gauss-Seidel smoother */
#define SMOOTHER_CG             4  /**< CG as a smoother */
#define SMOOTHER_SOR            5  /**< SOR smoother */
#define SMOOTHER_SSOR           6  /**< SSOR smoother */
#define SMOOTHER_GSOR           7  /**< GS + SOR smoother */
#define SMOOTHER_SGSOR          8  /**< SGS + SSOR smoother */
#define SMOOTHER_POLY           9  /**< Polynomial smoother */
#define SMOOTHER_L1DIAG        10  /**< L1 norm diagonal scaling smoother */
#define SMOOTHER_FJACOBI       11  /**< Fractional Jacobi smoother */
#define SMOOTHER_FGS           12  /**< Fractional Gauss-Seidel smoother */
#define SMOOTHER_FSGS          13  /**< Fractional Symmetric Gauss-Seidel smoother */
#define SMOOTHER_USERDEF       20  /**< User defined smoother (NB! requires fptr to smoother mxv */
#define SMOOTHER_JACOBI_GS     31  /**< Jacobi GS smoother for block dCSRmat (Block Jacobi + GS for each block)  */
#define SMOOTHER_JACOBI_SGS    32  /**< Jacobi SGS smoother for block dCSRmat (Block Jacobi + SGS for each block)  */
#define SMOOTHER_METRIC_ADD    40  /**< Additive metric smoother for block dCSRmat  */
#define SMOOTHER_METRIC_MUL    41  /**< Multiplicative metric smoother for block dCSRmat */

/**
 * \brief Type of vertices (DOFs) for coarsening
 */
#define UNPT                   -1  /**< Undetermined points */

/**
 * \brief Type of ordering for smoothers
 */
#define USERDEFINED             0  /**< User defined order */
#define ASCEND                 12  /**< Ascending order */
#define DESCEND                21  /**< Descending order */


/**
 * \brief Definition of max, min, abs
 */
#define MAX(a,b) (((a)>(b))?(a):(b)) /**< max of a and b */
#define MIN(a,b) (((a)<(b))?(a):(b)) /**< min of a and b */
#define ABS(a) (((a)>=0.0)?(a):-(a)) /**< absolute value of a (abs or fabs or fabsl shoould be preferred) */

/**
 * \brief Definition of >, >=, <, <=, and isnan
 */
#define GT(a,b) (((a)>(b))?(TRUE):(FALSE))   /**< is a > b? */
#define GE(a,b) (((a)>=(b))?(TRUE):(FALSE))  /**< is a >= b? */
#define LS(a,b) (((a)<(b))?(TRUE):(FALSE))   /**< is a < b? */
#define LE(a,b) (((a)<=(b))?(TRUE):(FALSE))  /**< is a <= b? */
#define ISNAN(a) (((a)!=(a))?(TRUE):(FALSE)) /**< is a == NAN? */

/**
 * \brief Definition of Norms
 */
#define ONENORM                 1  /**< 1 norm */
#define TWONORM                 2  /**< 2 norm */

#endif
