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
 * \brief integer and floating point numbers
 */
#define SHORT            short      /**< short integer type */
#define INT              int        /**< regular integer type: int or long */
#define LONG             long       /**< long integer type */
#define LONGLONG         long long  /**< long integer type */
#define REAL             double     /**< float type */

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
