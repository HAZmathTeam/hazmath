//
//  sparse.h
//  
//
//  Created by Hu, Xiaozhe on 1/9/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _sparse_h
#define _sparse_h


/**
 * \struct dCSRmat
 * \brief Sparse matrix of REAL type in CSR format
 *
 * CSR Format (IA,JA,A) in REAL
 *
 * \note The starting index of A is 0.
 */
typedef struct dCSRmat{
    
    //! row number of matrix A, m
    INT row;
    
    //! column of matrix A, n
    INT col;
    
    //! number of nonzero entries
    INT nnz;
    
    //! integer array of row pointers, the size is m+1
    INT *IA;
    
    //! integer array of column indexes, the size is nnz
    INT *JA;
    
    //! nonzero entries of A
    REAL *val;
    
} dCSRmat; /**< Sparse matrix of REAL type in CSR format */

/**
 * \struct iCSRmat
 * \brief Sparse matrix of INT type in CSR format
 *
 * CSR Format (IA,JA,A) in integer
 *
 * \note The starting index of A is 0.
 */
typedef struct iCSRmat{
    
    //! row number of matrix A, m
    INT row;
    
    //! column of matrix A, n
    INT col;
    
    //! number of nonzero entries
    INT nnz;
    
    //! integer array of row pointers, the size is m+1
    INT *IA;
    
    //! integer array of column indexes, the size is nnz
    INT *JA;
    
    //! nonzero entries of A
    INT *val;
    
} iCSRmat; /**< Sparse matrix of INT type in CSR format */


#endif
