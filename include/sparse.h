//
//  sparse.h
//
//
//  Created by Adler, James and Hu, Xiaozhe on 1/9/15.
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

    //! m: number of rows in the matrix A
    INT row;

    //! n: number of columns in the matrix A
    INT col;

    //! nnz: number of the nonzero entries in A
    INT nnz;

    //! integer array of row pointers, the size is m+1
    INT *IA;

    //! integer array of column indexes, the size is nnz
    INT *JA;

    //! nonzero entries of A
    INT *val;

} iCSRmat; /**< Sparse matrix of INT type in CSR format */

/**
 * \struct dCOOmat
 * \brief Sparse matrix of REAL type in COO (or IJ) format
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
 * \note Change I to rowind, J to colind. To avoid with complex.h confliction on I.
 */
typedef struct dCOOmat{

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row indices, the size is nnz
    INT *rowind;

    //! integer array of column indices, the size is nnz
    INT *colind;

    //! nonzero entries of A
    REAL *val;

} dCOOmat; /**< Sparse matrix of REAL type in COO format */

/**
 * \struct dCOOmat
 * \brief Sparse matrix of REAL type in COO (or IJ) format
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
 * \note Change I to rowind, J to colind. To avoid with complex.h confliction on I.
 */
typedef struct iCOOmat{

    //! row number of matrix A, m
    INT row;

    //! column of matrix A, n
    INT col;

    //! number of nonzero entries
    INT nnz;

    //! integer array of row indices, the size is nnz
    INT *rowind;

    //! integer array of column indices, the size is nnz
    INT *colind;

    //! nonzero entries of A
    INT *val;

} iCOOmat; /**< Sparse matrix of INT type in COO format */


/**
 * \struct block_dCSRmat
 * \brief Block REAL CSR matrix format
 *
 * \note The starting index of A is 0.
 */
typedef struct block_dCSRmat {

    //! row number of blocks in A, m
    INT brow;

    //! column number of blocks A, n
    INT bcol;

    //! blocks of dCSRmat, point to blocks[brow][bcol]
    dCSRmat **blocks;

} block_dCSRmat; /**< Matrix of REAL type in Block CSR format */

/**
 * \struct block_iCSRmat
 * \brief Block INT CSR matrix format
 *
 * \note The starting index of A is 0.
 */
typedef struct block_iCSRmat {

    //! row number of blocks in A, m
    INT brow;

    //! column number of blocks A, n
    INT bcol;

    //! blocks of iCSRmat, point to blocks[brow][bcol]
    iCSRmat **blocks;

} block_iCSRmat; /**< Matrix of INT type in Block CSR format */

#endif
