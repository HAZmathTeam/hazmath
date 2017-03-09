/*! \file src/utilities/format.c
 *
 *  Created by James Adler and Xiaozhe Hu on 02/17/16.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/28/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/28/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
{
    /*!
     * \fn dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
     *
     * \brief Transform a block_dCSRmat matrix to a dCSRmat matrix
     *
     * \param Ab   Pointer to a block_dCSRmat matrix
     *
     * \return A   dCSRmat matrix if succeed, NULL if fail
     *
     * \note Memory space for the dCSRmat matrix is allocated inside this functions! -- Xiaozhe Hu
     *
     */
    
    // local variables
    INT m=0,n=0,nnz=0;
    const INT mb=Ab->brow, nb=Ab->bcol, nbl=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
    INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrow1;
    INT *row, *col;

    // flag for errors
    SHORT status = SUCCESS;
    
    row = (INT *)calloc(mb+1,sizeof(INT));
    col = (INT *)calloc(nb+1,sizeof(INT));
    
    // get the size of A
    row[0]=0; col[0]=0;

    // count number of rows
    for (i=0;i<mb;++i) {
        
        status = ERROR_BLKMAT_ZERO;

        for (j=0; j<nb; ++j){
            if (blockptr[i*nb+j]) {
                m+=blockptr[i*nb+j]->row;
                row[i+1]=m;
                status = SUCCESS;
                break;
            }
        }

        // check error
        if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);

    }
    
    // count number of columns
    for (i=0;i<nb;++i) {
        
        status = ERROR_BLKMAT_ZERO;

        for (j=0;j<mb;++j){
            if (blockptr[j*mb+i]) {
                n+=blockptr[j*mb+i]->col;
                col[i+1]=n;
                status = SUCCESS;
                break;
            }
            
        }

        // check error
        if (status < SUCCESS) check_error(ERROR_BLKMAT_ZERO, __FUNCTION__);

    }
    
    for (i=0;i<nbl;++i) {
        if (blockptr[i]) {
            nnz+=blockptr[i]->nnz;
        }
    }
    
    // memory space allocation
    A = dcsr_create(m,n,nnz);
    
    // set dCSRmat for A
    A.IA[0]=0;
    for (i=0;i<mb;++i) {
        
        for (ir=row[i];ir<row[i+1];ir++) {
            
            for (length=j=0;j<nb;++j) {
                
                ij=i*nb+j;
                blockptrij=blockptr[ij];
                
                if (blockptrij && blockptrij->nnz>0) {

                    start=A.IA[ir]+length;
                    irmrow=ir-row[i]; irmrow1=irmrow+1;
                    ilength=blockptrij->IA[irmrow1]-blockptrij->IA[irmrow];

                    if (ilength>0) {
                        memcpy(&(A.val[start]),&(blockptrij->val[blockptrij->IA[irmrow]]),ilength*sizeof(REAL));
                        memcpy(&(A.JA[start]), &(blockptrij->JA[blockptrij->IA[irmrow]]), ilength*sizeof(INT));
                        for (i1=0;i1<ilength;i1++) A.JA[start+i1]+=col[j];
                        length+=ilength;
                    }

                 }
                
            } // end for j
            
            A.IA[ir+1]=A.IA[ir]+length;

        } // end for ir
        
    } // end for i
    
    free(row);
    free(col);
    
    return(A);
}

/***********************************************************************************************/
block_dCSRmat dcsr_2_bdcsr (dCSRmat *A, int bnum, int *bsize)
{
    /*!
     * \fn block_dCSRmat dcsr_2_bdcsr (dCSRmat *A)
     *
     * \brief
     *
     * \param A         Pointer to dCSRmat matrix
     * \param bnum      Block size of the block_dCSRmat matrix
     * \param bsize     Pointer to the size of each diagonal block
     *
     * \return Ab       block_dCSRmat matrix if succeed, NULL if fail
     *
     * \note Assume dCSRmat A has block form already and just need to explicitly form the blocks without reordering! -- Xiaozhe Hu
     * \note SUM(bsize) = A->row = A->col, i.e., size has to be consistent!!  -- Xiaozhe Hu
     *
     */
    
    // local variable
    int i, j;
    SHORT status = SUCCESS;
    
    // allocate block dCSRmat
    block_dCSRmat Ab;
    Ab.brow = bnum; Ab.bcol = bnum;
    
    Ab.blocks = (dCSRmat **)calloc(bnum*bnum, sizeof(dCSRmat));
    for (i=0; i<bnum*bnum; i++) Ab.blocks[i] = (dCSRmat *) calloc(1, sizeof(dCSRmat));
    
    // allocate
    int *idx_row = (int *)calloc(A->row, sizeof(INT));
    int *idx_col = (int *)calloc(A->col, sizeof(INT));

    for (i=0; i<A->row; i++) idx_row[i] = i;  // generate row index
    for (i=0; i<A->col; i++) idx_col[i] = i;  // generate col index
    
    int *idx_start = (int *)calloc(bnum, sizeof(INT));
    idx_start[0] = 0;
    for (i=0; i<bnum-1; i++) idx_start[i+1] = idx_start[i] + bsize[i];
    
    for (i=0; i<bnum; i++){
        
        for (j=0; j<bnum; j++) {   

            status = dcsr_getblk(A, idx_row+idx_start[i], idx_col+idx_start[j], bsize[i], bsize[j], Ab.blocks[i*bnum+j]);
            if (status < SUCCESS) check_error(status, __FUNCTION__);

        }
        
    }

    
    // free memory
    free(idx_row);
    free(idx_col);
    free(idx_start);
    
    return Ab;
}

/***********************************************************************************************/
SHORT dcoo_2_dcsr (dCOOmat *A,
                   dCSRmat *B)
{
    /*!
     * \fn SHORT dcoo_2_dcsr (dCOOmat *A, dCSRmat *B)
     *
     * \brief Transform a dCOOmat matrix to a dCSRmat format.
     *
     * \param A   Pointer to dCOOmat matrix
     * \param B   Pointer to dCSRmat matrix
     *
     * \return    SUCCESS if successed; otherwise, error information.
     *
     */
    
    // get size
    const INT m=A->row, n=A->col, nnz=A->nnz;

    // allocate
    dcsr_alloc(m,n,nnz,B);
    
    // local variable
    INT * ia = B->IA, i;
    INT   iind, jind;
    
    INT *ind = (INT *) calloc(m+1,sizeof(INT));
    
    //for (i=0; i<=m; ++i) ind[i]=0; // initialize
    memset(ind, 0, sizeof(INT)*(m+1));
    
    for (i=0; i<nnz; ++i) ind[A->rowind[i]+1]++; // count nnz in each row
    
    ia[0] = 0; // first index starting from zero
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i]; // set row pointer
        ind[i] = ia[i];
    }
    
    // loop over nnz and set column index and value
    for (i=0; i<nnz; ++i) {
        iind = A->rowind[i]; jind = ind[iind];
        B->JA [jind] = A->colind[i];
        B->val[jind] = A->val[i];
        ind[iind] = ++jind;
    }
    
    if (ind) free(ind);
    
    return SUCCESS;
}

/********************************  END  ********************************************************/

