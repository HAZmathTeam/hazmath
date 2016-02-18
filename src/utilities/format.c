/*
 *  format.c
 *
 *  Created by James Adler and Xiaozhe Hu on 02/17/16.
 *  Copyright 2016__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/***********************************************************************************************/
dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
{
    /**
     * \fn dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
     *
     * \brief Form A in dCSRmat format using blocks given by Ab which is in block_dCSRmat format
     *
     * \param Ab   Pointer to block_dCSRmat matrix
     *
     * \return     dCSRmat matrix if succeed, NULL if fail
     *
     * \author Xiaozhe Hu
     * \date   02/17/2016
     */
    
    // local variables
    INT m=0,n=0,nnz=0;
    const INT mb=Ab->brow, nb=Ab->bcol, nbl=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
    INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrow1;
    INT *row, *col;
    
    row = (INT *)calloc(mb+1,sizeof(INT));
    col = (INT *)calloc(nb+1,sizeof(INT));
    
    // count the size of A
    row[0]=0; col[0]=0;
    
    for (i=0;i<mb;++i) {
        
        for (j=0; j<nb; ++j){
        
            if (blockptr[i*nb+j]) {
                m+=blockptr[i*nb]->row;
                row[i+1]=m;
                break;
            }
            
        }
    }
    
    for (i=0;i<nb;++i) {
        
        for (j=0;j<mb;++j){
        
            if (blockptr[j*mb+i]) {
                n+=blockptr[j*mb+i]->col;
                col[i+1]=n;
                break;
            }
            
        }
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
                
                ij=i*nb+j; blockptrij=blockptr[ij];
                
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
    /**
     * \fn block_dCSRmat dcsr_2_bdcsr (dCSRmat *A)
     *
     * \brief Form Ab in block_dCSRmat format using A which is in dCSRmat format
     *
     * \param A   Pointer to dCSRmat matrix
     *
     * \return     block_dCSRmat matrix if succeed, NULL if fail
     *
     * \author Xiaozhe Hu
     * \date   02/17/2016
     */
    
    // local variable
    int i, j;
    
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
            
            dcsr_getblk(A, idx_row+idx_start[i], idx_col+idx_start[j], bsize[i], bsize[j], Ab.blocks[i*bnum+j]);
            
        }
        
    }
    
    
    // free memory
    free(idx_row);
    free(idx_col);
    free(idx_start);
    
    return Ab;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
