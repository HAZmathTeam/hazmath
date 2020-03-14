/*! \file src/utilities/format.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 02/17/16.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/28/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/28/2016
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
 *
 * \brief   Transform a block_dCSRmat matrix to a dCSRmat matrix
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \return  A    dCSRmat matrix if succeed, NULL if fail
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this functions! -- Xiaozhe Hu
 *
 */
dCSRmat bdcsr_2_dcsr (block_dCSRmat *Ab)
{
    // local variables
    INT m=0,n=0,nnz=0;
    const INT mb=Ab->brow, nb=Ab->bcol, n_blocks=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
    INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrowp1;
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

    // count number of nonzeros
    for (i=0;i<n_blocks;++i) {
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
                    irmrow=ir-row[i];irmrowp1=irmrow+1;
                    ilength=blockptrij->IA[irmrowp1]-blockptrij->IA[irmrow];

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
/*!
 * \fn dCSRmat bdcsr_subblk_2_dcsr (block_dCSRmat *Ab)
 *
 * \brief   Transform a block_dCSRmat matrix to a dCSRmat matrix
 *
 * \param   Ab   Pointer to a block_dCSRmat matrix
 *
 * \return  A    dCSRmat matrix if succeed, NULL if fail
 *
 * \note    Memory space for the dCSRmat matrix is allocated inside this functions! -- Xiaozhe Hu
 *
 */
dCSRmat bdcsr_subblk_2_dcsr (block_dCSRmat *Ab, INT brow_start, INT brow_end, INT bcol_start, INT bcol_end)
{
    INT i, j;
    INT ii, jj;
    block_dCSRmat Asub;
    bdcsr_alloc_minimal( brow_end-brow_start+1, bcol_end-bcol_start+1, &Asub);

    for(j = bcol_start; j<bcol_end+1; j++){
      jj = j-bcol_start;
      for(i = brow_start; i<brow_end+1; i++){
        ii = i-brow_start;
        Asub.blocks[ii*Asub.bcol + jj] = Ab->blocks[i*Ab->bcol + j];
      }
    }

    dCSRmat A = bdcsr_2_dcsr( &Asub );
//    bdcsr_free_minimal( &Asub );
    return(A);
}

/***********************************************************************************************/
/*!
 * \fn block_dCSRmat dcsr_2_bdcsr (dCSRmat *A, int bnum, int *bsize)
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
block_dCSRmat dcsr_2_bdcsr (dCSRmat *A,
                            int bnum,
                            int *bsize)
{
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
SHORT dcoo_2_dcsr (dCOOmat *A,
                   dCSRmat *B)
{
    // get size
    const INT m=A->row, n=A->col, nnz=A->nnz;

    // allocate
    dcsr_alloc(m,n,nnz,B);

    // local variable
    INT *ia = B->IA;
    INT *ja = B->JA;
    REAL *Bval = B->val;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    REAL *Aval = A->val;
    INT i, iind, jind;

    INT *ind = (INT *) calloc(m+1,sizeof(INT));

    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));

    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;

    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }

    // set column index and values
    for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        Bval[jind] = Aval[i];
        ind[iind] = ++jind;
    }

    if (ind) free(ind);

    return SUCCESS;
}

/***********************************************************************************************/
/*!
 * \fn coordinates* array_2_coord (REAL *xyz, INT ndof, INT dim)
 *
 * \brief Transform an array of coordinates from [x0,y0,z0,x1,y1,z1,...] to a coordinates structure
 *
 * \param xyz  Array with coordinate values
 * \param ndof Number of degrees of freedom
 * \param dim  Dimension
 *
 * \return cv Coordinate struct
 *
 */
coordinates* array_2_coord ( REAL* xyz, INT ndof, INT dim)
{
  // Local Variables
  INT dof;//,i,j;
  // Allocate
  coordinates* cv = allocatecoords(ndof, dim);
  // Fill
  switch (dim)
  {
  case 1:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
    }
    break;
  case 2:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
      cv->y[dof] = xyz[dof*dim+1];
    }
    break;
  case 3:
    for( dof=0; dof<ndof; dof++){
      cv->x[dof] = xyz[dof*dim+0];
      cv->y[dof] = xyz[dof*dim+1];
      cv->z[dof] = xyz[dof*dim+2];
    }
    break;
  default:
    check_error( ERROR_DIM, __FUNCTION__);
  }
  return cv;
}

/********************************  END  ********************************************************/
