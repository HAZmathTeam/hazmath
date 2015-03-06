//
//  sparse.c
//  
//
//  Created by Hu, Xiaozhe on 1/9/15.
//
//

#include <math.h>
#include <time.h>

#include "macro.h"
#include "sparse.h"
#include "vec.h"

/***********************************************************************************************/
dCSRmat dcsr_create (const INT m,
                     const INT n,
                     const INT nnz)
{
    /**
     * \fn dCSRmat dcsr_create (const INT m, const INT n, const INT nnz)
     *
     * \brief Create CSR sparse matrix data memory space
     *
     * \param m    Number of rows
     * \param n    Number of columns
     * \param nnz  Number of nonzeros
     *
     * \return A   the new dCSRmat matrix
     *
     */
    
    dCSRmat A;
    
    if ( m > 0 ) {
        A.IA = (INT *)calloc(m+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( n > 0 ) {
        A.JA = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.JA = NULL;
    }
    
    if ( nnz > 0 ) {
        A.val = (REAL *)calloc(nnz, sizeof(REAL));
    }
    else {
        A.val = NULL;
    }
    
    A.row = m; A.col = n; A.nnz = nnz;
    
    return A;
}

/***********************************************************************************************/
iCSRmat icsr_create (const INT m,
                     const INT n,
                     const INT nnz)
{
    
    /**
     * \fn iCSRmat icsr_create (const INT m, const INT n, const INT nnz)
     *
     * \brief Create CSR sparse matrix data memory space
     *
     * \param m    Number of rows
     * \param n    Number of columns
     * \param nnz  Number of nonzeros
     *
     * \return A   the new iCSRmat matrix
     *
     */
    
    iCSRmat A;
    
    if ( m > 0 ) {
        A.IA = (INT *)calloc(m+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( n > 0 ) {
        A.JA = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.JA = NULL;
    }
    
    if ( nnz > 0 ) {
        A.val = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.val = NULL;
    }
    
    A.row = m; A.col = n; A.nnz = nnz;
    
    return A;
}

/***********************************************************************************************/
void dcsr_free (dCSRmat *A)
{
    /**
     * \fn void dcsr_free (dCSRmat *A)
     *
     * \brief Free CSR sparse matrix data memory space
     *
     * \param A   Pointer to the dCSRmat matrix
     *
     */
    
    if ( A == NULL ) return;
    
    if (A->IA) {
       free(A->IA);
        A->IA  = NULL;
    }
    
    if (A->JA) {
        free(A->JA);
        A->JA  = NULL;
    }
    
    if (A->val) {
        free(A->val);
        A->val = NULL;
    }
}

/***********************************************************************************************/
void icsr_free (iCSRmat *A)
{
    
    /**
     * \fn void icsr_free (iCSRmat *A)
     *
     * \brief Free CSR sparse matrix data memory space
     *
     * \param A   Pointer to the iCSRmat matrix
     *
     */
    
    if ( A == NULL ) return;
    
    if (A->IA) {
        free(A->IA);
        A->IA  = NULL;
    }
    
    
    if (A->JA) {
        free(A->JA);
        A->JA  = NULL;
    }
    
    if (A->val) {
        free(A->val);
        A->val = NULL;
    }
}

/***********************************************************************************************/
INT dcsr_trans (dCSRmat *A,
                dCSRmat *AT)
{
    
    /**
     * \fn void dcsr_trans (dCSRmat *A, dCSRmat *AT)
     *
     * \brief Find transpose of dCSRmat matrix A (index of the array starts with 0!!)
     *
     * \param A   Pointer to the dCSRmat matrix
     * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
     *
     */
    
    const INT n=A->row, m=A->col, nnz=A->nnz;
    
    // Local variables
    INT i,j,k,p;
    
    AT->row=m;
    AT->col=n;
    AT->nnz=nnz;
    
    AT->IA=(INT*)calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)calloc(nnz,sizeof(INT));
    
    if (A->val) {
        AT->val=(REAL*)calloc(nnz,sizeof(REAL));
        
    }
    else { AT->val=NULL; }
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
    
    // fasp_iarray_set(m+1, AT->IA, 0);
    memset(AT->IA, 0, sizeof(INT)*(m+1));
    
    for (j=0;j<nnz;++j) {
        i=A->JA[j]; // column Number of A = row Number of A'
        if (i<m-1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val) {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->val[k]=A->val[p];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    }
    else {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for (p=ibegin;p<iend1;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end of i
    } // end if
    
    return 0;
}

/***********************************************************************************************/
void dcsr_trans_1 (dCSRmat *A,
                   dCSRmat *AT)
{
 
    /**
     * \fn void dcsr_trans_1 (dCSRmat *A, dCSRmat *AT)
     *
     * \brief Find transpose of dCSRmat matrix A (index of the array starts with 1!!)
     *
     * \param A   Pointer to the dCSRmat matrix
     * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
     *
     */
    
    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
    
    dcsr_trans(A, AT);
    
    INT *AT_IA = AT->IA;
    INT *AT_JA = AT->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<AT->row+1; i++) AT_IA[i] = AT_IA[i]+1;
    for (i=0; i<A->nnz; i++) {
        A_JA[i] = A_JA[i]+1;
        AT_JA[i] = AT_JA[i]+1;
    }
    
}

/***********************************************************************************************/
void icsr_trans (iCSRmat *A,
                 iCSRmat *AT)
{

    /**
     * \fn void icsr_trans (iCSRmat *A, iCSRmat *AT)
     *
     * \brief Find transpose of iCSRmat matrix A (index of the array starts with 0!!)
     *
     * \param A   Pointer to the iCSRmat matrix
     * \param AT  Pointer to the transpose of iCSRmat matrix A (output)
     *
     */
    
    const INT n=A->row, m=A->col, nnz=A->nnz, m1=m-1;
    
    // Local variables
    INT i,j,k,p;
    INT ibegin, iend;
    
#if DEBUG_MODE
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",m,n,nnz);
#endif
    
    AT->row=m; AT->col=n; AT->nnz=nnz;
    
    AT->IA=(INT*)calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)calloc(nnz,sizeof(INT));
    
    if (A->val) {
        AT->val=(INT*)calloc(nnz,sizeof(INT));
    }
    else {
        AT->val=NULL;
    }
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
    memset(AT->IA, 0, sizeof(INT)*(m+1));
    
    for (j=0;j<nnz;++j) {
        i=A->JA[j]; // column Number of A = row Number of A'
        if (i<m1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val != NULL) {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->val[k]=A->val[p];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    }
    else {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    } // end if
}

/***********************************************************************************************/
void icsr_trans_1 (iCSRmat *A,
                 iCSRmat *AT)
{
    /**
     * \fn void icsr_trans_1 (iCSRmat *A, iCSRmat *AT)
     *
     * \brief Find transpose of iCSRmat matrix A (index of the array starts with 1!!)
     *
     * \param A   Pointer to the iCSRmat matrix
     * \param AT  Pointer to the transpose of iCSRmat matrix A (output)
     *
     */
    
    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
    
    icsr_trans(A, AT);
    
    INT *AT_IA = AT->IA;
    INT *AT_JA = AT->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<AT->row+1; i++) AT_IA[i] = AT_IA[i]+1;
    for (i=0; i<A->nnz; i++) {
        A_JA[i] = A_JA[i]+1;
        AT_JA[i] = AT_JA[i]+1;
    }
    
}


/***********************************************************************************************/
void dcsr_mxv (dCSRmat *A,
                REAL *x,
                REAL *y)
{
    
    /**
     * \fn void dcsr_mxv (dCSRmat *A, REAL *x, REAL *y)
     *
     * \brief Matrix-vector multiplication y = A*x (index starts with 0!!)
     *
     * \param A   Pointer to dCSRmat matrix A
     * \param x   Pointer to array x
     * \param y   Pointer to array y
     *
     */
    
    const INT m=A->row;
    const INT *ia=A->IA, *ja=A->JA;
    const REAL *aj=A->val;
    INT i, k, begin_row, end_row, nnz_num_row;
    register REAL temp;
    
    for (i=0;i<m;++i) {
        temp=0.0;
        begin_row=ia[i];
        end_row=ia[i+1];
        nnz_num_row = end_row - begin_row;
        switch(nnz_num_row) {
            case 3:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 4:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 5:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 6:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 7:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            default:
                for (k=begin_row; k<end_row; ++k) {
                    temp+=aj[k]*x[ja[k]];
                }
                break;
        }
        
        y[i]=temp;
        
    }
}



/***********************************************************************************************/
void dcsr_mxm (dCSRmat *A,
               dCSRmat *B,
               dCSRmat *C)
{
    
    /**
     * \fn void dcsr_mxm (dCSRmat *A, dCSRmat *B, dCSRmat *C)
     *
     * \brief Sparse matrix multiplication C=A*B (index starts with 0!!)
     *
     * \param A   Pointer to the dCSRmat matrix A
     * \param B   Pointer to the dCSRmat matrix B
     * \param C   Pointer to dCSRmat matrix equal to A*B
     *
     */
    
    INT i,j,k,l,count;
    
    INT *JD = (INT *)calloc(B->col,sizeof(INT));
    
    C->row=A->row;
    C->col=B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) JD[i]=-1;
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==count) {
                    JD[count]=B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i+1]=count;
        for (j=0;j<count;++j) {
            JD[j]=-1;
        }
    }
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==countJD) {
                    C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }
        
        //for (j=0;j<countJD;++j) JD[j]=-1;
        memset(JD, -1, sizeof(INT)*countJD);
    }
    
    free(JD);
    
    // step 3: Find the structure A of C
    C->val=(REAL*)calloc(C->IA[C->row],sizeof(REAL));
    
    for (i=0;i<C->row;++i) {
        for (j=C->IA[i];j<C->IA[i+1];++j) {
            C->val[j]=0;
            for (k=A->IA[i];k<A->IA[i+1];++k) {
                for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
                    if (B->JA[l]==C->JA[j]) {
                        C->val[j]+=A->val[k]*B->val[l];
                    } // end if
                } // end for l
            } // end for k
        } // end for j
    }    // end for i
    
    C->nnz = C->IA[C->row]-C->IA[0];
    
}

/***********************************************************************************************/
void dcsr_mxm_1 (dCSRmat *A,
                 dCSRmat *B,
                 dCSRmat *C)
{
    
    /**
     * \fn void dcsr_mxm_1 (dCSRmat *A, dCSRmat *B, dCSRmat *C)
     *
     * \brief Sparse matrix multiplication C=A*B (index starts with 1!!)
     *
     * \param A   Pointer to the dCSRmat matrix A
     * \param B   Pointer to the dCSRmat matrix B
     * \param C   Pointer to dCSRmat matrix equal to A*B
     *
     */
    
    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    INT *B_IA = B->IA;
    INT *B_JA = B->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]-1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]-1;
    
    dcsr_mxm(A, B, C);
    
    INT *C_IA = C->IA;
    INT *C_JA = C->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]+1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]+1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]+1;
    
    for (i=0; i<C->row+1; i++) C_IA[i] = C_IA[i]+1;
    for (i=0; i<C->nnz; i++) C_JA[i] = C_JA[i]+1;
    
}

/***********************************************************************************************/
void icsr_mxm (iCSRmat *A,
               iCSRmat *B,
               iCSRmat *C)
{
    
    /**
     * \fn void icsr_mxm (iCSRmat *A, iCSRmat *B, iCSRmat *C)
     *
     * \brief Sparse matrix multiplication C=A*B (index starts with 0!!)
     *
     * \param A   Pointer to the iCSRmat matrix A
     * \param B   Pointer to the iCSRmat matrix B
     * \param C   Pointer to iCSRmat matrix equal to A*B
     *
     */
    
    INT i,j,k,l,count;
    
    INT *JD = (INT *)calloc(B->col,sizeof(INT));
    
    C->row = A->row;
    C->col = B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) JD[i]=-1;
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==count) {
                    JD[count]=B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i+1]=count;
        for (j=0;j<count;++j) {
            JD[j]=-1;
        }
    }
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==countJD) {
                    C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }
        
        //for (j=0;j<countJD;++j) JD[j]=-1;
        memset(JD, -1, sizeof(INT)*countJD);
    }
    
    free(JD);
    
    // step 3: Find the structure A of C
    C->val=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        for (j=C->IA[i];j<C->IA[i+1];++j) {
            C->val[j]=0;
            for (k=A->IA[i];k<A->IA[i+1];++k) {
                for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
                    if (B->JA[l]==C->JA[j]) {
                        C->val[j]+=A->val[k]*B->val[l];
                    } // end if
                } // end for l
            } // end for k
        } // end for j
    }    // end for i
    
    C->nnz = C->IA[C->row]-C->IA[0];
    
    
}

/***********************************************************************************************/
void icsr_mxm_1 (iCSRmat *A,
               iCSRmat *B,
               iCSRmat *C)
{
    
    /**
     * \fn void icsr_mxm_1 (iCSRmat *A, iCSRmat *B, iCSRmat *C)
     *
     * \brief Sparse matrix multiplication C=A*B (index starts with 1!!)
     *
     * \param A   Pointer to the iCSRmat matrix A
     * \param B   Pointer to the iCSRmat matrix B
     * \param C   Pointer to iCSRmat matrix equal to A*B
     *
     */
    
    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    INT *B_IA = B->IA;
    INT *B_JA = B->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]-1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]-1;
    
    icsr_mxm(A, B, C);
    
    INT *C_IA = C->IA;
    INT *C_JA = C->JA;
            
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]+1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]+1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]+1;
    
    for (i=0; i<C->row+1; i++) C_IA[i] = C_IA[i]+1;
    for (i=0; i<C->nnz; i++) C_JA[i] = C_JA[i]+1;
    
}

/***********************************************************************************************/
void icsr_mxm_symb (iCSRmat *A,
                    iCSRmat *B,
                    iCSRmat *C)
{
    
    /**
     * \fn void icsr_mxm_symb (iCSRmat *A, iCSRmat *B, iCSRmat *C)
     *
     * \brief Symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
     *
     * \param A   Pointer to the iCSRmat matrix A
     * \param B   Pointer to the iCSRmat matrix B
     * \param C   Pointer to iCSRmat matrix equal to A*B
     *
     */
    
    INT i,j,k,l,count;
    
    INT *JD = (INT *)calloc(B->col,sizeof(INT));
    
    C->row = A->row;
    C->col = B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) JD[i]=-1;
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==count) {
                    JD[count]=B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i+1]=count;
        for (j=0;j<count;++j) {
            JD[j]=-1;
        }
    }
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==countJD) {
                    C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }
        
        //for (j=0;j<countJD;++j) JD[j]=-1;
        memset(JD, -1, sizeof(INT)*countJD);
    }
    
    free(JD);
    
    C->nnz = C->IA[C->row]-C->IA[0];
    
}

/***********************************************************************************************/
void icsr_mxm_symb_1 (iCSRmat *A,
                 iCSRmat *B,
                 iCSRmat *C)
{
    
    /**
     * \fn void icsr_mxm_symb_1 (iCSRmat *A, iCSRmat *B, iCSRmat *C)
     *
     * \brief Symbolic sparse matrix multiplication C=A*B (index starts with 1!!)
     *
     * \param A   Pointer to the iCSRmat matrix A
     * \param B   Pointer to the iCSRmat matrix B
     * \param C   Pointer to iCSRmat matrix equal to A*B
     *
     */
    
    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    INT *B_IA = B->IA;
    INT *B_JA = B->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]-1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]-1;
    
    icsr_mxm_symb(A, B, C);
    
    INT *C_IA = C->IA;
    INT *C_JA = C->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]+1;
    
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]+1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]+1;
    
    for (i=0; i<C->row+1; i++) C_IA[i] = C_IA[i]+1;
    for (i=0; i<C->nnz; i++) C_JA[i] = C_JA[i]+1;
    
}

/***********************************************************************************************/

/**
 * \fn void icsr_mxm_symb_max (iCSRmat *A, iCSRmat *B, iCSRmat *C, INT multmax)
 *
 * \brief symbolic sparse matrix multiplication C=A*B (index starts with 0!!)
 *
 * \param A         Pointer to the iCSRmat matrix A
 * \param B         Pointer to the iCSRmat matrix B
 * \param C         Pointer to iCSRmat matrix equal to A*B
 * \param multimax  value allowed in the iCSRmat matrix C, any entry that is not equal to multimax will be deleted
 *
 */
void icsr_mxm_symb_max (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax)
{
 
    INT i,j,k,l,count;
    
    INT *JD = (INT *)calloc(B->col,sizeof(INT));
    INT *entry_count = (INT *)calloc(B->col,sizeof(INT));
    
    C->row = A->row;
    C->col = B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) {
        JD[i] = -1;
        entry_count[i] = 0;
    }
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) {
                        entry_count[l] = entry_count[l]+1;
                        break;
                    }
                }
                
                if (l==count) {
                    JD[count]=B->JA[j];
                    entry_count[count] = 1;
                    count++;
                }
            }
        
                
        }
        
        
        C->IA[i+1]=count;
        
        for (j=0;j<count;++j) {
            
            JD[j]=-1;
            
            if (entry_count[j] != multmax) C->IA[i+1] = C->IA[i+1]-1;

            entry_count[j] = 0;
            
        }
        
        
    }
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) {
                        entry_count[l] = entry_count[l]+1;
                        break;
                    }
                }
                
                if (l==countJD) {
                    //C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    entry_count[countJD] = 1;
                    //count++;
                    countJD++;
                }
            }
            
        }
        
        
        for (j=0;j<countJD;++j) {
            
            
            if (entry_count[j] == multmax) {
                C->JA[count]=JD[j];
                count++;
            }
            
            JD[j]=-1;
            entry_count[j] = 0;
            
        }
        
        
        
    }
    
    // free
    free(JD);
    free(entry_count);
    
    C->nnz = C->IA[C->row]-C->IA[0];
}


/*
 void icsr_mxm_symb_max (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax)
{
    printf("max-1\n");
    
    INT *ia = A->IA;
    INT *ja = A->JA;
    INT na = A->row;
    //INT mab = A->col;
    INT mb = B->col;
    INT *ib = B->IA;
    INT *jb = B->JA;
    
    C->row = A->row;
    C->col = B->col;
    C->val = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    C->JA  = NULL;
    INT *ic = C->IA;
    //INT *jc = C->JA;
    
    INT i,jk,icp,icp_temp,iaa,iab,ibb,iba,jbk,j,if1; // Counters and such
    INT *ix=NULL;
    INT *col=NULL;
    INT jck=0;
    ix = calloc(mb,sizeof(INT));
    col = calloc(mb,sizeof(INT));
	
    printf("max-2\n");
    
    for(i = 0;i<mb;i++) {
        ix[i]=0;
        col[i]=0;
    }
    icp = 0;
    
    printf("max-3\n");
    
    // first loop to figure out the structure of C->JA
    for(i = 1;i<=na;i++) {
        ic[i-1]=icp;
        icp_temp = icp;
        iaa=ia[i-1];
        iab=ia[i]-1;
        if(iab>=iaa) {
            for(jk = iaa;jk<=iab;jk++){
                if1 = ja[jk-1];
                iba = ib[if1-1];
                ibb = ib[if1]-1;
                if(ibb>=iba) {
                    for (jbk = iba;jbk<=ibb;jbk++){
                        j = jb[jbk-1];
                        col[j-1]++;
                        if(ix[j-1] != i) {
                            ix[j-1]=i;
                            //jc[icp_temp-1] = j;
                            icp_temp++;
                        }
                    }
                }
            }
            for (jck=ic[i-1];jck<icp_temp;jck++){
                j = jc[jck-1];
                if(col[j-1] == multmax || (!multmax)) {
                    col[j-1]=0;
                    jc[icp-1] = j;
                    icp++;
                } else {
                    col[j-1]=0;
                }
            }
        }
    }
    
    printf("max-3\n");
    
    printf("ic[%d]=%d\n",0, ic[0]);
    
    ic[na] = icp;
    
    printf("max-4\n");

    if(ix) free(ix);
    
    printf("max-5\n");

    if(col) free(col);
    
    printf("max-6\n");

    
    return;
}
 */


/***********************************************************************************************/
/**
 * \fn void icsr_mxm_symb_max (iCSRmat *A, iCSRmat *B, iCSRmat *C INT multmax)
 *
 * \brief symbolic sparse matrix multiplication C=A*B (index starts with 1!!)
 *
 * \param A         Pointer to the iCSRmat matrix A
 * \param B         Pointer to the iCSRmat matrix B
 * \param C         Pointer to iCSRmat matrix equal to A*B
 * \param multimax  value allowed in the iCSRmat matrix C, any entry that is not equal to multimax will be deleted
 *
 */
void icsr_mxm_symb_max_1 (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax)
{

    INT i;
    
    INT *A_IA = A->IA;
    INT *A_JA = A->JA;
    INT *B_IA = B->IA;
    INT *B_JA = B->JA;
    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]-1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]-1;
            
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]-1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]-1;
    
    icsr_mxm_symb_max(A, B, C, multmax);

    
    INT *C_IA = C->IA;
    INT *C_JA = C->JA;
                    
    for (i=0; i<A->row+1; i++) A_IA[i] = A_IA[i]+1;
    for (i=0; i<A->nnz; i++) A_JA[i] = A_JA[i]+1;
                            
    for (i=0; i<B->row+1; i++) B_IA[i] = B_IA[i]+1;
    for (i=0; i<B->nnz; i++) B_JA[i] = B_JA[i]+1;
                                    
    for (i=0; i<C->row+1; i++) C_IA[i] = C_IA[i]+1;
    for (i=0; i<C->nnz; i++) C_JA[i] = C_JA[i]+1;
                                            
}

